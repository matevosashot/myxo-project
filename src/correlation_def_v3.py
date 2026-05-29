import argparse
import os
from termios import PARODD
import numpy as np
import numexpr as ne
from .utils import log
from .correlation_def import Gcalc, S as _S_field


class GcalcForQ(Gcalc):
    """
    Gcalc variant for the Q-side of FourierSolver. Replaces v2.2's Gaussian
    `G = exp(-|r2-r1|^2 / (2*lm^2))` and the P-product structure by a
    wrapped-Gaussian / nematic-axis form built from per-grid fields

        S(r)              = scalar nematic-strength field
                            (correlation_def.S, stashed by calcQ)
        sigma(r)          = sqrt(-2 * log max(S(r), eps))
                            -- wrapped-Gaussian standard deviation
                            (stored as `self.varphi_sigma`, since
                            `Gcalc.__init__` already binds `self.sigma`
                            to a scalar covariance knob)
        varphi(r)         = arctan2(y, x)
                            -- nematic-axis identifier, varphi = 2*theta

    and the exponential (NOT Gaussian) inter-point kernel

        G(r1, r2)         = exp(-|r2 - r1| / lm).

    All four (a == b, m == n) index cases collapse to one master expression

        c[r1,r2; a,b,m,n] = factor * factor_T * S(r1) * S(r2) * (
              T_ab(varphi(r1)) * T_mn(varphi(r2)) * cosh(arg)
              + sign * U_ab(varphi(r1)) * U_mn(varphi(r2)) * sinh(arg)
        )

    with arg = G(r1, r2) * sigma(r1) * sigma(r2) and the (a, b, m, n)-only
    constants

        T_ab(x)  = cos(x) if a == b else sin(x)
        U_ab(x)  = sin(x) if a == b else cos(x)
        sign     = +1 if (a == b) == (m == n) else -1
        factor   = -1 if (a, b) == (1, 1) else 1
        factor_T = -1 if (m, n) == (1, 1) else 1.

    Symmetry under r1 <-> r2: the S, exp(-|r2-r1|/lm), and sigma*sigma_T
    pieces are symmetric, and factor*factor_T does not depend on (r1, r2).
    Swapping r1 <-> r2 only swaps varphi <-> varphi_T inside the cos/sin
    pairings:

      - "coscos" (a == b, m == n) and "sinsin" (a != b, m != n), i.e. the
        sign = +1 cases, are symmetric: c[r1,r2] == c[r2,r1].
      - "cossin" (a == b, m != n) and "sincos" (a != b, m == n), i.e. the
        sign = -1 cases, swap into each other under r1 <-> r2; the
        symmetrised value collapses via cosh - sinh = exp(-arg) into a
        single sin(varphi + varphi_T) * exp(-arg) product (handled by
        `calc_C_Q_index_symmetrized` and `calc_C_Q` via the unified
        expression above).
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._el = self.lm

    def calcQ(self):
        """
        Populate self.Q (via the parent) and stash the v3 wrapped-Gaussian
        / nematic-axis fields on the (N, N) grid so calc_C_Q_index et al.
        can pick them up without recomputing per call:

            self.S            : S(r)       = correlation_def.S(r, S0, l, R)
            self.varphi_sigma : sigma(r)   = sqrt(-2 * log max(S, eps))
            self.varphi       : varphi(r)  = arctan2(y, x).

        The eps-clamp on `S` keeps `log` finite at r == 0 (where the
        nematic-strength field smoothly vanishes) and at r >> R (where the
        Gaussian cutoff drives `S` toward 0). All three fields are cast to
        `self.dtype` so the parent's `_validate_dtype` accepts them.
        """
        super().calcQ()
        x = self.grid[:, None]
        y = self.grid[None, :]
        r = np.sqrt(x * x + y * y)
        eps = self.dtype.type(1e-15)
        self.S = _S_field(r, self.S0, self.l, self.R).astype(self.dtype)
        self.varphi_sigma = np.sqrt(
            -2 * np.log(np.maximum(self.S, eps))
        ).astype(self.dtype)
        self.varphi = np.arctan2(y, x).astype(self.dtype)

        self._validate_dtype(self.S)
        self._validate_dtype(self.varphi_sigma)
        self._validate_dtype(self.varphi)

    def calc_phi(self):
        # v2 computes the exp kernel inline; there's no stored phi buffer.
        # Set a non-None sentinel so FourierSolver._ensure_phi treats this
        # Gcalc as "phi already materialised" and skips the no-op rebuild.
        self.phi = True

    def drop_phi(self):
        # No buffer to release; keep the sentinel so _ensure_phi continues
        # to see phi as "materialised".
        pass

    def _periodic_displacements(self):
        dx = self.x2 - self.x1
        dx = dx - self.L_grid * np.round(dx / self.L_grid)
        dy = self.y2 - self.y1
        dy = dy - self.L_grid * np.round(dy / self.L_grid)
        return dx, dy

    def calc_C_Q_index(self, a, b, m, n, out=None):
        """
        Single (a, b, m, n) slice of the C_Q tensor on the (r1, r2) grid;
        shape (N, N, N, N). Dispatches on (a == b, m == n) by selecting one
        of four (cos/sin, cos/sin) pair templates and substituting the
        exponential-kernel expression `c231 = exp(-|r2-r1|/lm)` for the
        inter-point coupling. See the class docstring for the closed form.
        """
        S = self.S[:, :, None, None]                      # (N, N, 1, 1)
        ST = S.transpose(2, 3, 0, 1)                      # (1, 1, N, N)
        sigma = self.varphi_sigma[:, :, None, None]       # (N, N, 1, 1)
        sigma_T = sigma.transpose(2, 3, 0, 1)             # (1, 1, N, N)
        varphi = self.varphi[:, :, None, None]            # (N, N, 1, 1)
        varphi_T = varphi.transpose(2, 3, 0, 1)           # (1, 1, N, N)
        el = self.dtype.type(self._el)
        dx, dy = self._periodic_displacements()

        exprs = {
            "coscos": "S*ST*(cos(j)*cos(jT)*cosh(c231*sigma*sigma_T) + sin(j)*sin(jT)*sinh(c231*sigma*sigma_T))",
            "cossin": "S*ST*(cos(j)*sin(jT)*cosh(c231*sigma*sigma_T) - sin(j)*cos(jT)*sinh(c231*sigma*sigma_T))",
            "sincos": "S*ST*(sin(j)*cos(jT)*cosh(c231*sigma*sigma_T) - cos(j)*sin(jT)*sinh(c231*sigma*sigma_T))",
            "sinsin": "S*ST*(sin(j)*sin(jT)*cosh(c231*sigma*sigma_T) + cos(j)*cos(jT)*sinh(c231*sigma*sigma_T))",
        }
        c231 = "exp(-sqrt(dx * dx + dy * dy) / ell)"

        factor = self.dtype.type(-1 if a == 1 and b == 1 else 1)
        factor_T = self.dtype.type(-1 if m == 1 and n == 1 else 1)
        expr = ("cos" if a == b else "sin") + ("cos" if m == n else "sin")

        to_evaluate = "factor * factor_T * " + exprs[expr].replace("c231", c231)

        out = ne.evaluate(
            to_evaluate,
            local_dict={
                'ell': el,
                'dx': dx, 'dy': dy,
                'S': S, 'ST': ST,
                'sigma': sigma, 'sigma_T': sigma_T,
                'j': varphi, 'jT': varphi_T,
                'factor': factor, 'factor_T': factor_T,
            },
            out=out)

        self._validate_dtype(out)
        return out

    def calc_C_Q_index_symmetrized(self, a, b, m, n, out=None):
        """
        c[r1,r2] + c[r2,r1] for the v3 form, computed in a single
        ne.evaluate without materialising the un-symmetrised intermediate.
        Equivalent to `c + c.transpose(2, 3, 0, 1)` where
        `c = calc_C_Q_index(a, b, m, n)`.

        Derivation: swapping r1 <-> r2 only swaps varphi <-> varphi_T
        inside the cos/sin pairings (the S, exp(-|r2-r1|/lm), and
        sigma*sigma_T pieces are themselves r1<->r2 symmetric, and
        factor*factor_T does not depend on (r1, r2)). Writing the
        unified form

            c[r1,r2] = factor*factor_T * S*ST * (
                  T_ab(j)  * T_mn(jT) * cosh(arg)
                + sign     * U_ab(j)  * U_mn(jT) * sinh(arg)
            )

        with arg = exp(-|r2-r1|/lm) * sigma * sigma_T and the (a,b,m,n)
        constants T, U, sign from the class docstring, adding the
        swapped expression yields a single closed form valid for all
        four (a == b, m == n) cases:

            c + c.T = factor*factor_T * S*ST * (
                  ( T_ab(j) * T_mn(jT) + T_ab(jT) * T_mn(j) ) * cosh(arg)
                + sign * ( U_ab(j) * U_mn(jT) + U_ab(jT) * U_mn(j) ) * sinh(arg)
            ).

        Sanity checks of the four cases:
          - "coscos" / "sinsin"  (sign = +1): the pairings collapse to
            2*T_ab(j)*T_mn(jT) and 2*U_ab(j)*U_mn(jT), recovering 2*c.
          - "cossin" / "sincos"  (sign = -1): both pairings collapse to
            sin(j + jT), and cosh - sinh = exp(-arg) gives the closed form
                c + c.T = factor*factor_T*S*ST * sin(j + jT)
                          * exp(-exp(-|r2-r1|/lm) * sigma * sigma_T).
        """
        S = self.S[:, :, None, None]
        ST = S.transpose(2, 3, 0, 1)
        sigma = self.varphi_sigma[:, :, None, None]
        sigma_T = sigma.transpose(2, 3, 0, 1)
        varphi = self.varphi[:, :, None, None]
        varphi_T = varphi.transpose(2, 3, 0, 1)
        el = self.dtype.type(self._el)
        dx, dy = self._periodic_displacements()

        factor = self.dtype.type(-1 if a == 1 and b == 1 else 1)
        factor_T = self.dtype.type(-1 if m == 1 and n == 1 else 1)
        # `sign` collides with numexpr's builtin sign(x); use `sgn` instead.
        sgn = self.dtype.type(1 if (a == b) == (m == n) else -1)
        T_ab = "cos" if a == b else "sin"
        T_mn = "cos" if m == n else "sin"
        U_ab = "sin" if a == b else "cos"
        U_mn = "sin" if m == n else "cos"

        c231 = "exp(-sqrt(dx * dx + dy * dy) / ell)"
        expr = (
            "factor * factor_T * S * ST * ("
            f"({T_ab}(j) * {T_mn}(jT) + {T_ab}(jT) * {T_mn}(j))"
            " * cosh(c231 * sigma * sigma_T)"
            f" + sgn * ({U_ab}(j) * {U_mn}(jT) + {U_ab}(jT) * {U_mn}(j))"
            " * sinh(c231 * sigma * sigma_T))"
        )

        out = ne.evaluate(
            expr.replace("c231", c231),
            local_dict={
                'ell': el,
                'dx': dx, 'dy': dy,
                'S': S, 'ST': ST,
                'sigma': sigma, 'sigma_T': sigma_T,
                'j': varphi, 'jT': varphi_T,
                'factor': factor, 'factor_T': factor_T,
                'sgn': sgn,
            },
            out=out)

        self._validate_dtype(out)
        return out

    def calc_C_Q(self, a=None, b=None, m=None, n=None):
        """
        Full rank-4 C_Q tensor of shape (N, N, N, N, 2, 2, 2, 2). The
        (a, b, m, n) keyword arguments are accepted for signature parity
        with v1 / v2.x but ignored -- this method always returns the full
        tensor. Trailing axes are (a, b, m, n) so that
        out[..., a, b, m, n] == calc_C_Q_index(a, b, m, n).

        Vectorised assembly of the master expression in the class
        docstring. The four (a == b, m == n) cases are unified through
        small (2, 2) and (2, 2, 2, 2) constant arrays:

            M_diag[a, b] = 1 if a == b else 0
            M_off[a, b]  = 1 - M_diag[a, b]
            s_arr[a, b]  = +1 if a == b else -1
            factor_arr[a, b] = -1 if (a, b) == (1, 1) else 1
            factor_full[a,b,m,n] = factor_arr[a,b] * factor_arr[m,n]
            sign_full[a,b,m,n]   = s_arr[a,b]      * s_arr[m,n]

        and the per-grid trig fields

            T_grid[i, j, a, b] = cos(varphi[i, j]) if a == b else sin(varphi[i, j])
            U_grid[i, j, a, b] = sin(varphi[i, j]) if a == b else cos(varphi[i, j]).

        Per-grid kernels (cosh, sinh of arg = G * sigma * sigma_T, and the
        S * ST prefactor) are built once on the (N, N, N, N) axes, then
        broadcast against T_grid / U_grid (placed on r1 and r2 grid axes
        respectively) in a single numexpr call.
        """
        el = self.dtype.type(self._el)
        dx, dy = self._periodic_displacements()

        sigma = self.varphi_sigma[:, :, None, None]    # (N, N, 1, 1)
        sigma_T = self.varphi_sigma[None, None, :, :]  # (1, 1, N, N)
        S = self.S[:, :, None, None]                    # (N, N, 1, 1)
        ST = self.S[None, None, :, :]                   # (1, 1, N, N)

        # arg(r1, r2) = exp(-|r2 - r1| / lm) * sigma(r1) * sigma(r2)
        arg = ne.evaluate(
            "exp(-sqrt(dx * dx + dy * dy) / ell) * sigma * sigma_T",
            local_dict={
                'ell': el,
                'dx': dx, 'dy': dy,
                'sigma': sigma, 'sigma_T': sigma_T,
            })  # (N, N, N, N)
        cosh_term = ne.evaluate("cosh(arg)", local_dict={'arg': arg})
        sinh_term = ne.evaluate("sinh(arg)", local_dict={'arg': arg})
        pref_grid = ne.evaluate(
            "S * ST", local_dict={'S': S, 'ST': ST})    # (N, N, N, N)

        # 2x2 helpers for the (a, b) and (m, n) index axes.
        M_diag = np.eye(2, dtype=self.dtype)            # (2, 2): 1 on diag
        M_off = self.dtype.type(1) - M_diag             # (2, 2): 1 off-diag
        s_arr = self.dtype.type(2) * M_diag - self.dtype.type(1)  # +1 / -1
        factor_arr = np.ones((2, 2), dtype=self.dtype)
        factor_arr[1, 1] = self.dtype.type(-1)
        factor_full = (factor_arr[:, :, None, None]
                       * factor_arr[None, None, :, :])   # (2, 2, 2, 2)
        sign_full = (s_arr[:, :, None, None]
                     * s_arr[None, None, :, :])          # (2, 2, 2, 2)

        # Per-grid cos / sin of varphi (shared between the r1 and r2 sides).
        cos_j = np.cos(self.varphi).astype(self.dtype)   # (N, N)
        sin_j = np.sin(self.varphi).astype(self.dtype)   # (N, N)
        # T_grid[i, j, a, b] = cos(varphi[i, j]) if a == b else sin(varphi[i, j])
        T_grid = (M_diag[None, None, :, :] * cos_j[:, :, None, None]
                  + M_off[None, None, :, :] * sin_j[:, :, None, None])
        # U_grid is T_grid with cos/sin swapped (the "other" pick).
        U_grid = (M_off[None, None, :, :] * cos_j[:, :, None, None]
                  + M_diag[None, None, :, :] * sin_j[:, :, None, None])

        # Place T_grid / U_grid against the (r1, r2, a, b, m, n) layout:
        # axes 0,1 = r1; 2,3 = r2; 4,5 = a, b; 6,7 = m, n.
        T_ab_j = T_grid[:, :, None, None, :, :, None, None]   # (N,N,1,1,2,2,1,1)
        T_mn_jT = T_grid[None, None, :, :, None, None, :, :]  # (1,1,N,N,1,1,2,2)
        U_ab_j = U_grid[:, :, None, None, :, :, None, None]
        U_mn_jT = U_grid[None, None, :, :, None, None, :, :]

        value = ne.evaluate(
            "factor_full * pref_grid * "
            "(T_ab_j * T_mn_jT * cosh_term"
            " + sign_full * U_ab_j * U_mn_jT * sinh_term)",
            local_dict={
                'factor_full': factor_full[None, None, None, None, :, :, :, :],
                'sign_full': sign_full[None, None, None, None, :, :, :, :],
                'pref_grid': pref_grid[..., None, None, None, None],
                'T_ab_j': T_ab_j, 'T_mn_jT': T_mn_jT,
                'U_ab_j': U_ab_j, 'U_mn_jT': U_mn_jT,
                'cosh_term': cosh_term[..., None, None, None, None],
                'sinh_term': sinh_term[..., None, None, None, None],
            })

        self._validate_dtype(value)
        return value


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Compute G_{ij} on a 4D grid")
    parser.add_argument("--S0",       type=float, default=1.0)
    parser.add_argument("--l",        type=float, default=1.0)
    parser.add_argument("--lM",       type=float, default=7.0/3)
    parser.add_argument("--lm",       type=float, default=0.7/3)
    parser.add_argument("--R",        type=float, default=10, help="Cutoff radius")
    parser.add_argument("--N",        type=int,   default=50, help="Grid size per dimension")
    parser.add_argument("--parallel", action="store_true",     help="Compute G in parallel using threads")
    parser.add_argument("--workers",  type=int,   default=None, help="Number of worker threads (default: CPU count)")
    parser.add_argument("--dtype",    type=str,   default="float64", help="NumPy dtype for computation (e.g. float16, float32, float64)")
    args = parser.parse_args()

    dtype = np.dtype(args.dtype)

    print("O(N^4) memory is {:0.1f}GB".format(args.N**4 / 1024**3 * dtype.itemsize))

    grid_1d = np.linspace(-5, 5, args.N, dtype=dtype)
    gcalc = GcalcForQ(grid_1d, S0=args.S0, l=args.l, lM=args.lM, lm=args.lm, R=args.R, dtype=dtype)
    gcalc.precompute()
    print(gcalc.Sigma.mean())

    value = gcalc.calc_C_Q()
    value_2 = np.empty_like(value)
    for a in range(2):
        for b in range(2):
            for m in range(2):
                for n in range(2):
                    value_2[..., a, b, m, n] = gcalc.calc_C_Q_index(a, b, m, n)

    print("absdiff", np.abs(value - value_2).mean())
    exit(0)
