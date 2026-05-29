import argparse
import os
from termios import PARODD
import numpy as np
import numexpr as ne
from .utils import log
from .correlation_def import Gcalc, S as _S_field


class GcalcForQ(Gcalc):
    """
    Gcalc variant for the Q-side of FourierSolver. Extends v2.1 by an
    additional scalar prefactor S(r1) * S(r2):

        c[r1,r2; a,b,m,n] = G(r1,r2) * S(r1) * S(r2) * (
              P_ab*P_mn   + P_am*P_bn   + P_an*P_bm   - P_ab*PT_mn
            + PT_ab*PT_mn + PT_am*PT_bn + PT_an*PT_bm - P_ab*PT_mn
        )

    where G(r1,r2) = exp(-|r2-r1|^2 / (2 * lm^2)), S(r) is the scalar
    nematic-strength field from correlation_def.S (precomputed on the
    grid into self.S by the overridden calcQ), P_xy = P(r1)_{x,y}, and
    PT_xy = P(r2)_{x,y}. The G and S*ST prefactors are symmetric under
    r1<->r2, but the two `-P_ab*PT_mn` cross-terms break r1<->r2 symmetry
    of the full expression: in general c[r1,r2] != c[r2,r1].
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._el = self.lm

    def calcQ(self):
        # Parent populates self.Q; we additionally stash the scalar
        # field S(r) on the (N, N) grid so calc_C_Q_index et al. can pick
        # it up as `self.S` without recomputing it on every call.
        super().calcQ()
        x = self.grid[:, None]
        y = self.grid[None, :]
        r = np.sqrt(x * x + y * y)
        self.S = _S_field(r, self.S0, self.l, self.R).astype(self.dtype)
        self._validate_dtype(self.S)

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
        P = self.P[:, :, None, None]                      # N x N x 1 x 1 x 2 x 2
        PT = P.transpose(2, 3, 0, 1, 4, 5)                # 1 x 1 x N x N x 2 x 2
        S = self.S[:, :, None, None]                      # N x N x 1 x 1   (scalar field)
        ST = S.transpose(2, 3, 0, 1)                      # 1 x 1 x N x N   (scalar field)
        el = self.dtype.type(self._el)
        dx, dy = self._periodic_displacements()

        out = ne.evaluate(
            "exp(-(dx * dx + dy * dy) / two_el_squared) * S * ST * "
            "( (P_ab  * P_mn  + P_am  * P_bn  + P_an  * P_bm  - P_ab * PT_mn) + "
            "(PT_ab * PT_mn + PT_am * PT_bn + PT_an * PT_bm - P_ab * PT_mn) )",
            local_dict={
                'two_el_squared': 2 * el * el,
                'dx': dx, 'dy': dy,
                'S': S, 'ST': ST,
                'P_ab': P[..., a, b], 'P_mn': P[..., m, n],
                'P_am': P[..., a, m], 'P_bn': P[..., b, n],
                'P_an': P[..., a, n], 'P_bm': P[..., b, m],
                'PT_ab': PT[..., a, b], 'PT_mn': PT[..., m, n],
                'PT_am': PT[..., a, m], 'PT_bn': PT[..., b, n],
                'PT_an': PT[..., a, n], 'PT_bm': PT[..., b, m],
            },
            out=out)

        self._validate_dtype(out)
        return out

    def calc_C_Q_index_symmetrized(self, a, b, m, n, out=None):
        """
        c[r1,r2] + c[r2,r1] for the v2.2 form. Swapping r1<->r2 swaps
        P<->PT (and leaves both G and S(r1)*S(r2) invariant). Adding the
        swapped expression to the original doubles the symmetric P-product
        terms and produces the symmetric cross-term `P_ab*PT_mn + PT_ab*P_mn`:

            c[r1,r2] + c[r2,r1] = 2 * G * S * ST * (
                  P_ab*P_mn  + P_am*P_bn  + P_an*P_bm
                + PT_ab*PT_mn + PT_am*PT_bn + PT_an*PT_bm
                - P_ab*PT_mn - PT_ab*P_mn
            )
        """
        P = self.P[:, :, None, None]
        PT = P.transpose(2, 3, 0, 1, 4, 5)
        S = self.S[:, :, None, None]
        ST = S.transpose(2, 3, 0, 1)
        el = self.dtype.type(self._el)
        dx, dy = self._periodic_displacements()

        out = ne.evaluate(
            "two * exp(-(dx * dx + dy * dy) / two_el_squared) * S * ST"
            " * (P_ab * P_mn + P_am * P_bn + P_an * P_bm"
            "    + PT_ab * PT_mn + PT_am * PT_bn + PT_an * PT_bm"
            "    - P_ab * PT_mn - PT_ab * P_mn)",
            local_dict={
                'two': self.dtype.type(2),
                'two_el_squared': 2 * el * el,
                'dx': dx, 'dy': dy,
                'S': S, 'ST': ST,
                'P_ab': P[..., a, b], 'P_mn': P[..., m, n],
                'P_am': P[..., a, m], 'P_bn': P[..., b, n],
                'P_an': P[..., a, n], 'P_bm': P[..., b, m],
                'PT_ab': PT[..., a, b], 'PT_mn': PT[..., m, n],
                'PT_am': PT[..., a, m], 'PT_bn': PT[..., b, n],
                'PT_an': PT[..., a, n], 'PT_bm': PT[..., b, m],
            },
            out=out)

        self._validate_dtype(out)
        return out

    def calc_C_Q(self, a=None, b=None, m=None, n=None):
        """
        Full rank-4 C_Q tensor of shape (N, N, N, N, 2, 2, 2, 2). The
        (a, b, m, n) keyword arguments are accepted for signature parity
        with v1 but ignored -- this method always returns the full tensor.
        Trailing axes are (a, b, m, n) so that out[..., a, b, m, n] equals
        calc_C_Q_index(a, b, m, n) element-wise.

        Assembled via the same outer-product/transpose split as v1's
        calc_C_Q: build the r1-only `part_left` (which carries the
        `(P_ab*P_mn + P_am*P_bn + P_an*P_bm - P_ab*PT_mn)` half), then
        get the r2-only `part_right` by index-transposing it; sum.
        The leading `half` of v1 is dropped because here phi == phi_T = G.
        The S(r1)*S(r2) scalar prefactor is baked into the kernel
        H = G * S * ST; H is symmetric under r1<->r2 so the transpose
        identity still gives the correct r2-half.
        """
        el = self.dtype.type(self._el)
        dx, dy = self._periodic_displacements()

        # H(r1, r2) = G(r1, r2) * S(r1) * S(r2), symmetric under r1<->r2.
        S = self.S[:, :, None, None]                       # (N, N, 1, 1)
        ST = self.S[None, None, :, :]                      # (1, 1, N, N)
        H = ne.evaluate(
            "exp(-(dx * dx + dy * dy) / two_el_squared) * S * ST",
            local_dict={
                'two_el_squared': 2 * el * el,
                'dx': dx, 'dy': dy,
                'S': S, 'ST': ST,
            })  # (N, N, N, N)

        P = self.P[:, :, None, None]                       # (N, N, 1, 1, 2, 2)
        P_transpose = P.transpose(2, 3, 0, 1, 4, 5)         # (1, 1, N, N, 2, 2)

        # P_outer_rr_abmn[..., a, b, m, n] = P(r1)_{a,b} * P(r1)_{m,n}
        P_outer_rr_abmn = ne.evaluate(
            "P_left * P_right",
            local_dict={
                'P_left': P[..., :, :, None, None],
                'P_right': P[..., None, None, :, :],
            })  # (N, N, 1, 1, 2, 2, 2, 2)
        P_outer_rr_ambn = P_outer_rr_abmn.transpose(0, 1, 2, 3, 4, 6, 5, 7)
        P_outer_rr_anbm = P_outer_rr_abmn.transpose(0, 1, 2, 3, 4, 6, 7, 5)

        # Cross term P(r1)_{a,b} * P(r2)_{m,n}, full (N,N,N,N) on the grid axes.
        P_outer_rr1_abmn = ne.evaluate(
            "P_left * P_right",
            local_dict={
                'P_left': P[..., :, :, None, None],
                'P_right': P_transpose[..., None, None, :, :],
            })  # (N, N, N, N, 2, 2, 2, 2)

        part_left = ne.evaluate(
            "H * (P_outer_rr_abmn + P_outer_rr_ambn + P_outer_rr_anbm - P_outer_rr1_abmn)",
            local_dict={
                'H': H[..., None, None, None, None],
                'P_outer_rr_abmn': P_outer_rr_abmn,
                'P_outer_rr_ambn': P_outer_rr_ambn,
                'P_outer_rr_anbm': P_outer_rr_anbm,
                'P_outer_rr1_abmn': P_outer_rr1_abmn,
            })  # (N, N, N, N, 2, 2, 2, 2)
        # part_right = part_left evaluated at swapped (r1, r2) and swapped
        # (a, b) <-> (m, n): grid axes 0,1,2,3 -> 2,3,0,1 ; index axes
        # 4,5,6,7 -> 6,7,4,5. P is symmetric so the index swap matches the
        # r2-half of calc_C_Q_index term-for-term.
        part_right = part_left.transpose(2, 3, 0, 1, 6, 7, 4, 5)

        value = ne.evaluate(
            "part_left + part_right",
            local_dict={'part_left': part_left, 'part_right': part_right})

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
