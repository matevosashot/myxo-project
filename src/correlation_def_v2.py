import argparse
import os
from termios import PARODD
import numpy as np
import numexpr as ne
from .utils import log
from .correlation_def import Gcalc


class GcalcForQ(Gcalc):
    """
    Gcalc variant for the Q-side of FourierSolver. Replaces v1's Gaussian
    phi(r1, r2) by an exponential kernel  phi_v2(r1, r2) = 2 * exp(-|r2-r1|/lM)
    computed inline in calc_C_Q_index, calc_C_Q_index_symmetrized, and
    calc_C_Q. Since phi_v2 is symmetric under r1<->r2 and the P-product
    structure of the v2 expression is too, c[r1,r2] == c[r2,r1].
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._el = self.lm

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
        el = self.dtype.type(self._el)
        dx, dy = self._periodic_displacements()

        out = ne.evaluate(
            "two * exp(-sqrt(dx * dx + dy * dy) / el)"
            " * (P_am * P_bn + P_an * P_bm + PT_am * PT_bn + PT_an * PT_bm)",
            local_dict={
                'two': self.dtype.type(2),
                'el': el,
                'dx': dx, 'dy': dy,
                'P_am': P[..., a, m], 'P_bn': P[..., b, n],
                'P_an': P[..., a, n], 'P_bm': P[..., b, m],
                'PT_am': PT[..., a, m], 'PT_bn': PT[..., b, n],
                'PT_an': PT[..., a, n], 'PT_bm': PT[..., b, m],
            },
            out=out)

        self._validate_dtype(out)
        return out

    def calc_C_Q_index_symmetrized(self, a, b, m, n, out=None):
        """
        c[r1,r2] + c[r2,r1] for the v2 form. Since phi_v2 is symmetric
        under r1<->r2 and so is the P-product piece, c[r1,r2] == c[r2,r1]
        and the symmetrised value reduces to 2*c[r1,r2] -- i.e. the same
        expression with the front factor doubled from 2 to 4.
        """
        P = self.P[:, :, None, None]
        PT = P.transpose(2, 3, 0, 1, 4, 5)
        el = self.dtype.type(self._el)
        dx, dy = self._periodic_displacements()

        out = ne.evaluate(
            "four * exp(-sqrt(dx * dx + dy * dy) / el)"
            " * (P_am * P_bn + P_an * P_bm + PT_am * PT_bn + PT_an * PT_bm)",
            local_dict={
                'four': self.dtype.type(4),
                'el': el,
                'dx': dx, 'dy': dy,
                'P_am': P[..., a, m], 'P_bn': P[..., b, n],
                'P_an': P[..., a, n], 'P_bm': P[..., b, m],
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
        """
        el = self.dtype.type(self._el)
        dx, dy = self._periodic_displacements()

        # Add four singleton axes to dx/dy so they broadcast against the
        # 2x2x2x2 outer-product structure in axes 4..7.
        dx_b = dx[..., None, None, None, None]
        dy_b = dy[..., None, None, None, None]

        # P with r2 broadcasting axes baked in, then PT as its swap.
        P = self.P[:, :, None, None, :, :]                 # (N, N, 1, 1, 2, 2)
        PT = P.transpose(2, 3, 0, 1, 4, 5)                  # (1, 1, N, N, 2, 2)

        # Outer products over the trailing (i, j) axes of P / PT, placing
        # the resulting (a, b, m, n) indices in axes 4..7.
        # P_outer_r1_ambn[..., a, b, m, n] = P(r1)_{a,m} * P(r1)_{b,n}
        P_outer_r1_ambn = ne.evaluate(
            "P_am * P_bn",
            local_dict={
                'P_am': P[..., :, None, :, None],   # axes 4=a, 6=m
                'P_bn': P[..., None, :, None, :],   # axes 5=b, 7=n
            })
        P_outer_r1_anbm = ne.evaluate(
            "P_an * P_bm",
            local_dict={
                'P_an': P[..., :, None, None, :],   # axes 4=a, 7=n
                'P_bm': P[..., None, :, :, None],   # axes 5=b, 6=m
            })
        PT_outer_r2_ambn = ne.evaluate(
            "PT_am * PT_bn",
            local_dict={
                'PT_am': PT[..., :, None, :, None],
                'PT_bn': PT[..., None, :, None, :],
            })
        PT_outer_r2_anbm = ne.evaluate(
            "PT_an * PT_bm",
            local_dict={
                'PT_an': PT[..., :, None, None, :],
                'PT_bm': PT[..., None, :, :, None],
            })

        value = ne.evaluate(
            "two * exp(-sqrt(dx_b * dx_b + dy_b * dy_b) / el)"
            " * (P_r1_ambn + P_r1_anbm + PT_r2_ambn + PT_r2_anbm)",
            local_dict={
                'two': self.dtype.type(2),
                'el': el,
                'dx_b': dx_b, 'dy_b': dy_b,
                'P_r1_ambn': P_outer_r1_ambn,
                'P_r1_anbm': P_outer_r1_anbm,
                'PT_r2_ambn': PT_outer_r2_ambn,
                'PT_r2_anbm': PT_outer_r2_anbm,
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
