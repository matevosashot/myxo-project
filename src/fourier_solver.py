import argparse
import logging
import os
import sys

import numexpr as ne
import numpy as np
from scipy.fft import irfftn, rfftn, set_workers

from . import _NWORKERS
from .correlation_def import Gcalc
from .saving_utils import HDF5Saver
from .utils import log



def _build_wavevectors_1d(n: int, dx: float) -> np.ndarray:
    return np.fft.fftfreq(n, d=dx / (2.0 * np.pi))


def _build_rfft_wavevectors_1d(n: int, dx: float) -> np.ndarray:
    return np.fft.rfftfreq(n, d=dx / (2.0 * np.pi))


def apply_fft(x_rank4):
    # log(f"apply_fft: C-contig={x_rank4.flags['C_CONTIGUOUS']}, "
    #     f"F-contig={x_rank4.flags['F_CONTIGUOUS']}, "
    #     f"aligned={x_rank4.flags['ALIGNED']}, "
    #     f"strides={x_rank4.strides}, dtype={x_rank4.dtype}")
    with set_workers(_NWORKERS):
        return rfftn(x_rank4, axes=(0, 1, 2, 3), overwrite_x=True)

def apply_ifft(q_rank4):
    # log(f"apply_ifft: C-contig={q_rank4.flags['C_CONTIGUOUS']}, "
    #     f"F-contig={q_rank4.flags['F_CONTIGUOUS']}, "
    #     f"aligned={q_rank4.flags['ALIGNED']}, "
    #     f"strides={q_rank4.strides}, dtype={q_rank4.dtype}")
    with set_workers(_NWORKERS):
        return irfftn(q_rank4, axes=(0, 1, 2, 3), overwrite_x=True)
    


class FourierSolver(HDF5Saver):
    def __init__(self, n, L,  model_params, dtype="float32", eps=None, verbose=False,
                 mmap=False, mmap_path="/scratch/phi_temp.mmap"):
        self.n = n
        self.L = L
        self.model_params = model_params
        self.dtype = np.dtype(dtype)

        self.dx = 2.0 * L / n
        if eps is None:
            self.eps = (np.pi / (self.n * self.L)) ** 4
        else:
            self.eps = eps

        self.grid_1d = np.linspace(-L, L - self.dx, n)

        self.gcalc = Gcalc(self.grid_1d, **self.model_params, dtype=self.dtype,
                           mmap=mmap, mmap_path=mmap_path)

        self.logger = logging.getLogger("main")

        if verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)

        self.logger.info(f"Required memory: {self.n**4 * self.dtype.itemsize / 1024**3:.2f} GB for N^4 {self.dtype}")

    def precompute(self):
        log("Precomputing Gcalc", notime=True)
        self.gcalc.precompute()
        log("Building wavevectors")
        self.qfull = _build_wavevectors_1d(self.n, self.dx).astype(self.dtype)
        self.qhalf = _build_rfft_wavevectors_1d(self.n, self.dx).astype(self.dtype)

        
        self.qx  = self.qfull[:, None, None, None]
        self.qy  = self.qfull[None, :, None, None]
        self.qxp = self.qfull[None, None, :, None]
        self.qyp = self.qhalf[None, None, None, :]
     
        self.q1_inv_sq = 1.0 / (self.qx ** 2 + self.qy ** 2 + self.eps)
        self.q2_inv_sq = 1.0 / (self.qxp ** 2 + self.qyp ** 2 + self.eps)

        self.component = {0: self.qx, 1: self.qy}
        self.component_prime = {0: self.qxp, 1: self.qyp}

        log("Allocating reusable c buffer")
        self._c_buffer = np.empty((self.n,) * 4, dtype=self.dtype)
        log("Done precomputing")


    def get_C_P_fft(self, a, b):
        c = self.gcalc.calc_C_P_index(a, b)

    def get_C_Q_fft(self, a, b, m, n):
        c = self.gcalc.calc_C_Q_index(a, b, m, n)
        

    def _k2_kernel(self, a, b):
        k1, k2 = self.component[a], self.component_prime[b]
        return k1 * k2 * self.q1_inv_sq * self.q2_inv_sq

    def _k4_kernel(self, a, b, m, n):
        k1, k2 = self.component[a], self.component_prime[b]
        k3, k4 = self.component[m], self.component_prime[n]
        return k1 * k2 * k3 * k4 * self.q1_inv_sq * self.q2_inv_sq

    def solve_P_contribution(self):
        log("Solving P contribution", notime=True)
        s_pressure = None
        for a in range(2):
            for b in range(2):
                log(f"Calculating C_P_index({a}, {b})")
                c = self.gcalc.calc_C_P_index(a, b, out=self._c_buffer)
                log(f"calculating fft (c dtype {c.dtype})")
                s_hat = apply_fft(c)
                del c
                ne.evaluate(
                    "- s_hat * k1 * k2 * q1 * q2",
                    local_dict={
                        "s_hat": s_hat,
                        "k1": self.component[a],
                        "k2": self.component_prime[b],
                        "q1": self.q1_inv_sq,
                        "q2": self.q2_inv_sq,
                    },
                    out=s_hat,
                )
                if s_pressure is None:
                    s_pressure = s_hat
                else:
                    ne.evaluate("s_pressure + s_hat",
                                local_dict={"s_pressure": s_pressure, "s_hat": s_hat},
                                out=s_pressure)
                    del s_hat
                log(f"done (s_pressure dtype {s_pressure.dtype})")
        s_pressure[0, 0, 0, 0] = 0.0

        log("calculating ifft")
        c_pressure = apply_ifft(s_pressure)
        log(f"done (c_pressure dtype {c_pressure.dtype})")

        return c_pressure

    def solve_P_contribution__symmetry_optimized(self):
        """
        Same result as solve_P_contribution but exploits P_ab = P_ba so that
        C_P(a,b) = C_P(b,a). This collapses the 4 (a,b) tuples into 3
        canonical classes, so calc_C_P_index and the FFT each run 3 times
        instead of 4. The (0,1) class folds the kernels for (0,1) and (1,0).
        """
        log("Solving P contribution (symmetry-optimized)", notime=True)

        pair_factor_expr = {
            (0, 0): "(qx * qxp)",
            (0, 1): "(qx * qyp + qy * qxp)",
            (1, 1): "(qy * qyp)",
        }
        canonical_classes = [(0, 0), (0, 1), (1, 1)]

        s_pressure = None
        for (a, b) in canonical_classes:
            log(f"Calculating C_P_index({a},{b}) [class ({a}{b})]")
            c = self.gcalc.calc_C_P_index(a, b, out=self._c_buffer)
            log(f"calculating fft (c dtype {c.dtype})")
            s_hat = apply_fft(c)
            log(f"applying kernel (s_hat dtype {s_hat.dtype} {s_hat.shape})")
            del c

            expr = f"- s_hat * {pair_factor_expr[(a, b)]} * q1 * q2"
            ne.evaluate(
                expr,
                local_dict={
                    "s_hat": s_hat,
                    "qx": self.component[0],
                    "qy": self.component[1],
                    "qxp": self.component_prime[0],
                    "qyp": self.component_prime[1],
                    "q1": self.q1_inv_sq,
                    "q2": self.q2_inv_sq,
                },
                out=s_hat,
            )
            if s_pressure is None:
                s_pressure = s_hat
            else:
                ne.evaluate(
                    "s_pressure + s_hat",
                    local_dict={"s_pressure": s_pressure, "s_hat": s_hat},
                    out=s_pressure,
                )
                del s_hat
            log(f"done (s_pressure dtype {s_pressure.dtype})")

        s_pressure[0, 0, 0, 0] = 0.0

        log("calculating ifft")
        c_pressure = apply_ifft(s_pressure)
        log(f"done (c_pressure dtype {c_pressure.dtype})")
        log("Done solving P contribution (symmetry-optimized)", notime=True)
        return c_pressure

    def solve_Q_contribution(self):
        log("Solving Q contribution", notime=True)
        s_pressure = None
        for a in range(2):
            for b in range(2):
                for m in range(2):
                    for n in range(2):
                        log(f"Calculating C_Q_index({a}, {b}, {m}, {n})")
                        c = self.gcalc.calc_C_Q_index(a, b, m, n, out=self._c_buffer)
                        log(f"calculating fft (c dtype {c.dtype})")
                        s_hat = apply_fft(c)
                        del c
                        ne.evaluate(
                            "s_hat * k1 * k2 * k3 * k4 * q1 * q2",
                            local_dict={
                                "s_hat": s_hat,
                                "k1": self.component[a],
                                "k2": self.component_prime[b],
                                "k3": self.component[m],
                                "k4": self.component_prime[n],
                                "q1": self.q1_inv_sq,
                                "q2": self.q2_inv_sq,
                            },
                            out=s_hat,
                        )
                        if s_pressure is None:
                            s_pressure = s_hat
                        else:
                            ne.evaluate("s_pressure + s_hat",
                                        local_dict={"s_pressure": s_pressure, "s_hat": s_hat},
                                        out=s_pressure)
                            del s_hat
                        log(f"done (s_pressure dtype {s_pressure.dtype})")
        log("calculating ifft")
        c_pressure = apply_ifft(s_pressure)
        log(f"done (c_pressure dtype {c_pressure.dtype})")
        log("Done solving Q contribution", notime=True)
        return c_pressure

        
    def solve_Q_contribution__symmetry_optimized(self):
        """
        Same result as solve_Q_contribution but exploits the three index
        symmetries of C_Q:
            c[a,b,m,n]      = c[b,a,m,n]            (P symmetric)
            c[a,b,m,n]      = c[a,b,n,m]
            c[a,b,m,n][r1,r2] = c[m,n,a,b][r2,r1]
        This collapses the 16 (a,b,m,n) tuples into 6 canonical classes,
        so calc_C_Q_index and the FFT each run 6 times instead of 16.

        The third symmetry is folded in by adding c.transpose(2,3,0,1) to c
        in real space before the FFT (rfft makes the same swap in Fourier
        space awkward because axis 3 is half-sized).

        For each canonical class the per-class kernel factorises as
            ab_part * mn_part * q1_inv_sq * q2_inv_sq
        with
            ab_part = qx*qxp                   if (a,b) = (0,0)
                    = qx*qyp + qy*qxp          if (a,b) = (0,1)
                    = qy*qyp                   if (a,b) = (1,1)
        and similarly for mn_part, by commutativity of the kernel product.
        """
        log("Solving Q contribution (symmetry-optimized)", notime=True)

        pair_factor_expr = {
            (0, 0): "(qx * qxp)",
            (0, 1): "(qx * qyp + qy * qxp)",
            (1, 1): "(qy * qyp)",
        }
        # 6 canonical (ab, mn) classes: 3 diagonal + 3 unordered off-diagonal.
        canonical_classes = [
            ((0, 0), (0, 0)),
            ((0, 1), (0, 1)),
            ((1, 1), (1, 1)),
            ((0, 0), (0, 1)),
            ((0, 0), (1, 1)),
            ((0, 1), (1, 1)),
        ]

        s_pressure = None
        for ab, mn in canonical_classes:
            a, b = ab
            m, n = mn
            log(f"Calculating C_Q_index({a},{b},{m},{n}) "
                f"[class ({a}{b})x({m}{n})]")
            if ab == mn:
                c = self.gcalc.calc_C_Q_index(a, b, m, n, out=self._c_buffer)
            else:
                # Off-diagonal: must symmetrize via c.transpose(2,3,0,1),
                # which aliases c and so forbids writing the calc result
                # directly into the buffer. Allocate raw, symmetrize into
                # the buffer, then drop the raw alloc before the FFT.
                c_raw = self.gcalc.calc_C_Q_index(a, b, m, n)
                log("symmetrizing under r1<->r2")
                ne.evaluate(
                    "c_raw + c_raw_T",
                    local_dict={
                        "c_raw": c_raw,
                        "c_raw_T": c_raw.transpose(2, 3, 0, 1),
                    },
                    out=self._c_buffer,
                )
                del c_raw
                c = self._c_buffer

            log(f"calculating fft (c dtype {c.dtype})")
            s_hat = apply_fft(c)
            del c

            expr = (
                f"s_hat * {pair_factor_expr[ab]} "
                f"* {pair_factor_expr[mn]} * q1 * q2"
            )
            ne.evaluate(
                expr,
                local_dict={
                    "s_hat": s_hat,
                    "qx": self.component[0],
                    "qy": self.component[1],
                    "qxp": self.component_prime[0],
                    "qyp": self.component_prime[1],
                    "q1": self.q1_inv_sq,
                    "q2": self.q2_inv_sq,
                },
                out=s_hat,
            )
            if s_pressure is None:
                s_pressure = s_hat
            else:
                ne.evaluate(
                    "s_pressure + s_hat",
                    local_dict={"s_pressure": s_pressure, "s_hat": s_hat},
                    out=s_pressure,
                )
                del s_hat
            log(f"done (s_pressure dtype {s_pressure.dtype})")

        log("calculating ifft")
        c_pressure = apply_ifft(s_pressure)
        log(f"done (c_pressure dtype {c_pressure.dtype})")
        log("Done solving Q contribution (symmetry-optimized)", notime=True)
        return c_pressure


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

# Keep the model-param defaults next to the Gcalc defaults so the CLI surface
# stays in sync with the Python API.
_MODEL_PARAM_DEFAULTS = {
    "S0":    1.0,
    "l":     1.0,
    "lM":    7.0 / 3,
    "lm":    0.7 / 3,
    "R":     10.0,
    "sigma": 6.0,
}


def _build_arg_parser():
    parser = argparse.ArgumentParser(
        prog="python -m myxo.fourier_solver",
        description=(
            "Compute the requested correlation tensors and pressure-"
            "fluctuation contributions on a 4-D grid via the FFT-based "
            "biharmonic solver, and save them to an HDF5 file (model "
            "parameters are stored as file-level attributes). What gets "
            "computed is determined by which --save_* flags are set."
        ),
    )

    grid = parser.add_argument_group("grid")
    grid.add_argument("--N",       type=int,   required=True,
                      help="Grid size per dimension (N**4 total points).")
    grid.add_argument("--L",       type=float, required=True,
                      help="Half-side length of the periodic box (domain is [-L, L)).")

    model = parser.add_argument_group("model parameters")
    for name, default in _MODEL_PARAM_DEFAULTS.items():
        model.add_argument(f"--{name}", type=float, default=default,
                           help=f"Model parameter {name} (default: {default!r}).")

    solver = parser.add_argument_group("solver options")
    solver.add_argument("--dtype",     type=str,   default="float32",
                        choices=["float32", "float64"],
                        help="NumPy dtype used internally (default: float32).")
    solver.add_argument("--eps",       type=float, default=None,
                        help="Wavevector regulariser; defaults to (pi/(N*L))**4.")
    solver.add_argument("--mmap",      action="store_true",
                        help="Memory-map phi to disk instead of holding it in RAM.")
    solver.add_argument("--mmap-path", type=str,   default="/scratch/phi_temp.mmap",
                        help="Path used when --mmap is set (default: %(default)s).")

    save = parser.add_argument_group(
        "output fields",
        "At least one of these must be set; they decide what gets computed.",
    )
    save.add_argument("--save_C_P",    action="store_true",
                      help="Save all 4 components of C_P as C_P/{ab} datasets "
                           "(uses calc_C_P_index, ignores symmetry).")
    save.add_argument("--save_C_Q",    action="store_true",
                      help="Save all 16 components of C_Q as C_Q/{abmn} "
                           "datasets (uses calc_C_Q_index, ignores symmetry).")
    save.add_argument("--save_P",      action="store_true",
                      help="Save the full 4-D P-pressure contribution as 'P_contrib'.")
    save.add_argument("--save_P_diag", action="store_true",
                      help="Save the diagonal P[i,j,i,j] as 'P_contrib_diag' (shape N, N).")
    save.add_argument("--save_Q",      action="store_true",
                      help="Save the full 4-D Q-pressure contribution as 'Q_contrib'.")
    save.add_argument("--save_Q_diag", action="store_true",
                      help="Save the diagonal Q[i,j,i,j] as 'Q_contrib_diag' (shape N, N).")

    run = parser.add_argument_group("run")
    run.add_argument("--method", choices=["symmetry", "naive"], default="symmetry",
                     help=("Algorithm used for P_contrib / Q_contrib: "
                           "symmetry-optimized (default) or naive index loop."))
    run.add_argument("--output", "-o", required=True,
                     help=("Destination for the HDF5 result. May be either a "
                           "file path or a directory; if a directory (or a "
                           "path ending in '/'), the filename is "
                           "auto-generated from the model parameters."))
    run.add_argument("--overwrite", action="store_true",
                     help="Overwrite --output if it already exists.")
    run.add_argument("--verbose", action="store_true",
                     help="Enable debug-level logging.")
    run.add_argument("--draft", action="store_true",
                     help="Run in draft mode (no precompute, no save).")

    return parser


def main(argv=None):
    args = _build_arg_parser().parse_args(argv)

    logger = logging.getLogger("main")
    logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)

    save_flags = {
        flag: getattr(args, flag) for flag in HDF5Saver._SAVE_FLAG_TAGS
    }
    if not any(save_flags.values()):
        # Fail fast: no save flags means nothing to compute.
        print(
            "Error: at least one --save_* flag must be set "
            f"(choose from: {', '.join('--' + f for f in save_flags)}).",
            file=sys.stderr,
        )
        return 2

    model_params = {name: getattr(args, name) for name in _MODEL_PARAM_DEFAULTS}

    logger.info("Building solver: N=%d, L=%g, dtype=%s, method=%s, save=%s",
                args.N, args.L, args.dtype, args.method,
                ",".join(HDF5Saver._save_tags_from_flags(**save_flags)))
    logger.info("Model params: %s", model_params)

    solver = FourierSolver(
        n=args.N,
        L=args.L,
        model_params=model_params,
        dtype=args.dtype,
        eps=args.eps,
        mmap=args.mmap,
        mmap_path=args.mmap_path,
        verbose=args.verbose,
    )

    output_path = os.path.abspath(
        solver.resolve_output_path(args.output, method=args.method, **save_flags)
    )
    if output_path != os.path.abspath(args.output):
        logger.info("Auto-generated output filename: %s", output_path)

    if os.path.isfile(output_path) and not args.overwrite:
        # Fail fast before doing the expensive precompute / FFTs.
        parser_error = (
            f"Output file {output_path!r} already exists; pass --overwrite "
            f"to replace it."
        )
        print(parser_error, file=sys.stderr)
        return 2

    if args.draft:
        exit(0)

    solver.precompute()

    logger.info("Saving requested fields to %s", output_path)
    solver.save_results_h5(output_path, method=args.method, **save_flags)
    logger.info("Done.")
    # Emit the resolved output path as the sole stdout payload so callers
    # can capture it via $(python -m myxo.fourier_solver ...). All logging
    # is on stderr (see toolbox.setup_loggers), so this stays parseable.
    print(output_path, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
