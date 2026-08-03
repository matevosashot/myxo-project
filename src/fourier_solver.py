import argparse
import gc
import logging
import os
import sys
import time
import json 

import numexpr as ne
import numpy as np
from scipy.fft import irfftn, rfftn, set_workers, fftn, ifftn

from . import _NWORKERS
from .correlation_def import Gcalc
from .correlation_def_v2 import GcalcForQ
from .saving_utils import HDF5Saver
from .utils import log, diag_large_ndarrays, clear_numexpr_last_cache



def _build_wavevectors_1d(n: int, dx: float) -> np.ndarray:
    return np.fft.fftfreq(n, d=dx / (2.0 * np.pi))


def _build_rfft_wavevectors_1d(n: int, dx: float) -> np.ndarray:
    return np.fft.rfftfreq(n, d=dx / (2.0 * np.pi))


def apply_fft(x_rank4, fft_type="rfft", overwrite_x=True):
    # log(f"apply_fft: C-contig={x_rank4.flags['C_CONTIGUOUS']}, "
    #     f"F-contig={x_rank4.flags['F_CONTIGUOUS']}, "
    #     f"aligned={x_rank4.flags['ALIGNED']}, "
    #     f"strides={x_rank4.strides}, dtype={x_rank4.dtype}")
    fn = rfftn if fft_type == "rfft" else fftn
    with set_workers(_NWORKERS):
        return fn(x_rank4, axes=(0, 1, 2, 3), overwrite_x=overwrite_x)

def apply_ifft(q_rank4, fft_type="rfft", overwrite_x=True, s=None):
    # log(f"apply_ifft: C-contig={q_rank4.flags['C_CONTIGUOUS']}, "
    #     f"F-contig={q_rank4.flags['F_CONTIGUOUS']}, "
    #     f"aligned={q_rank4.flags['ALIGNED']}, "
    #     f"strides={q_rank4.strides}, dtype={q_rank4.dtype}")
    # `s` is the (full) output shape on the transformed axes. Required for
    # rfft mode when the original last-axis length is odd, because
    # irfftn defaults to 2*(M-1) which is off-by-one for odd n. fft mode
    # ignores `s`: ifftn never has a shape ambiguity since the input
    # already carries the full spectrum.
    with set_workers(_NWORKERS):
        if fft_type == "rfft":
            return irfftn(q_rank4, s=s, axes=(0, 1, 2, 3), overwrite_x=overwrite_x)
        # ifftn returns a complex array; copying .real gives a fresh
        # contiguous real buffer to match irfftn's dtype/shape semantics.
        out = ifftn(q_rank4, axes=(0, 1, 2, 3), overwrite_x=overwrite_x)
        return out.real.copy()


_SPATIAL_AXES_2D = (0, 1)


def apply_fft_spatial(x, fft_type="rfft", overwrite_x=True):
    fn = rfftn if fft_type == "rfft" else fftn
    with set_workers(_NWORKERS):
        return fn(x, axes=_SPATIAL_AXES_2D, overwrite_x=overwrite_x)


def apply_ifft_spatial(x, fft_type="rfft", overwrite_x=True, s=None):
    with set_workers(_NWORKERS):
        if fft_type == "rfft":
            return irfftn(x, s=s, axes=_SPATIAL_AXES_2D, overwrite_x=overwrite_x)
        out = ifftn(x, axes=_SPATIAL_AXES_2D, overwrite_x=overwrite_x)
        return out.real.copy()


class FourierSolver(HDF5Saver):
    def __init__(self, n, L,  model_params, dtype="float32", eps=None, verbose=False,
                 mmap=False, mmap_path="/scratch/phi_temp.mmap", drop_phi=False,
                 gcalc_Q_cls=None, fft_type="rfft", ):

        if n % 2 == 0 and not os.environ.get("MYXO_TESTING", "0") == "1":
            raise ValueError(
                f"n must be odd, got {n!r}"
            )

        if fft_type not in ("rfft", "fft"):
            raise ValueError(
                f"fft_type must be 'rfft' or 'fft', got {fft_type!r}"
            )
        self.fft_type = fft_type
        self.n = n
        self.L = L
        self.model_params = model_params
        self.dtype = np.dtype(dtype)
        self.drop_phi = drop_phi

        self.dx = 2.0 * L / n
        if eps is None:
            self.eps = (np.pi / (self.n * self.L)) ** 4
        else:
            self.eps = eps

        self.grid_1d = np.linspace(-L, L - self.dx, n)

        # P side: always the default Gcalc.
        self.gcalc = Gcalc(self.grid_1d, **self.model_params, dtype=self.dtype,
                           mmap=mmap, mmap_path=mmap_path)
        # Q side: optionally a different Gcalc subclass (e.g. correlation_def_v2)
        # so the C_Q computation can use a different phi formula. When not
        # provided, it aliases to self.gcalc so memory and behaviour are
        # unchanged. When provided, the second instance gets its own mmap
        # path so the two phi buffers don't collide.
        if gcalc_Q_cls is None or gcalc_Q_cls is Gcalc:
            self.gcalc_Q = self.gcalc
        else:
            self.gcalc_Q = gcalc_Q_cls(self.grid_1d, **self.model_params, dtype=self.dtype,
                                       mmap=mmap, mmap_path=mmap_path + ".Q")

        self.logger = logging.getLogger("main")

        if verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)

        self.logger.info(f"Required memory: {self.n**4 * self.dtype.itemsize / 1024**3:.2f} GB for N^4 {self.dtype}")

    def _ensure_phi(self, gcalc):
        if self.drop_phi and gcalc.phi is None:
            log("Re-materializing phi")
            gcalc.calc_phi()

    def _maybe_drop_phi(self, gcalc):
        if self.drop_phi:
            gcalc.drop_phi()

    def _drop_inactive_phi(self, active):
        # Each solve_* uses exactly one of (gcalc, gcalc_Q). If they're
        # distinct instances and the inactive one is holding an N^4 phi
        # (any Gcalc subclass that inherits v1's calc_phi will), free it
        # for the duration of this solve so we don't pay 4*N^4 peak.
        # No-op when the two are aliased, or when drop_phi=False.
        if not self.drop_phi:
            return
        other = self.gcalc_Q if active is self.gcalc else self.gcalc
        if other is not active:
            log(f"Dropping inactive {type(other).__name__}.phi")
            other.drop_phi()

    def precompute(self):
        log("Precomputing Gcalc", notime=True)
        self.gcalc.precompute()
        if not self.drop_phi:
            log("Calculating phi")
            self.gcalc.calc_phi()
        if self.gcalc_Q is not self.gcalc:
            log("Precomputing Gcalc (Q-side)", notime=True)
            self.gcalc_Q.precompute()
            if not self.drop_phi:
                log("Calculating phi (Q-side)")
                self.gcalc_Q.calc_phi()
                
        log("Building wavevectors")
        self.qfull = _build_wavevectors_1d(self.n, self.dx).astype(self.dtype)
        if self.fft_type == "rfft":
            self.qhalf = _build_rfft_wavevectors_1d(self.n, self.dx).astype(self.dtype)
            qyp_axis = self.qhalf
        else:
            qyp_axis = self.qfull

        self.qx  = self.qfull[:, None, None, None]
        self.qy  = self.qfull[None, :, None, None]
        self.qxp = self.qfull[None, None, :, None]
        self.qyp = qyp_axis[None, None, None, :]
     
        self.q1_inv_sq = 1.0 / (self.qx ** 2 + self.qy ** 2 + self.eps)
        self.q2_inv_sq = 1.0 / (self.qxp ** 2 + self.qyp ** 2 + self.eps)

        self.component = {0: self.qx, 1: self.qy}
        self.component_prime = {0: self.qxp, 1: self.qyp}
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
        self._drop_inactive_phi(self.gcalc)
        
        c_buffer = np.empty((self.n,) * 4, dtype=self.dtype)
        s_pressure = None
        for a in range(2):
            for b in range(2):
                log(f"Calculating C_P_index({a}, {b})")
                self._ensure_phi(self.gcalc)
                c = self.gcalc.calc_C_P_index(a, b, out=c_buffer)
                self._maybe_drop_phi(self.gcalc)

                
                log(f"calculating fft (c dtype {c.dtype})")
                s_hat = apply_fft(c, fft_type=self.fft_type)
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
                    out=s_hat                )
                if s_pressure is None:
                    s_pressure = s_hat
                else:
                    ne.evaluate("s_pressure + s_hat",
                                local_dict={"s_pressure": s_pressure, "s_hat": s_hat},
                                out=s_pressure)
                    del s_hat
                log(f"done (s_pressure dtype {s_pressure.dtype})")
        s_pressure[0, 0, 0, 0] = 0.0

        del c_buffer
        log("calculating ifft")
        c_pressure = apply_ifft(s_pressure, fft_type=self.fft_type, s=(self.n,) * 4)
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
        self._drop_inactive_phi(self.gcalc)

        pair_factor_expr = {
            (0, 0): "(qx * qxp)",
            (0, 1): "(qx * qyp + qy * qxp)",
            (1, 1): "(qy * qyp)",
        }
        canonical_classes = [(0, 0), (0, 1), (1, 1)]

        c_buffer = np.empty((self.n,) * 4, dtype=self.dtype)
        s_pressure = None
        for (a, b) in canonical_classes:
            log(f"Calculating C_P_index({a},{b}) [class ({a}{b})]")
            self._ensure_phi(self.gcalc)
            c = self.gcalc.calc_C_P_index(a, b, out=c_buffer)
            self._maybe_drop_phi(self.gcalc)


            log(f"calculating fft (c dtype {c.dtype})")
            s_hat = apply_fft(c, fft_type=self.fft_type)
            del c
            log(f"applying kernel (s_hat dtype {s_hat.dtype} {s_hat.shape})")

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
                out=s_hat)
            if s_pressure is None:
                s_pressure = s_hat
            else:
                
                ne.evaluate(
                    "s_pressure + s_hat",
                    local_dict={"s_pressure": s_pressure, "s_hat": s_hat},
                    out=s_pressure
                )
                
                del s_hat
            log(f"done (s_pressure dtype {s_pressure.dtype})")
        del c_buffer
        clear_numexpr_last_cache()
        s_pressure[0, 0, 0, 0] = 0.0

        log("calculating ifft")
        c_pressure = apply_ifft(s_pressure, fft_type=self.fft_type, s=(self.n,) * 4)
        log(f"done (c_pressure dtype {c_pressure.dtype})")
        log("Done solving P contribution (symmetry-optimized)", notime=True)
        return c_pressure

    def solve_Q_contribution(self):
        log("Solving Q contribution", notime=True)
        self._drop_inactive_phi(self.gcalc_Q)
        c_buffer = np.empty((self.n,) * 4, dtype=self.dtype)
        s_pressure = None
        for a in range(2):
            for b in range(2):
                for m in range(2):
                    for n in range(2):
                        log(f"Calculating C_Q_index({a}, {b}, {m}, {n})")
                        self._ensure_phi(self.gcalc_Q)
                        c = self.gcalc_Q.calc_C_Q_index(a, b, m, n, out=c_buffer)
                        self._maybe_drop_phi(self.gcalc_Q)
                        log(f"calculating fft (c dtype {c.dtype})")
                        s_hat = apply_fft(c, fft_type=self.fft_type)
                        del c
                        ne.evaluate(
                            "s_hat * k1 * k2 * k3 * k4 * q1 * q2",
                            local_dict={
                                "s_hat": s_hat,
                                "k1": self.component[a],
                                "k2": self.component[b],
                                "k3": self.component_prime[m],
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
        s_pressure[0, 0, 0, 0] = 0.0
        del c_buffer
        clear_numexpr_last_cache()
        log("calculating ifft")
        c_pressure = apply_ifft(s_pressure, fft_type=self.fft_type, s=(self.n,) * 4)
        log(f"done (c_pressure dtype {c_pressure.dtype})")
        log("Done solving Q contribution", notime=True)
        return c_pressure


    def solve_Q_contribution__symmetry_optimized(self):
        """
        Same result as solve_Q_contribution but exploits the three C_Q
        symmetries:
            c[a,b,m,n]        = c[b,a,m,n]               (P symmetric)
            c[a,b,m,n]        = c[a,b,n,m]
            c[a,b,m,n](r1,r2) = c[m,n,a,b](r2,r1)
        With kernel K(a,b,m,n)(q1,q2) = q1_a q1_b q2_m q2_n / (|q1|^2|q2|^2),
        the kernel respects the same symmetries in Fourier space (the last
        one with a q1<->q2 swap), so the 16 tuples collapse into 6 canonical
        classes -> calc_C_Q_index runs 6 times instead of 16.

        Each class's kernel-sum factorises as
            f_q1[ab](q1) * f_q2[mn](q2) / (|q1|^2|q2|^2)
        with
            f_q1[ab] in {qx^2,   2*qx*qy,    qy^2}   (q1 uses full grid)
            f_q2[mn] in {qxp^2,  2*qxp*qyp,  qyp^2}  (q2 uses rfft grid)

        Off-diagonal classes (ab != mn) come in (A, B=A-mirror) pairs related
        by the third symmetry. ĉ_B(q1,q2) = ĉ_A(q2,q1) but K_A != K_B (only
        K_A(q1,q2) = K_B(q2,q1)), so the symmetrised-c trick that worked
        with the old k1*k2*k3*k4 = q1_a q2_b q1_m q2_n kernel is no longer
        valid. Instead, one calc_C_Q_index call feeds two FFTs: the second
        FFT runs on c.transpose(2,3,0,1), which by the third symmetry is the
        mirror-class c in real space, and yields ĉ_B(q1,q2) in Fourier space.
        """
        log("Solving Q contribution (symmetry-optimized)", notime=True)
        self._drop_inactive_phi(self.gcalc_Q)

        f_q1 = {
            (0, 0): "(qx * qx)",
            (0, 1): "(2 * qx * qy)",
            (1, 1): "(qy * qy)",
        }
        f_q2 = {
            (0, 0): "(qxp * qxp)",
            (0, 1): "(2 * qxp * qyp)",
            (1, 1): "(qyp * qyp)",
        }
        diagonal_classes = [
            ((0, 0), (0, 0)),
            ((0, 1), (0, 1)),
            ((1, 1), (1, 1)),
        ]
        off_diagonal_pairs = [
            ((0, 0), (0, 1)),
            ((0, 0), (1, 1)),
            ((0, 1), (1, 1)),
        ]

        c_buffer = np.empty((self.n,) * 4, dtype=self.dtype)
        s_pressure = None

        kernel_locals = {
            "qx": self.component[0],
            "qy": self.component[1],
            "qxp": self.component_prime[0],
            "qyp": self.component_prime[1],
            "q1": self.q1_inv_sq,
            "q2": self.q2_inv_sq,
        }

        def _apply_kernel_and_accumulate(s_hat, ab_class, mn_class):
            nonlocal s_pressure
            ne.evaluate(
                f"s_hat * {f_q1[ab_class]} * {f_q2[mn_class]} * q1 * q2",
                local_dict={"s_hat": s_hat, **kernel_locals},
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

        for ab, mn in diagonal_classes:
            a, b = ab
            m, n = mn
            log(f"Calculating C_Q_index({a},{b},{m},{n}) "
                f"[diag class ({a}{b})x({m}{n})]")
            self._ensure_phi(self.gcalc_Q)
            c = self.gcalc_Q.calc_C_Q_index(a, b, m, n, out=c_buffer)
            self._maybe_drop_phi(self.gcalc_Q)
            log(f"calculating fft (c dtype {c.dtype})")
            s_hat = apply_fft(c, fft_type=self.fft_type)

            _apply_kernel_and_accumulate(s_hat, ab, mn)
            if s_pressure is not s_hat:
                del s_hat
            del c
            log(f"done (s_pressure dtype {s_pressure.dtype})")

        for ab, mn in off_diagonal_pairs:
            a, b = ab
            m, n = mn
            log(f"Calculating C_Q_index({a},{b},{m},{n}) "
                f"[off-diag pair ({a}{b}),({m}{n})]")
            self._ensure_phi(self.gcalc_Q)
            c = self.gcalc_Q.calc_C_Q_index(a, b, m, n, out=c_buffer)
            c_T = c.transpose(2, 3, 0, 1)

            self._maybe_drop_phi(self.gcalc_Q)

            # FFT of c -> ĉ_A. overwrite_x=False so c stays live for the
            # second FFT below.
            log(f"calculating fft for ({a}{b})x({m}{n})")
            s_hat = apply_fft(c, fft_type=self.fft_type, overwrite_x=False)
            del c
            _apply_kernel_and_accumulate(s_hat, ab, mn)
            if s_pressure is not s_hat:
                del s_hat

            # c.transpose(2,3,0,1) is, by the third C_Q symmetry, the
            # mirror class's c in real space; its FFT gives ĉ_A(q2,q1)
            # = ĉ_B(q1,q2). rfftn/fftn internally materialise a
            # contiguous copy of the non-contiguous view, so we leave
            # the underlying c_buffer untouched and reuse it on the
            # next iteration.
            log(f"calculating fft for mirror ({m}{n})x({a}{b})")
            s_hat = apply_fft(c_T, fft_type=self.fft_type, overwrite_x=False)
            del c_T
            _apply_kernel_and_accumulate(s_hat, mn, ab)
            if s_pressure is not s_hat:
                del s_hat
            log(f"done (s_pressure dtype {s_pressure.dtype})")

        s_pressure[0, 0, 0, 0] = 0.0
        del c_buffer
        clear_numexpr_last_cache()
        
        log("calculating ifft")
        c_pressure = apply_ifft(s_pressure, fft_type=self.fft_type, s=(self.n,) * 4)

        log(f"done (c_pressure dtype {c_pressure.dtype})")
        log("Done solving Q contribution (symmetry-optimized)", notime=True)
        return c_pressure



    def solve_pressure_field(self):
        r"""
        0 = grad p + \xi v + \zeta_c div Q

        p = - k_a k_b / |k|^2 * Q_ab  (\zeta_c = 1)  in fourier space
        """
        Q = self.gcalc_Q.Q  # shape (N, N, 2, 2)
        Q_fourier = apply_fft_spatial(Q, fft_type=self.fft_type)

        if self.fft_type == "rfft":
            k_comp = (self.qfull[:, None], self.qhalf[None, :])
        else:
            k_comp = (self.qfull[:, None], self.qfull[None, :])

        k_inv_squared = 1.0 / (k_comp[0] ** 2 + k_comp[1] ** 2 + self.eps)

        k_outer = np.empty(k_inv_squared.shape + (2, 2), dtype=self.dtype)
        for a in range(2):
            for b in range(2):
                k_outer[..., a, b] = k_comp[a] * k_comp[b]

        kernel = k_outer * k_inv_squared[..., None, None]
        p_fourier = -np.sum(kernel * Q_fourier, axis=(2, 3))
        p_fourier[0, 0] = 0.0

        return apply_ifft_spatial(
            p_fourier, fft_type=self.fft_type, s=(self.n, self.n),
        )








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
    solver.add_argument("--drop_phi",  action="store_true",
                        help="Recompute phi each iteration instead of caching it. "
                             "Cuts inner-loop peak from 4*N^4 to 3*N^4 at the cost "
                             "of recomputing phi N_iter times per solve.")

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

    import toolbox
    toolbox.setup_loggers(
        base_path="./myxo.log",
        debug=False,
        stdout=True,
        train_logger=False,
    )

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
        drop_phi=args.drop_phi,
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
