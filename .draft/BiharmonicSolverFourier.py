"""
Biharmonic Fourier Solver
=========================
Solves:  Lap_r Lap_{r'} f(r,r') = d_{r_i} d_{r'_j} G_{ij}(r,r')

where r, r' in R^2, G_{ij} is a given 2x2 matrix function,
f is regular and vanishes at infinity.

Fourier space solution:
    f_hat(q, q') = -q_i q'_j G_hat_{ij}(q, q') / (|q|^2 |q'|^2 + eps)

Output: HDF5 file readable by Mathematica via Import["file.h5"].

Usage:
    python biharmonic_solver.py --n 32 --L 8.0 --output solution.h5
    python biharmonic_solver.py --n 64 --L 10.0 --eps 1e-8 --workers 4
"""

import os
import numpy as np
from scipy.fft import rfftn, irfftn, set_workers
from concurrent.futures import ThreadPoolExecutor
import h5py
import time
import argparse
from typing import Callable, Optional

from .gdef_periodic import Gcalc  # Example G_{ij} function (not used in main code)

def _build_wavevectors_1d(n: int, dx: float) -> np.ndarray:
    return np.fft.fftfreq(n, d=dx / (2.0 * np.pi))


def _build_rfft_wavevectors_1d(n: int, dx: float) -> np.ndarray:
    return np.fft.rfftfreq(n, d=dx / (2.0 * np.pi))


def solve_biharmonic_fourier(
    n: int,
    L: float,
    output_file: str = "biharmonic_solution.h5",
    eps: Optional[float] = None,
    workers: int = -1,
    compression: int = 4,
    dtype: np.dtype = np.float64,
) -> None:
    """
    Solve the biharmonic equation and write the result to HDF5.

    Parameters
    ----------
    g_func : callable
        Vectorised function with signature
            g_func(x1, y1, x2, y2) -> array of shape (..., 2, 2)
        where x1, y1, x2, y2 are broadcastable arrays.
    n : int
        Number of grid points per dimension.
    L : float
        Half-width of the domain [-L, L] in each coordinate.
    output_file : str
        Path to the output HDF5 file.
    eps : float, optional
        Regularisation parameter. Default: (pi / (n * L))^4.
    workers : int
        Number of FFT threads. -1 = all cores (default).
    dtype : numpy dtype
        Floating-point dtype for grid and computation (e.g. np.float32).
    """
    if workers <= 0:
        workers = os.cpu_count()

    dtype = np.dtype(dtype)
    dx = 2.0 * L / n
    pts = np.linspace(-L, L - dx, n, dtype=dtype)

    if eps is None:
        eps = (np.pi / (n * L)) ** 4

    nrfft = n // 2 + 1

    print("=" * 55)
    print("  Biharmonic Fourier Solver (Python)")
    print("=" * 55)
    print(f"  dtype: {dtype}")
    print(f"  Grid size: {n}^4 = {n**4:,} points")
    print(f"     Memory: {n**4 * dtype.itemsize / 1024**3:.2f} GB")
    print(f"  Domain: [-{L}, {L}]^2 x [-{L}, {L}]^2")
    print(f"  Regularisation: eps = {eps:.3e}")
    print(f"  FFT workers: {workers}")
    print(f"  Output: {output_file}")
    print("-" * 55)

    # ------------------------------------------------------------------
    # Step 1: Evaluate G_{ij} on the 4D grid (fully vectorised)
    # ------------------------------------------------------------------
    print("[Step 1/3] Evaluating G_{ij} on 4D grid ...", flush=True)
    t0 = time.perf_counter()

    G = Gcalc(grid=pts, workers=workers, dtype=dtype)
    G.precompute()

    # x1 = pts[:, None, None, None]
    # y1 = pts[None, :, None, None]
    # x2 = pts[None, None, :, None]
    # y2 = pts[None, None, None, :]

    # g_all = g_func(x1, y1, x2, y2)
    # del x1, y1, x2, y2

    # mem_mb = g_all[..., 0, 0].nbytes / 1024**2
    dt = time.perf_counter() - t0
    print(f"  Done. ({dt:.2f} s)")
    

    # ------------------------------------------------------------------
    # Steps 2+3: Stream rFFT + pointwise multiply, one component at a time.
    #   Peak memory: g_all (4×) + 1 real component + 1 complex fhat ≈ 6×
    #   vs. old approach: g_all (4×) + 4 real + 4 complex ≈ 12×
    # ------------------------------------------------------------------
    print("[Step 2/3] rFFT + Fourier-space accumulation (streaming) ...", flush=True)
    t0 = time.perf_counter()

    fft_axes = (0, 1, 2, 3)

    qfull = _build_wavevectors_1d(n, dx).astype(dtype)
    qhalf = _build_rfft_wavevectors_1d(n, dx).astype(dtype)

    qx  = qfull[:, None, None, None]
    qy  = qfull[None, :, None, None]
    qxp = qfull[None, None, :, None]
    qyp = qhalf[None, None, None, :]

    q1sq = qfull[:, None] ** 2 + qfull[None, :] ** 2
    q2sq = qfull[:, None] ** 2 + qhalf[None, :] ** 2
    # den = q1sq[:, :, None, None] * q2sq[None, None, :, :] + eps
    # del q1sq, q2sq

    _components = [
        (0, 0, qx, qxp),
        (1, 1, qy, qyp),
        (0, 1, qx, qyp),
        (1, 0, qy, qxp),
        
    ]
    # del qx, qy, qxp, qyp


    fhat = None
    g_ij = np.empty((n, n, n, n), dtype=dtype)  
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for i, j, qr, qc in _components:
            print(f"  Computing G_{i}{j} ...")

            if (i, j) != (1, 0):
                G.on_grid(i,j, out=g_ij)
            else:
                # g_10=g_01
                pass

            # g_ij = G.on_grid(i, j)
            print("  Computing rFFT ...")
            with set_workers(workers):
                comp_hat = rfftn(g_ij, axes=fft_axes)
            # del g_ij
        
            def compute_chunk(idx):
                if qr.shape[0] == 1:
                    comp_hat[idx] *= qr
                else:
                    comp_hat[idx] *= qr[idx]

                comp_hat[idx] *= qc
                if fhat is not None:
                    fhat[idx] += comp_hat[idx]

            futures = [pool.submit(compute_chunk, idx) 
                    for idx in np.array_split(np.arange(n), n)]
            for fut in futures:
                fut.result()
            if fhat is None:
                fhat = comp_hat

        del g_ij, comp_hat



        # den = q1sq[:, :, None, None] * q2sq[None, None, :, :] + eps
        # del q1sq, q2sq
        def compute_chunk(idx):
            den_1 = q1sq[idx, :, None, None] * q2sq[None, None, :, :] + eps

            fhat[idx] = -fhat[idx] / den_1

        futures = [pool.submit(compute_chunk, idx) 
                for idx in np.array_split(np.arange(n), n)]
        for fut in futures:
            fut.result()

    # np.negative(fhat, out=fhat)
    # fhat /= den
    # del den

    fhat[0, 0, 0, 0] = 0.0

    dt = time.perf_counter() - t0
    print(f"  Done. ({dt:.2f} s)")

    # ------------------------------------------------------------------
    # Step 4: Inverse 4D rFFT
    # ------------------------------------------------------------------
    print("[Step 3/3] Computing inverse 4D rFFT ...", flush=True)
    t0 = time.perf_counter()

    with set_workers(workers):
        f_grid = irfftn(fhat, s=(n, n, n, n), axes=fft_axes)
    del fhat

    dt = time.perf_counter() - t0
    print(f"  Done. ({dt:.2f} s)")

    print("-" * 55)
    # print(f"  Max |f| = {np.max(np.abs(f_grid)):.6e}")

    # ------------------------------------------------------------------
    # Step 5: Write to HDF5
    # ------------------------------------------------------------------
    print(f"[Export] Writing to {output_file} ...", flush=True)
    t0 = time.perf_counter()

    f_diag = f_grid[np.arange(n)[:, None], np.arange(n)[None, :], 
               np.arange(n)[:, None], np.arange(n)[None, :]]
    with h5py.File(output_file, "w") as hf:
        kw = {"compression": "gzip", "compression_opts": compression} \
             if compression > 0 else {}
        hf.create_dataset("f", data=f_grid, **kw)
        hf.create_dataset("fdiag", data=f_diag, **kw)
        hf.create_dataset("grid", data=pts)
        hf.attrs["n"] = n
        hf.attrs["L"] = L
        hf.attrs["dx"] = dx
        hf.attrs["eps"] = eps

    file_size_mb = os.path.getsize(output_file) / 1024**2
    dt = time.perf_counter() - t0
    print(f"  Done. ({dt:.2f} s, {file_size_mb:.1f} MB on disk)")
    print(f"  Solution complete.")
    print("=" * 55)


# ==================================================================
# TEST EXAMPLE
# ==================================================================



def g_test(x1, y1, x2, y2):
    """Isotropic Gaussian: G_{ij} = delta_{ij} exp(-|r-r'|^2 / 2)."""
    sigma2 = 1.0
    r2 = (x1 - x2) ** 2 + (y1 - y2) ** 2
    scalar = np.exp(-r2 / (2.0 * sigma2))
    out = np.zeros(scalar.shape + (2, 2), dtype=np.float64)
    out[..., 0, 0] = scalar
    out[..., 1, 1] = scalar
    return out


def parse_args():
    parser = argparse.ArgumentParser(
        description="Biharmonic Fourier Solver: "
                    "Lap_r Lap_{r'} f = d_{r_i} d_{r'_j} G_{ij}",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--n", type=int, default=32,
        help="Number of grid points per dimension",
    )
    parser.add_argument(
        "--L", type=float, default=8.0,
        help="Half-width of domain [-L, L]",
    )
    parser.add_argument(
        "--eps", type=float, default=None,
        help="Regularisation parameter (default: (pi/(n*L))^4)",
    )
    parser.add_argument(
        "--workers", type=int, default=-1,
        help="Number of FFT threads (-1 = all cores)",
    )
    parser.add_argument(
        "--output", "-o", type=str, default="biharmonic_solution.h5",
        help="Output HDF5 file path",
    )
    parser.add_argument(
        "--compress", type=int, default=0, choices=range(0, 10),
        help="Gzip compression level (0 = none, 9 = max)",
    )
    parser.add_argument(
        "--dtype", type=str, default="float64",
        choices=["float16", "float32", "float64"],
        help="Floating-point dtype for computation",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    output_folder = os.path.dirname(args.output)
    os.makedirs(output_folder, exist_ok=True)

    solve_biharmonic_fourier(
        n=args.n,
        L=args.L,
        output_file=args.output,
        eps=args.eps,
        workers=args.workers,
        compression=args.compress,
        dtype=np.dtype(args.dtype),
    )

    # ---- Quick read-back check ----
    with h5py.File(args.output, "r") as hf:
        f = hf["f"][:]
        grid = hf["grid"][:]
        print(f"\nRead-back check:")
        print(f"  f shape: {f.shape}")
        print(f"  grid: [{grid[0]:.2f}, ..., {grid[-1]:.2f}], {len(grid)} pts")
        print(f"  f range: [{f.min():.4e}, {f.max():.4e}]")

    print(f"""
To load in Mathematica:
    f = Import["{args.output}", {{"Datasets", "/f"}}];
    grid = Import["{args.output}", {{"Datasets", "/grid"}}];
""")