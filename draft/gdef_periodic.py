"""
Vectorised G_{ij}(r, r') for a half-charge nematic defect.

G_{ij}(x1,y1,x2,y2) = (1/2)[ A(r) * N(r'; r, Sigma(r))
                              + A(r') * N(r; r', Sigma(r')) ]_{ij}

where:
  S(r)     = S0 * f(r/l),  f(x) = x * sqrt((0.34 + 0.07 x^2) / (1 + 0.41 x^2 + 0.07 x^4))
  theta    = phi / 2       (half-vortex)
  Q        = S * [[cos 2theta, sin 2theta], [sin 2theta, -cos 2theta]]
  A        = (I + Q) / 2
  Sigma    = (lM^2 + lm^2)/2 * I  +  (lM^2 - lm^2)/2 * Q
  N(x; mu, Sigma) = multinormal PDF

All functions are fully vectorised over arbitrary broadcast shapes.
"""

import numpy as np
import argparse
import os
import time
from concurrent.futures import ThreadPoolExecutor
from functools import partial

def f(x):
    """f[x_] := x Sqrt[(0.34 + 0.07 x^2)/(1 + 0.41 x^2 + 0.07 x^4)]"""
    x2 = x**2
    return x * np.sqrt((0.34 + 0.07 * x2) / (1.0 + 0.41 * x2 + 0.07 * x2**2))


# def S(r, S0, l):
#     """S[r_] := S0 f[r/l]"""
#     return S0 * f(r / l)

def S(r, S0, l, R):
    """S[r_] := S0 f[r/l]"""
    return S0 * f(r / l) * np.exp(-(r / R))  # add Gaussian cutoff to ensure finite integrals


# def Q(x, y, S0, l):
#     """
#     Q tensor (2x2) at positions (x, y).
#     x, y can be arrays of shape (...).
#     Returns shape (..., 2, 2).
#     """
#     x = np.asarray(x)
#     y = np.asarray(y)
#     dtype = np.result_type(x, y)
#     r = np.sqrt(x**2 + y**2)
#     phi = np.arctan2(y, x)          # Mathematica ArcTan[x, y] = atan2(y, x)
#     theta = phi / dtype.type(2)     # Theta[phi] := phi/2
#     s = S(r, S0, l)                 # scalar field, shape (...)
#     c2 = np.cos(dtype.type(2) * theta)
#     s2 = np.sin(dtype.type(2) * theta)

#     # Build (..., 2, 2) tensor:  [[S cos2θ, S sin2θ], [S sin2θ, -S cos2θ]]
#     out = np.empty(np.broadcast_shapes(np.shape(x), np.shape(y)) + (2, 2), dtype=dtype)
#     out[..., 0, 0] = s * c2
#     out[..., 0, 1] = s * s2
#     out[..., 1, 0] = s * s2
#     out[..., 1, 1] = -s * c2
#     return out


# def A(x, y, S0, l):
#     """A = 1/2 (I + Q)  — shape (..., 2, 2)"""
#     x = np.asarray(x)
#     y = np.asarray(y)
#     dtype = np.result_type(x, y)
#     eye = np.eye(2, dtype=dtype)
#     return dtype.type(0.5) * (eye + Q(x, y, S0, l))


# def Sigma(x, y, S0, l, lM, lm):
#     """
#     Σ = 1/2 I (lM² + lm²) + 1/2 Q (lM² - lm²)
#     Returns shape (..., 2, 2) covariance matrix.
#     """
#     x = np.asarray(x)
#     y = np.asarray(y)
#     dtype = np.result_type(x, y)
#     eye = np.eye(2, dtype=dtype)
#     q = Q(x, y, S0, l)
#     return dtype.type(0.5) * eye * dtype.type(lM**2 + lm**2) + dtype.type(0.5) * q * dtype.type(lM**2 - lm**2)


def _mvn_pdf_2d(dx, dy, cov):
    """
    Vectorised 2-D multivariate normal PDF from displacement components.
    dx, dy : broadcastable scalar arrays  (displacement r2 - r1)
    cov    : (..., 2, 2)  — may have fewer broadcast dims than dx/dy
    returns: broadcast shape of (dx, dy, cov[..., 0, 0])
    """
    dtype = cov.dtype
    det  = cov[..., 0, 0] * cov[..., 1, 1] - cov[..., 0, 1] * cov[..., 1, 0]
    # Mahalanobis: uses only 3 unique inverse entries (cov is symmetric)
    maha2 = (  cov[..., 1, 1] / det * dx**2
             - dtype.type(2) * cov[..., 0, 1] / det * dx * dy
             +  cov[..., 0, 0] / det * dy**2 )
    return np.exp(dtype.type(-0.5) * maha2) / (dtype.type(2 * np.pi) * np.sqrt(det))


# def Gcomponents(x1, y1, x2, y2, S0, l, lM, lm):
#     """
#     G = 1/2 [ A(r1) · N(r2 | r1, Σ(r1))  +  A(r2) · N(r1 | r2, Σ(r2)) ]

#     All inputs can be broadcastable arrays.
#     Returns shape (..., 2, 2).
#     """
#     # Anisotropy tensors — sparse shapes, e.g. (N,N,1,1,2,2) and (1,1,N,N,2,2)
#     A1 = A(x1, y1, S0, l)
#     A2 = A(x2, y2, S0, l)
#     print("A1", A1.shape)
#     print("A2",A2.shape)


#     # Covariance matrices — same sparse shapes
#     S1 = Sigma(x1, y1, S0, l, lM, lm)
#     S2 = Sigma(x2, y2, S0, l, lM, lm)

#     # Displacement components — sparse, e.g. (N,1,N,1) and (1,N,1,N)
#     dx = x2 - x1
#     dy = y2 - y1

#     # PDF scalars expand to full broadcast shape only here
#     pdf1 = _mvn_pdf_2d( dx,  dy, S1)   # N(r2 | r1, Σ1)
#     pdf2 = _mvn_pdf_2d(-dx, -dy, S2)   # N(r1 | r2, Σ2)

#     return A1,pdf1, A2, pdf2

# def G1(x1, y1, x2, y2, S0, l, lM, lm):
#     """
#     G = 1/2 [ A(r1) · N(r2 | r1, Σ(r1))  +  A(r2) · N(r1 | r2, Σ(r2)) ]

#     All inputs can be broadcastable arrays.
#     Returns shape (..., 2, 2).
#     """
#     A1, pdf1, A2, pdf2 = Gcomponents(x1, y1, x2, y2, S0, l, lM, lm)
#     return 0.5 * (A1 * pdf1[..., None, None] + A2 * pdf2[..., None, None])

# ==================================================================
# Usage example
# ==================================================================


class Gcalc:
    def __init__(self, grid, S0=1.0, l=1.0, lM=7.0/3, lm=0.7/3, R=10, workers=0, 
                 dtype=np.float64):   
        self.S0 = S0
        self.l = l
        self.lM = lM
        self.lm = lm
        self.R = R
        self.workers = os.cpu_count() if workers <= 0 else workers

        self.grid = grid.astype(dtype)

        dx = grid[1] - grid[0]
        self.L_grid = grid[-1] - grid[0] + dx  # total length of the periodic box
    
        self.dtype = dtype

        self.pool = ThreadPoolExecutor(max_workers=self.workers) 


    def calcQ(self, x, y):
        """
        Q tensor (2x2) at positions (x, y).
        x, y can be arrays of shape (...).
        Returns shape (..., 2, 2).
        """
        x = np.asarray(x)
        y = np.asarray(y)
        dtype = np.result_type(x, y)
        r = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)          # Mathematica ArcTan[x, y] = atan2(y, x)
        theta = phi / dtype.type(2)     # Theta[phi] := phi/2
        s = S(r, self.S0, self.l, self.R)                 # scalar field, shape (...)
        c2 = np.cos(dtype.type(2) * theta)
        s2 = np.sin(dtype.type(2) * theta)

        # Build (..., 2, 2) tensor:  [[S cos2θ, S sin2θ], [S sin2θ, -S cos2θ]]
        out = np.empty(np.broadcast_shapes(np.shape(x), np.shape(y)) + (2, 2), dtype=dtype)
        out[..., 0, 0] = s * c2
        out[..., 0, 1] = s * s2
        out[..., 1, 0] = s * s2
        out[..., 1, 1] = -s * c2
        return out

    def calcA(self, x, y, Q):
        """A = 1/2 (I + Q)  — shape (..., 2, 2)"""
        x = np.asarray(x)
        y = np.asarray(y)
        dtype = np.result_type(x, y)
        eye = np.eye(2, dtype=dtype)
        return dtype.type(0.5) * (eye + Q)
    
    def calcSigma(self, x, y, Q):
        """
        Σ = 1/2 I (lM² + lm²) + 1/2 Q (lM² - lm²)
        Returns shape (..., 2, 2) covariance matrix.
        """
        x = np.asarray(x)
        y = np.asarray(y)
        dtype = np.result_type(x, y)
        eye = np.eye(2, dtype=dtype)
        q = Q
        return dtype.type(0.5) * eye * dtype.type(self.lM**2 + self.lm**2) + dtype.type(0.5) * q * dtype.type(self.lM**2 - self.lm**2)



    def precompute(self):
        self.x1 = self.grid[:, None, None, None]
        self.y1 = self.grid[None, :, None, None]
        self.x2 = self.grid[None, None, :, None]
        self.y2 = self.grid[None, None, None, :]

        Q = self.calcQ(self.x1[..., 0, 0], self.y1[..., 0, 0])
        self.A = self.calcA(self.x1[..., 0, 0], self.y1[..., 0, 0], Q)
        self.Sigma = self.calcSigma(self.x1[..., 0, 0], self.y1[..., 0, 0], Q)

        
        n_workers = os.cpu_count() if self.workers <= 0 else self.workers
        print(f"Precomputing PDFs in parallel with {n_workers} threads")
        # Split along x1 axis; numpy releases the GIL so threads run concurrently
        N = len(self.grid)
        chunk_indices = np.array_split(np.arange(N), N)
        pdf1 = np.empty((N, N, N, N), dtype=self.dtype)
        # pdf2 = np.empty((N, N, N, N), dtype=self.dtype)
        
        S1 = self.Sigma[:,:,None,None]
        # S2 = self.Sigma[None,None,:,:]

        dy = self.y2 - self.y1
        
        dy = dy - self.L_grid * np.round(dy / self.L_grid)  # periodic boundary conditions

        def compute_chunk(idx):
            x1c = self.x1[idx]
            dxc = self.x2 - x1c

            dxc = dxc - self.L_grid * np.round(dxc / self.L_grid)  # periodic boundary conditions

            S1c = self.Sigma[idx,:, None, None]

            pdf1[idx] = _mvn_pdf_2d(dxc, dy, S1c)
            # pdf2[idx] = _mvn_pdf_2d(-dxc, -dy, S2)

        # with ThreadPoolExecutor(max_workers=n_workers) as pool:
        futures = [self.pool.submit(compute_chunk, idx) for idx in chunk_indices]
        for fut in futures:
            fut.result()
            
        self.pdf1 = pdf1
        self.pdf2 = pdf1.transpose(2, 3, 0, 1)  # N(r1 | r2, Σ2) = N(r2 | r1, Σ1) with x1,x2 swapped

    def on_grid(self, i=None, j=None, out=None):
        N = len(self.grid)
        chunk_indices = np.array_split(np.arange(N), N)

        if out is not None:
            values = out
        elif i is None or j is None:
            values = np.empty((N, N, N, N, 2, 2), dtype=self.dtype)
        else:
            values = np.empty((N, N, N, N), dtype=self.dtype)

        def compute_chunk(idx):
            # if len(idx) == 0:
            #     return
            if i is None or j is None:
                # return N, N, N, N, 2, 2 matrix
                # tmp = 0.5 * (self.A[idx, :, None, None] * self.pdf1[idx,..., None, None])  # shape (K, N, N, N, 2, 2)
                # # print("tmp", tmp.shape)
                # values[idx, ...] += tmp 
                # values[:, :, idx, :] += tmp.transpose(2,3,0,1,4,5)  # shape (N, N, K, N, 2, 2)

                values[idx, ...] = 0.5 * (self.A[idx, :, None, None] * self.pdf1[idx,..., None, None] + 
                    self.A[None, None, :, :] * self.pdf2[idx,..., None, None])
            else:
                values[idx, ...] = 0.5 * (self.A[idx,:, None, None, i, j] * self.pdf1[idx] + 
                    self.A[None, None, :, :, i, j] * self.pdf2[idx])

        # with ThreadPoolExecutor(max_workers=self.workers) as pool:
        futures = [self.pool.submit(compute_chunk, idx) for idx in chunk_indices]
        for fut in futures:
            fut.result()
        
        return values
    

    def __call__(self, x1, y1, x2, y2):
        
        Q1 = self.calcQ(x1, y1)
        Q2 = self.calcQ(x2, y2)
        A1 = self.calcA(x1, y1, Q1)
        A2 = self.calcA(x2, y2, Q2)
        S1 = self.calcSigma(x1, y1, Q1)
        S2 = self.calcSigma(x2, y2, Q2)

        dx = x2 - x1
        dy = y2 - y1

        dx = dx - self.L_grid * np.round(dx / self.L_grid)  # periodic boundary conditions
        dy = dy - self.L_grid * np.round(dy / self.L_grid)

        pdf1 = _mvn_pdf_2d( dx,  dy, S1)   # N(r2 | r1, Σ1)
        pdf2 = _mvn_pdf_2d( dx,  dy, S2)   # N(r1 | r2, Σ2)

        return 0.5 * (A1 * pdf1[..., None, None] + A2 * pdf2[..., None, None])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute G_{ij} on a 4D grid")
    parser.add_argument("--S0",       type=float, default=1.0)
    parser.add_argument("--l",        type=float, default=1.0)
    parser.add_argument("--lM",       type=float, default=7.0/3)
    parser.add_argument("--lm",       type=float, default=0.7/3)
    parser.add_argument("--R",        type=float, default=10, help="Cutoff radius")
    parser.add_argument("--N",        type=int,   default=50, help="Grid size per dimension")
    parser.add_argument("--parallel", action="store_true",     help="Compute G in parallel using threads")
    parser.add_argument("--workers",  type=int,   default=-1, help="Number of worker threads (default: CPU count)")
    parser.add_argument("--dtype",    type=str,   default="float64", help="NumPy dtype for computation (e.g. float16, float32, float64)")
    args = parser.parse_args()

    dtype = np.dtype(args.dtype)

    print("O(N^4) memory is {:0.1f}GB".format(args.N**4 / 1024**3 * dtype.itemsize))

    G = Gcalc(grid=np.linspace(-5, 5, args.N, dtype=dtype), S0=args.S0, l=args.l, lM=args.lM, lm=args.lm, workers=args.workers, dtype=dtype)
    G.precompute()
    print(G.Sigma.mean())
    print(G.pdf1[4,10:12,3,5:7])
    print(G.on_grid(i=0, j=0).mean())
    print(G.on_grid().mean())

    exit(0)
    # Quick test: single point
    Gval = G(1.0, 0.0, 0.5, 0.3)
    Gval = G(G.grid[6], G.grid[7], G.grid[8], G.grid[9])
    print("G(1.0, 0.0, 0.5, 0.3) =")
    print(Gval)
    
    t0 = time.perf_counter()
    G.precompute()
    print("Precomputation done")
    time.sleep(4)
    G_all = np.empty((1, args.N, args.N, args.N, args.N), dtype=dtype)
    G.on_grid(i=0, j=0, out=G_all[0])
    # G(i=0, j=1, out=G_all[1])
    # G(i=1, j=1, out=G_all[2])

    # G_all = G.on_grid()

    print("Full G computed")
    
    dt = time.perf_counter() - t0
    print(f"\nG on {args.N}^4 grid: shape {1}, {dt:.2f} s")
    # print(f"Memory: {G_all.nbytes / 1024**3:.2f} GB")

    print("Value from grid calculation:")
    # print(G_all[6, 7, 8, 9])
    exit(0)

    # Create a partial function with default arguments
    g_func = partial(G, S0=args.S0, l=args.l, lM=args.lM, lm=args.lm)

    # Quick test: single point
    Gval = g_func(1.0, 0.0, 0.5, 0.3)
    print("G(1,0, 0.5,0.3) =")
    print(Gval)

    # Grid computation
    pts = np.linspace(-5, 5, args.N, dtype=dtype)
    N = args.N

    x1 = pts[:, None, None, None]
    y1 = pts[None, :, None, None]
    x2 = pts[None, None, :, None]
    y2 = pts[None, None, None, :]

    t0 = time.perf_counter()

    if args.parallel:
        n_workers = args.workers or os.cpu_count()
        if n_workers == 0:
            n_workers = N
        G_all = np.empty((N, N, N, N, 2, 2), dtype=dtype)
        A_all = np.empty((2, N, N, 1, 1), dtype=dtype)

        # Split along x1 axis; numpy releases the GIL so threads run concurrently
        chunk_indices = np.array_split(np.arange(N), n_workers)

        def compute_chunk(idx):
            x1c = pts[idx, None, None, None]
            G_all[idx] = g_func(x1c, y1, x2, y2)

        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            futures = [pool.submit(compute_chunk, idx) for idx in chunk_indices]
            for fut in futures:
                fut.result()

        print(f"Parallel computation with {n_workers} threads")
    else:
        
        G_all = g_func(x1, y1, x2, y2)

    dt = time.perf_counter() - t0
    print(f"\nG on {N}^4 grid: shape {G_all.shape}, {dt:.2f} s")
    print(f"Memory: {G_all.nbytes / 1024**3:.2f} GB")