#!/usr/bin/env python3
"""Lyapunov solve with the noise split by source, for lyapunov_full_solver.m.

    A C + C A^T + W N W = 0,   N = sum_k N_k,   W = blockdiag(E (x) E) x 3

One real Schur factorisation A = Z T Z^T serves every right-hand side.  The
smoothing is moved onto the Schur vectors, Z^T W N_k W Z = Zw^T N_k Zw with
Zw = W Z, so the sparse N_k are never densified.

  total     C from the summed N, refined once, residual checked; diagonals of the
            six blocks and full rows at selected points
  sources   for every N_k: the six block diagonals of C_k (no refinement); the
            sum over k is compared with the total
  pathway   drho_d obeys d/dt drho_d = A11 drho_d + xi_rho (density with dQ frozen).
            Var drho_d by the sine basis of A11 = gamma0 Lap; Cov(drho_d, psi) by
            the Sylvester equation A11 X + X A^T + (W N W)_{rho,.} = 0, solved
            with the same Schur factor
  modes     the slowest eigenpairs of A (shift-invert at 0)

Inputs (scratch dir, written by lfWriteScratch): A_*.bin, N_<src>_*.bin (CSR,
0-based), E.bin, Wsym.bin, meta.json.  Outputs: out.npz, result.json and, if
meta["saveC"], C.bin.

Reuses _lyap_tri, _split, _symmetrize_inplace, _eigs_from_schur from
../lyapunov_solver/lyap_solve.py (imported, not modified).
"""

import gc
import json
import os
import sys
import time

_d = sys.argv[1] if len(sys.argv) > 1 else None
if _d and os.path.exists(os.path.join(_d, "meta.json")):
    with open(os.path.join(_d, "meta.json")) as _fh:
        _t = str(int(json.load(_fh).get("threads", 0)) or os.cpu_count() or 1)
    for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[_v] = _t

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lyapunov_solver"))
from lyap_solve import (_lyap_tri, _split, _symmetrize_inplace,  # noqa: E402
                        _eigs_from_schur)

COLBLK = 1024
ROWBLK = 2048
LEAF = 64


def load_csr(d, name, n):
    indptr = np.fromfile(os.path.join(d, name + "_indptr.bin"), dtype=np.int64)
    indices = np.fromfile(os.path.join(d, name + "_indices.bin"), dtype=np.int64)
    data = np.fromfile(os.path.join(d, name + "_data.bin"), dtype=np.float64)
    return sp.csr_matrix((data, indices, indptr), shape=(n, n))


# ---------------------------------------------------------------------------
# W = blockdiag(E (x) E): grid index k = i + j N (x fastest), so a block row
# vector reshapes to (N_j, N_i) and W acts as E X E.
# ---------------------------------------------------------------------------
def apply_W_rows(M, E, N):
    """M W for M (b x n): every row is smoothed, block by block (in place)."""
    b, m = M.shape[0], N * N
    for a in range(3):
        blk = np.array(M[:, a * m:(a + 1) * m], order="C").reshape(b, N, N)   # always a copy
        M[:, a * m:(a + 1) * m] = np.matmul(np.matmul(E, blk), E).reshape(b, m)
    return M


def W_times(Z, E, N, out=None):
    """W Z for Z (n x m), in column strips."""
    n, m = Z.shape
    if out is None:
        out = np.empty_like(Z)
    for j0 in range(0, m, COLBLK):
        j1 = min(j0 + COLBLK, m)
        # a COPY: for a Fortran-ordered Z (the Schur vectors) Z[:, j0:j1].T is
        # already C-contiguous and ascontiguousarray would return a view of Z
        strip = np.array(Z[:, j0:j1].T, order="C")
        apply_W_rows(strip, E, N)
        out[:, j0:j1] = strip.T
        del strip
    return out


def transform(Nk, Zw):
    """F = -(Zw^T Nk Zw), column-blocked; Nk sparse symmetric."""
    n = Zw.shape[0]
    F = np.empty((n, n))
    for j0 in range(0, n, COLBLK):
        j1 = min(j0 + COLBLK, n)
        tmp = Nk @ Zw[:, j0:j1]
        np.matmul(Zw.T, tmp, out=F[:, j0:j1])
        del tmp
    np.negative(F, out=F)
    _symmetrize_inplace(F)
    return F


def block_diags(U, Z, N):
    """Diagonals of C = U Z^T for the six blocks (rr, 11, 22, r1, r2, 12)."""
    m = N * N
    out = {}
    pairs = {"rr": (0, 0), "11": (1, 1), "22": (2, 2),
             "r1": (0, 1), "r2": (0, 2), "12": (1, 2)}
    for key, (a, b) in pairs.items():
        out[key] = np.einsum("ij,ij->i", U[a * m:(a + 1) * m], Z[b * m:(b + 1) * m])
    return out


def resid_strips(A, Ntot, E, N, C, R=None):
    """max|A C + C A^T + W N W|; fills R with the residual when given."""
    n = C.shape[0]
    mx = 0.0
    for i0 in range(0, n, ROWBLK):
        i1 = min(i0 + ROWBLK, n)
        # rows of W: smooth the identity strip
        Wr = np.zeros((i1 - i0, n))
        Wr[np.arange(i1 - i0), np.arange(i0, i1)] = 1.0
        apply_W_rows(Wr, E, N)                      # rows i0:i1 of W
        M = np.asarray((Ntot @ Wr.T).T)             # (W N)[i0:i1]
        apply_W_rows(M, E, N)                       # (W N W)[i0:i1]
        M += np.asarray(A[i0:i1] @ C)
        M += np.asarray(A @ C[i0:i1].T).T
        mx = max(mx, float(np.max(np.abs(M))))
        if R is not None:
            R[i0:i1] = M
        del Wr, M
    return mx


def sylv_diag(lam, T, F):
    """Solve diag(lam) X + X T^T = F in place (T upper quasi-triangular)."""
    q = T.shape[0]
    if q > LEAF:
        k = _split(T, q)
        sylv_diag(lam, T[k:, k:], F[:, k:])
        F[:, :k] -= F[:, k:] @ T[:k, k:].T
        sylv_diag(lam, T[:k, :k], F[:, :k])
        return F
    j = q - 1
    while j >= 0:
        if j > 0 and T[j, j - 1] != 0.0:              # 2x2 block (j-1, j)
            i = j - 1
            if j + 1 < q:
                ri = F[:, i] - F[:, j + 1:] @ T[i, j + 1:]
                rj = F[:, j] - F[:, j + 1:] @ T[j, j + 1:]
            else:
                ri, rj = F[:, i].copy(), F[:, j].copy()
            a11 = lam + T[i, i]; a12 = T[i, j]
            a21 = T[j, i];       a22 = lam + T[j, j]
            det = a11 * a22 - a12 * a21
            F[:, i] = (a22 * ri - a12 * rj) / det
            F[:, j] = (a11 * rj - a21 * ri) / det
            j -= 2
        else:
            r = F[:, j] - (F[:, j + 1:] @ T[j, j + 1:] if j + 1 < q else 0.0)
            F[:, j] = r / (lam + T[j, j])
            j -= 1
    return F


def sine_basis(N):
    i = np.arange(1, N + 1)
    return np.sqrt(2.0 / (N + 1)) * np.sin(np.pi * np.outer(i, i) / (N + 1))


def kron_rows(M, S, N, left=True):
    """For M (m x c) with m = N^2: (S (x) S)^T M  (left) on the row index."""
    c = M.shape[1]
    X = M.reshape(N, N, c)                          # [j, i, c]
    X = np.tensordot(S.T, X, axes=(1, 0))           # [q, i, c]
    X = np.tensordot(S.T, X, axes=(1, 1))           # [p, q, c]
    return np.ascontiguousarray(X.transpose(1, 0, 2)).reshape(N * N, c)


def main():
    d = sys.argv[1]
    res = {"ok": False}
    tm = {}
    try:
        with open(os.path.join(d, "meta.json")) as fh:
            meta = json.load(fh)
        n, N, h, g0 = int(meta["n"]), int(meta["nInt"]), float(meta["h"]), float(meta["gamma0"])
        m = N * N
        srcs = meta["sources"]
        E = np.fromfile(os.path.join(d, "E.bin")).reshape(N, N)
        wsym = np.fromfile(os.path.join(d, "Wsym.bin"))
        S = sine_basis(N)
        res["E_vs_sine"] = float(np.max(np.abs(E - (S * wsym) @ S.T)))

        t0 = time.perf_counter()
        A = load_csr(d, "A", n)
        Nk = {s: load_csr(d, "N_" + s, n) for s in srcs}
        Ntot = sum(Nk.values())
        tm["read"] = time.perf_counter() - t0

        # the density block must be gamma0 Lap for the pathway split
        mu = -(4.0 / h**2) * np.sin(np.pi * np.arange(1, N + 1) / (2 * (N + 1)))**2
        v = np.random.default_rng(1).standard_normal(m)
        A11v = A[:m, :m] @ v
        # (S (x) S) diag(lam) (S (x) S)^T v with the [j,i] layout of kron_rows
        lamPQ = g0 * (mu[None, :] + mu[:, None])                # [q(j), p(i)] symmetric
        vt = kron_rows(v[:, None], S, N)[:, 0] * lamPQ.ravel()
        A11v_s = kron_rows(vt[:, None], S, N)[:, 0]             # S symmetric
        res["A11_vs_sine"] = float(np.max(np.abs(A11v - A11v_s)) / np.max(np.abs(A11v)))
        lam11 = lamPQ.ravel()
        wPQ = (wsym[None, :] * wsym[:, None]).ravel()

        # slowest modes (cheap, sparse shift-invert)
        t0 = time.perf_counter()
        try:
            kmodes = int(meta.get("nModes", 8))
            ev, evec = spla.eigs(A.tocsc(), k=kmodes, sigma=0.0, which="LM")
            o = np.argsort(-ev.real)
            ev, evec = ev[o], evec[:, o]
        except Exception as exc:                                  # noqa: BLE001
            ev, evec = np.zeros(0, complex), np.zeros((n, 0), complex)
            res["modes_error"] = str(exc)
        tm["modes"] = time.perf_counter() - t0

        # Schur
        t0 = time.perf_counter()
        Ad = A.toarray()
        T, Z = sla.schur(Ad, output="real", overwrite_a=True, check_finite=False)
        del Ad
        gc.collect()
        tm["schur"] = time.perf_counter() - t0
        lamA = _eigs_from_schur(T)

        out = {"gridInt": np.array(meta["gridInt"]), "eig": lamA,
               "modes_val": ev, "modes_vec": evec.astype(np.complex64)}

        # ---- total -------------------------------------------------------
        t0 = time.perf_counter()
        Zw = W_times(Z, E, N)
        F = transform(Ntot, Zw)
        tm["total_transform"] = time.perf_counter() - t0
        t0 = time.perf_counter()
        _lyap_tri(T, F)
        _symmetrize_inplace(F)
        U = Z @ F
        C = np.dot(U, Z.T, out=F)
        del U
        gc.collect()
        _symmetrize_inplace(C)
        tm["total_solve"] = time.perf_counter() - t0
        t0 = time.perf_counter()
        R = np.empty((n, n))
        r1 = resid_strips(A, Ntot, E, N, C, R)
        res["residual"] = r1
        if meta.get("refine", True):
            W = Z.T @ R
            dF = np.dot(W, Z, out=R)
            del W
            np.negative(dF, out=dF)
            _symmetrize_inplace(dF)
            _lyap_tri(T, dF)
            _symmetrize_inplace(dF)
            W = Z @ dF
            dC = np.dot(W, Z.T, out=dF)
            del W
            _symmetrize_inplace(dC)
            C += dC
            del dC
        del R
        gc.collect()
        res["residual_refined"] = resid_strips(A, Ntot, E, N, C)
        res["scale_C"] = float(np.max(np.abs(np.diag(C))))
        tm["total_residual_refine"] = time.perf_counter() - t0

        I = np.arange(m)
        dC = np.diag(C).copy()
        out["total"] = np.stack([dC[:m], dC[m:2 * m], dC[2 * m:],
                                 C[I, I + m], C[I, I + 2 * m], C[I + m, I + 2 * m]])
        # rows at selected points: the nodes nearest to each requested point,
        # plus the extrema of Var drho within r < 5
        g = np.array(meta["gridInt"])
        X, Y = np.meshgrid(g, g)                       # [j, i]
        Xf, Yf = X.ravel(), Y.ravel()
        pts = [tuple(p) for p in meta.get("points", [])]
        idx = [int(np.argmin((Xf - px)**2 + (Yf - py)**2)) for px, py in pts]
        near = (Xf**2 + Yf**2) < 25.0
        vr = out["total"][0]
        imax = int(np.flatnonzero(near)[np.argmax(vr[near])])
        imin = int(np.flatnonzero(near)[np.argmin(vr[near])])
        idx += [imax, imin]
        out["row_index"] = np.array(idx)
        out["row_xy"] = np.array([[Xf[k], Yf[k]] for k in idx])
        out["rows"] = np.stack([np.stack([C[k], C[k + m], C[k + 2 * m]]) for k in idx])
        if meta.get("saveC", False):
            t0 = time.perf_counter()
            C.tofile(os.path.join(d, "C.bin"))
            tm["write_C"] = time.perf_counter() - t0
        del C
        gc.collect()

        # ---- sources -----------------------------------------------------
        sumd = np.zeros_like(out["total"])
        for s in srcs:
            t0 = time.perf_counter()
            F = transform(Nk[s], Zw)
            _lyap_tri(T, F)
            _symmetrize_inplace(F)
            U = Z @ F
            del F
            gc.collect()
            bd = block_diags(U, Z, N)
            del U
            gc.collect()
            arr = np.stack([bd[k] for k in ("rr", "11", "22", "r1", "r2", "12")])
            out["src_" + s] = arr
            sumd += arr
            tm["source_" + s] = time.perf_counter() - t0
        res["sources_sum_vs_total"] = float(np.max(np.abs(sumd - out["total"])) /
                                            np.max(np.abs(out["total"][0])))

        # ---- pathway -----------------------------------------------------
        t0 = time.perf_counter()
        # Var drho_d: C_d = V Ct V^T, Ct_pq = -w_p w_q (V^T Nrr V)_pq / (lam_p + lam_q)
        Nrr = Ntot[:m, :m].toarray()
        Mt = kron_rows(Nrr, S, N)                     # V^T Nrr
        Mt = kron_rows(np.ascontiguousarray(Mt.T), S, N).T   # V^T Nrr V
        Ct = -(wPQ[:, None] * Mt * wPQ[None, :]) / (lam11[:, None] + lam11[None, :])
        del Mt, Nrr
        Vs = kron_rows(np.eye(m), S, N)               # V (m x m), symmetric
        Cd_half = kron_rows(Ct, S, N)                 # V Ct
        out["var_d"] = np.einsum("ij,ij->i", Cd_half, Vs)
        del Cd_half, Ct
        gc.collect()
        # X = Cov(drho_d, psi) = V Xt Z^T,  Lam Xt + Xt T^T = -Dw V^T Nrr. Zw
        G = np.asarray(Ntot[:m, :] @ Zw)              # m x n
        G = kron_rows(G, S, N)                        # V^T (N Zw)_rho
        G *= -wPQ[:, None]
        sylv_diag(lam11, T, G)                        # G <- Xt
        VX = kron_rows(G, S, N)                       # V Xt  (m x n)
        del G
        out["cov_d"] = np.stack([np.einsum("ij,ij->i", VX, Z[a * m:(a + 1) * m])
                                 for a in range(3)])  # Cov(drho_d, drho / dQ1 / dQ2)
        del VX, Vs
        gc.collect()
        tm["pathway"] = time.perf_counter() - t0

        out["timings_keys"] = np.array(list(tm.keys()))
        np.savez(os.path.join(d, "out.npz"), **out)
        res.update({"ok": True, "n": n, "max_re_lambda": float(lamA.real.max()),
                    "timings": tm, "peak_rss_gib": _peak_rss_gib(),
                    "threads": os.environ.get("OPENBLAS_NUM_THREADS"),
                    "numpy": np.__version__})
    except Exception as exc:                                      # noqa: BLE001
        import traceback
        res.update({"ok": False, "error": "%s: %s" % (type(exc).__name__, exc),
                    "traceback": traceback.format_exc(), "timings": tm})
    with open(os.path.join(d, "result.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    if not res["ok"]:
        sys.stderr.write(res.get("traceback", res.get("error", "")) + "\n")
        return 1
    return 0


def _peak_rss_gib():
    try:
        import resource
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20
    except Exception:                                             # noqa: BLE001
        return None


if __name__ == "__main__":
    sys.exit(main())
