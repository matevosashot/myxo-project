#!/usr/bin/env python3
"""Dense continuous-time Lyapunov solve for lyapunov_solver_optimized.m.

Solves  A C + C A^T + Q = 0  for symmetric C, by real Schur decomposition plus a
recursive (GEMM-bound) Bartels-Stewart back-substitution.

Why not the obvious routes
--------------------------
* scipy.linalg.solve_continuous_lyapunov calls LAPACK ?trsyl, which is serial and
  unblocked.  Measured on an EPYC 9655: 8.9 s at n=2000 against 2.9 s for the
  eigendecomposition route -- i.e. slower than what it would replace.
* LAPACK's blocked ?trsyl3 would fix that, but scipy 1.18 does not expose it.
The recursion below splits the triangular Sylvester equation until the leaves are
small enough for ?trsyl, leaving the bulk of the flops in GEMM.  Measured 3.9x
faster than the Mathematica eigendecomposition route at n=4000, and it stays in
Real64 throughout, where the eigen route needs six Complex128 n x n temporaries.

Contract with the caller (all paths under one scratch directory)
----------------------------------------------------------------
  in   A.bin, Q.bin   row-major Real64, n x n (Q symmetric)
       meta.json      {"n": int, "neumann": false, "refine": bool, "threads": int}
  out  C.bin          row-major Real64, n x n (symmetric)
       eig.bin        2n Real64, interleaved (Re, Im) of the eigenvalues of A
       result.json    diagnostics; {"ok": false, "error": ...} on failure

Exit status is 0 only if result.json says ok.  Neumann is rejected: A is then
singular by design and the deflation is not implemented here (the caller falls
back to the reference kernel).
"""

import gc
import json
import os
import sys
import time

# Thread count must be set before numpy is imported, or OpenBLAS has already
# sized its pool.  The caller passes it because SLURM's allocation, not the
# node's core count, is what we are entitled to.
_meta_path = os.path.join(sys.argv[1], "meta.json") if len(sys.argv) > 1 else None
if _meta_path and os.path.exists(_meta_path):
    with open(_meta_path) as _fh:
        _meta_early = json.load(_fh)
    _t = str(int(_meta_early.get("threads", 0)) or os.cpu_count() or 1)
    for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[_v] = _t

import numpy as np
import scipy.linalg as sla
from scipy.linalg.lapack import dtrsyl

# Leaf size for the recursion.  Below this the serial ?trsyl is cheaper than the
# extra GEMM traffic from splitting further; 192 measured best on Zen5, but the
# curve is flat between ~128 and ~384 so this is not a tuned constant.
BLK = 192


def _split(T, m):
    """Midpoint of [0, m) that does not cut a 2x2 block of the real Schur form.

    Splitting inside a conjugate pair would hand ?trsyl a matrix that is not
    quasi-upper-triangular and silently corrupt that block.
    """
    k = m // 2
    if 0 < k < m and T[k, k - 1] != 0.0:
        k += 1
    return k


def _sylv_t(T, U, F):
    """Solve  T X + X U^T = F  with T, U upper quasi-triangular."""
    p, q = F.shape
    if p <= BLK and q <= BLK:
        X, scale, info = dtrsyl(T, U, F, trana='N', tranb='T', isgn=1)
        if info < 0:
            raise RuntimeError("dtrsyl: illegal argument %d" % -info)
        return X / scale
    if p >= q:
        k = _split(T, p)
        X2 = _sylv_t(T[k:, k:], U, F[k:, :])
        X1 = _sylv_t(T[:k, :k], U, F[:k, :] - T[:k, k:] @ X2)
        return np.vstack([X1, X2])
    k = _split(U, q)
    X2 = _sylv_t(T, U[k:, k:], F[:, k:])
    X1 = _sylv_t(T, U[:k, :k], F[:, :k] - X2 @ U[:k, k:].T)
    return np.hstack([X1, X2])


def _lyap_tri(T, F, out=None):
    """Solve  T Y + Y T^T = F  with T upper quasi-triangular and F symmetric.

    Writes the answer into `out` (default: into F itself).  The obvious version
    returns np.block([[Y11, Y12], [Y12.T, Y22]]), which allocates a fresh n x n
    at EVERY level of the recursion -- the single largest avoidable allocation in
    this file.  Recursing into views of one preallocated array removes all of it.
    """
    n = F.shape[0]
    if out is None:
        out = F
    if n <= BLK:
        Y, scale, info = dtrsyl(T, T, F, trana='N', tranb='T', isgn=1)
        if info < 0:
            raise RuntimeError("dtrsyl: illegal argument %d" % -info)
        np.divide(Y, scale, out=out)
        return out
    k = _split(T, n)
    T11, T12, T22 = T[:k, :k], T[:k, k:], T[k:, k:]

    # (2,2) is independent; (1,2) needs it; (1,1) needs both.  Order matters:
    # each sub-block of F is consumed before the corresponding block of out is
    # written, so aliasing out onto F is safe.
    Y22 = _lyap_tri(T22, F[k:, k:], out[k:, k:])

    G12 = F[:k, k:] - T12 @ Y22
    Y12 = _sylv_t(T11, T22, G12)
    del G12

    G = F[:k, :k]
    G = G - T12 @ Y12.T
    G -= Y12 @ T12.T
    G += G.T.copy()
    G *= 0.5
    _lyap_tri(T11, G, out[:k, :k])
    del G

    out[:k, k:] = Y12
    out[k:, :k] = Y12.T
    return out


def _eigs_from_schur(T):
    """Eigenvalues off the real Schur form: free, no second O(n^3) pass."""
    n = T.shape[0]
    lam = np.zeros(n, dtype=np.complex128)
    i = 0
    while i < n:
        if i + 1 < n and T[i + 1, i] != 0.0:
            a, b, c, d = T[i, i], T[i, i + 1], T[i + 1, i], T[i + 1, i + 1]
            tr, det = a + d, a * d - b * c
            disc = np.sqrt(complex((tr / 2) ** 2 - det))
            lam[i], lam[i + 1] = tr / 2 + disc, tr / 2 - disc
            i += 2
        else:
            lam[i] = T[i, i]
            i += 1
    return lam


BLK_SYM = 2048          # block size for the out-of-place-free symmetrizer
RESID_BLK = 2048        # row strip for the streamed residual
TRANSFORM_BLK = 1024    # column block for the streamed Z^T Q Z


def _symmetrize_inplace(M, blk=BLK_SYM):
    """M <- (M + M^T)/2 using only blk x blk temporaries.

    M += M.T.copy() would allocate a second full n x n at the worst possible
    moment (34 GiB at Nint = 150).
    """
    n = M.shape[0]
    for i in range(0, n, blk):
        ii = min(i + blk, n)
        d = M[i:ii, i:ii]
        t = d + d.T
        t *= 0.5
        d[:] = t
        del t
        for j in range(ii, n, blk):
            jj = min(j + blk, n)
            a, b = M[i:ii, j:jj], M[j:jj, i:ii]
            t = a + b.T
            t *= 0.5
            a[:] = t
            b[:] = t.T
            del t
    return M


def _resid_max(Apath, Qpath, C, n):
    """max |A C + C A^T + Q|, streaming A and Q from scratch in row strips.

    A and Q are deleted from RAM before this runs, so the residual costs one
    strip rather than the 3 x 8n^2 that A @ C + C @ A.T + Q would allocate.
    """
    A = np.memmap(Apath, dtype=np.float64, mode="r", shape=(n, n))
    Q = np.memmap(Qpath, dtype=np.float64, mode="r", shape=(n, n))
    m = 0.0
    try:
        for i in range(0, n, RESID_BLK):
            j = min(i + RESID_BLK, n)
            # Both products read the mapped files through dgemm (A.T is handled
            # by transB, not by a transposed copy), so the only anonymous
            # allocation here is the s x n strip r.  Pages touched land in the
            # page cache, which is reclaimable and not anonymous RSS.
            r = np.matmul(A[i:j], C)
            r += np.matmul(C[i:j], A.T)
            r += Q[i:j]
            m = max(m, float(np.max(np.abs(r))))
            del r
    finally:
        del A, Q
    return m


def solve_from_files(d, n, refine=True):
    """Solve out of core.  Returns (C, eigenvalues, timings, resid, resid_ref).

    Peak live set is ~5 x 8n^2.  The ordering below is deliberate: every array is
    dropped at the first point it is provably dead, and A/Q are re-read from
    scratch for the residual rather than held.
    """
    tm = {}
    Apath = os.path.join(d, "A.bin")
    Qpath = os.path.join(d, "Q.bin")

    t0 = time.perf_counter()
    A = _load_A(d, n)
    tm["read_A"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    # overwrite_a lets LAPACK reduce in A's own buffer instead of copying it,
    # saving one full n x n exactly where the peak sits.  Safe because A is
    # re-read from scratch for the residual and never used from RAM again.
    # check_finite=False skips a full O(n^2) scan; the residual would catch a
    # NaN anyway, and the caller aborts on a non-finite one.
    T, Z = sla.schur(A, output="real", overwrite_a=True, check_finite=False)
    tm["schur"] = time.perf_counter() - t0
    del A                                             # now garbage as well as dead
    gc.collect()

    lam = _eigs_from_schur(T)

    t0 = time.perf_counter()
    # F = -(Z^T Q Z), one COLUMN BLOCK at a time with Q streamed from scratch.
    # Holding Q and the Z^T Q product together would put five full n x n arrays
    # live at once (T, Z, Q, W, F); this way Q never enters RAM and the only
    # extra is an n x blk strip.
    F = np.empty((n, n), dtype=np.float64)
    Qmm = np.memmap(Qpath, dtype=np.float64, mode="r", shape=(n, n))
    for j0 in range(0, n, TRANSFORM_BLK):
        j1 = min(j0 + TRANSFORM_BLK, n)
        tmp = np.matmul(Qmm, Z[:, j0:j1])             # n x blk
        np.matmul(Z.T, tmp, out=F[:, j0:j1])
        del tmp
    del Qmm
    gc.collect()
    np.negative(F, out=F)
    _symmetrize_inplace(F)
    tm["transform"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    _lyap_tri(T, F)                                   # in place: Y is F
    tm["back_substitute"] = time.perf_counter() - t0

    if not refine:
        del T                                         # dead unless we refine
        T = None
        gc.collect()

    t0 = time.perf_counter()
    _symmetrize_inplace(F)
    W = Z @ F
    C = np.dot(W, Z.T, out=F)                         # reuse F's storage
    del W
    gc.collect()
    _symmetrize_inplace(C)
    tm["congruence"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    resid = _resid_max(Apath, Qpath, C, n)
    tm["residual_1"] = time.perf_counter() - t0
    if not np.isfinite(resid):
        raise RuntimeError("residual is not finite (%r): the solve produced "
                           "NaN or Inf" % resid)

    # One refinement pass.  Bartels-Stewart is backward stable, but the FORWARD
    # error scales with the conditioning of the Lyapunov operator, and on the
    # real drift one pass alone leaves a residual of ~5e-11 against ~1e-14
    # refined.  A random stable A shows ~2e-14 unrefined and is NOT evidence
    # that this can be skipped.
    resid_ref = resid
    if refine:
        t0 = time.perf_counter()
        Amm = np.memmap(Apath, dtype=np.float64, mode="r", shape=(n, n))
        Qmm = np.memmap(Qpath, dtype=np.float64, mode="r", shape=(n, n))
        Aa = np.asarray(Amm)
        R = Aa @ C
        R += C @ Aa.T
        R += np.asarray(Qmm)
        del Amm, Qmm, Aa
        gc.collect()
        W = Z.T @ R
        dF = np.dot(W, Z, out=R)
        del W
        gc.collect()
        np.negative(dF, out=dF)
        _symmetrize_inplace(dF)
        _lyap_tri(T, dF)
        _symmetrize_inplace(dF)
        W = Z @ dF
        dC = np.dot(W, Z.T, out=dF)
        del W
        gc.collect()
        _symmetrize_inplace(dC)
        C += dC
        del dC, R
        gc.collect()
        tm["refine"] = time.perf_counter() - t0
        t0 = time.perf_counter()
        resid_ref = _resid_max(Apath, Qpath, C, n)
        tm["residual_2"] = time.perf_counter() - t0

    del T, Z
    gc.collect()
    return C, lam, tm, resid, resid_ref


def _load_A(d, n):
    """A is written either dense (A.bin) or as CSR (A_indptr/indices/data.bin).

    The sparse form is what lyapunov_solver_optimized.m writes: the assembled
    drift has ~13 nonzeros per row, so CSR is ~10 MB where the dense form is
    34 GiB at Nint = 150, and it spares Mathematica the Normal[] copy entirely.
    The dense form is still accepted so the script stays usable on its own.
    """
    ip = os.path.join(d, "A_indptr.bin")
    if os.path.exists(ip):
        import scipy.sparse as sp
        indptr = np.fromfile(ip, dtype=np.int64)
        indices = np.fromfile(os.path.join(d, "A_indices.bin"), dtype=np.int64)
        data = np.fromfile(os.path.join(d, "A_data.bin"), dtype=np.float64)
        Asp = sp.csr_matrix((data, indices, indptr), shape=(n, n))
        del data, indices, indptr
        A = Asp.toarray()
        del Asp
        gc.collect()
        # the residual streams A back off disk, so it must exist densely there
        if not os.path.exists(os.path.join(d, "A.bin")):
            A.tofile(os.path.join(d, "A.bin"))
        return A
    return np.fromfile(os.path.join(d, "A.bin"), dtype=np.float64).reshape(n, n)


def solve(A, Q, refine=True):
    """In-memory entry point, kept for verify_optimized.wls and direct use."""
    import tempfile
    n = A.shape[0]
    d = tempfile.mkdtemp(prefix="lyap_mem_")
    try:
        A.tofile(os.path.join(d, "A.bin"))
        Q.tofile(os.path.join(d, "Q.bin"))
        return solve_from_files(d, n, refine=refine)
    finally:
        import shutil
        shutil.rmtree(d, ignore_errors=True)


def main():
    if len(sys.argv) != 2:
        sys.stderr.write("usage: lyap_solve.py <scratch-dir>\n")
        return 2
    d = sys.argv[1]
    res = {"ok": False}
    try:
        with open(os.path.join(d, "meta.json")) as fh:
            meta = json.load(fh)
        n = int(meta["n"])
        if meta.get("neumann"):
            raise RuntimeError("Neumann is not supported by this backend: A is "
                               "singular by design and the null-mode deflation "
                               "is not implemented")

        C, lam, tm, resid, resid_ref = solve_from_files(
            d, n, refine=bool(meta.get("refine", True)))

        t0 = time.perf_counter()
        C.tofile(os.path.join(d, "C.bin"))
        np.stack([lam.real, lam.imag], axis=1).ravel().tofile(
            os.path.join(d, "eig.bin"))
        tm["write"] = time.perf_counter() - t0

        res = {
            "ok": True,
            "n": n,
            "residual": resid,
            "residual_refined": resid_ref,
            "max_re_lambda": float(lam.real.max()),
            "timings": tm,
            "peak_rss_gib": _peak_rss_gib(),
            "threads": os.environ.get("OPENBLAS_NUM_THREADS"),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "backend": "python/schur+recursive-bartels-stewart",
        }
    except Exception as exc:                                  # noqa: BLE001
        import traceback
        res = {"ok": False, "error": "%s: %s" % (type(exc).__name__, exc),
               "traceback": traceback.format_exc()}
    with open(os.path.join(d, "result.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    if not res["ok"]:
        sys.stderr.write(res.get("traceback", res["error"]) + "\n")
        return 1
    return 0


def _peak_rss_gib():
    try:
        import resource
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20
    except Exception:
        return None


if __name__ == "__main__":
    sys.exit(main())
