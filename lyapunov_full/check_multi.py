#!/usr/bin/env python3
"""check_multi.py <scratch-dir> -- the source split and the pathway split of
lyap_multi.py against independent dense solves (small grids only).

  T3  every source: Lyapunov solve with W N_k W alone (lyap_solve.solve) against
      the six block diagonals in out.npz; and sum over sources = total
  T4  pathway: the augmented system (drho_d, psi) with drift diag(A11, A) and
      noise built from the same W N W, solved densely, against var_d and cov_d
"""
import json
import os
import sys

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lyapunov_solver"))
from lyap_solve import solve  # noqa: E402

d = sys.argv[1]
meta = json.load(open(os.path.join(d, "meta.json")))
n, N = meta["n"], meta["nInt"]
m = N * N


def csr(name):
    ip = np.fromfile(os.path.join(d, name + "_indptr.bin"), dtype=np.int64)
    ix = np.fromfile(os.path.join(d, name + "_indices.bin"), dtype=np.int64)
    da = np.fromfile(os.path.join(d, name + "_data.bin"))
    return sp.csr_matrix((da, ix, ip), shape=(n, n)).toarray()


A = csr("A")
E = np.fromfile(os.path.join(d, "E.bin")).reshape(N, N)
W2 = np.kron(E, E)
Wb = np.zeros((n, n))
for a in range(3):
    Wb[a * m:(a + 1) * m, a * m:(a + 1) * m] = W2
out = np.load(os.path.join(d, "out.npz"))
I = np.arange(m)


def diags(C):
    dC = np.diag(C)
    return np.stack([dC[:m], dC[m:2 * m], dC[2 * m:],
                     C[I, I + m], C[I, I + 2 * m], C[I + m, I + 2 * m]])


scale = np.max(np.abs(out["total"][0]))
Ntot = np.zeros((n, n))
worst = 0.0
for s in meta["sources"]:
    Nk = csr("N_" + s)
    Ntot += Nk
    C, *_ = solve(A, Wb @ Nk @ Wb, refine=True)
    e = np.max(np.abs(diags(C) - out["src_" + s])) / scale
    worst = max(worst, e)
    print("T3  %-9s max|diag - dense| / max Var_rho = %.2e" % (s, e))
Ct, *_ = solve(A, Wb @ Ntot @ Wb, refine=True)
print("T3  total     max|diag - dense| / max Var_rho = %.2e"
      % (np.max(np.abs(diags(Ct) - out["total"])) / scale))

# augmented system (drho_d, drho, dQ1, dQ2)
Nw = Wb @ Ntot @ Wb
Aaug = np.zeros((n + m, n + m))
Aaug[:m, :m] = A[:m, :m]
Aaug[m:, m:] = A
Naug = np.zeros((n + m, n + m))
Naug[:m, :m] = Nw[:m, :m]
Naug[:m, m:] = Nw[:m, :]
Naug[m:, :m] = Nw[:, :m]
Naug[m:, m:] = Nw
Ca, *_ = solve(Aaug, Naug, refine=True)
var_d = np.diag(Ca)[:m]
cov_d = np.stack([Ca[I, m + a * m + I] for a in range(3)])
print("T4  var_d     rel.err = %.2e" % (np.max(np.abs(var_d - out["var_d"])) / scale))
print("T4  cov_d     rel.err = %.2e" % (np.max(np.abs(cov_d - out["cov_d"])) / scale))
print("result.json:", {k: v for k, v in json.load(open(os.path.join(d, "result.json"))).items()
                       if k not in ("timings",)})
