"""Flexible cells: the spine is a free chain of M-1 segment angles.

Outer parts of the cell follow the director (strong alignment + tilt cap).
The middle part (|u| < MID_HALF) may bend freely, so a cell can break and
switch to a neighbouring streamline.

Parameters per cell: v = [cx, cy, phi_0 .. phi_{M-2}].
"""
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree
import pack
from pack import (W, L, w, r, ELL, M, A, U, R_EFF, director_angle, dense_spines,
                  min_spine_distances, largest_gaps, free_length, in_box_mask)

K = M - 1
MID_HALF = 0.15
MIDMASK = np.abs(U) < MID_HALF              # segments in the breakable middle
KAPPA_MAX = 1.0 / r                          # sharpest allowed bend: radius r

P = dict(k_o=1e5, k_w=1e5, lam_a=10.0, lam_mid=0.5, lam_b=0.005, k_c=10.0,
         g=0.3, k_t=300.0, tilt_max=np.radians(10))


def geom(V, ell):
    phi = V[:, 2:]
    h = ell / K
    seg = h[:, None, None] * np.stack([np.cos(phi), np.sin(phi)], -1)
    X = V[:, None, :2] + np.einsum('jk,nkd->njd', A, seg)
    return X, phi, h


def from_q(q, ell):
    _, phi, _, _ = pack.geom(q, ell)
    return np.hstack([q[:, :2], phi])


def energy(V, ell, fixed, p=P, want_grad=True):
    N = len(V)
    X, phi, h = geom(V, ell)
    GX = np.zeros_like(X)
    E = 0.0
    # overlap
    flat = X.reshape(-1, 2)
    cid = np.repeat(np.arange(N), M)
    D = 2 * R_EFF
    pairs = cKDTree(flat).query_pairs(D, output_type='ndarray')
    if len(pairs):
        pairs = pairs[cid[pairs[:, 0]] != cid[pairs[:, 1]]]
    if len(pairs):
        dv = flat[pairs[:, 0]] - flat[pairs[:, 1]]
        d = np.sqrt((dv ** 2).sum(1)) + 1e-12
        ov = D - d
        E += p['k_o'] * (ov ** 2).sum()
        g = (-2 * p['k_o'] * ov / d)[:, None] * dv
        gf = np.zeros_like(flat)
        np.add.at(gf, pairs[:, 0], g)
        np.add.at(gf, pairs[:, 1], -g)
        GX += gf.reshape(X.shape)
    # wall
    ex = np.abs(X) - (pack.EXT - r)
    m = ex > 0
    E += p['k_w'] * (ex[m] ** 2).sum()
    GX += np.where(m, 2 * p['k_w'] * ex * np.sign(X), 0.0)
    # pull toward the core
    if p['g'] > 0:
        rn = np.sqrt((X ** 2).sum(-1)) + 1e-9
        E += p['g'] * (h[:, None] * rn).sum()
        GX += p['g'] * h[:, None, None] * X / rn[..., None]
    # alignment (weak in the middle) and tilt cap (outer parts only)
    mids = 0.5 * (X[:, :-1] + X[:, 1:])
    psi = director_angle(mids[..., 0], mids[..., 1])
    dl = phi - psi
    hh = h[:, None]
    lam = np.where(MIDMASK, p['lam_mid'], p['lam_a'])[None, :]
    E += (lam * hh * np.sin(dl) ** 2).sum()
    ddl = lam * hh * np.sin(2 * dl)
    sd = np.sin(dl)
    exc = np.abs(sd) - np.sin(p['tilt_max'])
    mt = (exc > 0) & ~MIDMASK[None, :]
    E += p['k_t'] * (hh * np.where(mt, exc, 0) ** 2).sum()
    ddl = ddl + np.where(mt, 2 * p['k_t'] * hh * exc * np.sign(sd) * np.cos(dl), 0)
    rho2 = (mids ** 2).sum(-1) + 1e-4
    gm = (-ddl)[..., None] * np.stack([-mids[..., 1], mids[..., 0]], -1) / (2 * rho2[..., None])
    GX[:, :-1] += 0.5 * gm
    GX[:, 1:] += 0.5 * gm
    dphi = ddl.copy()
    # bending and curvature cap
    dp = np.diff(phi, axis=1)
    E += p['lam_b'] * (dp ** 2 / hh).sum()
    gb = 2 * p['lam_b'] * dp / hh
    kap = np.abs(dp) / hh
    ek = kap - KAPPA_MAX
    mk = ek > 0
    E += p['k_c'] * (hh * np.where(mk, ek, 0) ** 2).sum()
    gb = gb + np.where(mk, 2 * p['k_c'] * ek * np.sign(dp), 0)
    dphi[:, 1:] += gb
    dphi[:, :-1] -= gb
    if not want_grad:
        return E
    perp = np.stack([-np.sin(phi), np.cos(phi)], -1)
    ATG = np.einsum('jk,njd->nkd', A, GX)
    dphi = dphi + hh * (ATG * perp).sum(-1)
    gV = np.zeros_like(V)
    gV[:, :2] = GX.sum(1)
    gV[:, 2:] = dphi
    gV[fixed] = 0.0
    return E, gV


def relax(V, ell, fixed, maxiter=400, p=P):
    free = ~fixed
    nv = V.shape[1]

    def f(x):
        VV = V.copy()
        VV[free] = x.reshape(-1, nv)
        E, g = energy(VV, ell, fixed, p)
        return E, g[free].ravel()

    res = minimize(f, V[free].ravel(), jac=True, method='L-BFGS-B',
                   options=dict(maxiter=maxiter, maxcor=30, gtol=1e-9, ftol=1e-13))
    V = V.copy()
    V[free] = res.x.reshape(-1, nv)
    return V


def valid(V, ell):
    X = geom(V, ell)[0]
    ds = [d for _, _, d in min_spine_distances(X)]
    dmin = min(ds) if ds else np.inf
    return dmin >= w * (1 - 1e-9), dmin


def cover(V, ell, h=0.01):
    X = geom(V, ell)[0]
    g = np.arange(-W + h / 2, W, h)
    gx, gy = np.meshgrid(g, g)
    d, _ = cKDTree(dense_spines(X).reshape(-1, 2)).query(np.stack([gx.ravel(), gy.ravel()], 1))
    return (d <= r).mean()


def try_insert(V, ell, fixed, pc, grow_steps=12, max_da=np.radians(10)):
    X = geom(V, ell)[0]
    psi = director_angle(*pc)
    best = None
    for da in np.linspace(-max_da, max_da, 9):
        a = psi + da
        fl = free_length(pc, a, X)
        score = min(fl) + 0.5 * max(fl)
        if best is None or score > best[0]:
            best = (score, a, fl)
    _, a, fl = best
    c0 = pc + 0.5 * (fl[0] - fl[1]) * np.array([np.cos(a), np.sin(a)])
    Vn = np.vstack([V, np.concatenate([c0, np.full(K, a)])])
    fx = np.append(fixed, False)
    ell0 = max(0.0, min(ELL, fl[0] + fl[1]))
    lengths = np.linspace(ell0, ELL, grow_steps + 1)[1:] if ell0 < ELL else [ELL]
    for Ln in lengths:
        Vn = relax(Vn, np.append(ell, Ln), fx, maxiter=150)
    en = np.append(ell, ELL)
    Vn = relax(Vn, en, fx, maxiter=800)
    ok, _ = valid(Vn, en)
    return ok, Vn, en, fx


def prune_outside(V, ell, fixed):
    X = geom(V, ell)[0]
    keep = fixed.copy()
    for i in range(len(V)):
        dx = np.maximum(np.abs(X[i]) - W, 0)
        if (np.sqrt((dx ** 2).sum(-1)) < r).any():
            keep[i] = True
    return V[keep], ell[keep], fixed[keep]


def metrics(V, ell):
    X, phi, _ = geom(V, ell)
    mids = 0.5 * (X[:, :-1] + X[:, 1:])
    inb = in_box_mask(mids)
    dl = phi - director_angle(mids[..., 0], mids[..., 1])
    dev = np.degrees(np.abs((dl + np.pi / 2) % np.pi - np.pi / 2))
    outer = inb & ~MIDMASK[None, :]
    mid = inb & MIDMASK[None, :]
    gaps = largest_gaps(X, n=5, h=0.004)
    pc, cmax = gaps[0] if gaps else (None, 0.0)
    ok, dmin = valid(V, ell)
    return dict(n=len(V), cover=cover(V, ell, 0.004), max_gap_over_w=2 * max(cmax, 0) / w,
                gap_at=pc, dev_outer_mean=dev[outer].mean(), dev_outer_max=dev[outer].max(),
                dev_mid_max=dev[mid].max() if mid.any() else 0.0,
                n_broken=int((dev * mid).max(1).__gt__(15).sum()), valid=ok, dmin_over_w=dmin / w)


def plot(V, ell, fname, title=None, gap=None):
    X = geom(V, ell)[0]
    pack.plot(None, None, fname, title=title, gap=gap, X=X)
