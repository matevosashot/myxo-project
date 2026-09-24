"""Seed cells exactly along streamlines, random-sequential in order of distance from core."""
import numpy as np
from scipy.spatial import cKDTree
from pack import *

def parabola_arc(p, ymax=4.0, n=8001):
    y = np.linspace(-ymax, ymax, n)
    if p == 0:
        x = np.linspace(0, 4, n); y = 0 * x
    else:
        x = (y ** 2 - p ** 2) / (2 * p)
    pts = np.stack([x, y], 1)
    ds = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    s = np.concatenate([[0], np.cumsum(ds)])
    t = np.gradient(pts, axis=0); phi = np.unwrap(np.arctan2(t[:, 1], t[:, 0]))
    return pts, s, phi

def fit_cell(pts, s, phi, sc):
    ss = np.linspace(sc - ELL / 2, sc + ELL / 2, 41)
    ph = np.interp(ss, s, phi)
    c2, c1, c0 = np.polyfit(ss - sc, ph, 2)
    cen = np.array([np.interp(sc, s, pts[:, 0]), np.interp(sc, s, pts[:, 1])])
    return np.array([cen[0], cen[1], c0, c1, 2 * c2])

def seed(dp=0.02, ds=0.02, lim=EXT, order='dist', rng=None, init=None):
    cands = []
    for p in np.arange(0, 4.0, dp):
        pts, s, phi = parabola_arc(p)
        for sc in np.arange(s[0] + ELL, s[-1] - ELL, ds):
            cen = np.array([np.interp(sc, s, pts[:, 0]), np.interp(sc, s, pts[:, 1])])
            if np.abs(cen).max() > W + ELL / 2 + r:
                continue
            cands.append((p, sc, cen))
    key = [np.linalg.norm(c[2]) for c in cands]
    if rng is not None:
        lane_off = {p: rng.uniform(0, 0.6) for p in np.unique([c[0] for c in cands])}
        key = np.array(key) + np.array([lane_off[c[0]] for c in cands])
    idx = np.argsort(key)
    q = [np.array([ELL / 2, 0, 0, 0, 0])] if init is None else list(init)
    arcs = {}
    Xall = list(geom(np.array(q), np.full(len(q), ELL))[0])
    tree = cKDTree(dense_spines(np.array(Xall), 3).reshape(-1, 2))
    for k in idx:
        p, sc, cen = cands[k]
        if p not in arcs:
            arcs[p] = parabola_arc(p)
        qi = fit_cell(*arcs[p], sc)
        Xi = geom(qi[None], np.array([ELL]))[0][0]
        if np.abs(Xi).max() > lim - r or np.abs(Xi).max(1).min() > W + r:
            continue
        d, _ = tree.query(dense_spines(Xi[None], 3)[0])
        if d.min() >= w * 1.003:
            q.append(qi); Xall.append(Xi)
            tree = cKDTree(dense_spines(np.array(Xall), 3).reshape(-1, 2))
    q = np.array(q)
    ell = np.full(len(q), ELL); fixed = np.zeros(len(q), bool); fixed[0] = True
    return q, ell, fixed

def brick_seed(lim=2.0, gap=1.004, spacing=1.01, rng=None):
    q = [np.array([ELL / 2, 0, 0, 0, 0])]
    Xall = [geom(q[0][None], np.array([ELL]))[0][0]]
    k = 0
    while True:
        p = 2 * w * (k + 1) * gap
        if p / 2 > lim:
            break
        pts, s, phi = parabola_arc(p, ymax=6.0)
        sv = np.interp(0.0, pts[:, 1], s)
        off = 0.0 if k % 2 == 0 else 0.5
        for j in range(-8, 9):
            sc = sv + (j + off) * L * spacing
            if sc - ELL < s[0] or sc + ELL > s[-1]:
                continue
            qi = fit_cell(pts, s, phi, sc)
            Xi = geom(qi[None], np.array([ELL]))[0][0]
            if np.abs(Xi).max() > lim - r or np.abs(Xi).max(1).min() > W + r:
                continue
            Xs = np.array(Xall)
            tree = cKDTree(dense_spines(Xs, 3).reshape(-1, 2))
            d, _ = tree.query(dense_spines(Xi[None], 3)[0])
            if d.min() >= w * 1.003:
                q.append(qi); Xall.append(Xi)
        k += 1
    return seed(lim=lim, rng=rng, init=q)

if __name__ == '__main__':
    import pickle, time
    t0 = time.time()
    q, ell, fixed = seed()
    print(len(q), time.time() - t0, valid(q, ell))
    pickle.dump((q, ell, fixed), open('seed.pkl', 'wb'))
    m = metrics(*prune_outside(q, ell, fixed)[:2]); print(m)
    plot(q, ell, 'seed.png', gap=(m['gap_at'], m['max_gap_diam_over_w'] * w / 2))
