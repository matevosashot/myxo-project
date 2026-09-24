"""Seed flexible cells: straight streamline cells plus 'broken' cells whose front half
sits on a neighbouring streamline (lane change in the middle)."""
import numpy as np
from scipy.spatial import cKDTree
from pack import W, L, w, r, ELL, M, dense_spines
from lanes import parabola_arc
import flex

_arcs = {}


def arc(p):
    key = round(p, 4)
    if key not in _arcs:
        _arcs[key] = parabola_arc(max(key, 0.0), ymax=6.0, n=6001)
    return _arcs[key]


def pts_on(p, s0, s1, n):
    pts, s, _ = arc(p)
    ss = np.linspace(s0, s1, n)
    return np.stack([np.interp(ss, s, pts[:, 0]), np.interp(ss, s, pts[:, 1])], 1)


def resample(poly, ell=ELL, m=M):
    seg = np.linalg.norm(np.diff(poly, axis=0), axis=1)
    s = np.concatenate([[0], np.cumsum(seg)])
    # centre the cell on the polyline midpoint (by arclength)
    mid = s[-1] / 2
    t = np.linspace(mid - ell / 2, mid + ell / 2, m)
    X = np.stack([np.interp(t, s, poly[:, 0]), np.interp(t, s, poly[:, 1])], 1)
    # re-project to exact segment length ell/(m-1)
    phi = np.arctan2(*(np.diff(X, axis=0)[:, ::-1].T))
    phi = np.unwrap(phi)
    return np.concatenate([X[(m - 1) // 2], phi])


def cell(p1, sc, dn):
    """Vector V for a cell centred at arclength sc on lane p1; front half shifted by dn."""
    pts1, s1, _ = arc(p1)
    m = np.array([np.interp(sc, s1, pts1[:, 0]), np.interp(sc, s1, pts1[:, 1])])
    if dn == 0:
        return resample(pts_on(p1, sc - ELL / 2, sc + ELL / 2, 41))
    th = np.arctan2(m[1], m[0]) % (2 * np.pi)
    p2 = p1 + 2 * np.sin(th / 2) * dn
    if p2 < 0:
        return None
    pts2, s2, _ = arc(p2)
    j = np.argmin(((pts2 - m) ** 2).sum(1))
    h = flex.MID_HALF * ELL
    rear = pts_on(p1, sc - ELL / 2 - 0.02, sc - h, 20)
    front = pts_on(p2, s2[j] + h, s2[j] + ELL / 2 + 0.05, 20)
    return resample(np.vstack([rear, front]))


def seed(init_V, rng, lim=2.0, dp=0.02, ds=0.02, dns=(0, 0.25, -0.25, 0.5, -0.5, 0.75, -0.75),
         pen=0.04):
    cands = []
    for p in np.arange(0, 3.2, dp):
        pts, s, _ = arc(p)
        for sc in np.arange(s[0] + ELL, s[-1] - ELL, ds):
            cen = np.array([np.interp(sc, s, pts[:, 0]), np.interp(sc, s, pts[:, 1])])
            if np.abs(cen).max() > W + ELL / 2 + r:
                continue
            for dn in dns:
                cands.append((p, sc, dn * w, np.linalg.norm(cen)))
    lane_off = {p: rng.uniform(0, 0.6) for p in np.unique([c[0] for c in cands])}
    key = np.array([c[3] + lane_off[c[0]] + pen * abs(c[2]) / w for c in cands])
    V = list(init_V)
    Xall = list(flex.geom(np.array(V), np.full(len(V), ELL))[0])
    tree = cKDTree(dense_spines(np.array(Xall), 3).reshape(-1, 2))
    for k in np.argsort(key):
        p, sc, dn, _ = cands[k]
        v = cell(p, sc, dn)
        if v is None:
            continue
        Xi = flex.geom(v[None], np.array([ELL]))[0][0]
        if np.abs(Xi).max() > lim - r or np.abs(Xi).max(1).min() > W + r:
            continue
        d, _ = tree.query(dense_spines(Xi[None], 3)[0])
        if d.min() >= w * 1.003:
            V.append(v); Xall.append(Xi)
            tree = cKDTree(dense_spines(np.array(Xall), 3).reshape(-1, 2))
    V = np.array(V)
    fixed = np.zeros(len(V), bool); fixed[0] = True
    return V, np.full(len(V), ELL), fixed


def brick_part(lim=2.0, gap=1.004, spacing=1.01):
    """Only the left-side brick lanes (p >= 2w), as V vectors."""
    from lanes import brick_seed
    import lanes
    orig = lanes.seed
    lanes.seed = lambda lim, rng, init: (np.array(init), None, None)
    try:
        q, _, _ = brick_seed(lim=lim, gap=gap, spacing=spacing)
    finally:
        lanes.seed = orig
    return flex.from_q(q, np.full(len(q), ELL))
