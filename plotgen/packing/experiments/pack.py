"""Pack bendable capsules around a +1/2 nematic defect.

Units: W = 1, box [-1,1]^2. Cell: capsule of total length L = 0.9, width w = L/5.
Spine length ell = L - w, spine radius r = w/2.

Each cell is a spine of fixed length whose tangent angle varies quadratically along
arclength: phi(s) = theta + k0 s + k1 s^2/2 (linear curvature). The spine is sampled
as a polyline with M nodes and exact segment length ell/(M-1).

Energy = overlap penalty (node disks, slightly inflated radius)
       + wall penalty (extended box)
       + alignment  lam_a * int sin^2(phi - psi) ds,  psi = polar/2
       + bending    lam_b * int kappa^2 ds  (+ curvature cap)

Cells are added one by one at the largest empty circle inside the box and grown from a
disk to full length while everything relaxes (cells push and tilt together).
"""
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree

W = 1.0
L = 0.9 * W
w = L / 5
r = w / 2
ELL = L - w
M = 21                      # nodes per spine
EXT = 2.0                   # confining wall (cells may stick out of the true box)
R_EFF = r * 1.012           # inflated radius for the node-disk overlap term
KAPPA_MAX = 1.0 / w         # min radius of curvature = w

U = (np.arange(M - 1) + 0.5) / (M - 1) - 0.5          # segment midpoints (fraction)
MID = (M - 1) // 2
A = (np.arange(M - 1)[None, :] < np.arange(M)[:, None]).astype(float) \
    - (np.arange(M - 1) < MID).astype(float)[None, :]   # (M, M-1)

P = dict(k_o=1e5, k_w=1e5, lam_a=1.0, lam_b=0.002, k_c=10.0, g=0.0, k_t=0.0, tilt_max=np.radians(10))


def director_angle(x, y):
    return 0.5 * np.arctan2(y, x)


def geom(q, ell):
    """q: (N,5) [cx,cy,theta,k0,k1]; ell: (N,). Returns nodes X (N,M,2), phi, s, h."""
    s = ell[:, None] * U[None, :]
    phi = q[:, 2:3] + q[:, 3:4] * s + 0.5 * q[:, 4:5] * s ** 2
    h = ell / (M - 1)
    seg = h[:, None, None] * np.stack([np.cos(phi), np.sin(phi)], -1)
    X = q[:, None, :2] + np.einsum('jk,nkd->njd', A, seg)
    return X, phi, s, h


def energy(q, ell, fixed, p=P, want_grad=True):
    N = len(q)
    X, phi, s, h = geom(q, ell)
    GX = np.zeros_like(X)
    E = 0.0
    # --- overlap between node disks of different cells
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
    # --- wall
    ex = np.abs(X) - (EXT - r)
    m = ex > 0
    E += p['k_w'] * (ex[m] ** 2).sum()
    GX += np.where(m, 2 * p['k_w'] * ex * np.sign(X), 0.0)
    # --- alignment
    mids = 0.5 * (X[:, :-1] + X[:, 1:])
    psi = director_angle(mids[..., 0], mids[..., 1])
    dl = phi - psi
    hh = h[:, None]
    E += p['lam_a'] * (hh * np.sin(dl) ** 2).sum()
    dphi = p['lam_a'] * hh * np.sin(2 * dl)
    if p['k_t'] > 0:                      # hard-ish cap on tilt away from the director
        sd = np.sin(dl)
        exc = np.abs(sd) - np.sin(p['tilt_max'])
        mt = exc > 0
        E += p['k_t'] * (hh * np.where(mt, exc, 0) ** 2).sum()
        dphi = dphi + np.where(mt, 2 * p['k_t'] * hh * exc * np.sign(sd) * np.cos(dl), 0)
    rho2 = (mids ** 2).sum(-1) + 1e-4
    dpsi = -dphi
    gm = dpsi[..., None] * np.stack([-mids[..., 1], mids[..., 0]], -1) / (2 * rho2[..., None])
    GX[:, :-1] += 0.5 * gm
    GX[:, 1:] += 0.5 * gm
    # --- weak pull toward the core (compaction)
    if p['g'] > 0:
        rn = np.sqrt((X ** 2).sum(-1)) + 1e-9
        E += p['g'] * (h[:, None] * rn).sum()
        GX += p['g'] * h[:, None, None] * X / rn[..., None]
    # --- bending
    k0, k1 = q[:, 3], q[:, 4]
    E += p['lam_b'] * (ell * k0 ** 2 + k1 ** 2 * ell ** 3 / 12).sum()
    gq = np.zeros_like(q)
    gq[:, 3] += p['lam_b'] * 2 * ell * k0
    gq[:, 4] += p['lam_b'] * 2 * k1 * ell ** 3 / 12
    for sgn in (1, -1):
        ke = k0 + sgn * k1 * ell / 2
        exc = np.abs(ke) - KAPPA_MAX
        mm = exc > 0
        E += p['k_c'] * (exc[mm] ** 2).sum()
        gk = np.where(mm, 2 * p['k_c'] * exc * np.sign(ke), 0.0)
        gq[:, 3] += gk
        gq[:, 4] += gk * sgn * ell / 2
    if not want_grad:
        return E
    # --- chain rule
    perp = np.stack([-np.sin(phi), np.cos(phi)], -1)
    ATG = np.einsum('jk,njd->nkd', A, GX)
    dphi = dphi + hh * (ATG * perp).sum(-1)
    gq[:, :2] += GX.sum(1)
    gq[:, 2] += dphi.sum(1)
    gq[:, 3] += (s * dphi).sum(1)
    gq[:, 4] += (0.5 * s ** 2 * dphi).sum(1)
    gq[fixed] = 0.0
    return E, gq


def relax(q, ell, fixed, maxiter=400, p=P):
    free = ~fixed

    def f(x):
        qq = q.copy()
        qq[free] = x.reshape(-1, 5)
        E, g = energy(qq, ell, fixed, p)
        return E, g[free].ravel()

    res = minimize(f, q[free].ravel(), jac=True, method='L-BFGS-B',
                   options=dict(maxiter=maxiter, maxcor=30, gtol=1e-9, ftol=1e-13))
    q = q.copy()
    q[free] = res.x.reshape(-1, 5)
    return q


# ---------------------------------------------------------------- exact geometry
def seg_dist(p1, p2, q1, q2):
    """Distance between segments p1p2 and q1q2 (arrays (...,2)), 2D."""
    def pt_seg(a, b, c):
        ab = c - b
        t = np.clip(((a - b) * ab).sum(-1) / np.maximum((ab ** 2).sum(-1), 1e-18), 0, 1)
        return np.sqrt(((b + t[..., None] * ab - a) ** 2).sum(-1))

    def cross(u, v):
        return u[..., 0] * v[..., 1] - u[..., 1] * v[..., 0]
    d = np.minimum(np.minimum(pt_seg(p1, q1, q2), pt_seg(p2, q1, q2)),
                   np.minimum(pt_seg(q1, p1, p2), pt_seg(q2, p1, p2)))
    r1, r2 = p2 - p1, q2 - q1
    o1 = cross(r1, q1 - p1) * cross(r1, q2 - p1)
    o2 = cross(r2, p1 - q1) * cross(r2, p2 - q1)
    return np.where((o1 < 0) & (o2 < 0), 0.0, d)


def min_spine_distances(X):
    """Exact min spine-spine distance for all cell pairs that are close. Returns (i,j,d) list."""
    N = len(X)
    cen = X.mean(1)
    out = []
    for i in range(N):
        for j in range(i + 1, N):
            if np.linalg.norm(cen[i] - cen[j]) > ELL + 2 * w:
                continue
            a1, a2 = X[i, :-1][:, None], X[i, 1:][:, None]
            b1, b2 = X[j, :-1][None], X[j, 1:][None]
            out.append((i, j, seg_dist(a1, a2, b1, b2).min()))
    return out


def dense_spines(X, refine=6):
    t = np.linspace(0, 1, refine, endpoint=False)
    pts = X[:, :-1, None, :] * (1 - t)[None, None, :, None] + X[:, 1:, None, :] * t[None, None, :, None]
    pts = pts.reshape(len(X), -1, 2)
    return np.concatenate([pts, X[:, -1:]], 1)


def clearance(pts, X):
    """Radius of largest circle centred at pts that stays in the box and misses all cells."""
    D = dense_spines(X).reshape(-1, 2)
    d, _ = cKDTree(D).query(pts)
    wall = np.minimum(W - np.abs(pts[:, 0]), W - np.abs(pts[:, 1]))
    return np.minimum(d - r, wall)


def largest_gaps(X, n=10, h=0.01, excl=None):
    g = np.arange(-W + h / 2, W, h)
    gx, gy = np.meshgrid(g, g)
    pts = np.stack([gx.ravel(), gy.ravel()], 1)
    c = clearance(pts, X)
    order = np.argsort(-c)
    chosen = []
    for k in order:
        if c[k] <= 0:
            break
        if all(np.linalg.norm(pts[k] - pc) > 2 * w for pc, _ in chosen):
            if excl is None or all(np.linalg.norm(pts[k] - e) > 0.02 for e in excl):
                chosen.append((pts[k], c[k]))
        if len(chosen) >= n:
            break
    # local refinement
    out = []
    for pc, cc in chosen:
        best = (pc, cc)
        step = h
        for _ in range(6):
            gg = np.linspace(-step, step, 9)
            ox, oy = np.meshgrid(gg, gg)
            cand = best[0] + np.stack([ox.ravel(), oy.ravel()], 1)
            cv = clearance(cand, X)
            k = np.argmax(cv)
            if cv[k] > best[1]:
                best = (cand[k], cv[k])
            step /= 3
        out.append(best)
    out.sort(key=lambda t: -t[1])
    return out


def max_gap(X, h=0.004):
    g = largest_gaps(X, n=5, h=h)
    return g[0] if g else (None, 0.0)


# ---------------------------------------------------------------- insertion
def free_length(pc, ang, X, maxlen=ELL):
    """Half-length (each way) a straight spine through pc at angle ang can have before hitting."""
    D = cKDTree(dense_spines(X).reshape(-1, 2))
    t = np.linspace(0, maxlen / 2, 30)
    res = []
    for sgn in (1, -1):
        pts = pc + sgn * t[:, None] * np.array([np.cos(ang), np.sin(ang)])
        d, _ = D.query(pts)
        bad = np.nonzero(d < 2 * r)[0]
        res.append(t[bad[0] - 1] if len(bad) and bad[0] > 0 else (0 if len(bad) else maxlen / 2))
    return res


def try_insert(q, ell, fixed, pc, grow_steps=12, verbose=False, max_da=np.pi / 2):
    X, _, _, _ = geom(q, ell)
    psi = director_angle(*pc)
    best = None
    for da in np.linspace(-max_da, max_da, 19):
        a = psi + da
        fl = free_length(pc, a, X)
        score = min(fl) + 0.5 * max(fl) - 0.15 * np.sin(da) ** 2
        if best is None or score > best[0]:
            best = (score, a, fl)
    _, a, fl = best
    shift = 0.5 * (fl[0] - fl[1])
    c0 = pc + shift * np.array([np.cos(a), np.sin(a)])
    qn = np.vstack([q, [c0[0], c0[1], a, 0.0, 0.0]])
    fx = np.append(fixed, False)
    ell0 = max(0.0, min(ELL, fl[0] + fl[1]))
    lengths = np.linspace(ell0, ELL, grow_steps + 1)[1:] if ell0 < ELL else [ELL]
    for Ln in lengths:
        en = np.append(ell, Ln)
        qn = relax(qn, en, fx, maxiter=150)
    en = np.append(ell, ELL)
    qn = relax(qn, en, fx, maxiter=800)
    ok, dmin = valid(qn, en)
    if verbose:
        print(f'   insert at {pc.round(3)} ell0={ell0:.2f} ok={ok} dmin/w={dmin / w:.4f}')
    return ok, qn, en, fx


def valid(q, ell):
    X, _, _, _ = geom(q, ell)
    ds = [d for _, _, d in min_spine_distances(X)]
    dmin = min(ds) if ds else np.inf
    return dmin >= w * (1 - 1e-9), dmin


# ---------------------------------------------------------------- metrics / plot
def in_box_mask(pts):
    return (np.abs(pts[..., 0]) <= W) & (np.abs(pts[..., 1]) <= W)


def prune_outside(q, ell, fixed):
    X, _, _, _ = geom(q, ell)
    keep = fixed.copy()
    for i in range(len(q)):
        # keep if any part of the capsule is in the box
        dx = np.maximum(np.abs(X[i]) - W, 0)
        if (np.sqrt((dx ** 2).sum(-1)) < r).any():
            keep[i] = True
    return q[keep], ell[keep], fixed[keep]


def metrics(q, ell, h=0.004):
    X, phi, _, _ = geom(q, ell)
    g = np.arange(-W + h / 2, W, h)
    gx, gy = np.meshgrid(g, g)
    pts = np.stack([gx.ravel(), gy.ravel()], 1)
    D = dense_spines(X)
    d, _ = cKDTree(D.reshape(-1, 2)).query(pts)
    phi_frac = (d <= r).mean()
    mids = 0.5 * (X[:, :-1] + X[:, 1:])
    inb = in_box_mask(mids)
    dl = phi - director_angle(mids[..., 0], mids[..., 1])
    dev = np.abs((dl + np.pi / 2) % np.pi - np.pi / 2)
    pc, cmax = max_gap(X)
    return dict(n=len(q), phi=phi_frac, max_dev_deg=np.degrees(dev[inb].max()), mean_dev_deg=np.degrees(dev[inb].mean()),
                rms_dev_deg=np.degrees(np.sqrt((dev[inb] ** 2).mean())),
                max_gap_diam_over_w=2 * max(cmax, 0) / w, gap_at=pc)


def capsule_polygon(Xi, nc=16):
    t = Xi[1:] - Xi[:-1]
    t = t / np.linalg.norm(t, axis=1, keepdims=True)
    tn = np.vstack([t[:1], 0.5 * (t[:-1] + t[1:]), t[-1:]])
    tn /= np.linalg.norm(tn, axis=1, keepdims=True)
    nrm = np.stack([-tn[:, 1], tn[:, 0]], 1)
    left = Xi + r * nrm
    right = Xi - r * nrm
    a_end = np.arctan2(tn[-1, 1], tn[-1, 0])
    a_st = np.arctan2(tn[0, 1], tn[0, 0])
    ce = np.linspace(a_end + np.pi / 2, a_end - np.pi / 2, nc)
    cs = np.linspace(a_st - np.pi / 2, a_st - 3 * np.pi / 2, nc)
    cap_e = Xi[-1] + r * np.stack([np.cos(ce), np.sin(ce)], 1)
    cap_s = Xi[0] + r * np.stack([np.cos(cs), np.sin(cs)], 1)
    return np.vstack([left, cap_e, right[::-1], cap_s])


def streamline_curves(ps=None):
    if ps is None:
        ps = np.concatenate([[0], np.arange(0.36, 6, 0.36)])
    out = []
    for p_ in ps:
        y = np.linspace(-3, 3, 2001)
        x = (y ** 2 - p_ ** 2) / (2 * p_) if p_ > 0 else None
        if p_ == 0:
            out.append(np.array([[0, 0], [3, 0]]))
        else:
            out.append(np.stack([x, y], 1))
    return out


def plot(q, ell, fname, title=None, gap=None, X=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon, Rectangle, Circle
    if X is None:
        X, _, _, _ = geom(q, ell)
    fig, ax = plt.subplots(figsize=(6, 6))
    clip = Rectangle((-W, -W), 2 * W, 2 * W, transform=ax.transData)
    for c in streamline_curves():
        ln, = ax.plot(c[:, 0], c[:, 1], color='0.5', lw=1, zorder=1)
        ln.set_clip_path(clip)
    for i in range(len(X)):
        pg = Polygon(capsule_polygon(X[i]), closed=True, fc='#bcd8e8', ec='#1b3a5c', lw=1, zorder=2)
        ax.add_patch(pg)
        pg.set_clip_path(clip)
    if gap is not None and gap[0] is not None:
        ax.add_patch(Circle(gap[0], max(gap[1], 0), fc='none', ec='red', lw=1.2, zorder=3))
    ax.add_patch(Rectangle((-W, -W), 2 * W, 2 * W, fc='none', ec='k', lw=2, zorder=4))
    ax.set_xlim(-W * 1.02, W * 1.02)
    ax.set_ylim(-W * 1.02, W * 1.02)
    ax.set_aspect('equal')
    ax.axis('off')
    if title:
        ax.set_title(title, fontsize=9)
    fig.savefig(fname, bbox_inches='tight')
    plt.close(fig)
