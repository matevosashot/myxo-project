import numpy as np, time, pickle, sys
import pack
from pack import *
init, tag = sys.argv[1], sys.argv[2]
pack.P.update(lam_a=float(sys.argv[3]), g=float(sys.argv[4]), k_t=float(sys.argv[5]))
thr = 0.5
MAXDA = np.radians(10)

def cover(q, ell, h=0.01):
    X = geom(q, ell)[0]
    g = np.arange(-W + h / 2, W, h); gx, gy = np.meshgrid(g, g)
    d, _ = cKDTree(dense_spines(X).reshape(-1, 2)).query(np.stack([gx.ravel(), gy.ravel()], 1))
    return (d <= r).mean()

q, ell, fixed = pickle.load(open(init, 'rb'))
q = relax(q, ell, fixed, maxiter=1500)
ok, dmin = valid(q, ell); print('relaxed seed valid', ok, dmin / w, 'cover', cover(q, ell), flush=True)
tabu = []
t0 = time.time()
for it in range(300):
    X = geom(q, ell)[0]
    gaps = [g for g in largest_gaps(X, n=8, excl=tabu) if g[1] > thr * r]
    if not gaps: break
    c_old = cover(q, ell)
    done = False
    for pc, c in gaps:
        ok, qn, en, fx = try_insert(q, ell, fixed, pc, max_da=MAXDA)
        c_new = cover(qn, en) if ok else 0
        if ok and c_new > c_old + 0.005:
            q, ell, fixed = qn, en, fx; done = True
            print(f'{it:3d} N={len(q)} gap {c/r:.3f} at {pc.round(2)} cover {c_old:.3f}->{c_new:.3f} t={time.time()-t0:.0f}s', flush=True)
            break
        tabu.append(pc)
        print(f'    rejected at {pc.round(2)} (gap {c/r:.3f}) ok={ok} cover {c_new:.3f}', flush=True)
    if not done: break
    pickle.dump((q, ell, fixed), open(f'state_{tag}.pkl', 'wb'))
q2, e2, f2 = prune_outside(q, ell, fixed)
m = metrics(q2, e2); print(m, valid(q2, e2), flush=True)
pickle.dump((q2, e2, f2), open(f'state_{tag}.pkl', 'wb'))
plot(q2, e2, f'pack_{tag}.png', title=f"N={m['n']}  cover={m['phi']:.3f}  maxgap={m['max_gap_diam_over_w']:.3f}w  dev mean {m['mean_dev_deg']:.1f} max {m['max_dev_deg']:.1f} deg", gap=(m['gap_at'], m['max_gap_diam_over_w']*w/2))
