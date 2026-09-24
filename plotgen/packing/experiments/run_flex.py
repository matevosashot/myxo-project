import numpy as np, time, pickle, sys
import flex
from flex import *
init, tag = sys.argv[1], sys.argv[2]
flex.P['lam_mid'] = float(sys.argv[3]) if len(sys.argv) > 3 else flex.P['lam_mid']
thr = 0.4
q, ell, fixed = pickle.load(open(init, 'rb'))
V = q if q.shape[1] == 2 + flex.K else from_q(q, ell)
V = relax(V, ell, fixed, maxiter=2000)
print('relaxed', valid(V, ell), 'cover', cover(V, ell), flush=True)
tabu = []; t0 = time.time()
for it in range(300):
    X = geom(V, ell)[0]
    gaps = [g for g in largest_gaps(X, n=8, excl=tabu) if g[1] > thr * r]
    if not gaps: break
    c_old = cover(V, ell); done = False
    for pc, c in gaps:
        ok, Vn, en, fx = try_insert(V, ell, fixed, pc)
        c_new = cover(Vn, en) if ok else 0
        if ok and c_new > c_old + 0.005:
            V, ell, fixed = Vn, en, fx; done = True
            print(f'{it:3d} N={len(V)} gap {c/r:.3f} at {pc.round(2)} cover {c_old:.3f}->{c_new:.3f} t={time.time()-t0:.0f}s', flush=True)
            break
        tabu.append(pc)
        print(f'    rejected at {pc.round(2)} (gap {c/r:.3f}) ok={ok} cover {c_new:.3f}', flush=True)
    if not done: break
    pickle.dump((V, ell, fixed), open(f'flex_{tag}.pkl', 'wb'))
V, ell, fixed = prune_outside(V, ell, fixed)
m = metrics(V, ell); print(m, flush=True)
pickle.dump((V, ell, fixed), open(f'flex_{tag}.pkl', 'wb'))
plot(V, ell, f'flex_{tag}.png', title=f"N={m['n']} cover={m['cover']:.3f} maxgap={m['max_gap_over_w']:.2f}w  outer dev {m['dev_outer_mean']:.1f}/{m['dev_outer_max']:.1f}deg  broken={m['n_broken']}", gap=(m['gap_at'], m['max_gap_over_w'] * w / 2))
