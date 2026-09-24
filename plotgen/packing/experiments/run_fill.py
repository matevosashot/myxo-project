import numpy as np, time, pickle, sys
import pack
from pack import *
init, tag, lam = sys.argv[1], sys.argv[2], float(sys.argv[3])
thr = float(sys.argv[4]) if len(sys.argv) > 4 else 0.55
pack.P['g'] = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
pack.P['lam_a'] = lam
q, ell, fixed = pickle.load(open(init, 'rb'))
q = relax(q, ell, fixed, maxiter=1000)
tabu = []
t0 = time.time()
for it in range(200):
    X,_,_,_ = geom(q, ell)
    gaps = [g for g in largest_gaps(X, n=6, excl=tabu) if g[1] > thr*r]
    if not gaps: break
    done = False
    for pc, c in gaps:
        ok, qn, en, fx = try_insert(q, ell, fixed, pc)
        if ok:
            q, ell, fixed = qn, en, fx; done = True
            print(f'{it:3d} N={len(q)} inserted gap r={c/r:.3f} at {pc.round(2)} t={time.time()-t0:.0f}s', flush=True)
            break
        tabu.append(pc)
        print(f'    failed at {pc.round(2)} (gap {c/r:.3f})', flush=True)
    if not done: break
    pickle.dump((q, ell, fixed), open(f'state_{tag}.pkl', 'wb'))
q2, e2, f2 = prune_outside(q, ell, fixed)
m = metrics(q2, e2); print(m, valid(q2, e2))
pickle.dump((q2, e2, f2), open(f'state_{tag}.pkl', 'wb'))
plot(q2, e2, f'pack_{tag}.png', title=f"N={m['n']}  phi={m['phi']:.3f}  maxgap={m['max_gap_diam_over_w']:.3f}w  dev={m['mean_dev_deg']:.1f}deg", gap=(m['gap_at'], m['max_gap_diam_over_w']*w/2))
