import numpy as np, time, pickle, sys
from pack import *
tag = sys.argv[1] if len(sys.argv) > 1 else 'a'
q = np.array([[ELL/2, 0, 0, 0, 0]]); ell = np.array([ELL]); fixed = np.array([True])
tabu = []
t0 = time.time()
for it in range(200):
    X,_,_,_ = geom(q, ell)
    gaps = largest_gaps(X, n=6, excl=tabu)
    gaps = [g for g in gaps if g[1] > 0.55*r]
    if not gaps:
        break
    done = False
    for pc, c in gaps:
        ok, qn, en, fx = try_insert(q, ell, fixed, pc)
        if ok:
            q, ell, fixed = qn, en, fx; done = True
            print(f'{it:3d} N={len(q)} inserted gap r={c/r:.3f} at {pc.round(2)} t={time.time()-t0:.0f}s', flush=True)
            break
        tabu.append(pc)
        print(f'    failed at {pc.round(2)} (gap {c/r:.3f})', flush=True)
    if not done:
        break
    pickle.dump((q, ell, fixed), open(f'state_{tag}.pkl', 'wb'))
q2, e2, f2 = prune_outside(q, ell, fixed)
m = metrics(q2, e2); print(m)
pickle.dump((q2, e2, f2), open(f'state_{tag}.pkl', 'wb'))
plot(q2, e2, f'pack_{tag}.pdf', title=f"N={m['n']}  phi={m['phi']:.3f}  maxgap={m['max_gap_diam_over_w']:.3f}w  dev={m['mean_dev_deg']:.1f}deg", gap=(m['gap_at'], m['max_gap_diam_over_w']*w/2))
