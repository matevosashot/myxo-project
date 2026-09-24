import pickle, sys, numpy as np
import flex
from flex import *
V0, ell, fixed = pickle.load(open('polish_p0.pkl', 'rb'))
for lb in [0.02, 0.05, 0.1]:
    flex.P['lam_b'] = lb
    V = relax(V0, ell, fixed, maxiter=3000)
    m = metrics(V, ell)
    Xp = geom(V, ell)[1]
    wig = np.degrees(np.abs(np.diff(Xp, axis=1)).max())
    print(lb, f"cover={m['cover']:.4f} gap={m['max_gap_over_w']:.3f} valid={m['valid']} outer={m['dev_outer_mean']:.2f}/{m['dev_outer_max']:.2f} mid={m['dev_mid_max']:.1f} maxturn/seg={wig:.1f}deg", flush=True)
    pickle.dump((V, ell, fixed), open(f'tidy_{lb}.pkl', 'wb'))
    plot(V, ell, f'tidy_{lb}.png', title=f"lam_b={lb} cover={m['cover']:.3f} gap={m['max_gap_over_w']:.2f}w", gap=(m['gap_at'], m['max_gap_over_w'] * w / 2))
