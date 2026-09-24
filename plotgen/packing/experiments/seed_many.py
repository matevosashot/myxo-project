import numpy as np, pickle, sys
from lanes import *
for sd in range(6):
    q, ell, fixed = seed(lim=2.0, rng=np.random.default_rng(sd))
    q, ell, fixed = prune_outside(q, ell, fixed)
    m = metrics(q, ell)
    print(sd, len(q), round(m['phi'],3), round(m['max_gap_diam_over_w'],3), flush=True)
    pickle.dump((q, ell, fixed), open(f'seed{sd}.pkl', 'wb'))
    plot(q, ell, f'seed{sd}.png', title=f'seed {sd}')
