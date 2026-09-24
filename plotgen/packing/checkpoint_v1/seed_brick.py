import numpy as np, pickle
from lanes import *
for sd in range(4):
    q, ell, fixed = brick_seed(rng=np.random.default_rng(sd))
    q, ell, fixed = prune_outside(q, ell, fixed)
    m = metrics(q, ell)
    print(sd, len(q), round(m['phi'],3), round(m['max_gap_diam_over_w'],3), valid(q,ell)[0], flush=True)
    pickle.dump((q, ell, fixed), open(f'brick{sd}.pkl', 'wb'))
    plot(q, ell, f'brick{sd}.png', title=f'brick seed {sd}')
