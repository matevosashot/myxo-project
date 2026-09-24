import numpy as np, pickle, sys, time
from lanes2 import *
import flex
Vb = brick_part()
for sd in range(4):
    t0 = time.time()
    V, ell, fixed = seed(Vb, np.random.default_rng(sd))
    V, ell, fixed = flex.prune_outside(V, ell, fixed)
    m = flex.metrics(V, ell)
    print(sd, f"N={m['n']} cover={m['cover']:.3f} gap={m['max_gap_over_w']:.3f} broken={m['n_broken']} valid={m['valid']} t={time.time()-t0:.0f}s", flush=True)
    pickle.dump((V, ell, fixed), open(f'fseed{sd}.pkl', 'wb'))
    flex.plot(V, ell, f'fseed{sd}.png', title=f"seed {sd} N={m['n']} cover={m['cover']:.3f} gap={m['max_gap_over_w']:.2f}w broken={m['n_broken']}", gap=(m['gap_at'], m['max_gap_over_w']*w/2))
