"""Local search around the largest gaps: insert / shift / bend moves, each followed by
full relaxation. Accept moves that raise coverage or shrink the largest gap."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np, pickle, sys, time
from multiprocessing import Pool
import flex
from flex import *

MAXGAP_H = 0.006


def evaluate(V, ell):
    X = geom(V, ell)[0]
    g = largest_gaps(X, n=4, h=MAXGAP_H)
    return cover(V, ell, 0.006), (2 * g[0][1] / w if g else 0.0), g


def run_move(args):
    kind, V, ell, fixed, extra = args
    try:
        if kind == 'insert':
            ok, V, ell, fixed = try_insert(V, ell, fixed, extra)
            if not ok:
                return None
        else:
            V = relax(V, ell, fixed, maxiter=600)
            if not valid(V, ell)[0]:
                return None
        c, gmax, _ = evaluate(V, ell)
        return kind, V, ell, fixed, c, gmax
    except Exception as e:  # keep the pool alive
        return None


def moves(V, ell, fixed, gaps):
    X = geom(V, ell)[0]
    out = []
    for pc, cl in gaps:
        out.append(('insert', V, ell, fixed, pc))
        # nearest non-fixed cells to the gap
        d = np.sqrt(((X - pc) ** 2).sum(-1)).min(1)
        d[fixed] = np.inf
        for i in np.argsort(d)[:3]:
            j = np.argmin(((X[i] - pc) ** 2).sum(-1))
            dirv = pc - X[i, j]
            dirv /= np.linalg.norm(dirv) + 1e-12
            # shift whole cell toward the gap
            for f in (0.5, 1.0):
                Vn = V.copy(); Vn[i, :2] += f * cl * dirv
                out.append(('shift', Vn, ell, fixed, None))
            # bend the front or the rear half, both ways
            for half in (np.arange(K) >= K // 2, np.arange(K) < K // 2):
                for ang in (-0.3, -0.15, 0.15, 0.3):
                    Vn = V.copy(); Vn[i, 2:][half] += ang
                    out.append(('bend', Vn, ell, fixed, None))
    return out


def main(init, tag, iters=40, nproc=22):
    V, ell, fixed = pickle.load(open(init, 'rb'))
    V = relax(V, ell, fixed, maxiter=1000)
    c0, g0, gaps = evaluate(V, ell)
    print(f'start N={len(V)} cover={c0:.4f} maxgap={g0:.3f}w', flush=True)
    tabu = []
    t0 = time.time()
    with Pool(nproc) as pool:
        for it in range(iters):
            X = geom(V, ell)[0]
            gaps = [g for g in largest_gaps(X, n=6, h=MAXGAP_H, excl=tabu) if g[1] > 0.3 * r]
            if not gaps:
                break
            res = [x for x in pool.map(run_move, moves(V, ell, fixed, gaps)) if x is not None]
            best = None
            for kind, Vn, en, fn, c, gm in res:
                gain = (c - c0) + 0.02 * (g0 - gm)       # coverage first, gap second
                if c < c0 - 0.002 or gm > max(g0, 0.8) + 1e-3:
                    continue
                if gain > 1e-3 and (best is None or gain > best[0]):
                    best = (gain, kind, Vn, en, fn, c, gm)
            if best is None:
                tabu += [g[0] for g in gaps[:2]]
                print(f'{it:3d} no improving move ({len(res)} valid), tabu {len(tabu)}', flush=True)
                if len(tabu) > 24:
                    break
                continue
            _, kind, V, ell, fixed, c0, g0 = best
            V, ell, fixed = prune_outside(V, ell, fixed)
            print(f'{it:3d} {kind:6s} N={len(V)} cover={c0:.4f} maxgap={g0:.3f}w t={time.time()-t0:.0f}s', flush=True)
            pickle.dump((V, ell, fixed), open(f'polish_{tag}.pkl', 'wb'))
    m = metrics(V, ell)
    print(m, flush=True)
    pickle.dump((V, ell, fixed), open(f'polish_{tag}.pkl', 'wb'))
    plot(V, ell, f'polish_{tag}.png', title=f"N={m['n']} cover={m['cover']:.3f} maxgap={m['max_gap_over_w']:.2f}w outer dev {m['dev_outer_mean']:.1f}/{m['dev_outer_max']:.1f}deg", gap=(m['gap_at'], m['max_gap_over_w'] * w / 2))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
