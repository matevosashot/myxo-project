# Resource measurements — timing, threading, memory

All numbers measured 2026-09-18/19 on this cluster. Read this before choosing
`--mem`, `--time` or `--cpus-per-task` for a sweep.

**Hardware.** Compute nodes are AMD **EPYC 9655** (Zen 5), 192 logical cores,
2.3 TB RAM, 28 TB node-local NVMe at `/scratch`. `horus*` is Intel Cascade Lake
(96 cores, 4.1 TB). The login node is an Intel Core Ultra 9 285K and is **not**
representative — see [Measurement hazards](#measurement-hazards).

---

## 1. Where the time goes

$n = 3N_{\text{int}}^2$. Profiled at $n=3000$ on an exclusive node:

| step | time | share |
|---|---|---|
| `Eigensystem` | 7.04 s | **72%** |
| 6 × `LinearSolve` | 2.49 s | 25% |
| 12 × gemm | 0.30 s | 3% |

The solve is one irreducible $O(n^3)$ eigendecomposition plus change. This is a
healthy profile, not a defect: there is no hidden pathology to fix.

### Backend comparison, same node, both with a working BLAS

| $n$ | Mathematica eigen route | Python Schur + recursive B–S | speedup |
|---|---|---|---|
| 1000 | 1.92 s | 0.44 s | 4.3× |
| 2000 | 6.41 s | 1.71 s | 3.7× |
| 3000 | 12.10 s | 3.46 s | 3.5× |
| 4000 | 23.53 s | 5.99 s | 3.9× |

The gap widens with $n$: Mathematica's `Eigensystem` grew 7.6× from $n$ = 2000→4000
where numpy's `eig` grew 4.3×.

### Routes that look obvious and are slower

| route | $n=4000$ | verdict |
|---|---|---|
| Mathematica `LyapunovSolve` | 74.6 s | 3.2× **slower** than the eigen route, and scales worse than cubic. Not an option. |
| `scipy.linalg.solve_continuous_lyapunov` | 8.9 s @ $n$=2000 | slower than the eigen route it would replace — calls serial, unblocked LAPACK `?trsyl`. scipy 1.18 does not expose the blocked `?trsyl3`. |

This is why `lyap_solve.py` carries its own recursive Bartels–Stewart rather than
calling a library.

### Extrapolation to production sizes

Cubic in $n$ from the $n=4000$ anchor. **Treat as an estimate**: it spans 17× in $n$.

| `Nint` | $n$ | Mathematica route | Python backend |
|---|---|---|---|
| 50 | 7 500 | ~2.6 min | ~40 s |
| 100 | 30 000 | ~2.8 h | ~42 min |
| 130 | 50 700 | ~13.4 h | ~3.4 h |
| 140 | 58 800 | ~21 h | ~5.3 h |
| 150 | 67 500 | ~32 h | ~8.1 h |

---

## 2. Threading — the important one

**BLAS threads are spawned and are busy, but above ~24 they spin rather than
work.** `cores_busy` below is CPU-time ÷ wall-time, i.e. average cores actually
occupied.

### $n = 9000$, `scipy.linalg.schur` (the dominant step, 81–83% of the solve)

| threads | wall | CPU | cores busy |
|---|---|---|---|
| 8 | 35.1 s | 278 s | 7.9 |
| **24** | **32.7 s** | 777 s | 23.7 |
| 48 | 37.0 s | 1757 s | 47.5 |
| 64 | 36.0 s | 2279 s | 63.3 |
| 96 | 37.1 s | 2349 s | 63.3 |

Wall time is **best at 24 threads and degrades above it**. At 96 threads, 63 cores
burn 2349 CPU-seconds to deliver 37 s of wall — worse than 8 threads using 278.
Same shape at $n=3000$ (schur 2.92 s @ 8 threads, 2.81 s @ 64 — 4% for 8× the cores).

`_lyap_tri` (the recursive Bartels–Stewart) is flat above 24 threads: 5.71 / 5.01 /
5.04 / 4.98 / 4.99 s. Its GEMMs are leaf-sized (`BLK = 192`) and threading overhead
cancels the gain.

Pure GEMM, by contrast, scales properly — 1.63 s → 0.257 s from 8 → 64 threads, and
4025 GFLOPS at $n=3000$. The problem is specific to the eigensolve/Schur
reduction, whose QR sweep is inherently sequential.

### Hard cap: `MAX_THREADS = 64`

numpy's bundled OpenBLAS is built with `MAX_THREADS=64`. Requesting 96 or 192
reports `num_threads=64` and gives byte-identical timings. **Anything above 64 is
silently ignored.**

### Recommendation

Use **`--cpus-per-task=24`**. It is at the wall-time optimum, wastes ~3× less CPU
than 64, and schedules far faster because 24-core slots are abundant. The existing
`set2_2_comove` scripts ask for 72, which is both over-provisioned and past the
point where more cores make it slower.

---

## 3. Memory

Peak RAM is what sets `--mem`, hence the queue class, hence the wall-clock you
actually experience. One dense $n \times n$ Real64 is $8n^2$ bytes.

### Measured peak, in multiples of $8n^2$ (Mathematica kernel only)

| route | peak | at $n$=1000 / 1500 / 2000 |
|---|---|---|
| reference (eigendecomposition) | **15.5×** | 16.63 / 15.52 / 15.42 |
| `LyapunovSolve` | **7.00×** | 7.004 / 7.003 / 7.003 |

The reference route holds `eigvecs`, `Pev`, `lamSum`, `num`, `num/lamSum` and
`Transpose[Pev]` simultaneously, all **Complex128** (16 bytes/entry). The Python
backend never forms a complex matrix.

### Whole-job peak, measured

`sacct` MaxRSS over the batch step — Mathematica **and** its Python child together
— at `Nint = 70` ($n = 14700$, $8n^2 = 1.61$ GiB):

| backend | whole-job peak | ×$8n^2$ |
|---|---|---|
| reference | 29.1 GiB (Mathematica alone) | **18.1** |
| Python | 10.8 GiB | **6.70** |

So the Python path is ~2.7× smaller end to end. An earlier note in this repo
quoted ~3×; that counted only the Mathematica side and ignored the backend
process, which holds `A`, `T`, `Z` and the solution buffer itself. Size `--mem`
from **8×** $8n^2$: `sacct` samples at intervals and can undercount a short peak.

Inside the Python process alone, sampled per phase at $n=8000$ (one $8n^2$ =
0.477 GiB), peak 4.70×:

| phase | ×$8n^2$ |
|---|---|
| schur | 4.08 (the floor: `A` + `T` + `Z` + LAPACK workspace) |
| transform | 4.43 |
| back-substitute | 4.16 |
| congruence | 3.60 |
| residual | 4.64 (mostly reclaimable page cache from the memmaps) |

What got it there: streaming $Z^{\mathsf T} Q Z$ a column block at a time so `Q`
never enters RAM; recursing into views of one preallocated buffer instead of
`np.block`, which allocated a fresh $n \times n$ at every level; blockwise
symmetrisation instead of `M += M.T.copy()`; re-reading `A` and `Q` from scratch
for the residual instead of holding them; and dropping `T` once the
back-substitution is done.

**Do not turn off the refinement pass to save memory.** It costs two more
$n\times n$ buffers, and it is not optional: Bartels–Stewart is backward stable,
but the *forward* error scales with the conditioning of the Lyapunov operator,
which is poor for this drift. Measured on the real problem, one pass leaves a
residual of 5e-11 against 1e-14 refined — and on a production run at B = 2·10⁵,
refinement improved the residual by a factor of **39910**. A random stable `A`
shows 2e-14 unrefined and must not be used to justify skipping it.

### The printed estimate understates peak by ~8×

Block [6] prints `2*8*n^2` GiB, counting only `Adense` + `Qdense`. Real peak for the
reference route is ~15.5×, not 2×. Do not size `--mem` from that line.

### One dense $n \times n$ at production sizes

| `Nint` | $n$ | $8n^2$ | 3× | 9× | 15.5× |
|---|---|---|---|---|---|
| 50 | 7 500 | 0.42 GiB | 1.3 | 3.8 | 6.5 |
| 100 | 30 000 | 6.71 GiB | 20 | 60 | 104 |
| 130 | 50 700 | 19.15 GiB | 57 | 172 | 297 |
| 140 | 58 800 | 25.76 GiB | 77 | 232 | 399 |
| 150 | 67 500 | 33.95 GiB | 102 | 306 | 526 |

### Scratch

The backend stages `A`, `Q`, `C` as row-major Real64, so it needs $3 \cdot 8n^2$
bytes — 102 GiB at `Nint = 150`. Checked before the solve. Node-local `/scratch` has
**28 TB**; `$TMPDIR` points at shared GPFS, so set `scratchDir` to
`/scratch/$SLURM_JOB_ID` explicitly rather than relying on the `Automatic` default.

---

## 4. BLAS: get this wrong and nothing else matters

| environment | GEMM $n$=3000 | `eig` | `schur` | BLAS |
|---|---|---|---|---|
| `python/3.14.7` | 3865 GFLOPS | 3.55 s | 2.72 s | OpenBLAS **SkylakeX** |
| `python/3.13.9` | 4014 | 3.62 | 2.74 | OpenBLAS SkylakeX |
| `python/3.12.11` | 4077 | 3.52 | 2.74 | OpenBLAS |
| pip venv (isolated) | 4135 | 3.48 | 2.69 | OpenBLAS SkylakeX |
| `anaconda/2022.10` | 2854 | **6.64** | 5.52 | **MKL** |
| `intelpython3` | 589 | 7.59 | 6.25 | MKL 2018 |
| `anaconda3` → `/usr/bin/python3` | **6.4** | — | — | **reference netlib** |

Three things worth knowing:

* **`anaconda3` resolves to `/usr/bin/python3`, whose numpy has no optimized BLAS**
  — 6.4 GFLOPS, flat in thread count, ~600× slower than the module builds. A venv
  created with `--system-site-packages` inherits it on the compute nodes. This is
  the single largest performance effect measured in this whole exercise.
* **MKL is the wrong choice on these AMD nodes** — 1.9× slower on the eigensolve
  than OpenBLAS. Do not "upgrade" to an Intel stack here.
* OpenBLAS dispatches to its **SkylakeX** (AVX-512) kernel on Zen 5, which is
  correct. `OPENBLAS_CORETYPE=Zen` **segfaults** on this build — leave dispatch on
  auto.

The `python/*` module interpreters cannot be run by absolute path alone: they die
with `error while loading shared libraries: libpython3.14.so.1.0`, because
`PYTHONHOME`, `LD_LIBRARY_PATH` and `PYTHONPATH` come from the modulefile and
`libpython` lives in `lib64`, not `lib`. `lyapunov_solver_optimized.m` reconstructs
that environment itself, so jobscripts need no `module load`.

---

## Measurement hazards

Every one of these produced a wrong number during this work, and each wrong number
was believed for a while.

1. **The login node is shared and was at load average 98 on 24 cores.** It made
   `LinearSolve` look 43 s where it is 0.03 s, and produced an entirely fictitious
   "46× speedup from `LyapunovSolve`" that reversed on an exclusive node. **Never
   benchmark on the login node.** Use `sbatch --exclusive`.
2. **`--system-site-packages` silently swaps the BLAS** — same numpy version string,
   600× difference. Always print `threadpoolctl.threadpool_info()`; an empty list
   means no BLAS is loaded at all.
3. **`AbsoluteTiming[e]` with `e` an ordinary function argument times nothing** —
   the argument evaluates before `AbsoluteTiming` sees it. Needs `HoldFirst`/`HoldRest`.
   Symptom: a table of `0.` timings.
4. **Do not use `N` as a loop variable in Mathematica.** It shadows the built-in and
   every internal `N[...]` becomes `40[...]`; the failure surfaces somewhere else
   entirely (here: "steadyState carries no u").
5. **Exporting `PYTHONHOME`/`LD_LIBRARY_PATH` into a job shell corrupts a later
   Mathematica run in the same job** — measured `Eigensystem` 11.4 s against 7.14 s
   clean, and GEMM 359 GFLOPS against 4025. Scope those variables to the `python`
   invocation.
6. **`CForm`, not raw printing, for numbers in logs.** A `Real` printed with a
   superscript exponent is shredded by line-based log tools.
7. **The cell Péclet number does NOT scale with B.** Raising B from 10^4 to
   2·10^5 *lowered* the measured Péclet (0.0256 -> 0.0129), because the steady
   state flattens as B rises, so `(B/xi_0) grad rho_ss` self-regulates rather
   than growing linearly. Predicting it from a fixed steady state overestimates
   badly. (kappa_Q also rose from 4.147 to 6.48 with ell_d = 1.0, worth 0.64x of
   the change.)
8. **Beware printing a large expression into a log.** One `Rasterize::type` echo of a
   `Graphics3D` produced 2.9 MB at `Nint=50`; an unguarded condition-number call
   produced a 103 MB log here.

---

## Summary for sizing a job

| knob | value | why |
|---|---|---|
| `--cpus-per-task` | **24** | wall-time optimum; degrades above, and >64 is ignored |
| `--mem` | **8×** $8n^2$ | measured whole-job peak is 6.70×; 8× covers sacct's sampling |
| `scratchDir` | `/scratch/$SLURM_JOB_ID` | node-local NVMe, 28 TB; `$TMPDIR` is GPFS |
| python | leave to the backend | it finds and configures a module interpreter itself |
