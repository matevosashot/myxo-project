# Timing analysis of `fluctuations_dirichlet.wls`

Parses the batch log `../out_many.txt` into CSVs and plots how the run time
scales with the grid resolution `Nint`.

## Run it

```bash
PY=/data/biophys/ashmat/miniconda/envs/ntn/bin/python   # env with matplotlib + pandas
$PY parse_log.py          # ../out_many.txt -> runs.csv, timings_long.csv, timings_wide.csv
$PY plot_timings.py       # CSVs -> fig1..fig3 (.png and .pdf)
```

`parse_log.py` takes an optional log path as argv[1] if you want to point it at
a different batch output.

## What was in the log

17 attempted runs: `fbox=10, Nint=21` (the leading, headerless block) and
`fbox=20` for `Nint` = 20, 21, 30, 31, 40, 41, 60, 61, 80, 81, 100, 101, 120,
121, 140, 141.

* **13 completed**, each with a full 16-row timing table.
* **4 failed** — `Nint` = 120, 121, 140, 141, all `Quit (core dumped)`, exit 131.
  `Nint=120` got as far as `[0] densify A,Q  20.3 s  (27.8 GiB)` and died there.
  That is exactly `2 * 8 * n^2` bytes for `n = 43200`, i.e. dense `A` and `Q` as
  real doubles — the runs hit a memory wall, not a time limit.

## Files

| file | contents |
|---|---|
| `runs.csv` | one row per attempted run: `fbox`, `Nint`, `n`, `h`, `status`, `exit_code`, `total_tracked_s`, `lyapunov_solve_s`, `max_re_lambda`, `cov_size_MiB`, `cov_save_s`, `cov_path`, `timestamp` |
| `timings_long.csv` | tidy: one row per (run, stage) — `run_id`, `fbox`, `Nint`, `n`, `stage_group`, `stage`, `seconds`, `pct_tracked` |
| `timings_wide.csv` | one row per run, one column per stage (seconds), plus totals and covariance size |
| `fig1_stacked_absolute.png` | wall time vs `Nint`, stacked by stage; second panel zooms on `Nint <= 61` |
| `fig2_stacked_fraction.png` | the same as a 100% stack — where the time goes |
| `fig3_scaling_loglog.png` | log-log scaling in hours, with each stage's power law extrapolated past the last completed run |

`cov_size_MiB` / `cov_save_s` come from the `Saved: ... covariance_*.bin
(N MiB, T s)` lines. `T` duplicates the `[7] save covariance (binary)` row of
the timing table (they agree to the printed precision); the file size is the
part the table does not carry.

## Reading the plots

Two families of runs are parsed into the CSVs but excluded from the plots:

* the single `fbox=10` run — different box size (`h = 0.909` rather than
  `1.905`), so it is not on the same scaling curve as the rest;
* **all odd `Nint`** (21, 31, 41, 61, 81, 101) — the batch ran each resolution
  twice, with the origin on the grid (odd) and off it (even). Only the even
  family is plotted. The pairs agree to a few percent everywhere, so nothing is
  lost; set `DROP_ODD_NINT = False` in `plot_timings.py` to put them back.

The stacks draw one band per stage whose share of tracked time reaches 2.5% in
at least one run (7 stages); everything else is summed into `other`. The two
standalone residual evaluations (`#1 unrefined`, `#3 verification`) are merged
into one band — same operation, done twice. `[7] save covariance` falls into
`other` in the stacks but gets its own curve in fig3, since it is the only I/O
term.

## What the numbers say

| Nint | n | total [s] | eigensolve % | FEM + plot % |
|---|---|---|---|---|
| 20 | 1 200 | 44.8 | 5.4 | 89.8 |
| 30 | 2 700 | 54.2 | 9.5 | 82.0 |
| 40 | 4 800 | 74.8 | 23.7 | 59.8 |
| 60 | 10 800 | 284.8 | 60.7 | 15.9 |
| 80 | 19 200 | 1 134.8 | 72.0 | 4.5 |
| 100 | 30 000 | 4 266.3 | 79.3 | 1.0 |

Two regimes, crossing over around `Nint ≈ 40`:

* **Below it**, the run is fixed overhead. `[1] steady state` (13 s, genuinely
  N-independent) and `[7] plot + export` (~30 s, dominated by front-end startup
  and PNG rasterisation) are 90% of a `Nint=20` run.
* **Above it**, everything is the dense Lyapunov solve. Fitted over `N >= 60`:

  | stage | fit | expected |
  |---|---|---|
  | `[6] Eigensystem` | `N^5.8` | `n^3 = N^6` |
  | `[6] back-substitute` | `N^5.5` | `n^3` |
  | `[6] refinement step` | `N^5.4` | `n^3` |
  | `[6] residual evals` | `N^4.5` | `n^2..n^3` |
  | `[7] save covariance` | `N^4.0` | `n^2 = N^4` |
  | `[3] kappa(S) inversion` | `N^2.0` | `N^2` (one solve per grid point) |
  | **TOTAL tracked** | **`N^5.3`** | approaching `N^6` |

  The total is still below `N^6` only because the O(1) overhead has not fully
  washed out yet by `N=100`.

Extrapolating the total fit (`t = 1.14e-7 * N^5.28` s) and the `8 n^2`-byte
footprints:

| Nint | n | projected wall time | covariance file | dense A+Q |
|---|---|---|---|---|
| 120 | 43 200 | ~2.9 h | 13.9 GiB | 27.8 GiB (this is where it died) |
| 140 | 58 800 | ~6.6 h | 25.8 GiB | 51.5 GiB |
| 200 | 120 000 | ~44 h | 107 GiB | 215 GiB |

Memory, not time, is the binding constraint: `Nint=140` needs ~50 GiB just to
hold `A` and `Q` densely, before `Eigensystem` allocates a complex eigenvector
matrix of the same order again.

## Caveats

* The exponents are fitted on `N >= 60`, which after dropping the odd runs is
  only three points (60, 80, 100). They are the asymptotic part of the range,
  but treat the second decimal as noise. Lower `FIT_FROM_N` in
  `plot_timings.py` to widen the fit at the cost of contaminating it with the
  small-N overhead.
* `Nint=21, fbox=20` (not plotted) is an outlier in `[6] refinement step`
  (18 s against ~2 s at neighbouring N). That run is the badly-conditioned one:
  `max Re(λ) = 0.554` against `~1e-2` elsewhere, so the refinement had far more
  work to do.
* `[7] plot + export` fits as `N^-0.3` — that is a flat, noisy ~30 s of
  front-end overhead, not a real negative exponent.
* "TOTAL tracked" is the log's own total, i.e. the sum of the instrumented
  stages, not the wall time of the whole process.
