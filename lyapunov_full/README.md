# Fluctuations of the full model about the comoving +1/2 defect

Equal-time Gaussian fluctuations ψ = (δρ, δQ1, δQ2) about the full-model steady state of
`../comoving_defect_full/`, with the density variance split by noise source and by pathway.
The results and equations are in Sec. 7 of `../report_full_solution/report.pdf`.

- **Drift:** the linearization of the full model, with the same weights {wT, wR, wH, wE} as
  the steady solver. All weights 0 gives the paper's drift exactly.
- **Noise:** the paper's Dean noise with exact moments, split into five sources: ρρ isotropic,
  ρρ anisotropic, ρ–Q cross, QQ reversal, QQ rotational.
- **Grid, wall and cutoff:** those of `../lyapunov_solver/lyapunov_solver.m`, which is loaded
  read-only. The grid is uniform and node-centred; the wall is the reservoir wall ψ = 0 (and so
  δH = 0); the cutoff is the calibrated ℓ_W.

## Files

| file | content |
|---|---|
| `lyapunov_full_solver.m` | `lfSteadyFromRun`, `lfAssemble` (A + five N_k, sparse), `lfWriteScratch`, `lfRunBackend` |
| `lyap_multi.py` | backend: one real Schur of A; total refined; per-source block diagonals; pathway split (δρ with δQ frozen, via the sine basis and a Sylvester solve); slowest modes |
| `run_fluct.wls` | driver: `<Nint> <milestone> [runTag=A1] [saveC=0]` → `../data/fluctuations-full/<tag>/` |
| `jobs/run_fluct.sh` | SLURM wrapper (24 cpus; memory ≈ 5 × 8n², n = 3 Nint²) |
| `check_small.wls`, `check_multi.py`, `check_linearization.wls` | regression against the paper's module, the backend against dense solves, and the drift against the symbolic linearization; logs in `logs/` |

## Runs

| Nint | time | peak memory |
|---|---|---|
| 150 | 12–18 h | 205 GB |
| 180 | ~1.5 d | 424 GB |

The SLURM log of a 180 run is truncated after the backend. The jobscript therefore checks
`result.json` rather than the log.

## Pitfalls

- **Sparse-array backgrounds.** Mathematica `SparseArray` blocks with different backgrounds
  (0 vs 0.) make `ArrayFlatten` return a dense array. Here that made the Q–Q noise block nearly
  dense and crashed the kernel, which Mathematica reported as a "license error". Every block is
  built with background 0.
- **Schur vectors are a view.** `scipy.linalg.schur` returns Fortran-ordered Z, so
  `np.ascontiguousarray(Z[:, a:b].T)` is a view; copy it before modifying it.
