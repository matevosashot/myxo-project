# `lyapunov_graded_solver`

A graded-grid version of `lyapunov_solver/lyapunov_solver.m`. Same physics, same
scheme, same exact von Mises closure — the difference is that every operator
carries the metric of a graded tensor grid, and the UV cutoff is imposed by
a construction that survives a varying cell size.

**Scope.** Dirichlet (reservoir wall) only. The steady state comes from
`comoving_defect/comoving_steady_solver.m` unchanged.

---

## Why

Two constraints fight on a uniform mesh.

**Stability.** `maxReLambda` collapses onto `h` alone — box size and `Nint` do
not enter independently — and changes sign between `h = 0.294` and `h = 0.229`:

| h | 0.227 | 0.229 | 0.294 | 0.382 | 0.577 |
|---|---|---|---|---|---|
| maxReLambda | −0.00111 | −0.00104 | +0.00199 | +0.0116 | +0.0266 |

Above the threshold the drift operator has growing modes, so the Lyapunov
solution is not a covariance at all and can be indefinite — which is how
`readActiveStress`'s `fxfx` came out negative on most of the `set2_2_comove`
set.

**Cost.** Bartels–Stewart is `O((3N²)³)` in the *point count*. Holding
`h ≈ 0.23` uniformly over box40 needs `N ≈ 348`, i.e. `n ≈ 363k` and a ~1.1 TB
dump. Not possible.

The unstable mode is a **core** object — 85 % of its power is in the `Q2`
sector, 59 % inside `r < 3`, 90 % inside `r < 8`, and its density per unit area
falls ~85× from the core to `r ≈ 7.5`
(`claude_experiments/step0_unstable_mode.wls`). So only the core needs the fine
spacing, and a graded grid is the right answer.

---

## Result 1 — the adjoint pair on a graded grid

`lyapunov_solver.m` rests on one identity (its `:745` gate, `dean.tex`
`eq:G-lap-identity`):

```
Lap + (Gx^T Gx + Gy^T Gy) == 0
```

with a **plain** transpose. That is the adjoint in the unweighted ℓ² inner
product, and it is the correct adjoint only because a uniform mesh has mass
matrix `h² I` — a scalar multiple of the identity, which cancels from both
sides. On a graded grid the inner product is `⟨u,v⟩ = Σ w_k u_k v_k` with `w`
varying, the adjoint of `G` becomes `M⁻¹ Gᵀ S_f`, and a plain transpose silently
breaks the discrete fluctuation–dissipation balance. `dean.tex:1526` prices that
failure at *"exactly one half, at every h, for a scheme that is second-order
accurate and passes every classical consistency check"*, and `:1530` notes that
refinement cannot repair it.

### The fix is one similarity transform

With `M = diag(w)` the node dual cells (the trapezoid weights — exactly what
`cdGridPade` already returns), `S_f = diag(d)` the face spacings, and
`(G u)_m = (u_m − u_{m−1})/d_m`:

```
L    = −M⁻¹ Gᵀ S_f G                 self-adjoint in the M-inner product
Ghat = S_f^(1/2) G M^(−1/2)          psiHat = M^(1/2) psi
```

Then `Ghatᵀ Ghat = −M^(1/2) L M^(−1/2)`, so

```
Lhat = −Ghatᵀ Ghat        with a PLAIN transpose
```

— the uniform identity verbatim. The Kronecker lift survives because the
`M_y^(1/2)` factors cancel out of `Ghat_x`, so each 2-D operator is still a lift
of a purely 1-D hatted factor. **When `w` is constant `Ghat == G` exactly**, so
the uniform solver is the special case and the regression test is meaningful.

Three consequences:

- **The collocated derivative.** `D = M⁻¹ Pᵀ S_f G`. Then `M D` is the
  antisymmetric shift matrix `(δ_{j,i+1} − δ_{j,i−1})/2` *independent of the
  spacing*, so `Dhatᵀ = −Dhat` holds exactly — the identity `divBlk`'s
  transposed-gradient form depends on (`lyapunov_solver.m:1038-1039`). And it
  evaluates to `(u_{i+1} − u_{i−1})/(x_{i+1} − x_{i−1})`, the natural graded
  central difference.
- **`pf` and `locPf` disappear.** `Qhat_div = Ghatᵀ diag(M_face) Ghat`; the local
  `Λ` term becomes `diag(2Λ/ρ₀)` with no metric factor at all. The `1/h²` at
  `:1036` was always `1/(cell area)` (`dean.tex:1561`).
- **Second order, not fourth.** `L` is *defined* as `−M⁻¹GᵀS_f G`. A 4th-order
  `FiniteDifferenceDerivative` Laplacian — what `comoving_core.m:139` builds for
  the deterministic steady solve — does **not** satisfy the Gram identity, so it
  would break `R ≡ 1` and the flat pedestal. Deliberate trade, the same one the
  uniform Neumann branch already makes.

---

## Result 2 — the flat far-field pedestal

`eq:lattice-bz` (`dean.tex:1201`) is a **Brillouin-zone average**. It assumes
translation invariance and has no meaning on a graded grid, so `lEff` cannot be
calibrated the way `lyapunov_solver.m:449` does it.

It does not need to be. The bare discrete pedestal is `Σ_ii = T/(k_i w_i)` —
mesh-dependent *by construction*, because `δ(r−r') → δ_ij / w_i`. Build the
smoother from the analytic Gaussian sampled on the grid and weighted by the
quadrature:

```
W = Ghauss . M,   Ghauss[i,j] = g(x_i − x_j),   g of width lNoise/Sqrt[2]
```

(the `/√2` because `W` carries symbol `exp(−k²ℓ²/4)` and is applied twice). Then

```
Σ_ii = T Σ_j g(x_i − x_j)² w_j / k_j  →  (T/k_i) · 1/(2π ℓ²)
```

**The `1/w_i` of the bare pedestal is cancelled exactly by the `w_j` inside
`W`** — that is the whole mechanism, and it is what makes a varying `h` legal.
The result is `eq:iso-gauss`, mesh-independent.

### The per-node normalisation — easy to miss

Trapezoid quadrature of a Gaussian is spectrally accurate only at **equal**
spacing (Poisson summation). With `h` varying it reverts to `O(h²)`: measured on
a 6:1 graded grid the raw coincidence limit is off by **3.6 %**, against the
0.2 % flatness the uniform N131 reference achieves.

So each 1-D factor is row-scaled by

```
c_i = Sqrt[ Pcont_i / Pdisc_i ]
Pdisc_i = Σ_j gD(x_i,x_j)² w_j        Pcont_i = ∫ gD(x_i,y)² dy
```

the continuum integral evaluated once on a fine **uniform** subgrid where the
trapezoid rule *is* spectral. This makes the coincidence limit exact at every
node by construction, at no cost, and without breaking separability or the
congruence (`W → diag(c) W` is still a congruence; `W` need not be symmetric,
only `Q → W Q Wᵀ`).

Normalising against the **continuum** rather than against a constant is what
keeps the physics: near the wall the Dirichlet images genuinely suppress `gD`,
and that suppression appears in `Pcont` and `Pdisc` alike, so it survives. Only
the quadrature error is removed. Measured:

| | raw | corrected |
|---|---|---|
| `\|x\|,\|y\| < 12`, 6:1 graded | 3.6e−2 | **3.6e−4** |
| at the wall (x = 13.8) | 0.6158 | 0.6158 (physical droop kept) |
| uniform grid, `max\|c−1\|` | — | **2.2e−16** (inert) |

This is the graded-grid replacement for the scalar `lEff` calibration, and like
it, it is a **pedestal-channel** statement: it makes the ρ-sector coincidence
limit exact and does not claim to remove lattice corrections from the `δQ`
sector.

---

## Files

| file | contents |
|---|---|
| `graded_grid.m` | the grading law and the equidistribution that derives the point count from it |
| `graded_ops.m` | `G, P, L, D` graded, the `Ghat` transform, the Kronecker lift, the adjoint-pair gate |
| `graded_noise.m` | exact von Mises closure, metric-aware `divBlk`, the Gaussian-quadrature smoother with Dirichlet images and the per-node normalisation |
| `graded_lyapunov.m` | sampling, `A`, `Q`, the solve, the metric undo, export |
| `verify.wls` | the acceptance tests |

## Usage

```mathematica
Get["lyapunov_graded_solver/graded_grid.m"];
Get["lyapunov_graded_solver/graded_ops.m"];
Get["lyapunov_graded_solver/graded_noise.m"];
Get["lyapunov_graded_solver/graded_lyapunov.m"];

gg  = nuGrid[15, 0.20, 1.2, 4., 0.35];      (* fbox, hCore, hMax, rFine, growth *)
ss  = nuSteadyFromFile["data/.../steady_box15_N131_Dir_lN1.73_z8_ld0.8.m"];
res = nuSolve[gg, modelParams, ss];
nuExport[res, outDir, "box15_graded_lN1.73"];
```

`nuExport` writes the same `.bin` + `_meta.m` layout `binaryReadTools.m` already
reads. Two differences, both deliberate:

- `h` is `Missing[...]` — there is no scalar spacing. `gridInt`, `gridX`,
  `gridY`, `weights`, `hMin`, `hMax` carry the real geometry, and `covFn`'s
  `ListInterpolation` handles a non-uniform grid unchanged, so `readCovRho` /
  `readCovQ` / `readCovRhoQ` / `readCorrQ` work as-is.
- `maxReLambda`, `minEigQ`, `residual`, `lapIdentity` go **into the metadata**,
  not into a side file. In the uniform solver they live only in
  `sigma_rho_*.m`, which is why every dump could be unstable without anyone
  noticing.

`readActiveStress.m` is the one reader that does **not** carry over: its
`gradStencil[h]` is a scalar-`h` central difference and needs per-node stencil
coefficients.

## A caveat that is physics, not numerics: sigma_Q is not converged at box15

Head to head on an **identical** uniform grid (fbox 5, Nint 43) this solver and
`lyapunov_solver.m` agree to **0.12–0.14 %** on all three variance fields, and
their `maxReLambda` agree to nine digits. The residual is the reference's own
lattice correction (`lNoiseEff = 1.73375` against the physical `1.73`); the
Gaussian smoother here has no lattice error, so it uses `lNoise` directly.

But the graded box15 run and the stored uniform N131 run differ by a factor ~2
in `sigma_Q`, while agreeing to 0.6 % RMS in `sigma_rho`. That is not a
discretisation error. The arithmetic:

```
maxReLambda    N131 -0.00111    graded -0.00215     ratio 0.5172
sigma_Q2(0,0)  N131  85.84      graded  43.98       ratio 0.5124
```

— agreement to 0.9 %. The defect steady state at these parameters sits *just*
inside its stability threshold, and the least-damped mode is the same
85 %-`Q2`, core-localised mode that Step 0 found. A nearly-marginal mode
contributes `~ Q/(2|lambda|)` to the variance, so `sigma_Q` scales as
`1/|maxReLambda|` — and the residual damping is grid-dependent.

**Consequence.** `sigma_rho` and the pedestal are converged and trustworthy;
they do not couple to that mode. `sigma_Q` at box15 is controlled by how
strongly a given grid happens to damp a barely-stable mode, and is therefore not
converged on *either* grid. Reporting it requires a larger box, or an account of
why the steady state is this close to marginal — see the `\corr{}` note at
`dean.tex:782` about the comoving advection term that no longer cancels.

## Verification

```
wolframscript -file lyapunov_graded_solver/verify.wls
```

The sharpest test is that the graded machinery reproduces the uniform solver:
operators, `A` blocks and `Q` match to ~1e−15, and `maxReLambda` at box15/N51
reproduces the recorded `0.0265508793108064` to 9.7e−9.
