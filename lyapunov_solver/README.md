# Lyapunov solver — comoving drift, calibrated cutoff, two walls

Equal-time Gaussian fluctuations about the **comoving** $\rho$–$Q$ steady state, by direct
solution of the discrete Lyapunov equation

$$\tilde A\,\mathbf C + \mathbf C\,\tilde A^{\mathsf T} + \boldsymbol Q = 0,
\qquad \psi = (\delta\rho,\,\delta Q_1,\,\delta Q_2)^{\mathsf T}$$

This replaces the five `../lyapunov_solver_module*.m` files. Three things are different, and the
first is a correctness fix, not a refactor:

1. **The drift is the comoving one.** The old modules drop the advection $\nabla\cdot(\delta Q\,\mathbf w_{\mathrm{ss}})$
   entirely. See [The advection term](#the-advection-term).
2. **`lNoise` is calibrated** against the continuum theory by default, so a run at a given
   `lNoise` is directly comparable with `dean.tex` `eq:iso-gauss` without needing a fine grid.
   See [Cutoff calibration](#cutoff-calibration).
3. **Both** consistent boundary conditions of `eq:bc-triple` live here, chosen by an argument.

The steady state **must** come from `../comoving_defect/comoving_steady_solver.m`; a lab-frame
profile is refused.

## Files

| file | what it is |
|---|---|
| `lyapunov_solver.m` | the module — `SolveFluctuationsLyapunov`, `visualizeFluctuationsComoving`, `exportDataComoving` |
| `lyapunov_solver_optimized.m` | optional Python backend for the solve, ~3.5× faster and ~5× smaller in RAM ([Part 5](#part-5--the-python-backend)) |
| `lyap_solve.py` | the backend itself — real Schur + recursive Bartels–Stewart |
| `verify.wls` | small-grid verification; appends to `verify_results.csv` |
| `verify_optimized.wls` | proves the Python backend reproduces the reference kernel |
| `README.md` | this file |

Self-contained: it does **not** load `../iterative_solver_module.m`, and it opens its own
private context `LyapunovSolver`Private``, so it can share a kernel with an old module (which
`verify.wls` needs for its regression test). The one name to watch is `SolveFluctuationsLyapunov`
itself, which the superseded `../lyapunov_solver_module.m` also defines — don't load that file.

## Requirements

Mathematica 14.1. Rasterising plots needs the front end, which crashes headless unless `DISPLAY`
is cleared, so `unset DISPLAY` first. Every script locates the module relative to itself.

---

# Part 1 — Using it

```mathematica
Get[FileNameJoin[{repoRoot, "comoving_defect", "comoving_steady_solver.m"}]];
Get[FileNameJoin[{repoRoot, "lyapunov_solver", "lyapunov_solver.m"}]];

modelParams = {B -> 10000., a -> -0.172, b -> 0.1, L -> 0.05368,
               \[Zeta] -> 10., \[Xi]0 -> 5., \[Xi]r -> 0.1,
               \[Rho]0 -> 0.3, \[CapitalLambda] -> 0.1,
               lNoise -> 1.73};            (* the PHYSICAL cutoff *)

solverParams = {\[Delta]mesh -> 0.04, hMax -> 2., box -> 15.,
                maxSteps -> 25, relativeDiff -> 1.*^-8};

ss  = SolveActiveNematicComovingSteady[modelParams, solverParams, "Dirichlet"];
res = SolveFluctuationsLyapunov[{Nint -> 101, fbox -> 15},
                                modelParams, ss, "Dirichlet"];

res["sigmaRho"]     (* Var(drho) as a matrix, [[j,i]] = Sigma(x_i, y_j) *)
res["lNoiseEff"]    (* the length W was actually built from *)
res["boxFactor"]    (* predicted pedestal / eq:iso-gauss *)
```

`bcType` is the fourth argument and defaults to `"Dirichlet"`. `fbox` and `Nint` must be
**integers** (`h` is kept exact so odd `Nint` puts the defect core on a node), and `fbox` must not
exceed the steady state's own `box` or the steady interpolant is extrapolated at the wall.

### Optional `lyapunovParams` keys

| key | default | meaning |
|---|---|---|
| `rescaleCutoff` | `True` | calibrate `lNoise` against `eq:lattice-bz` |
| `advection` | `True` | include $\nabla\cdot(\delta Q\,\mathbf w_{\mathrm{ss}})$; `False` reproduces the superseded modules |
| `scratchDir` | `Automatic` | where the Python backend stages `A`/`Q`/`C`. `Automatic` picks `$TMPDIR`, then `/scratch/$SLURM_JOB_ID`, then `/scratch`. Ignored by the default kernel. |

### What comes back

The usual `"C"`, `"sigmaRho"`, `"sigmaRhoFn"`, `"sigmaQ1"`, `"sigmaQ2"`, `"gridInt"`, `"h"`,
`"Nint"`, `"fbox"`, `"eigenvalues"`, `"maxReLambda"`, `"residual"`, `"residualRefined"`,
`"timings"`, plus:

| key | meaning |
|---|---|
| `"u"` | the frame velocity taken from the steady state |
| `"peclet"` | cell Péclet number $\lvert w\rvert h/2\kappa$ — the advection's resolution |
| `"advIdentity"` | the $\mathrm{Adv}\cdot\mathbf 1$ tripwire, relative, interior nodes |
| `"lNoise"` / `"lNoiseEff"` | physical cutoff / the length `W` was built from |
| `"pedestalPredicted"`, `"pedestalContinuum"`, `"boxFactor"` | see [The `boxFactor`](#the-boxfactor) |
| `"lapIdentity"` | the adjoint-pair gate, exactly `0` |
| `"consA"`, `"consAT"`, `"consQ"`, `"nullOverlap"` | sealed-wall gates (`Missing` under Dirichlet) |

Output files are tagged by the **physical** `lNoise`, so a sweep stays keyed to the physics;
`lNoiseEff` is written alongside.

---

# Part 2 — How it works

## The advection term

`dean.tex` eq. `Q-comoving` is

$$\partial_t \mathbf{Q} = \nabla\cdot\Big(\mathbf{Q}\big(\mathbf{u} + \tfrac{B}{\xi_0}\nabla\rho\big)\Big) + \tfrac{4}{\xi_r}\mathbf{H}_\text{active}$$

with $\mathbf{u}$ a **constant** vector, obtained as the Lagrange multiplier conjugate to the pinning
constraints $Q_1(0,0)=Q_2(0,0)=0$. Linearising gives (`dean.tex` eq. `linearized-Q-equation`)

$$\partial_t \delta\mathbf{Q} = \tfrac{B}{\xi_0}\nabla\cdot(\mathbf{Q}_{\mathrm{ss}}\nabla\delta\rho)
  + \underbrace{\nabla\cdot(\delta\mathbf{Q}\,\mathbf w_{\mathrm{ss}})}_{\text{this}}
  + \tfrac{4}{\xi_r}\delta\mathbf{H}_\text{active} + \mathbf f^Q,
  \qquad \mathbf w_{\mathrm{ss}} = \mathbf{u} + \tfrac{B}{\xi_0}\nabla\rho_{\mathrm{ss}} .$$

**The old modules drop it.** Their `A22`/`A33` carry no first-derivative operator at all. The
reason is archaeological: earlier versions of `dean.tex` wrote the $\delta Q$ equation with
$+\frac{B}{\xi_0}(\nabla\rho_{\mathrm{ss}}\cdot\nabla)\delta Q$ on the *left*, which is $-(\mathbf{u}\cdot\nabla)\delta Q$
evaluated with the superseded pointwise estimate $\mathbf{u} = -\frac{B}{\xi_0}\nabla\rho$
(`../FRIDGE/wrong_comoving.tex:45`). With that substitution the advection cancels identically.
With $\mathbf{u}$ a constant it does not: the $\mathbf{u}$ of the frame cancels the $\mathbf{u}$ inside $\mathbf w_{\mathrm{ss}}$, and the
$\frac{B}{\xi_0}\nabla\rho_{\mathrm{ss}}$ advection **survives** — and it is the larger of the two pieces. At
the core, $\frac{B}{\xi_0}\lvert\nabla\rho_{\mathrm{ss}}\rvert \simeq 1.07$ against $\lvert\mathbf{u}\rvert = 0.28$ μm/min.

### The double-counting trap

Written in conservative form,

```mathematica
Adv = -(Transpose[Gx] . (wxE * Px) + Transpose[Gy] . (wyE * Py))
```

with `Px`, `Py` the node→face averaging operators (the $1/2$-valued twins of `Gx`, `Gy`). Because
$\mathbf{u}$ is constant,

$$\nabla\cdot(\delta\mathbf{Q}\,\mathbf w_{\mathrm{ss}}) = (\mathbf w_{\mathrm{ss}}\!\cdot\!\nabla)\delta\mathbf{Q} + (\nabla\cdot\mathbf w_{\mathrm{ss}})\,\delta\mathbf{Q},
\qquad \nabla\cdot\mathbf w_{\mathrm{ss}} = \tfrac{B}{\xi_0}\nabla^2\rho_{\mathrm{ss}}$$

so **`Adv` already contains the `(B/xi0) diag[rhoLap]` term the old modules have**, and that term
is deleted from `A22`/`A33` here. Keeping both double-counts it. The identity is asserted every
run in `[4b]` as `Adv.1` against `(B/xi0) rhoLap`, on interior nodes only — the constant field is
not representable under either wall condition, so the outermost rows differ legitimately.

Two things the advection does *not* disturb:

- **Mass conservation.** `Adv` enters `A22`/`A33` only. The $\rho$ row block and the $\delta\rho$
  column block are untouched, so the sealed-wall gates pass exactly as before. $\delta Q$ is not a
  conserved quantity, so nothing is owed.
- **The need for upwinding.** Cell Péclet is $\lvert w\rvert h/(2\kappa)$ with $\kappa = 4K'/\xi_r = 4.15$;
  measured, it is $0.06$–$0.11$ across every grid tried. Centred differencing is correct here by an
  order of magnitude. It is printed every run and warned on above 1.

At a sealed wall $\delta Q = 0$ by anchoring, so the advective flux through the wall vanishes and
the interior-face `Gx`/`Pav` pair is already the right operator.

## Cutoff calibration

The noise carries a Gaussian correlation length, imposed as a congruence $\boldsymbol Q\to W\boldsymbol QW$ with
$W = \exp[(\ell^2/4)\nabla^2]$, so the continuum equal-point variance is UV-finite:
$\Sigma^\text{iso} = \zeta\rho_{\mathrm{ss}}/(2\pi B\rho_0\ell_{\mathrm n}^2)$ (`eq:iso-gauss`).

But the **lattice** heat kernel decays more slowly than the continuum Gaussian, so it keeps too
much short-wavelength noise and the pedestal overshoots. `dean.tex` `eq:lattice-bz` gives the
factor in closed form: with $s=\ell/h$ and $F(s)=e^{-s^2}I_0(s^2)$,

$$\Sigma_\text{lattice}(\ell) = \frac{\zeta\rho_{\mathrm{ss}}}{B\rho_0h^2}\,F(s)^2 .$$

So by default, given the **physical** `lNoise`, the module solves

$$\boxed{\;F(\ell/h) = \frac{1}{\sqrt{2\pi}\,t}\;},\qquad t = \ell_{\mathrm n}/h$$

for the solver length $\ell$ (`lNoiseEff`) and builds `W` from that. Three properties say this is
the right rule rather than a fitted one:

- $F$ decreases monotonically from $F(0)=1$ to $F(s)\sim 1/(\sqrt{2\pi}s)$, so $\ell\to\ell_{\mathrm n}$ as
  $h\to0$: the correction vanishes on fine grids. Measured, `lNoiseEff/lNoise` at `fbox = 20`,
  `lNoise = 3` runs $1.0544,\,1.0130,\,1.0033,\,1.00085,\,1.00021$ for `Nint = 21…321`.
- It is solvable exactly when $t\ge 1/\sqrt{2\pi} = 0.3989$, i.e. $\ell_{\mathrm n} \ge h/\sqrt{2\pi}$ — which
  is `eq:ln-eff` recovered as the **solvability edge**, where $\ell=0$ and the bare grid already
  *is* the correct cutoff. Below it no smoothing can help (smoothing only removes variance) and
  the module emits `::uncalibratable` naming the minimum `Nint` and returns `$Failed`.
- Expanding, $\ell \simeq \ell_{\mathrm n} + h^2/(8\ell_{\mathrm n})$, which inverts the $+\frac14(h/\ell)^2$ bias the old
  modules documented. At $\ell_{\mathrm n} = 1.21h$ that is $+8.5\%$ in $\ell$, undoing the $17\%$ variance
  error tabulated in `../jobscripts/set2_2_comove/fluctuations.wls:58-69`.

**Numerics.** $x = s^2$ reaches a few thousand, where `BesselI[0, x]` overflows as a machine
number while `Exp[-x]` underflows. `scaledI0` evaluates at 50 digits above $x = 50$ — the two
factors are $\sim10^{\pm1085}$ there and their product is $O(10^{-2})$, which arbitrary precision
tracks exactly. Verified against a direct high-precision evaluation to machine epsilon out to
$x = 2500$. The inversion is by bisection on a smooth monotone function; it is called a few dozen
times per run, so robustness beats cleverness.

**Scope, honestly.** This is **one scalar, calibrated on the isotropic $\rho\rho$ pedestal**. It
cannot simultaneously be exact for the anisotropic term (`eq:aniso-kernel`) or for the $\delta Q$
sector, whose lattice correction factors are different functions. It buys direct comparability of
the pedestal, which is the quantity `dean.tex sec-cutoff-numerics` tests, and nothing more.

`rescaleCutoff -> False` restores the old behaviour.

## The `boxFactor`

In the $\rho$-only sector $W$ commutes with the drift, so $\mathbf C = W\mathbf C_0W$ and the
finite-box pedestal at the core follows in closed form from the very basis $W$ was built in. The
module computes it and reports the ratio to `eq:iso-gauss` as `boxFactor`.

Under **Dirichlet** it is $1$ to six digits at every grid and every `lNoise` — the BZ closed form
*is* the finite-box answer there, which is what makes it a legitimate calibration target.

Under **Neumann** it is not, and that is physics rather than error. The canonical projection
removes the uniform mode, so a correlation length comparable with the box leaves genuinely less
room for a density fluctuation. Measured, the factor depends on $\ell_{\mathrm n}/\text{fbox}$ **alone** —
identical at `fbox` 10, 20 and 40, and flat in `h`, which is what distinguishes geometry from
discretisation:

| $\ell_{\mathrm n}/\mathrm{fbox}$ | 0.05 | 0.1 | 0.2 | 0.3 | 0.4 |
|---|---|---|---|---|---|
| `boxFactor` | 0.997 | 0.985 | 0.938 | 0.860 | 0.750 |

So under a sealed wall, **do not expect the pedestal to equal `eq:iso-gauss`** — expect it to
equal `boxFactor` times it. Keep $\ell_{\mathrm n}/\text{fbox}\lesssim0.1$ if you want them within 1.5%.

## Boundary conditions

Posing the problem on a bounded domain is **three** choices, not one (`dean.tex` `eq:bc-triple`):
on the steady state, on the fluctuating field, and on the noise flux. Only two of the eight
combinations are coherent, and both are implemented:

| `bcType` | | |
|---|---|---|
| `"Dirichlet"` | **reservoir** | $\rho_{\mathrm{ss}}$ pinned, $\psi\rvert_{\partial\Omega}=0$, flux **free** — a window cut from a larger colony |
| `"Neumann"` | **sealed** | $\hat{\mathbf n}\cdot\nabla\rho_{\mathrm{ss}}=0$, $\hat{\mathbf n}\cdot\nabla\psi=0$, flux **blocked**, $\delta Q$ **anchored** |

The pairing is not a convention. Whether $\psi$ is pinned at the wall and whether $\mathbf F\cdot\hat{\mathbf n}$
vanishes there are the same question — can a bacterium cross — asked of the mean and of the
fluctuation; answering it differently for the two puts the wall out of detailed balance by fiat.

**Why $\delta Q$ is anchored at a sealed wall** and not left Neumann: sealing is a statement about
*mass*, and only $\rho$ is conserved, so the flux argument fixes the $\rho$ sector alone. Leaving
$\delta Q$ Neumann too leaves $\tilde A$ with exactly three unstable eigenvalues that **converge**
rather than refine away — the two defect translation modes and the global director rotation. A
free-anchoring wall images a $+1/2$ defect attractively, so the centred defect is a saddle; a wall
that blocks mass but exerts no torque cannot confine the defect and no equal-time covariance
exists. The anchored Laplacian `L1q` therefore enters `A22`/`A33` and nothing else.

What differs between the two, and only this:

| | `"Dirichlet"` | `"Neumann"` |
|---|---|---|
| grid | node-centred, $h=2\,\text{fbox}/(N{+}1)$ | cell-centred, $h=2\,\text{fbox}/N$ |
| `G1` | $(N{+}1)\times N$, edges incl. the wall | $(N{-}1)\times N$, interior faces only |
| `L1` | truncated `FiniteDifferenceDerivative` | **defined** as $-G_1^{\mathsf T}G_1$ |
| `Dxy` | `Kron[D1, D1]` | symmetrised double divergence |
| $\delta Q$ Laplacian | `Lap` | `Lapq` (anchored) |
| `W` basis | Dirichlet sine, vanishes at the wall | DCT-II, **reflects** off the wall |
| gates | `lapIdentity` | + `consA`, `consAT`, `consQ` |
| null mode | none | $(1,0,0)$ exact, canonically projected |

Neumann `Nint = M+1` has exactly the same $h$ as Dirichlet `Nint = M`, so the two are directly
comparable grid for grid.

## Inherited, and not to be "simplified"

Every item cost real debugging time in the modules this one replaces.

**Consistent stencils ($R\equiv1$).** The *diagonal* noise flux terms use `Gx`, `Gy` with
face-sampled coefficients; the *cross* terms keep the **collocated** differences with node-sampled
coefficients. That split is forced, not a compromise: a corner-interpolated cross term picks up a
phase $\cos((q-p)/2)$ from the half-cell offset between an $x$-face and a $y$-face, which destroys
the parity in $p$ behind the continuum angular cancellation for traceless $\mathbf{Q}_{\mathrm{ss}}$ and injects a
spurious anisotropic pedestal $0.087\,Q_2\,\zeta/(B\rho_0h^2)$ — an ~8% $h^{-2}$ modulation at
$S\sim0.86$.

**Exact rationals, then `ToPackedArray`.** `G1`, `L1`, `D1` are built from exact rationals so the
adjoint-pair identity can be tested with `== 0`. If any of that exactness survives into the dense
arrays they come back unpacked, every BLAS path is lost, and the back-substitute runs ~600×
slower.

**Every discrete divergence is minus a transposed gradient**, so $\mathbf 1^{\mathsf T}B_i = 0$
identically. Invisible under Dirichlet (where $D_1^{\mathsf T} = -D_1$); under Neumann it is the
difference between conserving mass and creating it in the four corner cells at $O(h^{-2})$.

**`Q -> W.Q.W` in place, one block pair at a time.** Building all nine blocks and
`ArrayFlatten`-ing them holds ~42 GiB of temporaries at `Nint = 131`.

**The sealed null mode.** $\lambda=0$ is exact, so the Lyapunov equation is solvable but not
unique. `eigLyapSolve` zeroes row and column `iNull` of the transformed right-hand side (the
canonical projection) and shifts that one denominator off zero — `0./0.` is `Indeterminate`, not
`0`, and a single such entry is smeared over all of `C` by the two matmuls.

---

# Part 3 — Verification

```bash
unset DISPLAY
wolframscript -file verify.wls        # tests 1-8,  ~50 s
wolframscript -file verify.wls all    # + end-to-end comoving solves, ~70 s
```

29/29 pass. Results appended to `verify_results.csv`.

| test | measured |
|---|---|
| adjoint-pair identity, both walls | exactly `0` |
| sealed-wall `consA`/`consAT`/`consQ` | $9\times10^{-16}$ |
| canonical: total $\rho\rho$ covariance | $-3\times10^{-18}$ |
| advection identity `Adv.1` vs `(B/xi0) rhoLap` | $3\times10^{-12}$ (D), $4\times10^{-12}$ (N) |
| `W` matches the exact finite-box pedestal | 0.07–0.08%, both walls |
| calibrated run matches it | 0.07–0.08%, both walls |
| BZ closed form = exact finite box, Dirichlet | $1\times10^{-15}$ |
| calibrated pedestal = `eq:iso-gauss`, Dirichlet | `boxFactor` $=1$ |
| `boxFactor` scale-invariant in $\ell_{\mathrm n}/\mathrm{fbox}$ | exactly `0` |
| $\ell_{\mathrm n}<h/\sqrt{2\pi}$ | clean `$Failed`, succeeds once refined |
| `lNoiseEff/lNoise` $\to1$ | $1.0544\to1.00021$ over `Nint` 21…321 |
| **$\mathbf w_{\mathrm{ss}}=0$ reproduces the old modules** | **exactly 0**, both walls |
| $\Sigma = $ pedestal $(I-\mathbf{11}^{\mathsf T}/N^2)$ | $9\times10^{-4}$ (the $\delta Q$ coupling) |
| end-to-end, real comoving steady state | converges, residual $\sim10^{-15}$ |

The 0.07–0.08% residual against the $\rho$-only prediction is the $\delta Q$ coupling that
`eq:sde` drops, which `dean.tex sec-cutoff-numerics` independently puts at 0.3–0.5%.

The decisive row is the regression: with $\rho_{\mathrm{ss}}$ constant and $u=0$, $\mathbf w_{\mathrm{ss}}$ vanishes
identically, so this module and the superseded ones must agree — and they do to *exactly* zero,
with $\mathbf{Q}_{\mathrm{ss}}\neq0$ exercising the closure, `A21`/`A31`, the cross noise terms and `W`. That
isolates the new advection from every other change and needs no compatibility flag.

---

# Part 4 — Caveats

**`maxReLambda > 0` on coarse grids is expected.** It is the defect core being under-resolved and
it refines away fast. Measured, reservoir wall, `fbox = 10`, $\ell_{\mathrm d} = 0.8$:

| `Nint` | 15 | 21 | 27 | 33 |
|---|---|---|---|---|
| $h$ | 1.25 | 0.909 | 0.714 | 0.588 |
| $\ell_{\mathrm d}/h$ | 0.64 | 0.88 | 1.12 | 1.36 |
| `maxReLambda` (advection on) | 0.498 | 0.0964 | 0.0342 | 0.0117 |
| `maxReLambda` (advection off) | 0.423 | 0.0789 | 0.0275 | 0.0056 |
| mean $\Sigma_\rho$ | 1.592e-4 | 1.540e-4 | 1.509e-4 | 1.488e-4 |

It falls faster than $h^3$ (closer to $h^5$ here) while $\Sigma_\rho$ converges, and the advection
changes it only marginally — so the advection is not the cause. Read the trend in `Nint`, never a
single grid.

**The pedestal calibration is the $\rho\rho$ channel only.** See
[Cutoff calibration](#cutoff-calibration). The anisotropic term and the $\delta Q$ sector still
carry uncorrected lattice factors.

**Under a sealed wall the pedestal is deliberately not `eq:iso-gauss`.** See
[The `boxFactor`](#the-boxfactor).

**The closure is not PSD above $S = 0.618$** with the leading-order moments, which is why the
exact max-entropy ones are used throughout; the defect steady state exceeds that bound over most
of the domain. `condition (ii)` is checked every run, on the *smoothed* $\boldsymbol Q$.

**Neumann costs several times more per solve** than Dirichlet at equal size, and the reason was
never established — three plausible explanations were measured and ruled out in the module this
one inherits from.

**Cost is the `Eigensystem`**, $O((3N^2)^3)$. `Nint = 33` is about a second; `Nint = 61` is
already minutes. That is what keeps `verify.wls` on small grids. The Python backend of
[Part 5](#part-5--the-python-backend) cuts it by ~3.5× but does not change the exponent.

**$u$ is not converged in box size** — that is a property of the steady state, not of this
module. Doubling the comoving solver's `box` from 8 to 16 μm moves $u$ by ~11%; the reservoir and
sealed walls bracket the answer from opposite sides. See `../comoving_defect/README.md`.

---

# Part 5 — The Python backend

Optional. It replaces **only** block [6] — the dense Lyapunov solve, which is ~98% of the
runtime. Everything else (stencils, closure, cutoff congruence, boundary conditions,
extraction, tags, export) stays in `lyapunov_solver.m` and is *shared, not copied*.

```mathematica
Get[FileNameJoin[{repoRoot, "lyapunov_solver", "lyapunov_solver.m"}]];            (* first  *)
Get[FileNameJoin[{repoRoot, "lyapunov_solver", "lyapunov_solver_optimized.m"}]];  (* second *)
```

**Load order matters**: loading the base module resets `$LyapunovSolveKernel` to the reference
implementation, so the override must come after. Nothing else in a driver changes — same call,
same returned association, same exports, same tags.

## Why it is faster

The reference kernel diagonalizes $\tilde A$. That is nonsymmetric, so the eigenvectors are
complex and the whole back-substitution runs in Complex128. This backend takes the **real Schur
form** and solves the triangular Lyapunov equation by a recursive Bartels–Stewart that bottoms
out in LAPACK `?trsyl` on small leaves, leaving the bulk of the flops in GEMM. Real64 throughout.

Measured on one exclusive EPYC 9655 node, both runtimes with a working BLAS:

| $n$ | Mathematica eigen route | this backend | speedup |
|---|---|---|---|
| 2000 | 6.41 s | 1.71 s | 3.7× |
| 3000 | 12.10 s | 3.46 s | 3.5× |
| 4000 | 23.53 s | 5.99 s | 3.9× |

Two routes that look obvious and are not: `scipy.linalg.solve_continuous_lyapunov` calls the
serial, unblocked LAPACK `?trsyl` and is *slower* than the eigen route it would replace (8.9 s
against 2.9 s at $n=2000$); and Mathematica's own `LyapunovSolve` is slower still (74 s at
$n = 4000$, 3.3× the eigen route) and scales worse than cubic.

## Why the memory matters more

Peak RAM in the Mathematica kernel, measured as multiples of one dense $n\times n$ Real64:

| route | peak |
|---|---|
| reference (eigendecomposition) | **15.5×** |
| `LyapunovSolve` | 7.0× |
| this backend | $A + Q + C$, i.e. **3×** |

The reference route has `eigvecs`, `Pev`, `lamSum`, `num`, `num/lamSum` and `Transpose[Pev]` all
live at once and all Complex128. This backend never forms a complex matrix. At `Nint = 151`
($n = 68403$, one dense $n\times n$ = 37.4 GiB) that is ~580 GiB against ~112 GiB — the
difference between the 2000G queue class, which reserves for days, and something that starts now.

Note the printed estimate in block [6] counts only `Adense` + `Qdense` and so understates the
reference route's peak by about 8×.

## Scope: Dirichlet only

Under Neumann the drift is singular **by design** — the conserved uniform mode — so
Bartels–Stewart would need an explicit null-mode deflation in the Schur basis that is not
implemented. Neumann therefore delegates to `LyapunovReferenceKernel`: unchanged, still correct,
just not accelerated. This is scope, not a failure path, and it is announced on every call.

## Failure policy: abort

Any backend failure — interpreter missing, numpy/scipy missing, scratch full, a non-finite
residual — **aborts the run**. It does not silently fall back, so every file in a sweep is known
to have come from one numerical path. The abort names the cause and leaves the scratch directory
in place for inspection.

## Choosing the interpreter

Handled automatically; you do **not** need `module load python/3.14.7` in a jobscript. Candidates
are tried in order and each is probed for numpy *and* scipy before use:

1. `$LYAPUNOV_PYTHON` (environment) or `$LyapunovPython` (a path you set)
2. the `python/3.14.7`, `python/3.13.9`, `python/3.12.11` module interpreters
3. a bare `python3` on `PATH`

Running a module interpreter by absolute path is **not** sufficient on its own — it dies with
`error while loading shared libraries: libpython3.14.so.1.0`, because `PYTHONHOME`,
`LD_LIBRARY_PATH` and `PYTHONPATH` come from the modulefile and `libpython` lives in `lib64`, not
`lib`. The backend reconstructs that environment itself. If a jobscript *has* done `module load`,
candidate 3 picks it up and inherits a correct environment anyway.

**Do not point it at `anaconda3`.** On this cluster that resolves to `/usr/bin/python3`, whose
numpy is linked against reference netlib BLAS: GEMM at 6.4 GFLOPS against ~4000 for the module
build, flat in thread count. The numpy/scipy probe is what rejects it. Nor is MKL the right
choice here — the nodes are AMD, and `anaconda/2022.10`'s MKL takes 6.6 s on the $n=3000$
eigensolve against OpenBLAS's 3.5 s.

## Scratch

`A`, `Q` and `C` are staged through `scratchDir` as row-major Real64, so the directory needs
$3\cdot 8n^2$ bytes free — 112 GiB at `Nint = 151`. This is checked *before* the solve, because
running out of scratch three hours in is the expensive way to discover it. Prefer the node-local
NVMe (`$TMPDIR` under SLURM, 28 TB on these nodes) over shared GPFS. `A.bin`/`Q.bin`/`C.bin` are
deleted on success unless `$LyapunovKeepScratch = True`, and always kept on failure.

`$LyapunovPythonThreads` defaults to `SLURM_CPUS_PER_TASK` when set, else `$ProcessorCount` — the
allocation, not the node, is what the job is entitled to.

## Verification

```
wolframscript -file lyapunov_solver/verify_optimized.wls
```

20 checks in two groups: kernel-level (same $A$, $Q$ through both kernels) and end-to-end (the
full `SolveFluctuationsLyapunov` on a synthetic comoving steady state). Agreement with the
reference is ~1e-15 on `C`, `sigmaRho`, `sigmaQ1`, `sigmaQ2` and every exported scalar, with an
identical `runTag`. `verify.wls` still passes 25/25 against the refactored base module, which is
what establishes that the kernel seam changed nothing.
