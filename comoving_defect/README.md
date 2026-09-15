# Comoving-frame $+1/2$ defect solver

Solves the coupled $\rho$–$Q$ active-nematic system in the frame that moves with a $+1/2$
nematic defect, with the defect core pinned at the origin, and returns the fields **and** the
frame velocity $\mathbf{u}$.

$$
\partial_t \rho = \frac{B}{\xi_0}\nabla^2\rho + \frac{\zeta}{\xi_0}\partial_a\partial_b Q_{ab}
\qquad
\partial_t Q = \nabla\cdot\!\Big(Q\big(\mathbf{u} + \tfrac{B}{\xi_0}\nabla\rho\big)\Big)
             + \frac{4}{\xi_r}H_\text{active}
$$

with $H_\text{active} = -a'Q - b\,(Q\!:\!Q)\,Q + K'\nabla^2 Q$, $\;a' = a + \Lambda\xi_r$,
$\;K' = K + \zeta\xi_r/(4\xi_0)$, and $Q = \begin{pmatrix} Q_1 & Q_2\\ Q_2 & -Q_1\end{pmatrix}$.

These are eqs. `rho-comoving` / `Q-comoving` of `../dean.tex`. The write-up is
Appendix `sec-comoving-solver` there; this file is the operating manual.

---

## Files

| file | what it is |
|---|---|
| `comoving_steady_solver.m` | **front door** — `SolveActiveNematicComovingSteady`, the repo's standard module interface |
| `comoving_transient_solver.m` | **front door** — `SolveActiveNematicComovingTransient`, ditto with time stepping |
| `comoving_core.m` | the numerics — everything the two modules wrap |
| `evolve_run.wls` | production driver: transient + steady polish, writes `evolve_<bc>.mx` |
| `evolve_plots.wls` | figures from a finished `.mx`; no solving |
| `run_case.wls` | verification driver (Jacobian, passive limit, lab-frame drift, single steady solves), appends to `results.csv` |
| `results.csv` | append-only log of every verification run |

Two ways in. The **two modules** take the same `modelParams` / `solverParams` rule lists as
`iterative_solver_module.m` and `transient_solver_module.m` and return the same
`<|"rho","Q1","Q2","mesh"|>` shape, so they drop straight into `visualizeSteadyState` and the
Lyapunov modules — that is what most callers want. The **`.wls` drivers** below go through
`comoving_core.m` directly and are what produced the figures in the paper.

Generated and **not** committed (see `.gitignore`): `evolve_<bc>.mx`, `out_<bc>/`.
The figures the paper actually includes are copied to `../figures/comoving/`.

## Requirements

Mathematica 14.1. On this cluster that is

```bash
WS=/usr/local/math/math1410/Executables/wolframscript
```

not `/usr/local/bin/wolframscript`, which is 13.3.1. Rasterising plots needs the front end,
which crashes on a headless node unless `DISPLAY` is cleared — `evolve_plots.wls` clears it
itself, so nothing extra is needed, but if you write your own plotting script, `unset DISPLAY`
first (as `../run.sh` does).

Every script locates the module relative to itself, so the directory can be moved as a whole
without editing anything.

---

# Part 1 — Using it

## The module interface

This is the normal way in, and it mirrors `iterative_solver_module.m` /
`transient_solver_module.m` exactly.

```mathematica
Import[FileNameJoin[{repoRoot, "comoving_defect", "comoving_steady_solver.m"}]];

modelParams = {B -> 10000., a -> -0.172, b -> 0.1, L -> 0.05368,
               \[Zeta] -> 10., \[Xi]0 -> 5., \[Xi]r -> 0.1,
               \[Rho]0 -> 0.3, \[CapitalLambda] -> 0.1};

solverParams = {\[Delta]mesh -> 0.05, hMax -> 1., box -> 8.,
                maxSteps -> 25, relativeDiff -> 1.*^-8};

res = SolveActiveNematicComovingSteady[modelParams, solverParams, "Dirichlet"];

res["u"]      (* {0.2842160764, 1.1e-5}  <- the frame velocity, um/min *)
res["rho"]    (* InterpolatingFunction of (x,y), as from the FEM solvers *)
```

Same parameter keys as the FEM modules, `L` being the Frank constant.
`\[Delta]mesh`, `hMax` and `box` mean the same thing here as there — the grid is graded by
the same $1/|S'(r)|$ criterion the FEM `targetCellSize` uses, and the point count is derived
from them rather than dialled in.

The transient module adds `timeInterval` and marches in blocks of exactly that length,
stopping when the relative change of all three fields across a block falls below
`relativeDiff`, or after `maxSteps` blocks:

```mathematica
Import[FileNameJoin[{repoRoot, "comoving_defect", "comoving_transient_solver.m"}]];

tr = SolveActiveNematicComovingTransient[modelParams,
       Join[solverParams, {timeInterval -> 2., returnTransient -> True}],
       "Dirichlet"];

tr["solver"]["history"]          (* per-block t, u, residual, substeps *)
tr["transient"]["u"][3.5]        (* u_x at t = 3.5 *)
tr["transient"]["rho"][3.5, {1., 2.}, {0., 1.}]   (* vectorised in space *)
```

### Optional `solverParams` keys

| key | default | meaning |
|---|---|---|
| `refine` | `1` | scales every cell together; grading preserved |
| `rhoGauge` | `"Point"` | sealed-wall gauge: `"Point"` / `"Mean"` / `"Replace"` |
| `order` | `4` | finite-difference order |
| `verbose` | `True` | per-step printing |
| `zetaSteps` | `5` | $\zeta$-continuation steps (steady only) |
| `stepTolerance` | `1.*^-10` | Newton tolerance at the target $\zeta$ (steady only) |
| `timeInterval` | *required* | block length (transient only) |
| `returnTransient` | `False` | keep the whole trajectory (transient only) |
| `dt0`, `dtMax` | `1.*^-5`, `timeInterval` | substepping (transient only) |
| `maxSubsteps` | `2000` | substep cap per block (transient only) |

`relativeDiff` means what it means in the FEM modules — the change across an outer step.
For the steady solver that is the Newton tolerance during $\zeta$-continuation; the final
solve at the target $\zeta$ uses `stepTolerance` instead, because a loose `relativeDiff`
would otherwise stop Newton one or two iterations short of the accuracy the published $u$
needs, at no real saving.

### What comes back

`"rho"`, `"Q1"`, `"Q2"`, `"mesh"` — the standard four — plus `"u" -> {ux, uy}` and
`"solver"` (convergence flag, residual, iteration count, pinned core, per-block history,
and the raw state vector under `"state"`). The fields are `ListInterpolation`s at
`InterpolationOrder -> 3` over the tensor grid, so they take vector arguments and survive
`Derivative[2,0]` — which is what the Lyapunov modules need of them.

`"mesh"` is **not** an `ElementMesh` — there is no FEM here. It is an association carrying
the grid (`"gx"`, `"gy"`, `"hMin"`, `"weights"`, ...) plus `"Coordinates"` in the same
$n\times2$ shape an `ElementMesh` would give, since that is the only mesh property anything
in the repo actually reads.

Two semantic warnings, both repeated in the `::usage` strings:

- **These are not the lab-frame steady states.** The core is pinned and the problem is posed
  in a frame translating at $\mathbf u$; results are not comparable with
  `SolveActiveNematicSteady` except in the passive limit $\zeta \to 0$, where $u \to 0$.
- **`"Neumann"` is not `iterative_solver_module.m`'s `"Neumann"`.** That one pins $\rho$
  along a whole ray. Here it is the sealed wall plus a single-node gauge — see
  [Boundary conditions](#boundary-conditions).

## Quick start (the `.wls` drivers)

```bash
cd comoving_defect

# 1. solve (a few minutes per wall at refine 1) -- prints every accepted step
$WS -file evolve_run.wls 8 1 both

# 2. figures from the saved state (no solving)
$WS -file evolve_plots.wls both
```

Arguments are `L refine bcType`, defaulting to `8 1 both`:

- `L` — half-width of the square box, in μm. The box is $[-L,L]^2$.
- `refine` — mesh refinement factor. `1` gives $75^2$ nodes; `2` gives $149^2$ and costs
  roughly $8\times$. Fractional values are allowed and are useful for smoke tests
  (`0.35` runs in seconds but is only good to a percent or so).
- `bcType` — `Dirichlet`, `Neumann`, or `both`.

`evolve_run.wls` prints one line per accepted step and checkpoints the `.mx` every 25 steps,
so it is safe to watch, and the output file is usable before the run finishes. A typical line:

```
   k =  147  t = 8.123e0      dt = 5.000e-1     ux = 0.2842161    uy = 1.097e-5
              |dX/dt| = 3.1e-8   nwt = 2   183.4s
```

`evolve_plots.wls <bc> [plotBox] [nFrames]` writes into `out_<bc>/`:

| output | content |
|---|---|
| `u_vs_t.png` / `.pdf` | $u_x(t)$ and $\lvert u_y(t)\rvert$, log in $t$ |
| `rho_slice`, `Q_slice` | $\rho(x,0)$, $Q_1(x,0)$, $Q_2(x,0)$ at log-spaced times |
| `rho_spacetime.png` | $\rho(x,0)$ as $x$ vs. $\log t$ |
| `slices.pdf` | the two-panel (a)/(b) figure used in the paper |
| `frames/frame_###.png` | defect shape over time — strength $Q_1^2+Q_2^2$ plus the headless director |
| `defect_steady.png` | the converged frame, under a stable name |
| `history.csv` | `t, dt, ux, uy, rate, detG, newton` per step |

## Verification runs

```bash
$WS -file run_case.wls jac                      #  ~2 s
$WS -file run_case.wls passive  6 1             #  zeta = 0: u must go to 0
$WS -file run_case.wls active   8 1 Dirichlet   #  ~9 s
$WS -file run_case.wls active   8 1 Neumann     #  ~110 s
$WS -file run_case.wls drift    8 1 Dirichlet   #  lab-frame cross-check
$WS -file run_case.wls sinh     8 51 Dirichlet  #  sinh grid instead of Pade-graded
```

A fifth argument selects the sealed-wall density gauge (`Point`, `Mean`, `Replace`; see
below). Each run appends one line to `results.csv`, so an interrupted sweep loses nothing.

## Calling `comoving_core.m` directly

Below the two modules. Use this when you want the raw state vector, a custom continuation,
or control the wrappers do not expose.

```mathematica
Get["comoving_core.m"];

(* 1. parameters *)
mp = <|"B" -> 10000., "a" -> -0.172, "b" -> 0.1, "K" -> 0.05368,
       "zeta" -> 10., "xi0" -> 5., "xir" -> 0.1, "rho0" -> 0.3,
       "Lambda" -> 0.1|>;
dp = cdDerived[mp];                       (* -> Bx, lam, kap, Aq, bet, ap, Kp, ld, S0 *)
dp = Append[dp, "rhoGauge" -> "Point"];   (* NB: cdDerived does not carry this across *)
pf = cdPadeFuns[dp];                      (* Pade +1/2 profile, used as seed and far field *)

(* 2. grid and operators *)
gg = cdGridPade[8., dp, 0.05, 1.0, 1];    (* L, dp, dmesh, hMax, refine *)
op = cdOps[gg, gg, 4];                    (* 4th-order difference operators *)
bc = cdBC[op, dp, pf, "Dirichlet"];

(* 3. initial state *)
X0 = cdInitial[op, pf];                            (* rho = 1, Pade core *)
X0 = cdInitialRho[X0, op, dp, bc, "Dirichlet"];    (* quasi-static rho; REQUIRED for Neumann *)
ex0 = ConstantArray[0., 2];                        (* {ux, uy}; 3 for Neumann, see below *)

(* 4. solve *)
ev = cdEvolve[X0, ex0, op, dp, bc, "Dirichlet", "Tolerance" -> 1.*^-7];
st = cdSteady[ev["X"], ev["ex"], op, dp, bc, "Dirichlet"];

st["ex"][[1]]        (* u_x  *)
st["X"]              (* Join[rho, Q1, Q2], length 3n *)
```

### State layout

```
X  = Join[rho, Q1, Q2]        length 3n,  n = nx*ny
ex = {ux, uy}                 Dirichlet, or Neumann under the "Replace" gauge
   = {ux, uy, mu}             Neumann under "Point" or "Mean"
```

Node ordering is `k = i + (j-1) nx` — **x fastest** — matching the rest of the repo.
The origin is node `op["k0"]`. Field blocks come back as matrices from
`cdFieldAt[X, "rho"|"Q1"|"Q2", op]` (rows = $y$) and as interpolants from `cdInterp`.

### Two gotchas

1. `cdDerived` returns a fresh association and **drops** `"rhoGauge"`. Append it afterwards
   or the default `"Point"` is used silently.
2. The length of `ex` must match the boundary condition and gauge:
   `ne = If[bcType === "Neumann" && rhoGauge =!= "Replace", 3, 2]`.

### Functions

**Setup** — `cdDerived[mp]`, `cdPadeFuns[dp]`, `cdGridPade[L, dp, dmesh, hMax, refine]`,
`cdGrid[L, n, s]` (uniform/sinh alternative), `cdOps[gx, gy, order]`,
`cdBC[op, dp, pf, bcType]`, `cdInitial[op, pf]`, `cdInitialRho[X, op, dp, bc, bcType]`.

**Solve** — `cdEvolve[...]` (implicit Euler; options `dt0`, `dtMax`, `dtGrow`, `MaxSteps`,
`Tolerance`, `StopTime`, `PrintEvery`, `SnapshotEvery`, `Checkpoint`, `CheckpointEvery`;
`StopTime` ends the run at exactly that $t$ and is reported separately from `Tolerance`
convergence via `"hitStopTime"`, with `"dtNext"` carrying the controller's step out so a
caller marching in blocks need not restart the ramp),
`cdSteady[...]` (bordered Newton; `MaxIterations`, `StepTolerance`),
`cdNewton[...]` (the shared inner solver), `cdRJ`/`cdRes` (Jacobian + residual, residual only).

**Diagnostics** — `cdRates` (residual split by field), `cdCoreG` (the $2\times2$ core gradient
matrix), `cdCoreEstimates` (the closed-form $\mathbf{u}$, as a check), `cdSymmetry` (mirror
symmetry, unimposed), `cdShiftRhoMean` (restore unit mean density),
`cdRadialBVP[ld, S0, R, m]` (independent 1D radial profile), `cdFmt` (compact single-line
numbers for progress printing).

### Parameters

`evolve_run.wls` and `run_case.wls` both solve the same constraint block for
$\xi_r, a, K$ from the measured quantities:

| fixed | | derived | |
|---|---|---|---|
| $\ell_d$ | 0.8 μm | $\xi_r$ | 0.1 |
| $S_0$ | 0.9 | $a$ | $-0.172$ |
| $\xi_0$ | 5 | $K$ | 0.05368 |
| $\zeta$ | 10 | $a'$ | $-0.162$ |
| $B$ | $10^4$ | $K'$ | 0.10368 |
| $b$ | 0.1 | $B/\xi_0$ | 2000 |
| $\gamma = 8b/\xi_r$ | 8 | $4K'/\xi_r$ | 4.147 |
| $\Lambda$ | 0.1 | $\alpha = \zeta\xi_r S_0/(4\xi_0 K')$ | 0.434 |

These match `sec-parameter-choice` of `../dean.tex`. The `.wls` drivers carry this block
verbatim, so **changing them there means editing both drivers**. The two `.m` modules take
the parameters as an argument instead, which is the better route if you are sweeping.

### Results at these parameters

$L = 8$ μm, refine 1 ($75^2$ nodes, $h_\text{min} = 0.075$, $\ell_d/h_\text{min} = 10.7$):

```
u_reservoir (Dirichlet) = 0.2842160765 um/min
u_sealed    (Neumann)   = 0.6980352321 um/min
```

with $u_y \sim 10^{-5}$ (not imposed), core exactly $0$, and sealed-wall mass conserved to
$2\times10^{-16}$.

---

# Part 2 — How it works

## The problem with the obvious approach

In the bulk the steady comoving equation is translation invariant: if $(\rho, Q)$ solves it,
so does any rigid translate. The discretised Jacobian therefore carries a **near-null mode**,
the translation $\partial_x(\rho, Q)$, broken only at $O(1/L)$ by the distant wall. Two
consequences, both of which look like bugs until you see the cause:

- integrated in the lab frame the defect simply drifts away — there is no stationary state to
  converge to;
- a plain steady Newton solve is ill-conditioned, because it is being asked to invert a matrix
  that is nearly singular in exactly the direction nobody cares about.

One can also get $\mathbf{u}$ in closed form. Differentiating the pinning condition
$Q_i(0,0,t) = 0$ in time, with $Q(0) = 0$, gives

$$
\mathbf{u} = -\frac{B}{\xi_0}\nabla\rho(0)
             - \frac{4K'}{\xi_r}G^{-1}\!\cdot\!\big(\nabla^2Q_1(0), \nabla^2Q_2(0)\big),
\qquad
G = \begin{pmatrix}\partial_x Q_1 & \partial_y Q_1\\ \partial_x Q_2 & \partial_y Q_2\end{pmatrix}_{\!0}
$$

but this is a bad way to *compute* it. With $G \approx c\,\mathbb{1}$, $c = S'(0) \approx 0.21$,
the prefactor $4K'/(\xi_r c) \approx 190$ multiplies a nodal second derivative that is itself
only a few percent of $S_0/\ell_d^2$: four clean digits in $\nabla^2Q_1(0)$ would be needed for
one percent in $u$. It is used here **only as a check** (`cdCoreEstimates`), where it agrees
with the converged answer to four digits. Note that the first term alone — the estimate quoted
in the body of `dean.tex` — is off by $3.75\times$ (reservoir) and $2.12\times$ (sealed).

## $\mathbf{u}$ as a Lagrange multiplier

The design decision. Treat the semi-discrete system as a **differential–algebraic** one and let
$\mathbf{u}$ be the multiplier conjugate to the pinning constraints. Unknowns are
$\rho, Q_1, Q_2$ at every node *plus* the two scalars $u_x, u_y$; equations are the discretised
PDE at every node *plus*

$$Q_1(0,0) = 0, \qquad Q_2(0,0) = 0 .$$

Two extra unknowns, two extra equations. The border replaces the near-null direction, so the
bordered matrix is well conditioned; the constraint holds **exactly at every step** — no drift,
no feedback gain to tune, no second derivative at the core ever needed; and $\mathbf{u}$ falls
out as a component of the solution vector.

The same bordered matrix serves both phases. An implicit Euler step solves
$(\mathbb{1}_d/\Delta t - J)\delta = \dots$ bordered; the steady solve is the same system with
$1/\Delta t \to 0$. Both are run, in that order, and for a reason: time stepping tells you
whether a fixed point exists *and is dynamically selected* — Newton alone converges happily to
unstable states and gives no diagnosis when there is none — and the steady solve then polishes
it to the round-off floor in a few iterations.

## Grid

A tensor-product grid graded by the criterion the FEM meshes in this repo already use
(`../iterative_solver_module.m:126`, `../transient_solver_module.m:182`): the cell size tracks
$1/|S'(r)|$, capped at $h_\text{max}$, so the mesh is fine where the defect amplitude turns
over and coarse in the far field where $S \to S_0$ and $S' \to 0$.

`targetCellSize` is a 2D *area* criterion, so the 1D version is built by counting cells:
$\xi(x) = \int_0^x \mathrm{d}x'/h(x')$ is the cell index, so sampling $x$ at uniform $\xi$
reproduces spacing $h(x)$ by construction. The point count is then **derived** from
`dmesh`/`hMax` rather than dialled in, and `refine` scales every cell together, so grading is
preserved under refinement and convergence orders stay meaningful. The grid is mirrored about
zero with an odd point count — the origin must be a node, since the pinning constraint lives
there.

Difference operators come from `NDSolve\`FiniteDifferenceDerivative` at 4th order (which
handles non-uniform grids with no change) and are lifted to 2D by Kronecker products. Boundary
nodes stay in the unknown vector and their rows are **overwritten** with the boundary
condition, which handles the inhomogeneous Padé far field correctly.

## Boundary conditions

$Q$ is held at the Padé $+1/2$ profile on the whole boundary, both cases. For $\rho$:

- **`"Dirichlet"`** — reservoir wall, $\rho = 1$.
- **`"Neumann"`** — sealed wall, $\hat n\cdot(B\nabla\rho + \zeta\nabla\cdot Q) = 0$.

Two things about the sealed wall are easy to get wrong.

**The normal must be signed.** At the corners where $\hat n \propto (+1,-1)$ or $(-1,+1)$ the
condition is $(\partial_x - \partial_y)\rho = 0$, not $(\partial_x + \partial_y)\rho = 0$. With
unsigned $0/1$ masks those two corner rows become exact linear combinations of their
neighbouring edge rows, and the bordered matrix loses rank by **exactly 2** at every grid size
and every stencil order. The symptom is that Neumann never takes a single step; the cause is
two rows out of 16875.

**The sealed problem is singular through the constant null vector** and needs one extra
condition. Three gauges are implemented (`dp["rhoGauge"]`):

| gauge | condition | cost |
|---|---|---|
| `"Mean"` | $\sum_k w_k\rho_k = \lvert\Omega\rvert$ — physical | dense row across the $\rho$ block |
| `"Point"` | $\rho(0,0) = 1$ — extra row + multiplier $\mu$ (**default**) | sparse row, dense $\mu$ column |
| `"Replace"` | $\rho = 1$ at one near-wall node, *replacing* that PDE row | no extra row, no multiplier |

They differ only by an additive constant in $\rho$, and $u$ is blind to that: the $Q$ equation
sees $\rho$ only through $\nabla\rho$ and the $\rho$ equation only through $\nabla^2\rho$.
Measured, the gauges agree on $u$ to ten significant figures. `cdShiftRhoMean` restores the
unit-mean normalisation afterwards.

This is *not* the ray pin of `../iterative_solver_module.m`, which fixes $\rho$ along a whole
ray — many nodes, over-constrained, hence its local distortion. One node with one multiplier is
exactly rank-1.

## Initial condition

$Q$ starts from the Padé profile. For $\rho$, uniform $\rho = 1$ is fine under the reservoir
wall — it satisfies it. Under the **sealed** wall it does not: $\zeta\,\hat n\cdot(\nabla\cdot Q)
\neq 0$ there, so the first implicit step has to resolve a violent boundary layer, and at
$B/\xi_0 = 2000$ the resulting $\nabla\rho$ feeds straight into the advection velocity
$\mathbf{w} = \mathbf{u} + (B/\xi_0)\nabla\rho$ and Newton diverges. `cdInitialRho` fixes this
with a quasi-static seed: hold $Q$ frozen and solve the $\rho$ equation alone, which is a single
**linear** solve.

## Residual and Jacobian

With $\lambda = \zeta/\xi_0$, $\kappa = 4K'/\xi_r$, $A = 4a'/\xi_r$, $\beta = 8b/\xi_r$ and
$\mathbf{w} = \mathbf{u} + (B/\xi_0)\nabla\rho$:

$$
\partial_t\rho = \frac{B}{\xi_0}\nabla^2\rho
  + \lambda\big[(\partial_x^2 - \partial_y^2)Q_1 + 2\partial_x\partial_y Q_2\big] + \mu
$$
$$
\partial_t Q_i = \nabla\cdot(Q_i\mathbf{w}) + \kappa\nabla^2 Q_i - A\,Q_i
  - \beta(Q_1^2 + Q_2^2)Q_i
$$

The Jacobian is assembled analytically and stays sparse. The only $Q_1 \leftrightarrow Q_2$
coupling is the cubic, since $\mathbf{w}$ depends on $\rho$ and $\mathbf{u}$ but never on $Q$.

The **full bordered matrix is factorised directly**. Eliminating through the unbordered block
via a Schur complement would defeat the entire purpose — that block is near-singular by
construction.

## Four things that are not optional

Each of these cost real debugging time, and none is visible from the equations.

**Row equilibration.** The $\rho$ rows carry $B/\xi_0 = 2000$ and the $Q$ rows
$4K'/\xi_r = 4.15$ — a factor 480 before $1/h^2$ multiplies both. At $h = 0.027$ that puts
matrix entries at $2.7\times10^6$ against $2.9\times10^3$ and the factorisation stops being
useful: measured, the passive solve simply *did not move*, with $\lvert R\rvert$ pinned at its
initial 0.174. Scaling each row by its own largest entry is an exact row scaling — the Newton
step is unchanged in exact arithmetic — and it is what makes fine grids solvable at all. This
did not matter at $B = 10^3$; it is essential at $10^4$.

**Continuation in $\zeta$.** At $\alpha = 0.43$ a cold Newton from the Padé seed does not
converge (measured: 40 wasted iterations, then a 30-step evolve, ~200 linear solves, 57 s).
Ramping $\zeta$ from 0 to its target in ~5 steps keeps every Newton a small perturbation of the
previous solution: ~3 iterations per step, ~20 solves, 14 s. The grid is built once at the
*target* $\zeta$ and held fixed.

**Constraint row scaling.** The pinning rows carry no $1/\Delta t$ while the differential rows
do, so the bordered Schur complement behaves like $\Delta t\,G$ and the matrix degenerates as
$\Delta t \to 0$ — the usual index-2 conditioning blow-up. Measured, the smallest singular value
is proportional to $\Delta t$ (2.2e-5 at $\Delta t = 10^{-4}$, 1.2e-1 at $\Delta t = 1$).
`cdBScale` scales the constraint rows by $1/\Delta t$, again exactly, making the conditioning
$\Delta t$-independent. Without it, cutting $\Delta t$ in response to a failed step makes matters
*worse*, which is a thoroughly misleading symptom.

**Bordering by index, not `ArrayFlatten`.** With ragged blocks ($3n\times3n$ beside
$3n\times n_e$) `ArrayFlatten` leaves its sparse path and materialises the result densely —
9.3 GB and 22.5 s of a 25 s assembly at $n = 81^2$. `cdBorder` works from the packed
`NonzeroPositions` arrays instead: assembly went 26.7 s → 0.078 s, a factor 340.

## Verification

All at the parameters above. `run_case.wls` reproduces every row.

| test | result |
|---|---|
| analytic Jacobian vs. central differences | $1.2\times10^{-10}$ relative, both walls |
| passive limit $\zeta = 0$: $u \to 0$ | order 4.3 in $h$ |
| passive profile vs. independent radial BVP | order 4.0, $1.3\times10^{-6}$ |
| pinning constraint at every step | exact |
| mass and multiplier, sealed wall | conserved to $2\times10^{-16}$, $\mu \to 6.6\times10^{-9}$ |
| mirror symmetry, $u_y$ | $\sim10^{-5}$, unimposed |
| closed-form core condition vs. multiplier | agrees to 4 digits |
| **lab-frame drift vs. bordered $\mathbf{u}$** | $6\times10^{-5}$ relative |

The last row is the decisive one. Take the converged comoving field, integrate it in the **lab**
frame with $\mathbf{u} \equiv 0$ and no border at all, track the core by root-finding on the
interpolated field, and compare $\dot x_\text{core}$ against the multiplier. Nothing in that
test uses the border, the constraint or the multiplier, so agreement at $6\times10^{-5}$ —
itself consistent with the first-order time stepping — is two unrelated numerical routes to the
same number.

The mirror symmetry is worth a note: nothing in the discretisation imposes $\rho, Q_1$ even and
$Q_2$ odd in $y$, so the measured asymmetry is a free check on every advection sign.

Two false alarms, recorded so they are not rediscovered. The passive profile error first came
out at $1.7\times10^{-2}$ with a nonsensical convergence order — an artifact of comparing
through a cubic interpolant laid across the corner that $S = |Q|$ has at $x = 0$. Comparing at
grid nodes gives $1.8\times10^{-5}$ and a clean 4th order. Separately, the Jacobian test broke
after row equilibration was added, at 0.66 relative error, because the equilibrated residual's
scale factors depend on the state — so the finite difference was not of a consistent function.
It compares the *unscaled* Jacobian against differences of `cdRes` now.

---

# Part 3 — Caveats

**$u$ is not converged in box size.** This is the real limitation. Doubling $L$ from 8 to 16 μm
moves $u$ by 10.8%. The two wall conditions bracket the answer from opposite sides with a gap
that closes like $1/L$, but the common extrapolate still drifts logarithmically, which is
suggestive of an unscreened far field and was never chased down (an Oseen-type hypothesis was
written up but not tested). The grid error, by contrast, is negligible — roughly 5th order in
$h$, five significant figures by the production resolution. **The box, not the mesh, is what
limits the quoted values.** Quoting $u$ for the infinite plane would need either a substantially
larger domain or an asymptotic far-field condition.

**The boundary conditions are imposed in the comoving frame,** so the far field is dragged along
with the defect. That is what defect stabilisation means here, but it is a physical modelling
choice and $u$ depends on it. The reservoir and sealed walls differ by a factor 2.5 in $u$ —
that spread is the honest uncertainty from the wall, not numerical error.

**Neumann costs ~7× more per solve** than Dirichlet at equal size — 110 s against 9 s at
refine 1 — and we do not know why. Three plausible explanations were measured and all three
ruled out: the dense mean-density row (the `"Point"` gauge removes it, no speedup), the dense
$\mu$ column (`"Replace"` removes that too, no speedup), and iteration count (identical, 59
either way, with nnz within 2%). It is stated as unresolved rather than explained away.

**The sealed wall is markedly more grid-sensitive** than the reservoir. Smoke-test resolutions
that give Dirichlet $u$ to 2% give Neumann $u$ to only ~30%. Do not trust a Neumann number from
`refine < 1`.

**The $\mathbf{u}\cdot\nabla\rho$ term is dropped by default.** `dean.tex` drops it as
$O((\nabla\rho)^2)$. It is implemented behind the switch `dp["ugr"] = True` so its size can be
measured rather than assumed — but that measurement was never actually made.

**Round-off floor.** With $B/\xi_0 = 2000$ and $1/h^2 \sim 180$, matrix entries reach $\sim4\times10^5$
and $\lvert R\rvert$ bottoms out near $10^{-10}$. A `StepTolerance` below that silently burns
every remaining iteration; `cdNewton` has a stagnation exit as the safety net, but do not set
`StepTolerance` below `1e-10` and expect it to mean anything.

**The parameter block is duplicated** in `evolve_run.wls` and `run_case.wls`, verbatim in both.
Changing parameters means changing both. It was left that way deliberately — each driver is then
self-contained and reproducible from its own source — but it is a trap if you forget. The two
`.m` modules take parameters as an argument and do not have this problem.

**The first transient block is much the most expensive.** `dt0` defaults to $10^{-5}$ and the
step grows by 1.6× only when Newton is comfortable, so reaching a block length of 2 takes
~175 substeps the first time and a handful thereafter. That ramp is not wasted — it is what
lets the run start from a flat density at $B/\xi_0 = 2000$ without diverging — but if you
are restarting from a state that is already close, raise `dt0`.

**The transient interpolates between block ends, it does not re-integrate.** The FEM module's
`"transient"` closures come from `NDSolveValue`, which is continuous in $t$ within a block.
Here `f[t,x,y]` is a `ListInterpolation` through the stored block-end states, so resolution
in $t$ is exactly `timeInterval`. Use a smaller `timeInterval` if you need a finer trajectory
— the cost is almost flat, since the substeps happen either way.

**`refine` changes the node count, not just the spacing.** Convergence studies are meaningful
because grading is preserved, but `refine 2` is ~$8\times$ the cost of `refine 1`, not $4\times$:
the node count quadruples and the factorisation is superlinear on top of that.
