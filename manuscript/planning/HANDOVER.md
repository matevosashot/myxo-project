# Handover — orienting yourself in this project

Written after a manuscript-planning session on 2026-09-15, for whoever picks this work up
next: a collaborator joining the paper, or a future assistant session starting cold. It is an
orientation map, not an audit. Read the coverage note at the end before trusting any claim
here about a part of `dean.tex` I did not open.

---

## What this project is

A theory project explaining an experimental observation: in *Myxococcus xanthus* monolayers,
new layers form preferentially near +1/2 topological defects (Copenhagen *et al.* 2020), and
the monolayer shows very large traction-force fluctuations (Han *et al.* 2024). Both papers
are in `manuscript/planning/`.

The hypothesis under test is that **microscopic fluctuations — specifically cell reversals —
drive layer formation**. The chain runs: microscopic Langevin dynamics of self-propelled
reversing rods → Dean's method → fluctuating active nematic equations *with the noise carried
through* → steady state around a +1/2 defect in the comoving frame → variance of the density
field → probability that the local free energy exceeds an extrusion threshold.

The research is finished. What remains is the manuscript, plus one missing figure.

---

## Where the physics actually lives

`dean.tex` is 2465 lines and is **a derivation log, not a draft**. It accumulated in working
order, so its section names are not a reliable guide. Two specific traps:

* `\section{Results}` at **line 1255** is an eight-line stub about active stress. It is not
  where the results are.
* The real results are at **895–1254** (the analytic variance) and **1269–1352** (the extrusion
  rate). The second of those — the step that actually connects the physics to the biology — sits
  immediately before `\appendix` and reads like an afterthought. It is not; it is the punchline.

The load-bearing chain, in order, with line numbers:

| Step | Lines | What happens |
|---|---|---|
| 1 | 109–173 | Microscopic Langevin equations. Reversals enter as dichotomous noise `k_i` at rate `frev`, Gaussian-approximated at 148–151. |
| 2 | 174–359 | Dean's method. The coarse-graining that keeps the noise. This is the methodological novelty. |
| 3 | 360–526 | Free-energy choice; integrating out the angle. |
| 4 | 681–747 | Noiseless dynamics; activity renormalises `a → a'` and `K → K'`; move to the comoving frame. |
| 5 | 748–818 | Linearise around the steady state, restore the noise. |
| 6 | 895–1254 | Set δQ = 0, solve for the density variance. Central result `eq:total-gauss` at ~1088. |
| 7 | 1269–1352 | Threshold → Gaussian tail → excursion duration → Poisson extrusion rate `eq:clump-rate`. |

Everything else — roughly 60% of the file — is numerical machinery and validation.

---

## Things that are not obvious and cost time to rediscover

**Reversals set two things, not one.** `zeta = xi0 v0^2 / (4 frev)` at line 223, and the noise
covariance carries `v0^2 / (frev rho0)` at line 311. So the same measured biological parameter
fixes the mean active stress *and* the fluctuation amplitude. They are not independent knobs.
This is the sharpest quantitative statement of the project's hypothesis and it is currently an
intermediate algebra step. The manuscript plan promotes it to the opening of §2A.

**The money plot does not exist.** Nothing in `dean.tex`, `figures/`, or `media/` compares the
predicted layer-formation probability against Fig. 3f of Copenhagen 2020, even though
`structure.md` describes the match as achieved. If you are looking for it, stop looking; it has
to be made.

**`ecrit` is fitted.** The extrusion threshold is calibrated against that same figure. Every
sentence about "prediction" has to be written with that in mind. What is genuinely free of
fitting is the *spatial structure* of the variance.

**Primed and unprimed parameters.** `a'` and `K'` are the activity-renormalised combinations
(`a' = a + Lambda xir`, `K' = K + zeta xir / (4 xi0)`), and it is `a'`, `K'` that set the defect
structure — `elld = sqrt(K'/-a')`. The bare `a`, `K` are what solve the renormalised relations
backwards. Lines 640–680 keep the primes visible deliberately. Mixing them silently gives a
defect of the wrong size.

**Two parameter sets exist in the literature and they disagree.** Copenhagen 2020 gives
`elld = 2.5 µm`; Han 2024 gives `elld = 0.8 µm`. This project follows Han. The set in `dean.tex`
§Parameter choice is explicitly described as mixed — some entries from one paper, some from the
other, some the project's own estimates. Do not assume any single source. Note also that job
scripts under `jobscripts/` may still carry an older parameter set than the one in the
`dean.tex` table; check before reusing a run.

**Two simplifications a referee will target.** The analytic line sets δQ = 0 ("assume rotational
diffusion is weak", line 896) — this has been checked numerically, and the check belongs in an
appendix. And the 4th-order Onsager-conjugate pair is dropped from both equations at ~line 700,
for numerical tractability. Both are defensible; neither is currently prominent.

**Two pleasant cancellations worth not losing in a rewrite.** In the anisotropic calculation,
piece (A) vanishes *exactly* by isotropy against the traceless `Q` (line ~1050), and the `e^x`
produced by the inner q-integral is cancelled by the external `e^{-x}` of the smoothing weight,
leaving `E_1` bare (line ~1065). Both make the final kernel much simpler than the intermediate
steps suggest.

**`elln` is physics, not a regulator.** The equal-point variance genuinely diverges; the cutoff
is the noise correlation length set by cell size, estimated as `l_rod/sqrt(12) ≈ sqrt(3) µm`. The
draft is careful about this and the manuscript should stay careful — it is the difference between
a controlled result and a fudge.

---

## Code layout, and one trap in it

| Path | What it is |
|---|---|
| `comoving_defect/` | The +1/2 defect steady-state solver in the comoving frame. Has a good README; it is the operating manual for appendix `sec-comoving-solver`. |
| `lyapunov_solver/` | **Current** equal-time fluctuation solver. Has a good README. |
| `lyapunov_solver_module*.m` (top level, five files) | **Superseded.** Per `lyapunov_solver/README.md` these drop the advection term `div(δQ w_ss)` entirely — a correctness bug, not a style difference. Use `lyapunov_solver/` instead. |
| `rho-only/` | Density-only calculations, matching §"Only density fluctuations". |
| `time_analysis/` | Solver timing study. Useful for justifying numerical choices; not paper content. |
| `jobscripts/` | Cluster runs. Check parameters before reuse. |
| `FRIDGE/` | Parked work. Assume nothing in it is current. |
| `claude_experiments/` | Sandbox for throwaway Wolfram scripts. |

Numerical results you can cite without re-running: the comoving frame velocities
`u_reservoir = 0.2842`, `u_sealed = 0.6980 µm/min` at L = 8 µm, and the seven-row verification
table at lines 2404–2424. Note the honest caveat attached to them — `u` is **not** converged in
box size (10.8% between L = 8 and 16 µm), and the two wall conditions bracket the true value
from opposite sides.

---

## Project conventions

* `CLAUDE.md` forbids editing Mathematica files (`.nb`, `.wls`, `.m`) unless explicitly asked —
  propose changes instead. Wolfram scripts may be created and run in `claude_experiments/`.
* `dean.before.tex` is a prior snapshot of the draft. `sealed_boundary.tex` and `section5.tex`
  are satellite files.
* Branch at time of writing: `dean-mathematica`.

---

## Where the planning landed

`structure.md` now carries Chapter Plan v2. The three decisions that shape everything:

1. **Main text is the analytic line, end to end. Appendices carry the numerics** — and their
   role is to show the main-text approximations are safe, not to stand as separate results.
2. **The claim is sufficiency, not uniqueness.** The author's own position is that the measured
   decay is simple and other models could fit it; what the work establishes is that the
   fluctuation hypothesis, carried through honestly, *suffices*. That caveat belongs in the
   Discussion in the paper's own voice.
3. **The `frev` prediction carries the falsifiability.** Because reversals set both the activity
   and the noise, and the threshold is already fixed, the model predicts the response to altered
   reversal frequency with no further parameters. A −1/2 defect test would have done the same job
   and was considered and declined as out of scope.

Blocking item: the Fig. 3f comparison figure.

---

## Coverage note — what I actually read

Read closely: lines 1–62, 527–760, 819–1120, 1255–1400, 2404–2465, plus the reversal-related
lines 143–151, 212, 223, 307–311 found by search, plus the four READMEs and `structure.md`.

**Read only as section headings:** the bulk of the Dean derivation (174–526), the real-space
cutoff and solver-comparison material (1121–1254), the whole OU / Lyapunov / finite-difference
appendix (1384–1952), the maximum-entropy closure and realizability analysis (1953–2313), and
the comoving discretisation (2349–2403).

So: statements above about the analytic chain and the extrusion argument rest on the text
itself. Statements about the numerical appendices rest on section titles and READMEs, and should
be verified before anyone relies on them. The claim that the closure is not load-bearing for the
analytic result is an inference from the fact that the main line sets δQ = 0 — it was put to the
author, who agreed the numerics can move to an appendix, but it was not verified line by line.
