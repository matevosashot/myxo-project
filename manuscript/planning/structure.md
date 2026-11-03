# Planning the manuscript

## High-level content

The paper "Topological defects promote layer formation in Myxococcus xanthus colonies" (./manuscript/planning/s41567-020-01056-4.pdf)
studies the monolayers of mycobacteria colonies. There is a experimental measurement on flow fields, orientation etc., An active nematic model is used to fit the data. 
Paper measures the layer formation events, they predominantly appear in regions with +1/2 defects.

With this paper we want to ultimatley explain this observation.

there is another key paper "Local polar order controls mechanical stress and triggers layer formation in Myxococcus xanthus colonies" (./manuscript/planning/s41467-024-55806-6.pdf), where more things are mesured, especially the fluctuations of the traction force. This reveals huge fluctuations in the monolayer.

We start with a microscopic model. We coarse grain and obtain an active nematic equations, including the noise term (this is a novelty). Then we solve it for the defect. And eventually we predict the probability of layer formation aroud the defect. It matches quite well with the figure 3f of the first paper above.

Preparation of the figures on the way.

We have two main results:

- coarse graining of the microscopic system and obtain fluctuating active nematic equations
- solving the equations (theoretically and numerically) for the system and obtaining consistent results with the experiments.

## Structure of the manuscript

*Chapter Plan v3 — v2 established 2026-09-15; revised 2026-09-21 once `dean.tex`'s
Extrusion-rate section converged on a single mechanism (classical nucleation theory,
barrier maximum descending to one cell) and produced the comparison-with-experiment
figure. Decisions are recorded as decided; anything still open is marked **TODO**.*

### Framing decisions

| Decision    | Choice                                                                                                                                                                        |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Venue shape | Long-form, **main text + appendices** (PRX / PRE shaped). Converts to Nat. Commun. Methods + Supplementary Notes if that becomes the target; the content split is identical. |
| Spine       | Phenomenon frames the paper; the equations are co-billed, not machinery.                                                                                                      |
| Claim type  | **Sufficiency, not uniqueness.**                                                                                                                                              |
| Hero figure | Two panels: mechanism, then experimental match. **Both now exist:** `figures/varrho.png`, `figures/final.png`.                                                              |
| Main text   | The analytic line, end to end — now including the nucleation-theory argument of §3.                                                                                          |
| Appendices  | The numerics and the parameter-sensitivity detail — *evidence that the main-text approximations are safe*, not independent results.                                          |

### The claim, stated precisely

We propose a mechanism to explain layer formation in these bacterial colonies:
fluctuations originating at the microscopic scale — in particular from cell reversals —
drive it. Starting from a microscopic model of self-propelled reversing rods,
coarse-grained with the noise retained, we carry this mechanism end to end and show
that it accounts for the observed pattern of layer formation around a +1/2 defect.

Uniqueness is **not** claimed. The measured decay is simple and other models could fit
it. This is written into the Discussion, not left implicit — see *What this does and does
not establish* below.

What is genuinely parameter-free and what is not — **settled, 2026-09-21** (the earlier
two-model fork is resolved; see *Resolved since v2* below):

- **Predicted, no freedom:** the spatial *shape* of the variance modulation around the
  defect — the decay over ~5 µm, the peak at r≈1.5 µm, and the core suppression by
  symmetry — all from the `E_1`-filtered `d_i d_j Q_ij` of `eq:aniso-kernel`. This is what
  a competing mechanism would also have to reproduce.
- **Fitted, and it is two parameters, not one:** `B` and `elln` both enter the single
  group `ln(ratio) = (p_crit^2/2)(2 pi rho_0 elln^2 / (B zeta)) nu_0/(1+nu_0)`, and `elln`
  is at least as consequential as `B` — `ln(ratio)` scales roughly as `elln^3`. dean.tex's
  own text calls this "a one-parameter fit (`B`, or equivalently the group)" — **worth
  prompting for a fix**, since project memory (built independently while deriving this
  section) is explicit that both are fitted. State it as two in the paper.
- **Result:** at `B = 2e5 Pa um`, `elln = 0.5 um`, the predicted enhancement is
  `ratio = 228`, against Copenhagen's first-bin value `~230` (error bar 40–500). This is
  the number that was `TODO` in v2 — it now exists, and it lands inside the error bar.

### Working title

> *Fluctuating active nematics account for layer formation in* Myxococcus xanthus *colonies*

From the author's own contribution sentence ("Fluctuating Active Nematic equations explain
Layer formation in Myxobacteria"), with *explain* softened to *account for* so the title does
not assert the uniqueness the paper deliberately does not claim.

