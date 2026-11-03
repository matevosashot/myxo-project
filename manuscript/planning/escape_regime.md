# Escape regime: does the threshold rule P(p > p_crit) hold?

2026-10-05. Written for Discussion topic N8 ("The threshold as a trigger"), from Challenge 1 of the
literature review (Doering & Gadoua 1992; Hänggi, Talkner & Borkovec 1990). Scripts:
`claude_experiments/escape_regime/escape_numbers.py` and `escape_visits.py`.

## Bottom line

The pressure stays above p_crit for about half a microsecond per visit. Any plausible extrusion
takes milliseconds, so strictly the model is in the fast-fluctuation regime, not the slow one that
eq:p-layer-def assumes. If a cell responds to the pressure averaged over its extrusion time T,
the layering probability keeps the same Gaussian form with a smaller variance. For T set by the
bottom layer relaxing over one cell, the reduction is a fixed factor of about 5, independent of B.
The fit then absorbs it, and B rises from 90 to about 460 kPa·µm. The **shape** of the profile is
not guaranteed to survive, because averaging over T also spreads the fluctuations over about
1.7 µm, comparable to the 1.2 µm offset of the peak. That last point needs a calculation.

## 1. Correlation time of the pressure

At one point, eq:rho-autocorrelation-isotropic gives

    C(t) = σ_p² / (1 + |t|/τ_c),   τ_c = ℓ_n²/(2γ₀),   γ₀ = B/ξ₀ = 1.8×10⁴ µm²/min.

With ℓ_n = 0.5 µm this gives **τ_c = 6.9×10⁻⁶ min = 0.42 ms**. At the same parameters,
σ_p,∞ = 1.69 kPa·µm and z∞ = p_crit/σ_p,∞ = 20.9, as in §V.

## 2. How long the pressure stays above p_crit

C(t) has a cusp at t = 0, C ≈ σ²(1 − |t|/τ_c), so near a high level the pressure moves like a
random walk with diffusivity σ²/τ_c. A Gaussian tail overshoots the level by about σ/z, and the
walk crosses back over that distance in a time of order τ_c/z².

Numerical check on an Ornstein–Uhlenbeck process, which has the same cusp: the time spent above
u per visit is 0.55–0.65 τ/u² for u = 2–3. It does not depend on the time step (dt = 10⁻³ and
2×10⁻⁴ τ). A visit counts as ended once the path falls below u − 1/u, because a path with a
cusp recrosses a level many times in quick succession.

At z = 20.9 the pressure therefore stays above p_crit for **about 0.6 τ_c/z² ≈ 0.6 µs** per visit.

## 3. How long leaving the plane takes

| Estimate | T | T/τ_c |
| --- | --- | --- |
| Drag-limited lift: one cell moves up by h under the force p_crit·L, with drag ξ₀/ρ₀ per cell | 2.8 ms | 6.8 |
| The bottom layer relaxes over one cell area, T = A_cell/γ₀ | 10 ms | 24 |
| Growth of a visible layer in Copenhagen et al. (frames every 8–12 s) | minutes | > 10⁵ |

The last row times the growth of a layer, not its onset, so it bounds nothing. Either of the
first two is 10⁴ times longer than one excursion above p_crit. Both the slow-limit rule and
"p > p_crit as a trigger" assume the cell completes its move while the pressure is still above
threshold, and that fails by four orders of magnitude.

## 4. What changes if the cell responds to the pressure averaged over T

Average p over a window T. The result is still Gaussian, with variance f(x)σ_p², where

    f(x) = (2/x) [ (1 + 1/x) ln(1+x) − 1 ],   x = T/τ_c.

The slow 1/t decay of C(t) gives the logarithm. The layering probability keeps the form of
eq:p-layer with z_eff = z/√f, and eq:R-value becomes ln R = (z²/2f)·ν₀/(1+ν₀), provided the
**relative** modulation ν₀ is unchanged by the averaging.

| T | x = T/τ_c | f | z_eff | B refitted to the same R (= B/f) |
| --- | --- | --- | --- | --- |
| 2.8 ms (drag-limited lift) | 6.8 | 0.40 | 33 | 230 kPa·µm |
| 10 ms (relaxation over one cell) | 24 | 0.20 | 47 | 460 kPa·µm |
| 1 s | 2.4×10³ | 5.7×10⁻³ | 280 | 1.6×10⁴ kPa·µm |

Points to note:
- **The relaxation estimate is B-independent.** x = (A_cell/γ₀)/(ℓ_n²/2γ₀) = 2A_cell/ℓ_n² = 24,
  so the refit is self-consistent: changing B changes τ_c and T together. A refitted
  B ≈ 0.46 MPa·µm still sits inside the App. A.2 ladder, above the *E. coli* effective modulus
  (190 kPa·µm) and below the cell wall (10⁴ kPa·µm). It is about 4× Pelling's M. xanthus envelope
  value (125 kPa·µm).
- **The shape may not survive.** Over T the fluctuations spread by √(γ₀T) = √A_cell ≈ 1.7 µm,
  the same size as the 1.2 µm offset of the variance peak and the few-µm decay of the
  modulation. Time averaging therefore also averages in space, and ν₀ of the averaged pressure
  is likely smaller and smoother than ν₀ of the instantaneous one. Whether the profile still
  peaks off the core, on the bend side, has to be computed. The App. G kernels allow it: the
  anisotropic piece integrated over a time window.
- **The attempt rate P/τ_c of §V is a convention.** With a cusp, the rate of crossings of p_crit
  is not defined, because it diverges with resolution. It cancels in the ratio either way.

## 5. Options for the text

1. **Keep N8 as it is.** State the slow-limit assumption and the trigger interpretation, which is
   your note. It is honest, but a referee who computes τ_c (0.4 ms) can show that the assumption
   fails by 10⁴.
2. **Add the averaged-pressure version to N8.** If the cell responds to the pressure averaged
   over its extrusion time, the variance drops by a geometric factor, 2A_cell/ℓ_n², independent
   of B. The fit absorbs it, and B rises about fivefold. This turns the assumption into a
   robustness statement, but only for the amplitude.
3. **Compute the time-averaged variance profile** (analytic, App. G kernels with a time window of
   T = A_cell/γ₀). This would settle whether the shape survives and would back option 2 for the
   shape too. It is roughly one working session; the Lyapunov solver cannot do it directly,
   because it gives only equal-time covariances.
