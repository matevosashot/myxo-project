# Dean–Kawasaki equation: what the review covers



---

## 2. Dean's equation is exact, with no coarse-graining (§2.2, comment i) ★★★

- The derivation makes **no approximation**. It is equivalent to the N coupled Langevin
equations and keeps "particle entity": each particle sits at one point at a time
(proved by Bothe et al. 2023, J. Phys. A 56, 175002).
- The review contrasts this with fluctuating hydrodynamics (§3.1): the DK equation "does
not rely on any coarse-graining approximation."
- Dean (1996) himself: the equation is "directly defined from an albeit simple microscopic
model" rather than "a phenomenological model … based around a known coarse grained
free energy."
- **Use in our paper:** supports the placeholder bullet "nothing is added at the
continuum level." Our only approximation up to eq:final-langevin is the Gaussian
white-noise limit of the reversal noise (see topic 11 for its extensions).

## 3. The density is a distribution, not a field (§2.2 ii, §4) ★★★

- ρ(x,t) = Σ δ(r_α(t) − x) is a sum of delta functions. It has a physical meaning only
after an ensemble average or a local spatial average. The equation lives in a space of
distributions.
- Mathematical regularization: replace δ by a Gaussian kernel w_σ of width σ (Cornalba,
Shardlow, Zimmer 2019, SIAM J. Math. Anal. 51, 1137; also 1D, higher-d, underdamped).
The original equation is the σ→0 limit. Alternatively, discretize space and check that
density fluctuations come out right (Cornalba et al. 2023, arXiv:2303.00429).
- **Use in our paper:** our short-length cutoff ℓ_n (noise with a finite correlation
length, Lyapunov appendix) plays the role of σ. We can cite the regularized DK
equation as the formal version of what the cutoff does. It also explains why results
depend on ℓ_n: the unregularized object has no scale of its own.

## 4. Well-posedness: the bare DK equation is almost ill-posed (§4) ★★

- Konarovskyi, Lehmann & von Renesse (2019, Electron. Commun. Probab. 24, 1): for
non-interacting particles, the DK equation has solutions only for a **discrete set of
diffusion coefficients D**, and then only the trivial (empirical-measure) ones.
Extended to smooth interaction potentials (2020, J. Stat. Phys. 178, 666), singular
kernels, specific initial conditions and underdamped dynamics.
- The review calls this "surprising and puzzling" given how well the equation works when
used less rigorously.
- **Use in our paper:** one sentence of honesty. Any continuum use of eq:final-langevin
(closure, linearization, Lyapunov solution) implicitly regularizes it. Pairs with
topic 3.

## 5. Two nonlinearities and multiplicative noise (§2.2 iii–iv, §6.1) ★★★

- The interaction term is ∝ ρ²; the noise is ∝ √ρ. So the DK equation is nonlinear
**even for non-interacting particles**.
- The noise is **multiplicative**: its amplitude depends on ρ itself. Dean's physical
reading: "the density cannot fluctuate in regions devoid of particles," so empty
regions stay empty.
- Even for V=0 the statistics are non-Gaussian. Velenich et al. (2008, J. Phys. A 41,
  1. solved this case exactly with a field theory: Poissonian statistics. The review
  says it is the only exactly solved case.
- Non-Gaussian corrections for interacting particles at high density and weak
interaction: Illien & Carof (2025, arXiv:2501.16206).
- **Use in our paper:** backs the placeholder bullet on multiplicative noise. Also
relevant to the caveat that we use the Gaussian tail 19σ out: the full DK equation has
non-Gaussian density statistics that linearization discards (topic 7).

## 6. Free-energy form and the equilibrium steady state (§2.2 v; Dean 1996) ★★★

- The DK equation can be written ∂_t ρ = ∇·(ρ ∇ δF/δρ) + ∇·(√ρ η), with
F = k_BT∫ρ ln ρ + ½∫∫ρVρ. The ideal-gas entropy appears with no Stirling formula and no
large-N assumption (Dean notes this as a surprise).
- The functional Fokker–Planck equation has exp(−F/T) as its steady state.
- Dean also rederives the equation Kawasaki's way: posit ∂_t ρ = −∇·j + noise, with
j = −ρ∇δF/δρ, and **fix the noise by the fluctuation–dissipation theorem**. This gives
the same noise.
- Dean: even expanded to fourth order about a background density, the equation does
**not** reduce to model A or model B, because third-order terms stay in the noise
prefactor.
- **Use in our paper:** a key contrast. In equilibrium, the FDT route and the microscopic
route agree, so one could skip the microscopic derivation. Our system is active: no FDT
fixes the noise, so the microscopic route is the **only** way to get it. This makes the
placeholder's point "neither can be repaired afterwards by an FDT argument" concrete,
with Dean's own paper as the reference case. The model A/B remark shows that a
phenomenological noise is not a small-fluctuation approximation of the true one.

## 8. Three densities that look alike: DK vs Kawasaki vs DDFT (§3.3) ★★★

- Averaging the DK equation gives an unclosed equation for ⟨ρ⟩ with ⟨ρρ⟩ inside: a
BBGKY hierarchy (Dean shows this in one line for the two-point function).
- **DDFT** (Evans; Marconi & Tarazona 1999) closes it with the adiabatic approximation,
using the equilibrium excess free energy. The result is **deterministic**: it gives the
mean evolution and nothing about fluctuations. The review calls fluctuations "the core
of SDFT and its main advantage."
- The three densities differ: Kawasaki's ρ̂ is coarse-grained in space, Dean's ρ is
microscopic (a distribution), and the DDFT ρ̄ is an ensemble average. Mixing them up
"has been a source of confusion in the literature" (Archer & Rauscher 2004, J. Phys. A
37, 9325; Frusawa & Hayakawa 2000, J. Phys. A 33, L155; Kawasaki 2006; Chavanis 2011,
2019).
- **Use in our paper:** we need to say which density our closed, linearized equations
describe. After the closure and the cutoff ℓ_n, our field is effectively Kawasaki-like
(locally averaged) while its noise is inherited from Dean's. Stating this preempts the
Archer–Rauscher objection. The DDFT point fits the placeholder bullet on mean-field
continuum limits: they give the mean and lose the noise.

## 10. Extensions of the original equation (§5) ★★

- Inertia (underdamped DK; also granular collisions).
- Hydrodynamic interactions (Donev & Vanden-Eijnden 2014, J. Chem. Phys. 140, 234115):
the diffusion tensor depends on all positions, so the noise becomes multiplicative in
a harder way.
- Orientational degrees of freedom: a "straightforward" extension (dipoles; active
matter).
- Multiple species (laning, electrolytes, mixtures at different temperatures,
non-reciprocal mixtures).
- Reflecting or partially absorbing boundaries; stochastic resetting (Bressloff 2024).
- **Use in our paper:** our equation is an orientational extension, with an active drift  
and a reversal noise. The review treats orientation as routine; our non-routine parts  
are the active and reversal terms (topic 11).

## References not yet in `literature.bib`

Only `dean1996` and `bertin2013` are present. If we cite anything above, we need to add
the entry: illien2025 (the review), kawasaki1994, archer2004, frusawa2000,
konarovskyi2019, cornalba2019, velenich2008, fodor2018, spinney2024, bressloff2024,
tailleur2008, lutsko2012, liu2024, frusawa2019, magaletti2022.