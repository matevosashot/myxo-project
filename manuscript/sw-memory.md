# Writing memory for this paper (sw-* skills)



## Paper

- Main file: manuscript.tex
- Type / venue / readership: journal article · Physical Review E (revtex4-2, aps/pre) · statistical and soft-matter physicists in active matter; biophysicists working on bacterial collectives (inferred)
- Central question(s): How does the first new layer form in a Myxococcus xanthus monolayer, and why does it form preferentially at +1/2 defects although the mean pressure there is far below the extrusion threshold? (inferred)
- Intended take-home message: two messages. (1) Layering is explained by nucleation theory over a fluctuating barrier. (2) The statistics of the pressure around a topological defect, obtained systematically from a microscopic model. (you said, 2026-11-21)
- Reader already knows: IK1 active nematics, the director and ±1/2 topological defects; IK2 continuum theory of dry active matter (active stress, friction with a substrate); IK3 Langevin equations, white noise, Itô calculus; IK4 classical nucleation theory; IK5 the fluctuation–dissipation theorem; IK6 that M. xanthus cells glide and reverse (inferred)
- Reader does not know: NK1 the cited M. xanthus measurements beyond what the paper reports; NK2 our microscopic model and its parameters; NK3 our coarse-graining; NK4 our numerical solver; NK5 unpublished observations of the experimental groups (inferred)
- Results (lab book): primary: the critical pressure of the nucleation model, three to four orders above the mean pressure (S2); fluctuating active-nematic equations in which one coefficient sets active stress and noise (S3); the density/pressure variance around a +1/2 defect, a few percent, peaked off the core at the head (S4.2, Eq. 44); the layering-rate enhancement of about 2×10² and its radial profile vs Copenhagen Fig. 3f, one fitted parameter B (S5); secondary: defect drift direction and closure (S4.1, App. E), the reversal-rate trend (S6.2), −1/2 defects (App. G), the Lyapunov solver check (App. H); SI candidates: the appendices (inferred)
- Results that speak against the conclusion: the white-noise limit fails at the cutoff scale (τc ≈ 0.4 ms vs 5 min); Gaussian tail with less than one cell per cutoff area; the lowest-order closure drifts the defect the wrong way; the core-centered profile misses the first measured bin (inferred)
- Native language of the author: <unknown> (you said, 2026-11-21: leave unknown)

## Intervention level

- full
- the main text (Abstract, S1–S6) at full; the appendices polish-only, run with --level polish-only (you said, 2026-11-21)

## Style and conventions

- Spelling: en-US (you said, before run 1) · Figure refs: \cref (prints "Fig.") · Equations: \cref · Dashes: --- unspaced, but see Rules adjusted (detected)
- Voice: "we" + active verbs; plain sentences LINKED by connectives (so, therefore, thus, because, since, however, but, instead, as a result) into a flow of thought; a paragraph must not read as a list of isolated facts (you said, before run 1)
- Motivate, then calculate: say why a quantity is needed before computing it; state the claim first, then the reason (you said, before run 1)
- No symbol before its definition; in summaries name the quantity in words before its equation number, say what a change is relative to, no metaphors (you said, before run 1)
- Angle brackets: use \langle … \rangle instead of \la … \ra when argument has \hat on it
- Caption panel labels bold: \textbf{(a)}, also for in-caption references to the same figure's panels (you said, before run 1)
- Prompt macros: \prompt (instruction, also when commented out); not addressed: \todo, % SW-TODO (inferred, sw-write 002)
- Highlighting: no \upd tags (this overrides CLAUDE.md); each round is shown as a latexdiff PDF next to the tex (you chose, run 9 and the 2026-10-26 rounds). Earlier \upd/\updb/\updc tags were all stripped on 2026-10-16 at your request
- Line width (CLAUDE.md): running text in manuscript.tex stays under 60 characters per line, broken at spaces; refill only the paragraph an edit touches; preamble, headings, equations and figure environments keep their layout. Note for sw runs: `sw.py apply` re-wraps changed lines at 60 and also re-wraps captions; run .sw/tools/fixwidth.py (≤59) and .sw/tools/restorecap.py after each apply (learned 2026-11-21)



## Terminology

- nematic field — the Q field around a defect — not: texture
- Dean–Kawasaki equation — the stochastic density equation — not: "Dean's method"; use "Dean–Kawasaki" only in the Introduction, "Dean's equation" elsewhere (you said, run 4)
- refractory period — a time after each reversal during which the cell cannot reverse again — not: refractory clock
- head of the defect — the side where the director bends around the core (Fig. 2(a)) — not: bend side
- time difference |t − t′| — not: lag
- charge selectivity (which defect) vs sign selectivity (layers vs holes) — keep apart
- v0/frev — distance glided between reversals — not: ℓ_run
- boundary — the edge of the domain (reservoir boundary, no-flux boundary, sealed boundary, boundary nodes) — not: wall (you said, 2026-10-16 prompt; "cell wall" of the bacterium stays)
- nematic force field — H = −δF/δQ — not: molecular field (you said, run 4)
- ξ0 — the friction coefficient per unit area; in the Fig. 3 caption explain ξ0 itself, not the ratio ξ0/ρ0 (Sec. III A keeps its ratio explanation) (you said, sw-write 002)
- nematic part (of the force) — Q̂(θ):∇H in eq:force-on-cells, Q:∇H in the density flux (= ∇·σ^n) — not: elastic force/term/part (both parts are elastic); compressibility part — −B∇ρ and the f_B part of the stress (not: isotropic part, because the nematic part has an isotropic part too; you said, review of stage two, 2026-11-21); sec:model-elastic is titled "Dropping the nematic part of the force" (you said, sw-write 006, 2026-11-20)
- F(r,θ) of eq:force-on-cells — the force on a rod; cite it instead of re-deriving it from eq:chain-rule (instruction, sw-write 002 prompts)
- f_passive = ⟨ρ0F⟩_ϱ ≃ −B∇ρ + Q:∇H is eq:passive-force (second line under eq:force-on-cells); Sec. III E cites it (sw-write 012)
- Sec. II rate notation: $\calP_\r(p)$ — probability density of the pressure at r (the old layering-probability $\calP$ is gone); $\calJ_p(\r)$ — layering rate at fixed pressure; $\calJ_{\Acell}$ — single-cell extrusion rate, always written with  (learned from your 2026-11 edits, sw-write 012)
- effective free energy — F with a′, K′ in place of a, K; H^a is minus its functional derivative (defined after Eq. (23)) (you said, run 4)
- lifting the white-noise limit — "keep the dichotomous noise of the reversals, with its finite correlation time, instead of its white-noise limit" [bressloff2024global, spinney2025dean] — not: "Poisson noise" (the model's reversals already form a Poisson process) (you said, sw-write 017)
- pressure as a slow variable / time-scale separation — what Eq. (layering-rate-approximation) requires — not: "slow limit" (only in the Sec. II remark; never rely on it elsewhere), not: "fast pressure" (you said, sw-write 019)
- J_Acell enhanced near defects — a possibility, never a fact: write "it is possible that ... is itself enhanced" (you chose, review of stage two, 2026-11-21; replaces "we speculate ... could")
- simplification — for the modeling choices (white-noise limit, closure, ...) — not: approximation ("in physics we have simplifications, not approximations"; you said, review of stage one, 2026-11-21; rephrase so it reads naturally, do not just swap the word; technical names such as Padé approximation and one-constant approximation stay; applies to the appendices too, review of round 3)
- Sec. III D: drop the nematic part of the force as a whole because both ∇H terms come from it; justify by the literature that drops them, not by an argument (you said, 2026-11-21)
- diffusive flux — $-B\nabla\rho/\xi_0$, the flux that balances the active force and with which the cross-coupling term convects $\Q$ — not: "pressure-driven flow", "the velocity $-B\nabla\rho/\xi_0$" (you said, review of round 3, 2026-11-21); "passive flux" stays as the general term, "flows" stays for the measured flows of cells
- active force — $\bfa$ / $\fa$ print as f_active / f^active (definitions.tex); figure labels must match ($|\bm f_{\!\text{active}}|$), not f^a (you changed the definition, 2026-11-21)
- anisotropic diffusion of nematic order — the transport by self-propulsion that the term Y carries (rods diffuse along their axis); the main text explains the drift reversal with it and never names Y (you said, review of round 3)
- cell division — not a cause of the layers; stated once, in the Introduction (copenhagen2021 Methods p. 6: 1-3 h of imaging, colony does not grow substantially, doubling time ~10 h at 22 °C); not repeated in the model section (you said, sw-write 2026-11-21)



## Rules adjusted for this paper

- L7 strict — no em dashes as parenthetical asides; use a comma, a colon or a new sentence (you said, before run 1)
- P13 relaxed — causal "since" is fine; keep it (you said, before run 1)
- B2 relaxed for the author's own sentences — flag a repetition in the author's text as a % SW-Q, never cut it (you said, before run 1)
- a caption of a figure may keep a clause that states the result its panel shows, because some readers read the captions first



## Do not touch

- \todo{...} commands: leave them exactly as they are, never add or remove one (you said, before run 1); count 2 active (+1 commented out) at sw-write 017; the \todo{this should be deleted:} you added was removed with its block at your answer (5 at run 10; all yours)
- the author's notation macros: \fa, \bEta_\rho, \bEta^\Q, \Ha, \sigmarho (you said, before run 1)



## Standing instructions

- Comments in manuscript.tex: only short pointers stay, i.e. "% Source:" lines to the scripts that make a figure or a number, "% dean.tex ..." line references, one-line literature page notes and one-line code notes. No dated edit history, no SW-TODO figure notes, no commented-out drafts (you cleaned them out on 2026-11-20)


## Learned from your edits

- Your review of round 5 (2026-11-21): text and caption must not repeat each other; the text gives the physical reason, the caption reads the mechanism off the figure's marks; no pronoun ("it") whose referent is a sentence away

- Your review of round 4 (2026-11-21): the anisotropic diffusion of Q comes from self-propulsion; the reversals only make it diffusive (never "with reversals, each rod diffuses"); rods and Q diffuse more along the director than across it. In the main text give the intuition only and point to the figure; the mechanism goes in the caption, which must stay complete and clear (do not shorten captions at the cost of information). Name arrows and curves by their mark ("red arrow", "top curve"); head/tail labels sit next to the diagram they refer to

- Your review of round 3 (2026-11-21): figures show the nematic field as continuous (gray) director lines, never as rod segments; a schematic that explains a main-text result belongs in the main text with the explanation in its caption

- Your review of stage two (2026-11-21): rewrite larger chunks when an answer needs it, not only local swaps; prefer "Padé profile" for Q_Padé (Fig. 17 caption)
- Your review of stage one (annotated PDF, 2026-11-21): you struck out added explanations you found unnecessary (J_Acell vs Boltzmann factor, "since the extruded cell then forms a patch…", "In the linearized theory the pressure is Gaussian", z∞≃21 in Sec. V, "the one parameter that we fit" in Sec. III C) and the nematic-stress sentence of footnote 6: prefer fewer bridging clauses

- Claim level of the comparison: say "we compare" the layering profile with the measurement; no "matches"/"reproduces the measured enhancement" in the abstract or the closing (learned sw-write 017: you deleted both, and chose a closing sentence without the fit claim)
- White-noise argument: Sec. IV B gives the numbers (tau_c, 5 ms, minutes) and the mismatch; Sec. VI A only reason and cost; the outlook the plan with the citations (you said, sw-write 017)



