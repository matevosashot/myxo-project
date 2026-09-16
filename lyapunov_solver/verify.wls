#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Verification for lyapunov_solver.m.  SMALL GRIDS ONLY -- the whole script is
   meant to run in well under a minute, so it is a correctness check and not a
   convergence study.  Most tests use SYNTHETIC steady states (analytic fields
   wrapped in a ListInterpolation) so that no FEM or comoving solve is needed;
   test 9 is the one end-to-end run and is opt-in.

       wolframscript -file verify.wls          # tests 1-8, seconds
       wolframscript -file verify.wls all      # + test 9 (real comoving solve)

   Every PASS/FAIL line is also appended to verify_results.csv. *)

thisDir = DirectoryName[$InputFileName];
repoRoot = ParentDirectory[thisDir];
Get[FileNameJoin[{thisDir, "lyapunov_solver.m"}]];

runAll = Length[$ScriptCommandLine] >= 2 && $ScriptCommandLine[[2]] === "all";

(* ------------------------------------------------------------------ *)
(* harness                                                            *)
(* ------------------------------------------------------------------ *)
results = {};
report[name_String, ok_, detail_] := (
  AppendTo[results, {name, TrueQ[ok], detail}];
  Print[If[TrueQ[ok], "  PASS  ", "  FAIL  "], StringPadRight[name, 46], detail];
  TrueQ[ok]);

(* close enough, relatively *)
rel[a_, b_] := Abs[a - b]/Max[1.*^-300, Abs[b]];

(* Single-line number formatting.  ToString on a small machine Real renders a
   SUPERSCRIPTED "m x 10^e", which OutputForm spreads over two lines and shreds
   the log -- the same trap comoving_steady_solver.m:245-248 records.  CForm
   never leaves one line. *)
fmt[x_] := If[NumericQ[x], ToString[CForm[N[x]]], ToString[x]];

(* ------------------------------------------------------------------ *)
(* synthetic steady states                                            *)
(* ------------------------------------------------------------------ *)
(* A steady state the module will accept: three 2-argument interpolations that
   survive Derivative[2,0], plus the frame velocity.  InterpolationOrder -> 3
   matches what comoving_steady_solver.m returns. *)
mkSS[rhoF_, q1F_, q2F_, uv_, L_, n_ : 120] := Module[{g, mk},
  g = N @ Subdivide[-L, L, n];
  mk[f_] := ListInterpolation[Outer[f, g, g], {g, g}, InterpolationOrder -> 3];
  <|"rho" -> mk[rhoF], "Q1" -> mk[q1F], "Q2" -> mk[q2F], "u" -> uv,
    "mesh" -> <|"box" -> L|>|>];

(* the dean.tex sec-parameter-choice block *)
mpOrdered[l_] := {B -> 10000., a -> -0.172, b -> 0.1, L -> 0.05368,
   \[Zeta] -> 10., \[Xi]0 -> 5., \[Xi]r -> 0.1, \[Rho]0 -> 0.3,
   \[CapitalLambda] -> 0.1, lNoise -> l};

(* Q_ss = 0 is only a fixed point in the DISORDERED phase, so tests that want a
   uniform isotropic base state flip the sign of a' = a + Lambda xir.  The
   pedestal zeta rho/(2 Pi B rho0 l^2) does not depend on a at all, so this
   costs the pedestal tests nothing. *)
mpIso[l_] := {B -> 10000., a -> 0.1, b -> 0.1, L -> 0.05368,
   \[Zeta] -> 10., \[Xi]0 -> 5., \[Xi]r -> 0.1, \[Rho]0 -> 0.3,
   \[CapitalLambda] -> 0.1, lNoise -> l};

ssFlat = mkSS[1. &, 0. &, 0. &, {0., 0.}, 26.];
quiet[e_] := Quiet[e, {SolveFluctuationsLyapunov::unstable,
                       SolveFluctuationsLyapunov::marginal,
                       SolveFluctuationsLyapunov::nullmode,
                       SolveFluctuationsLyapunov::notpsd}];

Print["================================================================"];
Print["lyapunov_solver.m verification   (small grids)"];
Print["================================================================"];

(* ------------------------------------------------------------------ *)
Print["\n[1] adjoint-pair identity  Lap + (Gx^T Gx + Gy^T Gy) == 0"];
(* ------------------------------------------------------------------ *)
Do[Module[{r},
   r = quiet @ SolveFluctuationsLyapunov[{Nint -> 11, fbox -> 8}, mpIso[3.], ssFlat, bc];
   report["lapIdentity exact, " <> bc, r["lapIdentity"] === 0,
          "lapIdentity = " <> fmt[r["lapIdentity"]]]],
  {bc, {"Dirichlet", "Neumann"}}];

(* ------------------------------------------------------------------ *)
Print["\n[2] sealed-wall conservation gates"];
(* ------------------------------------------------------------------ *)
Module[{r, sc},
  r = quiet @ SolveFluctuationsLyapunov[{Nint -> 11, fbox -> 8}, mpIso[3.], ssFlat, "Neumann"];
  (* consCheck aborts on failure, so reaching here at all is most of the test;
     assert the numbers are at round-off anyway. *)
  report["consA / consAT / consQ at round-off",
         Max[Abs[r["consA"]], Abs[r["consAT"]], Abs[r["consQ"]]] < 1.*^-9,
         "max = " <> fmt[Max[Abs[r["consA"]], Abs[r["consAT"]], Abs[r["consQ"]]]]];
  report["canonical: total rho-rho covariance = 0",
         Abs[Total @ Flatten @ r["C"][[1 ;; 121, 1 ;; 121]]] < 1.*^-12,
         "sum = " <> fmt[Total @ Flatten @ r["C"][[1 ;; 121, 1 ;; 121]]]]];

(* ------------------------------------------------------------------ *)
Print["\n[3] advection identity   Adv.1 == (B/xi0) Lap rho_ss"];
(* ------------------------------------------------------------------ *)
(* rho_ss = 1 + c(x^2+y^2)/2 has Lap rho = 2c and Grad rho = c(x,y), so the
   identity Div(w_ss) = (B/xi0) Lap rho_ss is non-trivially exercised.  c is
   kept small so rho stays near 1 over the box. *)
Do[Module[{r, c = 0.0004},
   r = quiet @ SolveFluctuationsLyapunov[{Nint -> 15, fbox -> 8}, mpIso[3.],
        mkSS[(1. + c (#1^2 + #2^2)/2) &, 0. &, 0. &, {0.3, -0.1}, 26.], bc];
   report["advection identity, " <> bc, NumericQ[r["advIdentity"]] && r["advIdentity"] < 0.02,
          "rel = " <> fmt[r["advIdentity"]]]],
  {bc, {"Dirichlet", "Neumann"}}];

(* ------------------------------------------------------------------ *)
Print["\n[4] cutoff calibration against dean.tex eq:lattice-bz / eq:iso-gauss"];
(* ------------------------------------------------------------------ *)
(* The reference is the EXACT finite-box rho-only pedestal, built from the same
   sine / DCT-II basis W is built in: in that sector W commutes with the drift,
   so C = W C0 W exactly.  The measured value carries in addition the dQ
   coupling that eq:sde drops, which dean.tex sec-cutoff-numerics puts at
   0.3-0.5%; that is the whole budget of this test. *)
Fbz[x_] := LyapunovSolver`Private`scaledI0[x];
exactPed[bc_, nI_, fb_, lEff_, zt_, Bb_, r0_] := Module[{h, W, v, ic, d},
  If[bc === "Neumann",
     h = 2.*fb/nI;
     W = Exp[-(lEff/h)^2 Sin[N[Pi] Range[0, nI - 1]/(2 nI)]^2];
     v = Table[Sqrt[If[p == 0, 1., 2.]/nI] Cos[N[Pi] p (i - 1/2)/nI], {i, nI}, {p, 0, nI - 1}],
     h = 2.*fb/(nI + 1);
     W = Exp[-(lEff/h)^2 Sin[N[Pi] Range[nI]/(2 (nI + 1))]^2];
     v = Sqrt[2./(nI + 1)] Table[Sin[N[Pi] p i/(nI + 1)], {i, nI}, {p, nI}]];
  ic = (nI + 1)/2; d = Total[W^2 v[[ic]]^2];
  (zt/(Bb r0 h^2)) (d^2 - If[bc === "Neumann", 1./nI^2, 0.])];

Print["    bc         Nint  l     t      meas/exact OFF  meas/exact ON   boxFactor"];
Do[Module[{nI = q[[1]], bc = q[[2]], l = q[[3]], rOff, rOn, ic, t, mOff, mOn,
           eOff, eOn, tag},
   tag = bc <> " N=" <> ToString[nI] <> " l=" <> ToString[l];
   rOff = quiet @ SolveFluctuationsLyapunov[
      {Nint -> nI, fbox -> 20, rescaleCutoff -> False}, mpIso[l], ssFlat, bc];
   rOn = quiet @ SolveFluctuationsLyapunov[
      {Nint -> nI, fbox -> 20, rescaleCutoff -> True}, mpIso[l], ssFlat, bc];
   ic = (nI + 1)/2; t = rOff["lNoiseOverH"];
   mOff = rOff["sigmaRho"][[ic, ic]]; mOn = rOn["sigmaRho"][[ic, ic]];
   eOff = exactPed[bc, nI, 20, rOff["lNoiseEff"], 10., 10000., 0.3];
   eOn  = exactPed[bc, nI, 20, rOn["lNoiseEff"], 10., 10000., 0.3];
   Print["    ", StringPadRight[bc, 11], StringPadRight[ToString@nI, 6],
      StringPadRight[ToString@l, 6], StringPadRight[ToString@NumberForm[t, 4], 7],
      StringPadRight[ToString@NumberForm[mOff/eOff, 6], 16],
      StringPadRight[ToString@NumberForm[mOn/eOn, 6], 15],
      ToString@NumberForm[rOn["boxFactor"], 6]];
   (* (a) W is built right: the run reproduces the exact finite-box rho-only
      pedestal.  The residual is the dQ coupling eq:sde drops, which
      dean.tex sec-cutoff-numerics puts at 0.3-0.5%. *)
   report["W matches exact finite box (uncal.), " <> tag, rel[mOff, eOff] < 0.01,
          "rel = " <> fmt[rel[mOff, eOff]]];
   report["calibrated matches exact finite box, " <> tag, rel[mOn, eOn] < 0.01,
          "rel = " <> fmt[rel[mOn, eOn]]];
   (* (b) Under Dirichlet the BZ closed form IS the finite-box answer, which is
      what makes it a legitimate calibration target, and the calibrated run then
      lands on eq:iso-gauss.  Under Neumann neither holds -- see below. *)
   If[bc === "Dirichlet",
      report["closed form = exact finite box, " <> tag,
             rel[rOff["boxFactor"], 2 Pi t^2 Fbz[t^2]^2] < 0.01,
             "rel = " <> fmt[rel[rOff["boxFactor"], 2 Pi t^2 Fbz[t^2]^2]]];
      report["calibrated pedestal = eq:iso-gauss, " <> tag,
             rel[rOn["boxFactor"], 1.] < 0.01,
             "boxFactor = " <> ToString@NumberForm[rOn["boxFactor"], 6]]]],
  {q, {{21, "Dirichlet", 2.}, {31, "Dirichlet", 4.}, {21, "Neumann", 2.}}}];

(* Under a SEALED wall the calibrated pedestal is deliberately NOT eq:iso-gauss:
   the canonical projection removes the uniform mode, so a correlation length
   comparable with the box leaves genuinely less room for a density fluctuation.
   The claim that this is geometry rather than discretisation is testable: the
   factor must depend on lNoise/fbox ALONE, i.e. be invariant under scaling the
   box and the cutoff together at fixed Nint. *)
Print["\n    sealed-box suppression is a function of lNoise/fbox alone:"];
Module[{bf},
  bf[fb_, ratio_] := quiet[SolveFluctuationsLyapunov[{Nint -> 21, fbox -> fb},
      mpIso[N[ratio fb]], ssFlat, "Neumann"]]["boxFactor"];
  Do[Module[{a = bf[10, ratio], b = bf[20, ratio]},
     Print["      lNoise/fbox = ", ratio, "   fbox=10: ",
           ToString@NumberForm[a, 6], "   fbox=20: ", ToString@NumberForm[b, 6]];
     report["boxFactor scale-invariant at l/fbox = " <> ToString[ratio],
            rel[a, b] < 1.*^-6, "rel = " <> fmt[rel[a, b]]]],
    {ratio, {0.1, 0.2}}]];

(* ------------------------------------------------------------------ *)
Print["\n[5] unreachable cutoff aborts cleanly (dean.tex eq:ln-eff)"];
(* ------------------------------------------------------------------ *)
Module[{hD = 2.*20/22, r, lBad},
  lBad = 0.9 hD/Sqrt[2 Pi];     (* just below the solvability edge *)
  r = Quiet @ SolveFluctuationsLyapunov[{Nint -> 21, fbox -> 20}, mpIso[lBad],
        ssFlat, "Dirichlet"];
  report["lNoise < h/Sqrt[2Pi] -> $Failed", r === $Failed,
         "lNoise = " <> ToString@NumberForm[lBad, 4] <>
         " < h/Sqrt[2Pi] = " <> ToString@NumberForm[hD/Sqrt[2 Pi], 4]];
  (* and the same lNoise is fine once the grid is fine enough.  Nint = 31 is
     already past the edge; anything larger would blow the runtime budget of
     this script on the O(n^3) Eigensystem for no extra information. *)
  r = quiet @ SolveFluctuationsLyapunov[{Nint -> 31, fbox -> 20}, mpIso[lBad],
        ssFlat, "Dirichlet"];
  report["...but succeeds on a finer grid", AssociationQ[r],
         "Nint = 31 gives h/Sqrt[2Pi] = " <>
         ToString@NumberForm[(2.*20/32)/Sqrt[2 Pi], 4]]];

(* ------------------------------------------------------------------ *)
Print["\n[6] calibration is a vanishing correction:  lNoiseEff/lNoise -> 1"];
(* ------------------------------------------------------------------ *)
Module[{rs, ratios},
  rs = Table[LyapunovSolver`Private`cutoffLength[3., 2.*20/(n + 1)]/3.,
             {n, {21, 41, 81, 161, 321}}];
  ratios = rs;
  report["lNoiseEff/lNoise decreasing to 1",
         And @@ Table[ratios[[i]] > ratios[[i + 1]], {i, Length[ratios] - 1}] &&
         Abs[Last[ratios] - 1] < 0.002,
         "Nint 21..321 -> " <> ToString[NumberForm[#, 6] & /@ ratios]];
  report["large-t asymptote lEff ~ l + h^2/(8 l)",
         rel[LyapunovSolver`Private`cutoffLength[3., 0.3], 3. + 0.3^2/(8*3.)] < 1.*^-3,
         "exact " <> ToString@NumberForm[LyapunovSolver`Private`cutoffLength[3., 0.3], 8] <>
         " vs asymptotic " <> ToString@NumberForm[3. + 0.3^2/(8*3.), 8]]];

(* ------------------------------------------------------------------ *)
Print["\n[7] REGRESSION: w_ss == 0 must reproduce the superseded modules"];
(* ------------------------------------------------------------------ *)
(* With rho_ss constant and u = 0 the advection velocity vanishes identically,
   so Adv == 0 and the old modules' (B/xi0) diag[Lap rho_ss] is 0 too: the two
   schemes must then agree to round-off.  Q_ss is NOT zero here, so the closure,
   A21/A31, the cross noise terms and W are all exercised.  This isolates the
   new advection from every other change and needs no compatibility flag.
   The old modules are loaded AFTER this one; they open a different private
   context and their public names differ, so nothing is clobbered. *)
Module[{ssQ, nI = 15, fb = 8, lN = 2.5, rNew, rOld, mp},
  ssQ = mkSS[1. &, 0.9 &, 0. &, {0., 0.}, 26.];
  mp = mpOrdered[lN];
  rNew = <||>; rOld = <||>;
  rNew["D"] = quiet @ SolveFluctuationsLyapunov[{Nint -> nI, fbox -> fb,
                 rescaleCutoff -> False}, mp, ssQ, "Dirichlet"];
  rNew["N"] = quiet @ SolveFluctuationsLyapunov[{Nint -> nI, fbox -> fb,
                 rescaleCutoff -> False}, mp, ssQ, "Neumann"];
  Print["    loading the superseded modules for comparison ..."];
  Block[{$Output = {}},
    Get[FileNameJoin[{repoRoot, "lyapunov_solver_module_cutoff.m"}]];
    rOld["D"] = Quiet @ SolveFluctuationsLyapunovCutoff[{Nint -> nI, fbox -> fb}, mp, ssQ];
    Get[FileNameJoin[{repoRoot, "lyapunov_solver_module_neumann_cutoff.m"}]];
    rOld["N"] = Quiet @ SolveFluctuationsLyapunovNeumannCutoff[{Nint -> nI, fbox -> fb}, mp, ssQ]];
  Do[Module[{d},
     d = Max @ Abs[rNew[k]["sigmaRho"] - rOld[k]["sigmaRho"]]/
         Max @ Abs[rOld[k]["sigmaRho"]];
     report["w_ss=0 reproduces old module, " <> If[k === "D", "Dirichlet", "Neumann"],
            d < 1.*^-10, "max rel diff in Sigma_rho = " <> fmt[d]]],
    {k, {"D", "N"}}]];

(* ------------------------------------------------------------------ *)
Print["\n[8] sealed-box closed form, Q_ss = 0 and lNoise = 0"];
(* ------------------------------------------------------------------ *)
(* Sigma = (zeta rho/(B rho0 h^2)) [I - 1 1^T/Nint^2] -- the flat pedestal
   (eq:fdt-split) and the canonical rank-one subtraction (eq:canonical) in one
   object.  Exact for the rho sector ALONE; here the full three-field system
   also carries the drho-dQ noise cross-blocks, which are non-zero even at
   Q_ss = 0 (QQ1r = pf divBlk[rho/2, -rho/2, 0, 0]).  So the right expectation
   is the ~1e-3 of that coupling, not round-off -- the old modules' 5e-16 was
   measured on a rho-only probe. *)
Module[{nI = 13, fb = 8, r, h, pred, C0, d},
  r = quiet @ SolveFluctuationsLyapunov[{Nint -> nI, fbox -> fb}, mpIso[0.],
        ssFlat, "Neumann"];
  h = 2.*fb/nI;
  pred = (10./(10000.*0.3 h^2));
  C0 = r["C"][[1 ;; nI^2, 1 ;; nI^2]];
  d = Max @ Abs[C0 - pred (IdentityMatrix[nI^2] - ConstantArray[1./nI^2, {nI^2, nI^2}])]/pred;
  report["Sigma = pedestal (I - 11^T/N^2) up to dQ", d < 5.*^-3,
         "max rel dev = " <> fmt[d] <> "  (dQ coupling)"]];

(* ------------------------------------------------------------------ *)
If[runAll,
  Print["\n[9] end-to-end with a real comoving steady state"];
  (* ---------------------------------------------------------------- *)
  Get[FileNameJoin[{repoRoot, "comoving_defect", "comoving_steady_solver.m"}]];
  Module[{mp, sp, ss, r},
    mp = {B -> 10000., a -> -0.172, b -> 0.1, L -> 0.05368, \[Zeta] -> 10.,
          \[Xi]0 -> 5., \[Xi]r -> 0.1, \[Rho]0 -> 0.3, \[CapitalLambda] -> 0.1,
          lNoise -> 2.0};
    sp = {\[Delta]mesh -> 0.12, hMax -> 2., box -> 10., maxSteps -> 25,
          relativeDiff -> 1.*^-8, verbose -> False};
    Do[Module[{},
       ss = SolveActiveNematicComovingSteady[mp, sp, bc];
       If[ss === $Failed, report["comoving steady state, " <> bc, False, "$Failed"],
         report["comoving steady state, " <> bc, TrueQ[ss["solver"]["converged"]],
                "u = " <> ToString[ss["u"]] <> ", |R| = " <>
                fmt[ss["solver"]["residual"]]];
         r = quiet @ SolveFluctuationsLyapunov[{Nint -> 15, fbox -> 8}, mp, ss, bc];
         report["end-to-end Lyapunov, " <> bc,
                AssociationQ[r] && r["residualRefined"] < 1.*^-9 && r["peclet"] < 1,
                "resid = " <> fmt[r["residualRefined"]] <>
                ", Peclet = " <> ToString@NumberForm[r["peclet"], 4] <>
                ", maxReLam = " <> ToString@NumberForm[r["maxReLambda"], 4]]]],
      {bc, {"Dirichlet", "Neumann"}}]],
  Print["\n[9] skipped (pass \"all\" to run the end-to-end comoving test)"]];

(* ------------------------------------------------------------------ *)
Print["\n================================================================"];
Module[{nPass = Count[results, {_, True, _}], nTot = Length[results], csv},
  Print[nPass, " / ", nTot, " passed"];
  If[nPass < nTot,
     Print["FAILED: ", StringRiffle[Cases[results, {n_, False, _} :> n], ", "]]];
  csv = FileNameJoin[{thisDir, "verify_results.csv"}];
  Export[csv, Prepend[
     Map[{DateString["ISODateTime"], #[[1]], If[#[[2]], "PASS", "FAIL"], #[[3]]} &, results],
     {"ISO8601", "test", "status", "detail"}]];
  Print["Wrote ", csv];
  If[nPass < nTot, Exit[1]]];
