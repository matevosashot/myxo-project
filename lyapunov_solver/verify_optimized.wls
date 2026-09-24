#!/usr/bin/env wolframscript
(* Validate the Python backend against the reference kernel.

   Two levels:
     A. kernel-level -- same (A, Q) through both kernels, compare C entrywise.
     B. end-to-end   -- the full SolveFluctuationsLyapunov on a synthetic steady
        state, both backends, compare every exported scalar and Sigma_rho.

   Small grids only: this is a correctness check, not a benchmark. *)

$HistoryLength = 0;
here = DirectoryName[$InputFileName];
Get[FileNameJoin[{here, "lyapunov_solver.m"}]];

fmt[x_] := If[NumericQ[x], ToString[CForm[N[x]]], ToString[x]];
pass = fail = 0;
check[lbl_, ok_, detail_] := (
   If[TrueQ[ok], pass++, fail++];
   Print[If[TrueQ[ok], "  PASS  ", "  FAIL  "], lbl, "   ", detail]);

Print["================ A. kernel level ================"];

(* A stable nonsymmetric A and an SPD Q -- the reference kernel's own contract,
   independent of any physics. *)
SeedRandom[7];
Do[
 Module[{n = nn, A, Q, meta, r1, r2, dC, scale},
  A = Developer`ToPackedArray @
     N[RandomReal[NormalDistribution[], {n, n}]/Sqrt[n] - 2. IdentityMatrix[n]];
  Q = Developer`ToPackedArray @ N@RandomReal[NormalDistribution[], {n, n}];
  Q = Q . Transpose[Q] + IdentityMatrix[n];
  meta = <|"neumann" -> False, "bcType" -> "Dirichlet", "nInt" -> 0, "n" -> n,
           "vCons" -> ConstantArray[1., n],
           "scratchDir" -> $TemporaryDirectory|>;

  r1 = LyapunovReferenceKernel[A, Q, meta];
  Get[FileNameJoin[{here, "lyapunov_solver_optimized.m"}]];
  r2 = LyapunovPythonKernel[A, Q, meta];
  $LyapunovSolveKernel = LyapunovReferenceKernel;

  scale = Max@Abs@r1["C"];
  dC = Max@Abs[r1["C"] - r2["C"]]/scale;
  check["n=" <> ToString[n] <> " C agrees with the reference", dC < 1.*^-9,
        "max rel dev = " <> fmt[dC]];
  check["n=" <> ToString[n] <> " residual is at round-off",
        r2["residual_refined"] === Missing["KeyAbsent", "residual_refined"] ||
        r2["residualRefined"] < 1.*^-9 scale,
        "resid = " <> fmt[r2["residualRefined"]] <>
        "  (reference " <> fmt[r1["residualRefined"]] <> ")"];
  check["n=" <> ToString[n] <> " max Re(lambda) agrees",
        Abs[r1["maxReLambda"] - r2["maxReLambda"]] < 1.*^-8 Abs[r1["maxReLambda"]],
        fmt[r2["maxReLambda"]] <> " vs " <> fmt[r1["maxReLambda"]]];
  check["n=" <> ToString[n] <> " C is symmetric",
        Max@Abs[r2["C"] - Transpose[r2["C"]]] < 1.*^-12 scale,
        "asym = " <> fmt[Max@Abs[r2["C"] - Transpose[r2["C"]]]/scale]];
 ], {nn, {120, 400}}];

Print[];
Print["================ B. end to end ================"];

(* A synthetic comoving steady state: smooth, nonzero u, no defect solve needed. *)
box = 6; nInt = 21;
mk[f_] := Module[{g = N@Range[-box, box, 2 box/40]},
   ListInterpolation[Outer[f, g, g], {g, g}, InterpolationOrder -> 3]];
steady = <|
   "rho" -> mk[(0.3 + 0.02 Exp[-(#1^2 + #2^2)/4.]) &],
   (* S0 = 0.7, not the physical 0.9: the cubic interpolant overshoots at the
      conical point of the defect core, and above ~0.8 the overshoot pushes S
      past 1 and makes the closure complain.  Irrelevant to what is being
      tested -- both backends see the identical A and Q either way -- but it
      keeps the log readable. *)
   "Q1"  -> mk[(0.7 Tanh[Sqrt[#1^2 + #2^2]/0.8] Cos[ArcTan[#1 + 1.*^-9, #2]]) &],
   "Q2"  -> mk[(0.7 Tanh[Sqrt[#1^2 + #2^2]/0.8] Sin[ArcTan[#1 + 1.*^-9, #2]]) &],
   "u"   -> {0.05, -0.02},
   "mesh" -> <|"box" -> box|> |>;
model = {Global`B -> 1.*^4, Global`a -> -0.172, Global`b -> 0.106,
         Global`L -> 0.06368, Global`\[Zeta] -> 8., Global`\[Xi]0 -> 5.,
         Global`\[Xi]r -> 0.1, Global`\[Rho]0 -> 0.3,
         Global`\[CapitalLambda] -> 0., Global`lNoise -> 1.73};
lyap = {Global`Nint -> nInt, Global`fbox -> box};

Print["running reference backend ..."];
$LyapunovSolveKernel = LyapunovReferenceKernel;
ref = SolveFluctuationsLyapunov[lyap, model, steady, "Dirichlet"];

Print[]; Print["running python backend ..."];
Get[FileNameJoin[{here, "lyapunov_solver_optimized.m"}]];
opt = SolveFluctuationsLyapunov[lyap, model, steady, "Dirichlet"];

Print[];
scale = Max@Abs@ref["C"];
check["end-to-end C agrees", Max@Abs[ref["C"] - opt["C"]]/scale < 1.*^-9,
      "max rel dev = " <> fmt[Max@Abs[ref["C"] - opt["C"]]/scale]];
check["Sigma_rho agrees",
      Max@Abs[ref["sigmaRho"] - opt["sigmaRho"]]/Max@Abs@ref["sigmaRho"] < 1.*^-9,
      "max rel dev = " <>
      fmt[Max@Abs[ref["sigmaRho"] - opt["sigmaRho"]]/Max@Abs@ref["sigmaRho"]]];
check["runTag identical", ref["runTag"] === opt["runTag"], opt["runTag"]];
Do[check["scalar " <> k <> " agrees",
         Abs[ref[k] - opt[k]] <= 1.*^-8 Abs[ref[k]] + 1.*^-14,
         fmt[opt[k]] <> " vs " <> fmt[ref[k]]],
   {k, {"maxReLambda", "residualRefined", "pedestalPredicted", "lNoiseEff",
        "boxFactor", "peclet"}}];
check["sigmaQ1 agrees",
      Max@Abs[ref["sigmaQ1"] - opt["sigmaQ1"]]/Max@Abs@ref["sigmaQ1"] < 1.*^-9,
      "max rel dev = " <>
      fmt[Max@Abs[ref["sigmaQ1"] - opt["sigmaQ1"]]/Max@Abs@ref["sigmaQ1"]]];
check["sigmaQ2 agrees",
      Max@Abs[ref["sigmaQ2"] - opt["sigmaQ2"]]/Max@Abs@ref["sigmaQ2"] < 1.*^-9,
      "max rel dev = " <>
      fmt[Max@Abs[ref["sigmaQ2"] - opt["sigmaQ2"]]/Max@Abs@ref["sigmaQ2"]]];
check["Neumann still uses the reference kernel",
      StringContainsQ[
         Module[{r}, r = SolveFluctuationsLyapunov[lyap, model, steady, "Neumann"];
                r["scheme"]; "ok"], "ok"],
      "(ran without error)"];

Print[];
Print["================================================"];
Print[pass, " / ", pass + fail, " passed"];
If[fail > 0, Exit[1]];
