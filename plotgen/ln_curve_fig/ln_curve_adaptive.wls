#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ln_curve_adaptive.wls -- ln_curve figure (fig:ln_curve), same style as
   ln_curve_draft.wls, with two changes (2026-10-02):

   1. nu0 is computed by ADAPTIVE integration of
         nu(x) = (lN^2/4) Int d^2w kappa(w) d_i d_j Q_ij(x - w),
         kappa(w) = (1 - Exp[-2 w^2/lN^2])/(2 Pi w^2),
      in coordinates centred on the evaluation point, with breakpoints at lN and 5 lN.
      The fixed log grid of ln_curve_draft.wls (radial spacing ~0.005 um near the peak)
      stops resolving the kernel for lN below ~0.01 um; there the logarithm is cut off by
      the grid and the curve turns into a spurious pure lN^2 law.
   2. The dashed reference line is the small-lN asymptote
         nu0 ~ lN^2 [ a Log[ld/lN] + c ],   a = max|d_i d_j Q_ij| / 4,
      with c fitted on the points lN <= 1e-3 um. It replaces the 0.79 lN^2 fit.

   Usage:  wolframscript -file ln_curve_adaptive.wls [outputDir] [recompute]
   Without an argument it exports to <repoRoot>/figures, as ln_curve_draft.wls does.
   Runtime: ~30 min (about 30 two-dimensional integrals per lN, 21 values of lN). The
   results are cached in <outputDir>/ln_curve_data.m and reused on the next run, so
   restyling the plot takes seconds; pass "recompute" as the second argument to redo them. *)

SetEnvironment["QT_QPA_PLATFORM" -> "offscreen"];   (* headless front end for Export *)

<< MaTeX`
w = 8.6*72/2.54;
tex[s_, fs_ : 11] :=
  MaTeX[s, FontSize -> fs,
   "Preamble" -> {"\\usepackage{amsmath,siunitx,newtxtext,newtxmath}"}];
colors = ColorData[97, "ColorList"];
magnify[m_ : 2] := Function[{fig}, Magnify[fig, m]];

repoRoot = If[$MachineName == "AMatevosyanL1M",
  "/Users/ashmat/Desktop/fluctuation-solver",
  "/home/ashmat/Projects/myxo-project-dean"];
figDir = FileNameJoin[{repoRoot, "figures"}];
If[Length[$ScriptCommandLine] >= 2, figDir = $ScriptCommandLine[[2]]];
If[! DirectoryQ[figDir], CreateDirectory[figDir]];

Get[FileNameJoin[{repoRoot, "binaryReadTools.m"}]];
Import[FileNameJoin[{repoRoot, "RemoteExport.m"}]];
export[name_, fig_] :=
 Block[{}, RemoteExport[FileNameJoin[{figDir, name <> ".pdf"}], fig];
  Export[FileNameJoin[{figDir, name <> "_high.png"}], fig, ImageResolution -> 900];
  Export[FileNameJoin[{figDir, name <> ".png"}], fig, ImageResolution -> 600]];

dataPath = FileNameJoin[{repoRoot,
   "data/fluctuations-set2_2_comove_optimized_2/covariance_box20_N180_Dir_lN0.5_z12.5_ld1_B90000_ss20.bin"}];
modelParams = readSteady[dataPath]["modelParams"];
padeParams = Table[x -> (x /. modelParams), {x, {a1, L1, b}}];
ld = Sqrt[-L1/a1] /. padeParams;

(* accurate Pade approximant, rational in z = r^2/ld^2 (within 1e-4 S0 of the exact profile) *)
fPade[z_] :=
  Sqrt[(0.34010999 + 0.04516965 z + 0.00233656 z^2)/(1 +
      0.38387160 z + 0.04750621 z^2 + 0.00233656 z^3)];
S0PadeOverR[r2_] := Sqrt[-a1/(2 b)] Sqrt[-a1/L1] fPade[-r2 a1/L1];
S0Pade[r_] := r S0PadeOverR[r^2];
gPade[r_] = D[D[S0Pade[r], r] + S0Pade[r]/r, r] // Simplify;  (* d_i d_j Q_ij = gPade(r) cos(phi) *)
gFn[r_] = gPade[r] /. padeParams // Evaluate;

(* kernel and source *)
phiS[t_] := If[t < 1.*^-6, 1 - t/2 + t^2/6, (1 - Exp[-t])/t];
kappa[wv_, lN_] := phiS[2 wv^2/lN^2]/(Pi lN^2);
src[xv_, yv_] := With[{r = Sqrt[xv^2 + yv^2]}, If[r < 1.*^-12, 0., gFn[r] xv/r]];

(* nu on the +x axis (the trough; the Pade profile is antisymmetric, so -nu there is nu0) *)
nuAt[x0_?NumericQ, lN_?NumericQ] := Module[{pts},
   pts = Union[{0., lN, 5. lN, x0, 3., 30., 400.}];
   (lN^2/4) NIntegrate[wv kappa[wv, lN] src[x0 - wv Cos[p], -wv Sin[p]],
     Evaluate@Prepend[pts, wv], {p, 0, 2 Pi},
     Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
     MaxRecursion -> 30, PrecisionGoal -> 7, AccuracyGoal -> 14]];

(* coarse scan, then two parabolic refinements (FindMaximum stalls on the
   integration noise); returns {lN, nu0, position of the extremum on the axis} *)
parabolaPeak[f_, xc_, h_] := Module[{fm = f[xc - h], f0 = f[xc], fp = f[xc + h], d},
   d = fm - 2 f0 + fp;
   If[d >= 0, {xc, f0}, With[{xv = xc + h (fm - fp)/(2 d)}, {xv, f[xv]}]]];
nu0At[lN_] := Module[{scan, x1, r1, r2, f},
   f = -nuAt[#, lN] &;
   scan = Table[{x, f[x]}, {x, 0.6, 3.2, 0.2}];
   x1 = First@MaximalBy[scan, Last];
   r1 = parabolaPeak[f, x1[[1]], 0.1];
   r2 = parabolaPeak[f, r1[[1]], 0.025];
   {lN, r2[[2]], r2[[1]]}];

cacheFile = FileNameJoin[{figDir, "ln_curve_data.m"}];
recompute = Length[$ScriptCommandLine] >= 3 && $ScriptCommandLine[[3]] === "recompute";
If[FileExistsQ[cacheFile] && ! recompute,
  Module[{c = Get[cacheFile]},
   {gPoints, g05, aCoef, cCoef} = Lookup[c, {"points", "nu05", "a", "c"}];
   Print["using cached results from ", cacheFile]],
  lNs = Exp /@ Subdivide[Log[0.0001], Log[3], 20];
  gPoints = Table[With[{r = nu0At[l]}, Print[r]; r], {l, lNs}];
  g05 = nu0At[0.5][[2]];
  (* small-lN asymptote: a = max|d_i d_j Q_ij|/4, c fitted on lN <= 1e-3 um *)
  aCoef = First@FindMaximum[{Abs[gFn[r]], 0.3 < r < 4}, {r, 1.2}]/4;
  cCoef = Mean[(#[[2]]/#[[1]]^2 - aCoef Log[ld/#[[1]]]) & /@ Select[gPoints, #[[1]] <= 0.001 &]];
  Export[cacheFile, <|"points" -> gPoints, "nu05" -> g05, "a" -> aCoef, "c" -> cCoef, "ld" -> ld|>]];
asym[lN_] := lN^2 (aCoef Log[ld/lN] + cCoef);
Print["nu0(0.5) = ", g05, "   asymptote: a = ", aCoef, "  c = ", cCoef];

nuLog = Interpolation[Log[gPoints[[All, ;; 2]]]];
nuFn[lN_] := Exp[nuLog[Log[lN]]];

(* log ticks without trailing zeros: decimals down to decLimit, powers of ten below it;
   unlabelled minor ticks at 2..9 per decade; the top and right frames carry the same
   ticks without labels *)
tickLabel[k_, decLimit_] := Which[k >= 0, ToString[10^k], k >= decLimit,
   "0." <> StringRepeat["0", -k - 1] <> "1", True, Superscript[10, k]];
logTicks[kmin_, kmax_, decLimit_, labelled_ : True] := Join[
   Table[{10.^k, If[labelled, tickLabel[k, decLimit], ""], {0.02, 0}}, {k, kmin, kmax}],
   Flatten[Table[{m 10.^k, "", {0.01, 0}}, {k, kmin - 1, kmax}, {m, 2, 9}], 1]];
xTicks = logTicks[-3, 0, -3];     (* 0.001, 0.01, 0.1, 1 *)
yTicks = logTicks[-6, -1, -2];    (* 10^-6 ... 10^-3, 0.01, 0.1 *)
frameTicks = {{yTicks, logTicks[-6, -1, -2, False]}, {xTicks, logTicks[-3, 0, -3, False]}};

Block[{lineStyle, aStr, cStr},
 lineStyle = {Directive[Black, Dashed, Thickness[0.01]],
   Directive[colors[[1]], Thickness[0.008]]};
 aStr = ToString@NumberForm[aCoef, {3, 3}];
 cStr = ToString@NumberForm[cCoef, {3, 3}];
 fig = Show[
   LogLogPlot[nuFn[lN], {lN, 0.001, 3},
    PlotStyle -> lineStyle[[2]],
    Frame -> True,
    FrameStyle -> Directive[Black, FontFamily -> "Times", FontSize -> 8],
    FrameLabel -> {tex["\\ell_{\\rm n}~(\\upmu \\rm{m})"],
      tex["\\text{maximum variation }~ \\nu_0"]},
    GridLines -> {{0.5}, {g05}},
    FrameTicks -> frameTicks,
    PlotRange -> All],
   LogLogPlot[asym[lN], {lN, 0.001, 1},
    PlotStyle -> lineStyle[[1]],
    PlotLegends -> Placed[LineLegend[
       {Directive[Black, Dashed, Thickness[0.005]]},
       {tex["\\ell_{\\rm n}^2\\left[" <> aStr <> "\\,\\ln(\\ell_{\\rm d}/\\ell_{\\rm n})+" <>
          cStr <> "\\right]", 8]},
       LabelStyle -> {FontSize -> 8, FontFamily -> "Times"},
       LegendMarkerSize -> 30,
       Spacings -> {0, 0.1}], {0.62, 0.2}]],
   ImageSize -> w
   ]];

export["ln_curve", fig];
Print["exported to ", figDir];
