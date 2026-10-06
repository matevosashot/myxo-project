#!/usr/bin/env wolframscript
(* ::Package:: *)

(* cutoff_calibration.wls -- fig:cutoff-calibration of the manuscript (subsection "The noise
   correlation length and its calibration", app:lyapunov-cutoff), 2026-10-06 version.
   Moved here from claude_experiments/cutoff_fig/ell_eff_ratio.wls. It replaces the pedestal
   plot of claude_experiments/cutoff_fig/cutoff_calibration_v2.wls, whose solver checks the
   text quotes (within 0.1%).

   The calibrated length solves (manuscript eq:ell-eff)
       I(l_eff/dx) = dx/(Sqrt[2 Pi] l_n),   I(s) = Exp[-s^2] BesselI[0, s^2].
   With t = l_eff/dx as the parameter the curve is explicit:
       dx/l_n = Sqrt[2 Pi] I(t),   l_eff/l_n = Sqrt[2 Pi] t I(t).
   Fine grids: l_eff/l_n ~ 1 + (dx/l_n)^2/8.  No solution for dx/l_n > Sqrt[2 Pi].
   Top axis: interior nodes per direction for the box [-20,20]^2 and l_n = 0.5 um,
   dx = 40/(N+1), i.e. dx/l_n = 80/(N+1).

   Output: cutoff_calibration.pdf and .png next to this script.
   Run:  ./run.sh, which also copies them to manuscript/figures/. *)

$HistoryLength = 0;
If[!StringQ[Environment["QT_QPA_PLATFORM"]],
   SetEnvironment["QT_QPA_PLATFORM" -> "offscreen"]];
thisDir  = DirectoryName[$InputFileName];
figDir   = thisDir;

scaledI0[x_?NumericQ] := If[x <= 50., Exp[-N[x]] BesselI[0, N[x]],
   With[{xp = SetPrecision[x, 40]}, N[Exp[-xp] BesselI[0, xp]]]];
cI[t_?NumericQ] := scaledI0[t^2];
xOf[t_] := Sqrt[2 Pi] cI[t];          (* dx/l_n *)
lamOf[t_] := Sqrt[2 Pi] t cI[t];      (* l_eff/l_n *)

pts = Table[{xOf[t], lamOf[t]}, {t, Join[Subdivide[0.02, 3., 300], Subdivide[3., 400., 600]]}];
pts = Select[pts, 0.005 <= #[[1]] <= 2.6 &];

(* production grid: box 20, N_int = 150, l_n = 0.5 *)
xProd = (40./151)/0.5;
tProd = t /. FindRoot[xOf[t] == xProd, {t, 1/xProd}];
lamProd = lamOf[tProd];
lamMax = Max[pts[[All, 2]]];
xAtMax = SelectFirst[pts, #[[2]] == lamMax &][[1]];
Print["production: dx/l_n = ", xProd, "  l_eff/l_n = ", lamProd, "  (pedestal overshoot without",
  " calibration 2 Pi s^2 I(s)^2 with s = l_n/dx: ", 2 Pi (1/xProd)^2 cI[1/xProd]^2, ")"];
tOne = t /. FindRoot[lamOf[t] == 1, {t, 0.45}];
Print["max l_eff/l_n = ", lamMax, " at dx/l_n = ", xAtMax, ";  l_eff/l_n = 1 again at dx/l_n = ",
  xOf[tOne]];
tFine = t /. FindRoot[xOf[t] == 0.05, {t, 20.}];
Print["fine-grid check: dx/l_n = 0.05 -> l_eff/l_n = ", lamOf[tFine], "  vs 1 + x^2/8 = ", 1 + 0.05^2/8];

base  = {FontFamily -> "Times", FontSize -> 11};
blue  = RGBColor[0.20, 0.40, 0.70];
orange = RGBColor[0.85, 0.37, 0.01];
italic[x_] := Style[x, Italic];
lblX = Row[{"\[CapitalDelta]", italic["x"], "/", Style["\[ScriptL]", Italic], Subscript["", "n"]}];
lblY = Row[{Style["\[ScriptL]", Italic], Subscript["", "eff"], "/", Style["\[ScriptL]", Italic],
   Subscript["", "n"]}];
lblA = Row[{"1 + (\[CapitalDelta]", italic["x"], "/", Style["\[ScriptL]", Italic],
   Subscript["", "n"], ")", Superscript["", 2], "/8"}];

nTicks = {50, 150, 500, 1500, 5000};
topTicks = {Log[80./(# + 1)], ToString[#]} & /@ nTicks;
botTicks = {Log[#], ToString[#]} & /@ {0.01, 0.1, 1};
minorBot = {Log[#], "", {0.006, 0}} & /@ Flatten[Table[k 10^e, {e, -2, 0}, {k, 2, 9}]];

Print["curve points: ", Length[pts], "  dx/l_n range ", MinMax[pts[[All, 1]]]];
curve = ListLogLinearPlot[SortBy[pts, First], Joined -> True, PlotRange -> All,
   PlotStyle -> Directive[blue, AbsoluteThickness[1.4]]];
asym  = LogLinearPlot[1 + x^2/8, {x, 0.005, 1.25},
   PlotStyle -> Directive[Black, Dashed, AbsoluteThickness[0.9]]];
one   = LogLinearPlot[1, {x, 0.005, 2.6}, PlotStyle -> Directive[GrayLevel[0.6], AbsoluteThickness[0.6]]];
edge  = Graphics[{GrayLevel[0.55], Dotted, AbsoluteThickness[0.9],
    Line[{{Log[Sqrt[2 Pi]], 0}, {Log[Sqrt[2 Pi]], 1.3}}]}];
prod  = Graphics[{orange, Disk[{Log[xProd], lamProd}, Offset[{3.5, 3.5}]]}];

fig = Show[one, edge, asym, curve, prod,
   Frame -> True, Axes -> False,
   PlotRange -> {{Log[0.005], Log[2.9]}, {0, 1.25}},
   FrameTicks -> {{Automatic, Automatic}, {Join[botTicks, minorBot], topTicks}},
   FrameLabel -> {{lblY, None}, {lblX, Row[{italic["N"], Subscript["", "int"]}]}},
   FrameStyle -> Directive[Black, AbsoluteThickness[0.8]],
   FrameTicksStyle -> base, LabelStyle -> base,
   ImageSize -> 300, ImagePadding -> {{48, 10}, {40, 34}}];
fig = Legended[fig, Placed[LineLegend[{Directive[blue, AbsoluteThickness[1.4]],
     Directive[Black, Dashed, AbsoluteThickness[0.9]]}, {"exact", lblA},
    LabelStyle -> base, LegendMarkerSize -> {22, 10}], {0.32, 0.35}]];

Export[FileNameJoin[{figDir, "cutoff_calibration.pdf"}], fig];
Export[FileNameJoin[{figDir, "cutoff_calibration.png"}], fig, ImageResolution -> 300];
Print["wrote ", FileNameJoin[{figDir, "cutoff_calibration.pdf"}]];
