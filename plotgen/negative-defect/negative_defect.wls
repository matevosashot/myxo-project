#!/usr/bin/env wolframscript
(* ::Package:: *)

(* negative_defect.wls -- analytic density variance and layering probability at a -1/2
   defect, compared with the +1/2 defect (appendix app:negative, negative_defect.tex).
   Written 2026-10-05.

   Both defects share the radial profile S(r) (eq:S-radial depends only on |m| = 1), and
   q = Q_xx + i Q_xy = S(r) Exp[i m phi], m = +1 (+1/2) or m = -1 (-1/2). Then
      d_i d_j Q_ij = g1(r) Cos[phi],    g1 = d/dr (S' + S/r)            (+1/2)
      d_i d_j Q_ij = g3(r) Cos[3 phi],  g3 = h' - 2 h/r, h = S' - S/r   (-1/2)
   The isotropic kernel kappa of eq:kappa-real keeps the angular harmonic, so
      Var/pedestal - 1 = (lN^2/4) H_n(r) Cos[n phi],
      H_n(r) = Int r' dr' g_n(r') Int dalpha kappa(|r - r'|) Cos[n alpha].
   H_n is computed on a log grid as in getRhoVariance of plotting_small_core+ln_2.nb, and
   the peaks are cross-checked with the adaptive integral of ln_curve_adaptive.wls.

   Usage:  module load mathematica/14.10; unset DISPLAY
           wolframscript -file negative-defect/negative_defect.wls [noFigures] [draft] [only=name]
   Outputs: negative-defect/numbers.txt, figures varrho3D_minus, varrho3D_pm, layering_pm
   (pdf, png, _high.png) in <repoRoot>/figures; "out=<dir>" writes all of them to <dir>.
   2026-11-01 (the author's request: the axis labels sit too far from the axes): in layering_pm the
   labels are placed by plotgen/frame_labels.wl, about 4 pt (printed) from the tick labels,
   with an explicit ImagePadding; the frame grows at the same width.
   2026-11-02 (the author's rule: no line thinner than about 1 px; a line meant to look
   thinner is gray): in layering_pm frame and
   ticks print at 0.75 pt (0.17 pt before), the grid line gray at 0.75 pt (frame_labels.wl);
   the legend, which touched the right frame, anchored in the upper right corner.
   2026-11-02, the author's later request: in layering_pm the grid line prints at 0.6 pt in a
   darker gray (gridStyle of frame_labels.wl). *)

SetEnvironment["QT_QPA_PLATFORM" -> "offscreen"];
SetEnvironment["GALLIUM_DRIVER" -> "softpipe"];   (* llvmpipe crashes the front end on 3D export (Illegal instruction) *)

<< MaTeX`
w = 8.6*72/2.54;
tex[s_, fs_ : 11] :=
  MaTeX[s, FontSize -> fs,
   "Preamble" -> {"\\usepackage{amsmath,siunitx,newtxtext,newtxmath}"}];
colors = ColorData[97, "ColorList"];

repoRoot = If[$MachineName == "AMatevosyanL1M",
  "/Users/ashmat/Desktop/fluctuation-solver",
  "/home/ashmat/Projects/myxo-project-dean"];
figDir = FileNameJoin[{repoRoot, "figures"}];
outDir = FileNameJoin[{repoRoot, "negative-defect"}];
Scan[If[StringStartsQ[#, "out="], figDir = outDir = StringDrop[#, 4]] &, $ScriptCommandLine];
Get[FileNameJoin[{repoRoot, "plotgen", "frame_labels.wl"}]];
doFigures = ! MemberQ[$ScriptCommandLine, "noFigures"];
(* smoke tests: "draft" = coarse plots, no numbers block, 150 dpi png only, into
   negative-defect/draft; "only=name" = build only that figure *)
draft = MemberQ[$ScriptCommandLine, "draft"];
onlyFig = StringDrop[#, 5] & /@ Select[$ScriptCommandLine, StringStartsQ[#, "only="] &];
want[name_] := onlyFig === {} || MemberQ[onlyFig, name];
{plotPts, maxRec} = If[draft, {30, 0}, {80, 2}];

Get[FileNameJoin[{repoRoot, "binaryReadTools.m"}]];
Import[FileNameJoin[{repoRoot, "RemoteExport.m"}]];
(* pdf = False for the 3D figures: png only *)
export[name_, fig_, pdf_ : True] := If[draft,
  Export[FileNameJoin[{outDir, "draft", name <> ".png"}], fig, ImageResolution -> 150],
  If[pdf, RemoteExport[FileNameJoin[{figDir, name <> ".pdf"}], fig]];
  Export[FileNameJoin[{figDir, name <> "_high.png"}], fig, ImageResolution -> 900];
  Export[FileNameJoin[{figDir, name <> ".png"}], fig, ImageResolution -> 600]];

(* 2026-10-26: parameters from the N200 ssDir production run (same values as the N180 run used
   before); pcrit from eq:pcrit, Ws - Wc + (gamma_w + Wc) h Sqrt[Pi/Acell] = 35258 Pa um, the
   value of the manuscript (it was hard-coded as 35000) *)
dataPath = FileNameJoin[{repoRoot,
   "data/fluctuations-set2_2_comove_optimized_2_ssDir/covariance_box20_N200_Dir_lN0.5_z12.5_ld1_B90000_ss20.bin"}];
modelParams = readSteady[dataPath]["modelParams"];
{lN, Bv, zetav, rho0v} = {lNoise, B, \[Zeta], \[Rho]0} /. modelParams;
pcrit = 2000. + 65000.*0.5*Sqrt[Pi/3.];

(* accurate Pade approximant (eq:pade-accurate), within 1e-4 S0 of the exact profile *)
fPade[z_] :=
  Sqrt[(0.34010999 + 0.04516965 z + 0.00233656 z^2)/(1 +
      0.38387160 z + 0.04750621 z^2 + 0.00233656 z^3)];
S0PadeOverR[r2_] := Sqrt[-a1/(2 b)] Sqrt[-a1/L1] fPade[-r2 a1/L1];
S0Pade[r_] := r S0PadeOverR[r^2];
Sfn[r_] = S0Pade[r] /. modelParams;
S0v = Sqrt[-a1/(2 b)] /. modelParams;
ldv = Sqrt[-L1/a1] /. modelParams;

g1[r_] = D[D[Sfn[r], r] + Sfn[r]/r, r] // Simplify;
hh[r_] = D[Sfn[r], r] - Sfn[r]/r;
g3[r_] = D[hh[r], r] - 2 hh[r]/r // Simplify;
gn[1] = g1; gn[3] = g3;

(* check the double divergence against the Cartesian components *)
Module[{Q1p, Q2p, Q1m, Q2m, x, y, ddp, ddm, pts},
  Q1p[x_, y_] := Sfn[Sqrt[x^2 + y^2]] x/Sqrt[x^2 + y^2];      (* S Cos[phi] *)
  Q2p[x_, y_] := Sfn[Sqrt[x^2 + y^2]] y/Sqrt[x^2 + y^2];      (* S Sin[phi] *)
  Q1m = Q1p; Q2m[x_, y_] := -Q2p[x, y];                         (* q = S Exp[-i phi] *)
  ddp[x_, y_] = D[Q1p[x, y], x, x] - D[Q1p[x, y], y, y] + 2 D[Q2p[x, y], x, y];
  ddm[x_, y_] = D[Q1m[x, y], x, x] - D[Q1m[x, y], y, y] + 2 D[Q2m[x, y], x, y];
  pts = {{0.7, 0.3}, {-1.2, 2.1}, {3., -0.4}};
  Print["check +1/2: ", Max@Abs[(ddp @@ # - g1[Norm[#]] Cos[ArcTan @@ #]) & /@ pts]];
  Print["check -1/2: ", Max@Abs[(ddm @@ # - g3[Norm[#]] Cos[3 ArcTan @@ #]) & /@ pts]]];

(* radial filter H_n(r) on a log grid (getRhoVariance of the notebook, harmonic n) *)
phiVec[t_] := Module[{s = UnitStep[1.*^-6 - t], td}, td = t + s;
   (1 - s) ((1 - Exp[-td])/td) + s (1 - t/2 + t^2/6)];
filterH[n_, rs_List, nr_ : 4000, nA_ : 720, rmax_ : 400.] :=
  Module[{rp, rp2, gp, ca, wgt, dw, il2},
   rp = Exp[Subdivide[N@Log[0.0001], N@Log[rmax], nr - 1]]; rp2 = rp^2;
   gp = gn[n] /@ rp;
   ca = Cos[(Range[0, nA - 1] + 0.5) 2. Pi/nA];
   wgt = Cos[n ArcCos[ca]] 2. Pi/nA;
   dw = Differences[rp]; il2 = 1./(Pi lN^2);
   Table[Module[{w2, k1, integ},
     w2 = (r^2 + rp2) - 2. r Outer[Times, rp, ca];
     k1 = il2 phiVec[2. w2/lN^2] . wgt;
     integ = rp gp k1;
     Total[(Most[integ] + Rest[integ])/2 dw]], {r, rs}]];

hRadii = Join[{0.}, Exp[Subdivide[N@Log[0.01], N@Log[48.], 119]]];
Hn[n_] := Hn[n] = Interpolation[Transpose[{hRadii, Join[{0.}, filterH[n, Rest[hRadii]]]}],
    InterpolationOrder -> 3];
(* nu_n(x, y) = Var/pedestal - 1 *)
nuPlus[x_, y_] := With[{r = Sqrt[x^2 + y^2]}, If[r == 0., 0., lN^2/4 Hn[1][r] x/r]];
nuMinus[x_, y_] := With[{r = Sqrt[x^2 + y^2]},
   If[r == 0., 0., lN^2/4 Hn[3][r] Cos[3 ArcTan[x, y]]]];

(* adaptive cross-check of the profile on an axis (ln_curve_adaptive.wls) *)
kappa[wv_] := phiVec[2 wv^2/lN^2]/(Pi lN^2);
src[n_][xv_, yv_] := With[{r = Sqrt[xv^2 + yv^2]},
   If[r < 1.*^-12, 0., gn[n][r] Cos[n ArcTan[xv, yv]]]];
nuAdaptive[n_, x0_?NumericQ] := Module[{pts = Union[{0., lN, 5. lN, Abs@x0, 3., 30., 400.}]},
   (lN^2/4) NIntegrate[wv kappa[wv] src[n][x0 - wv Cos[p], -wv Sin[p]],
     Evaluate@Prepend[pts, wv], {p, 0, 2 Pi},
     Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
     MaxRecursion -> 30, PrecisionGoal -> 7, AccuracyGoal -> 14]];

(* peaks: +1/2 on the -x axis, -1/2 on the +x axis (phi = 0, a raised lobe) *)
{nuP, xP} = Module[{s = FindMaximum[nuPlus[-x, 0], {x, 1.5, 0.5, 4}]}, {s[[1]], x /. s[[2]]}];
{nuM, xM} = Module[{s = FindMaximum[nuMinus[x, 0], {x, 3., 1, 6}]}, {s[[1]], x /. s[[2]]}];
If[! draft, nuPad = nuAdaptive[1, -xP]; nuMad = nuAdaptive[3, xM]];

(* layering: Gaussian tail of eq:p-layer, relative to the far field *)
pedestal = zetav/(2 Pi Bv rho0v lN^2);
zInf = pcrit/Sqrt[Bv^2 pedestal];
Qapprox[z_] := Exp[-z^2/2]/(z Sqrt[2 Pi]);
ratio[nu_] := Qapprox[zInf/Sqrt[1 + nu]]/Qapprox[zInf];
ratioPeak[nu_] := Exp[zInf^2/2 nu/(1 + nu)];          (* eq:R-value *)
circleAvg[f_, r_?NumericQ, nPts_ : 720] :=
  Mean@Table[ratio[f[r Cos[p], r Sin[p]]], {p, 2. Pi/nPts (Range[nPts] - 0.5)}];
annulusAvg[f_, {a_, b_}] := Module[{rs = a + (b - a) (Range[100] - 0.5)/100},
   Total[rs (circleAvg[f, #, 180] & /@ rs)] (b - a)/100/((b^2 - a^2)/2)];

rBins = {0.3, 1.3, 3., 6., 11., 18.};
edges = {{0, 0.5}, {0.5, 2}, {2, 4}, {4, 8}, {8, 14}};

If[! draft,
numbers = StringRiffle[{
    "params: B=" <> ToString[Bv] <> " zeta=" <> ToString[zetav] <> " rho0=" <>
     ToString[rho0v] <> " lN=" <> ToString[lN] <> " ld=" <> ToString[ldv] <> " S0=" <>
     ToString[S0v] <> " pcrit=" <> ToString[pcrit],
    "zInf = " <> ToString[zInf] <> ",  zInf^2/2 = " <> ToString[zInf^2/2],
    "+1/2: nu0 = " <> ToString[100 nuP] <> "% at x = -" <> ToString[xP] <>
     " um (adaptive " <> ToString[100 nuPad] <> "%), peak ratio " <> ToString[ratioPeak[nuP]],
    "-1/2: nu0 = " <> ToString[100 nuM] <> "% at r = " <> ToString[xM] <>
     " um, phi = 0, 120, 240 (adaptive " <> ToString[100 nuMad] <> "%), peak ratio " <>
     ToString[ratioPeak[nuM]],
    "-1/2 trough at phi = 60: " <> ToString[100 nuMinus[xM Cos[Pi/3], xM Sin[Pi/3]]] <> "%",
    "nu0 ratio (-)/(+) = " <> ToString[nuM/nuP] <> ",  log-enhancement ratio = " <>
     ToString[Log[ratioPeak[nuM]]/Log[ratioPeak[nuP]]],
    "source near core: g1(r)/r at r=0.05: " <> ToString[g1[0.05]/0.05] <>
     ",  g3(r)/r^3 at r=0.05: " <> ToString[g3[0.05]/0.05^3],
    "source far field: g1 r^2 at r=40: " <> ToString[g1[40.] 40.^2] <>
     " (-S0 = " <> ToString[-S0v] <> "),  g3 r^2 at r=40: " <> ToString[g3[40.] 40.^2] <>
     " (3 S0 = " <> ToString[3 S0v] <> ")",
    "max |g1| = " <> ToString[NMaxValue[{Abs@g1[r], 0.01 < r < 20}, r]] <>
     ",  max |g3| = " <> ToString[NMaxValue[{Abs@g3[r], 0.01 < r < 20}, r]],
    "nu at r = 1 um along peak direction: +1/2 " <> ToString[100 nuPlus[-1., 0]] <>
     "%,  -1/2 " <> ToString[100 nuMinus[1., 0]] <> "%",
    "core-centred circle average P(r)/P(inf) at the Copenhagen bin centres " <>
     ToString[rBins] <> ":",
    "  +1/2: " <> ToString[NumberForm[circleAvg[nuPlus, #], 3] & /@ rBins],
    "  -1/2: " <> ToString[NumberForm[circleAvg[nuMinus, #], 3] & /@ rBins],
    "annulus averages, edges " <> ToString[edges] <> ":",
    "  +1/2: " <> ToString[NumberForm[annulusAvg[nuPlus, #], 3] & /@ edges],
    "  -1/2: " <> ToString[NumberForm[annulusAvg[nuMinus, #], 3] & /@ edges],
    "max over r of circle average: +1/2 " <>
     ToString[NMaxValue[{circleAvg[nuPlus, r, 360], 0.2 < r < 8}, r]] <> ",  -1/2 " <>
     ToString[NMaxValue[{circleAvg[nuMinus, r, 360], 0.2 < r < 8}, r]]
    }, "\n"];
Print[numbers];
Export[FileNameJoin[{outDir, "numbers.txt"}], numbers, "Text"]];

If[! doFigures, Exit[]];

(* ---------------------------------------------------------------- figures *)
directorPlus[x_, y_] := With[{p = Pi - ArcTan[-x, y]}, {Cos[p/2], Sin[p/2]}];
directorMinus[x_, y_] := With[{p = Pi - ArcTan[-x, y]}, {Cos[-p/2], Sin[-p/2]}];
(* director streamlines, lifted onto the surface as 3D curves *)
streamLines[dir_, W_] := Cases[Normal@StreamPlot[dir[x, y], {x, -W, W}, {y, -W, W},
     StreamScale -> None, StreamMarkers -> "Line", StreamPoints -> Medium],
    (Line | Arrow)[pts_] :> pts, Infinity];
lifted[nu_, dir_, box_] := Module[{segs},
   segs = Select[#, Norm[#] < box &] & /@ streamLines[dir, box];
   segs = Select[segs, Length[#] > 2 &];
   {Directive[Black, AbsoluteThickness[1.2]],   (* 2026-10-26: was Opacity[0.6], 0.6 *)
    Line[Map[Append[#, 1 + 1.0005 nu @@ #] &, segs, {2}]]}];
zRange = {0.97, 1.03};   (* covers 1 +- nuP = 1 +- 0.0294 *)
surface[nu_, dir_, box_, zr_, zTicks_, opts___] :=
  Show[Plot3D[1 + nu[x, y], {x, y} \[Element] Disk[{0, 0}, box], opts, ImageSize -> w,
    AxesLabel -> {tex["x~(\\upmu \\rm{m})"], tex["y~(\\upmu \\rm{m})"],
      Rotate[tex["\\text{Var}[\\rho]~/~\\text{Var}[\\rho]_\\text{pedestal}", 10], 95 Degree]},
    LabelStyle -> Directive[Black, FontFamily -> "Times", FontSize -> 10],
    Mesh -> None, Ticks -> {Range[-20, 20, 4], Range[-20, 20, 4], zTicks},
    PlotRange -> {Automatic, Automatic, zr}, PlotPoints -> plotPts, MaxRecursion -> maxRec,
    ColorFunction -> (ColorData["TemperatureMap"][#3] &),
    Lighting -> "Neutral", BoxRatios -> {1, 1, 0.45},
    SphericalRegion -> False],   (* else Grid/Labeled pad the panel to a bounding sphere *)
   Graphics3D[lifted[nu, dir, box]]];

If[want["varrho3D_minus"],
 fig = surface[nuMinus, directorMinus, 8, All, {0.99, 1, 1.01}, ViewPoint -> {1.2, -2.0, 2.4}];
 export["varrho3D_minus", fig, False]];

(* both panels on one absolute colour scale over zRange, shared by the colorbar *)
jointColors = Sequence[ColorFunction -> (ColorData["TemperatureMap"][Rescale[#3, zRange]] &),
   ColorFunctionScaling -> False];
colorBar = BarLegend[{"TemperatureMap", zRange}, "Ticks" -> {0.97, 0.98, 0.99, 1, 1.01, 1.02, 1.03},
   LegendLayout -> "Column", LegendMarkerSize -> {8, 0.55 w},
   LabelStyle -> Directive[Black, FontFamily -> "Times", FontSize -> 10]];
If[want["varrho3D_pm"],
fig = Grid[{{
     Labeled[surface[nuPlus, directorPlus, 6, zRange, {0.98, 1, 1.02}, jointColors,
       ViewPoint -> {1.2, -2.0, 2.4}], tex["\\text{(a)}~+1/2", 13], Top],
     Labeled[surface[nuMinus, directorMinus, 6, zRange, {0.98, 1, 1.02}, jointColors,
       ViewPoint -> {1.2, -2.0, 2.4}], tex["\\text{(b)}~-1/2", 13], Top],
     colorBar}}, Spacings -> {0, 0}, Alignment -> Center];
export["varrho3D_pm", fig, False]];

If[want["layering_pm"],
dataExperiment = Apply[Function[{r, low, mean, high},
     {r, Around[mean, {mean - low, high - mean}]}],
    Rest@Import[FileNameJoin[{repoRoot, "top_plot.csv"}]], 1];
(* 2026-10-26: the core-centred +1/2 curve takes the orange of the core-centred curve of
   Fig. 2(c) (colors[[2]]), and the -1/2 curve a colour Fig. 2(c) does not use (colors[[3]]);
   labels 8 -> 10 pt, error bars 0.004 -> 0.006.
   2026-11-01, the author's request for one font and one set of sizes in all figures: all text
   of this figure in MaTeX Computer Modern (newtx and plain Times until then; tex, newtx, stays
   for the 3D figures above), ticks included; tick and axis labels print at 9 pt, the legend at
   8 pt, at 0.85 \linewidth (246 pt) *)
{fsT, fsL} = {9, 8} w/(0.85*246.);
texCM[str_, fs_ : fsT] := MaTeX[str, FontSize -> fs,
   "Preamble" -> {"\\usepackage{amsmath,siunitx,upgreek}"}];
yMajor = {1, 10, 100, 1000};
yMinor = Complement[Flatten@Table[m 10^k, {k, 0, 2}, {m, 2, 9}], yMajor];
xMinor = Complement[Range[0, 19], Range[0, 15, 5]];
layerTicks = {{Join[texTicks[yMajor, texCM], minorTicks[yMinor]],
    Join[{#, ""} & /@ yMajor, minorTicks[yMinor]]},
   {Join[texTicks[Range[0, 15, 5], texCM], minorTicks[xMinor]],
    Join[{#, ""} & /@ Range[0, 15, 5], minorTicks[xMinor]]}};
fig = Show[
   LogPlot[{circleAvg[nuPlus, r, 360], circleAvg[nuMinus, r, 360]}, {r, 0.05, 19},
    PlotRange -> {All, {0.7, 1000}}, PlotPoints -> 60, GridLines -> {None, {1}},
    GridLinesStyle -> gridStyle[w/(0.85*246.)],
    Frame -> True, FrameStyle -> frameStyle[w/(0.85*246.)], FrameTicks -> layerTicks,
    Epilog -> frameLabels[{texCM["\\text{Distance from the core}~~r~(\\upmu \\rm{m})"],
       texCM["\\text{Layering probability}~~\\mathcal P(r) / \\mathcal P(\\infty)"]}, {15.8, 25.4}],
    PlotStyle -> {Directive[colors[[2]], Thickness[0.006]],
      Directive[colors[[3]], Thickness[0.006]]},
    PlotLegends -> Placed[LineLegend[{Directive[colors[[2]], Thickness[0.004]],
        Directive[colors[[3]], Thickness[0.004]], Red},
       {texCM["+1/2~\\text{defect}", fsL], texCM["-1/2~\\text{defect}", fsL],
        texCM["\\text{Experiment},~+1/2", fsL]},
       LegendMarkers -> {None, None,
         Graphics[{Red, Disk[{0, 0}, 0.3]}, PlotRange -> {{-1, 1}, {-1, 1}}]},
       LegendMarkerSize -> {15, 9}, Spacings -> {0.3, 0.1}, LegendMargins -> 0],
      {{0.98, 0.97}, {1, 1}}]],   (* Scaled[{0.75, 0.8}] until 2026-11-02: it ran into the frame *)
   ListLogPlot[dataExperiment, PlotStyle -> Directive[Red, PointSize[0.015]],
    IntervalMarkersStyle -> Directive[Red, Thickness[0.006]]],
   ImageSize -> w, ImagePadding -> {{38.7, 6}, {29, 17}},   (* top: the y label is longer than the frame *)
   PlotRangeClipping -> False];   (* the labels lie outside the frame; the curves stay inside *)
export["layering_pm", fig]];
Print["figures exported to ", figDir];
