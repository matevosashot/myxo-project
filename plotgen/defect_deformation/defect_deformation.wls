#!/usr/bin/env wolframscript
(* defect_deformation.wls -- the three panels of fig:comoving-deformation (Fig. 12), written
   2026-11-02 on the author's request for a script that regenerates the figure in the style
   of the other figures:
     strength.pdf     (a) the nematic order S along +x, -x and +-y, with the Pade profile
     defectangle.pdf  (b) the director angle theta - phi/2 on circles r = 0.1, 1, 2 um
     activeforce.pdf  (c) (div Q)_x along the axis y = 0, numerical and Pade

   Source: the notebook remote_kernel/plotting_small_core+ln_2.nb (the copy of 2026-10-26;
   the older copy in negative-defect/ has no such cells), whose cell 33 (fig1, fig2, fig3)
   and cell 34 (export["strength.pdf", fig1], ...) made the pdfs in manuscript/figures; this
   script runs them, with the style changes listed below:
     tex, w, colors, ArcTan2           cell 4 (MaTeX in Computer Modern)
     modelParams, data                 readSteady of the ssDir N_int = 180 run (B = 9e4),
                                       solSteadyDirichletSaved: Dirichlet density wall, free
                                       Q, the steady state of Figs. 5, 6 and 11 (cells 16, 17)
     Pade profile                      S0Pade, forcePade (cell 12)
     the panels                        cell 33: S = |Q| (Q1, Q2 the interpolated components);
                                       theta = ArcTan2[Q1, Q2]/2 on circles; the numerical
                                       (div Q)_x = d_x Q1 + d_y Q2 at y = 0
   "repro" runs cell 33 as it stands (panel options, sizes, legends), to check the data
   against the pdfs in manuscript/figures (outputs <name>_repro.pdf).
   Check of 2026-11-02: repro reproduces the three pdfs of manuscript/figures pixel for pixel
   (rendered at 300 dpi). The numbers printed for the caption: the +-x curves of (a) split by
   at most 10.9% of S0 (at r = 1.9 um); each departs from the Pade profile by at most 6.2%
   (+x) and 4.7% (-x); in (c) the numerical value at the core is 1.4% below the Pade one.

   Style, the same as in the figures restyled on the author's requests of 2026-11-01 and
   2026-11-02 (plotgen/frame_labels.wl):
   - the axis labels about 4 pt (printed) from the tick labels, placed by frameLabels with an
     explicit ImagePadding (the notebook placed them by hand, 6-11 pt away);
   - all text MaTeX Computer Modern, ticks 9 pt and axis labels 9 pt, legends 8 pt as printed
     (7.1, 10.7 and 7.1 pt before), at the width 0.32 \textwidth = 163.2 pt of the figure*
     (revtex reprint, \textwidth 510 pt); the design width stays the notebook's 3 w/4;
   - frame and ticks 0.75 pt (0.53 before); the lines at zero, black 0.45 pt axes in (b) and
     (c) and a gray 0.18 pt grid line under the axis in (c), are grid lines (gridStyle: 0.6 pt,
     gray 0.5) under the curves; Axes off;
   - curves keep their printed widths (colours 1.2 pt; the Pade dashes 2.4 pt in (a), 1.8 pt
     in (c)), now as AbsoluteThickness, so that the legend samples are as thick as the curves
     (a relative Thickness in a legend refers to the sample: they printed at 2.3 pt, and the
     dashed Pade sample of (a) at 3.8 pt as a single block);
   - legends: rows 9.5 pt apart, no margin, anchored inside a corner of the frame (the one of
     (b) ran past the right frame);
   - the three panels as one row: the same height and the same top and bottom padding, so
     the frames align; the left padding of each holds its tick labels and y label (wider in
     (b), for the degree ticks), so the frames differ in width (148, 133 and 144 pt design).

   Usage:  plotgen/defect_deformation/run.sh [repro] [out=<dir>]
   Outputs: strength, defectangle, activeforce (.pdf, .png at 600 dpi) in out=<dir>, default
   next to this script. Copies nothing into manuscript/figures. *)

SetEnvironment["QT_QPA_PLATFORM" -> "offscreen"];
$HistoryLength = 0;

here = DirectoryName[$InputFileName];
repoRoot = ParentDirectory[ParentDirectory[here]];
args = Rest[$ScriptCommandLine];
repro = MemberQ[args, "repro"];
opts = Association[Rule @@@ (StringSplit[#, "=", 2] & /@ Select[args, StringContainsQ[#, "="] &])];
figDir = Lookup[opts, "out", here];
If[! DirectoryQ[figDir], CreateDirectory[figDir]];
t0 = AbsoluteTime[];
stamp[] := ToString[Round[AbsoluteTime[] - t0]] <> " s";

(* ---------------------------------------------------------------- notebook set-up *)
<< MaTeX`
w = 8.6*72/2.54;
tex[s_, fs_ : 11] := MaTeX[s, FontSize -> fs,
   "Preamble" -> {"\\usepackage{amsmath,siunitx,upgreek}"}];
colors = ColorData[97, "ColorList"];
ArcTan2[x_, y_] := Pi - ArcTan[-x, y];

Import[FileNameJoin[{repoRoot, "RemoteExport.m"}]];
export[name_, fig_] := Block[{base = FileNameJoin[{figDir, name}]},
   RemoteExport[base <> ".pdf", fig];
   RemoteExport[base <> ".png", fig, ImageResolution -> 600];
   Print["exported ", base, ".pdf  (", stamp[], ")"]];
Get[FileNameJoin[{repoRoot, "plotgen", "frame_labels.wl"}]];

Get[FileNameJoin[{repoRoot, "binaryReadTools.m"}]];
dataPath = FileNameJoin[{repoRoot,
    "data/fluctuations-set2_2_comove_optimized_2_ssDir/covariance_box20_N180_Dir_lN0.5_z12.5_ld1_B90000_ss20.bin"}];
solSteadyDirichletSaved = readSteady[dataPath];
modelParams = solSteadyDirichletSaved["modelParams"];
Print["modelParams ", modelParams];

(* the Pade profile *)
fPade[z_] := Sqrt[(0.34010999 + 0.04516965*z + 0.00233656*z^2)/(1 + 0.3838716*z +
      0.04750621*z^2 + 0.00233656*z^3)];
S0PadeOverR[r2_] := Sqrt[-a1/(2*b)]*Sqrt[-a1/L1]*fPade[(-r2)*(a1/L1)];
S0Pade[r_] := r*S0PadeOverR[r^2];
forcePade[r_] = D[S0Pade[r], r] + S0Pade[r]/r;

(* the three quantities of cell 33 *)
Q1 = solSteadyDirichletSaved["Q1"]; Q2 = solSteadyDirichletSaved["Q2"];
Snum[x_, y_] := Sqrt[Q1[x, y]^2 + Q2[x, y]^2];
thetaNum[r_, \[Phi]_] := With[{x = r Cos[\[Phi]], y = r Sin[\[Phi]]}, ArcTan2[Q1[x, y], Q2[x, y]]/2];
forceNum[x_] := Derivative[1, 0][Q1][x, 0] + Derivative[0, 1][Q2][x, 0];
rVals = {0.1, 1, 2};

(* numbers the caption quotes: the split of the +-x curves in (a), the force at the core (c) *)
Print["(a) max |S(r,0) - S(-r,0)| / S0 on 0 < r < 6: ",
  NMaxValue[{Abs[Snum[r, 0] - Snum[-r, 0]], 0 < r < 6}, r]/(Sqrt[-a1/(2 b)] /. modelParams)];
Print["(c) at x = 0: numerical ", forceNum[0.], ", Pade ", forcePade[10.^-6] /. modelParams,
  ", ratio - 1 = ", forceNum[0.]/(forcePade[10.^-6] /. modelParams) - 1];

If[repro,
 (* ---------------------------------------------------------------- cell 33 as it stands *)
 Module[{fs = 8, labelfs = 12, panelOpts, fmt, ticks, legend, axisLabels, fig1, fig2, fig3},
  panelOpts = {ImageSize -> 3*(w/4), AspectRatio -> 0.8, ImagePadding -> {{42, 6}, {30, 6}},
    Frame -> True, FrameStyle -> Directive[Black, AbsoluteThickness[0.6]],
    FrameTicksStyle -> Directive[Black, FontFamily -> "Times", FontSize -> fs],
    PlotRangeClipping -> False};
  fmt[x_] := Which[x == 0, "0", IntegerQ[Rationalize[x]], ToString[Round[x]], True,
    ToString[NumberForm[N[x], {3, 1}]]];
  ticks[vals_] := ({#1, tex[fmt[#1], fs]} &) /@ vals;
  legend[styles_, labels_, pos_] := Placed[LineLegend[styles, (tex[#1, fs] &) /@ labels,
     LabelStyle -> {FontSize -> fs, FontFamily -> "Times"}, LegendMarkerSize -> 14,
     Spacings -> {0.4, -0.3}], pos];
  axisLabels[xlab_, dx_, ylab_, dy_] := Epilog -> {Text[tex[xlab, labelfs],
      Offset[{0, -dx}, Scaled[{0.5, 0}]]],
     Text[Rotate[tex[ylab, labelfs], Pi/2], Offset[{-dy, 0}, Scaled[{0, 0.5}]]]};
  fig1 = Block[{tk1 = 0.01, lineStyle},
    lineStyle = {Directive[Dashed, Black, Thickness[0.02]], Directive[Thickness[tk1], colors[[1]]],
      Directive[Thickness[tk1], colors[[2]]], Directive[Thickness[tk1], colors[[3]]]};
    Plot[{S0Pade[r] /. modelParams, Snum[r, 0], Snum[-r, 0], Snum[0, r]}, {r, 0, 6},
     Evaluate[Sequence @@ panelOpts], PlotRange -> {{0, 6}, {0, 1}},
     PlotRangePadding -> {None, {0, Scaled[0.02]}}, FrameLabel -> None,
     Evaluate[axisLabels["r~(\\upmu\\rm{m})", 20, "S", 22]],
     FrameTicks -> {{ticks[Range[0, 0.8, 0.2]], None}, {ticks[Range[0, 6, 2]], None}},
     PlotStyle -> lineStyle, PlotLegends -> legend[lineStyle,
       {"\\text{Pad\\'e}", "S(r,0)", "S(-r,0)", "S(0,\\pm r)"}, {0.72, 0.3}]]];
  fig2 = Block[{tk1 = 0.01, lineStyle, frameXTicks},
    frameXTicks = Apply[{#1*Degree, tex[#2, fs]} &, {{-1, "-1^\\circ"}, {-0.5, "-0.5^\\circ"},
       {0, "0^\\circ"}, {0.5, "0.5^\\circ"}, {1, "1^\\circ"}}, 2];
    lineStyle = Table[Directive[Thickness[tk1], colors[[i]]], {i, Length[rVals]}];
    Plot[Evaluate[Table[thetaNum[r, \[Phi]] - \[Phi]/2, {r, rVals}]], {\[Phi], 0.001, 2*Pi - 0.001},
     Evaluate[Sequence @@ panelOpts], PlotRange -> {{0, 2*Pi}, {-1.05, 1.05}*Degree},
     PlotRangePadding -> None, FrameLabel -> None,
     Evaluate[axisLabels["\\phi", 20, "\\theta-\\phi/2", 30]],
     FrameTicks -> {{frameXTicks, None}, {{{0, tex["0", fs]}, {Pi, tex["\\pi", fs]},
         {2*Pi, tex["2\\pi", fs]}}, None}}, PlotStyle -> lineStyle,
     PlotLegends -> legend[lineStyle, StringTemplate["r=`1`~\\upmu\\rm{m}"] /@ rVals, {0.78, 0.22}]]];
  fig3 = Block[{tk1 = 0.01, lineStyle},
    lineStyle = {Directive[Dashed, Black], Directive[Thickness[tk1], colors[[1]]]};
    Plot[{forcePade[x] /. modelParams, forceNum[x]}, {x, -3, 3},
     Evaluate[Sequence @@ panelOpts], PlotRange -> {{-3, 3}, {0, 1.2}}, PlotRangePadding -> None,
     FrameLabel -> None, Evaluate[axisLabels["x~(\\upmu\\rm{m})", 20,
       "(\\nabla\\cdot\\mathbf{Q})_x~(\\upmu\\rm{m}^{-1})", 28]],
     FrameTicks -> {{ticks[{0, 0.5, 1}], None}, {ticks[{-2, 0, 2}], None}},
     GridLines -> {{0}, None}, GridLinesStyle -> GrayLevel[0.75], PlotStyle -> lineStyle,
     PlotLegends -> legend[lineStyle, {"\\text{Pad\\'e}", "\\text{numerical}"}, {0.75, 0.15}]]];
  export["strength_repro", fig1]; export["defectangle_repro", fig2];
  export["activeforce_repro", fig3]];
 Print["done (", stamp[], ")"]; Exit[]];

(* ---------------------------------------------------------------- figures *)
(* sizes: design pt that print at pt when a panel (width wP) is included at 0.32 \textwidth *)
wP = 3 w/4;                               (* the notebook's design width *)
printed[pt_] := pt wP/(0.32 510.);
fP = printed[1];                          (* design pt per printed pt, 1.12 *)
{fsT, fsL} = printed[{9, 8}];             (* tick and axis labels 9 pt, legends 8 pt *)
texT[s_] := tex[s, fsT];
texL[s_] := tex[s, fsL];
curve[c_] := Directive[c, AbsoluteThickness[1.2 fP]];   (* 1.2 pt printed, as Thickness[0.01] *)

(* one row: the same height and vertical padding in all three panels, so the frames align;
   the left padding holds the tick labels and the y label *)
hP = 143.9;                               (* the notebook's height: 0.8 of the frame + padding *)
{padB, padT, padR} = {27.6, 5., 7.};
(* label distances (design pt) measured to give 4 pt between tick-label and label ink; the
   left padding leaves about 1 pt of white beside the y label *)
padLeft = <|"a" -> 28.2, "b" -> 42.9, "c" -> 32.2|>;
dxLab = <|"a" -> 15.2, "b" -> 14.75, "c" -> 15.2|>;
dyLab = <|"a" -> 15.6, "b" -> 29.8, "c" -> 17.85|>;
panelOpts[k_] := Sequence[ImageSize -> {wP, hP}, AspectRatio -> Full,
   ImagePadding -> {{padLeft[k], padR}, {padB, padT}}, Frame -> True,
   FrameStyle -> frameStyle[fP], Axes -> False, PlotRangeClipping -> False];
labels[k_, xlab_, ylab_] := frameLabels[texT /@ {xlab, ylab}, {dxLab[k], dyLab[k]}];
(* legends: samples as thick as the curves, no margin beyond the rows, rows 9.5 pt apart *)
legendAt[styles_, labs_, pos_] := Placed[LineLegend[styles, texL /@ labs,
    LegendMarkerSize -> {16, 9}, Spacings -> {0.3, -0.2}, LegendMargins -> 0], pos];

(* (a) the order; the dashed Pade profile 2.4 pt printed, as Thickness[0.02] of the frame *)
styles1 = {Directive[Dashed, Black, AbsoluteThickness[2.4 fP]], curve[colors[[1]]],
   curve[colors[[2]]], curve[colors[[3]]]};
fig1 = Plot[{S0Pade[r] /. modelParams, Snum[r, 0], Snum[-r, 0], Snum[0, r]}, {r, 0, 6},
   Evaluate[panelOpts["a"]], PlotRange -> {{0, 6}, {0, 1}},
   PlotRangePadding -> {None, {0, Scaled[0.02]}},
   FrameTicks -> {{texTicks[Range[0, 0.8, 0.2], texT], None}, {texTicks[Range[0, 6, 2], texT], None}},
   PlotStyle -> styles1,
   PlotLegends -> legendAt[styles1, {"\\text{Pad\\'e}", "S(r,0)", "S(-r,0)", "S(0,\\pm r)"},
     {{0.98, 0.03}, {1, 0}}],
   Epilog -> labels["a", "r~(\\upmu\\rm{m})", "S"]];
export["strength", fig1];

(* (b) the director angle; the line at zero under the curves *)
styles2 = curve /@ colors[[;; Length[rVals]]];
fig2 = Plot[Evaluate[Table[thetaNum[r, \[Phi]] - \[Phi]/2, {r, rVals}]], {\[Phi], 0.001, 2 Pi - 0.001},
   Evaluate[panelOpts["b"]], PlotRange -> {{0, 2 Pi}, {-1.05, 1.05} Degree},
   PlotRangePadding -> None,
   GridLines -> {None, {0}}, GridLinesStyle -> gridStyle[fP],
   FrameTicks -> {{texTicks[{-1, -0.5, 0, 0.5, 1} Degree, texT,
       (tickString[#/Degree] <> "^\\circ") &], None},
     {{{0, texT["0"]}, {Pi, texT["\\pi"]}, {2 Pi, texT["2\\pi"]}}, None}},
   PlotStyle -> styles2,
   PlotLegends -> legendAt[styles2, StringTemplate["r=`1`~\\upmu\\rm{m}"] /@ rVals,
     {{0.98, 0.03}, {1, 0}}],
   Epilog -> labels["b", "\\phi", "\\theta-\\phi/2"]];
export["defectangle", fig2];

(* (c) the force; the dashed Pade curve at Plot's default width (1.8 pt printed); the core
   x = 0 as a grid line under the curves *)
styles3 = {Directive[Dashed, Black, AbsoluteThickness[2.]], curve[colors[[1]]]};
fig3 = Plot[{forcePade[x] /. modelParams, forceNum[x]}, {x, -3, 3},
   Evaluate[panelOpts["c"]], PlotRange -> {{-3, 3}, {0, 1.2}}, PlotRangePadding -> None,
   GridLines -> {{0}, None}, GridLinesStyle -> gridStyle[fP],
   FrameTicks -> {{texTicks[{0, 0.5, 1}, texT], None}, {texTicks[{-2, 0, 2}, texT], None}},
   PlotStyle -> styles3,
   PlotLegends -> legendAt[styles3, {"\\text{Pad\\'e}", "\\text{numerical}"},
     {{0.98, 0.03}, {1, 0}}],
   Epilog -> labels["c", "x~(\\upmu\\rm{m})", "(\\nabla\\cdot\\mathbf{Q})_x~(\\upmu\\rm{m}^{-1})"]];
export["activeforce", fig3];

Print["done (", stamp[], ")"];
