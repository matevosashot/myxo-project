#!/usr/bin/env wolframscript
(* drift_deformation.wls -- fig:drift-deformation of app:drift (2026-10-12b): how the closure
   term T of eq:drift-U-term rotates the director near the core of the comoving +1/2 defect.

   Director angle relative to the ideal profile, theta - phi/2, against the polar angle phi
   on circles of radius r = 0.1, 1 and 2 um about the core, as in defectangle.pdf
   (fig:comoving-deformation(b)), with theta = Arg[Q1 + i Q2]/2 taken mod pi.
     dashed: w = 0, the retained transport (the cross-coupling) alone, the paper's steady state;
     solid:  w = 1, with T at the exact moments.
   Colours ColorData[97] 1-3, one per radius, as in defectangle.pdf.

   Data: claude_experiments/drift_closure_2026_10_09/out/r1_exact_v2.mx (run_lambda.wls, box
   half-width 20 um, refine 1, Dirichlet reservoir wall, free Q; path w = 0, 0.25, 0.5, 0.75,
   1). Read as in analyze_deformation.wls: X = {rho, Q1, Q2} on the graded grid {gx, gy}.
   Expected extrema (printed): w = 0: 0.031, 0.341, 0.875 deg; w = 1: 0.826, 4.905, 8.730 deg.

   All text is MaTeX (Computer Modern): the headless front end garbles plain text in pdf.
   2026-11-01 (the author's request: smaller and inconsistent fonts, labels too far from the axes):
   ticks and axis labels at 9 pt, legends at 8 pt (printed), as in fig:drift-closure above it
   (matplotlib, font.size 9, legend 8); axis labels placed by frame_labels.wl, about 4 pt
   from the tick labels.
   2026-11-02 (the author's rule: no line thinner than about 1 px; a line meant to look
   thinner is gray): frame and ticks print at
   0.75 pt (0.51 pt before), the line at 0 gray 0.6 at 0.75 pt (gray 0.35 at 0.51 pt, over the
   black axis of ListLinePlot, now switched off), as a grid line under the curves (a Line
   over them before).
   2026-11-02, the author's later request: the line at 0 prints at 0.6 pt in a darker gray
   (gridStyle of frame_labels.wl). The x label names the angle: "Polar angle phi" (was "phi").
   Usage:  plotgen/drift_deformation/run.sh [out=<dir>]
   Outputs: drift_deformation.pdf and drift_deformation.png (600 dpi) next to this script;
   run.sh copies the pdf to manuscript/figures/. *)

SetEnvironment["QT_QPA_PLATFORM" -> "offscreen"];
$HistoryLength = 0;

here = DirectoryName[$InputFileName];
repoRoot = ParentDirectory[ParentDirectory[here]];
opts = Association[Rule @@@ (StringSplit[#, "=", 2] & /@
      Select[Rest@$ScriptCommandLine, StringContainsQ[#, "="] &])];
outDir = Lookup[opts, "out", here];

(* ---------------------------------------------------------------- data *)
Get[FileNameJoin[{repoRoot, "claude_experiments", "drift_closure_2026_10_09", "out",
    "r1_exact_v2.mx"}]];                                  (* defines out *)
gx = out["gx"]; gy = out["gy"]; n = out["n"]; nx = out["nx"];
interp[X_, k_] := ListInterpolation[Transpose@Partition[X[[(k - 1) n + 1 ;; k n]], nx],
   {gx, gy}, InterpolationOrder -> 3];
lamOf[p_] := If[ListQ[p["lam"]], p["lam"][[1]], p["lam"]];
state[w_] := First[Select[out["path"], Abs[lamOf[#] - w] < 10^-9 &]]["X"];

angWrap[a_] := Mod[a + Pi/2, Pi] - Pi/2;                  (* director angles are defined mod pi *)
phis = N@Subdivide[0., 2. Pi, 720];
rs = {0.1, 1., 2.};
dTheta[w_, r_] := With[{X = state[w]},
   angWrap[Arg[interp[X, 2][r Cos[phis], r Sin[phis]] +
         I interp[X, 3][r Cos[phis], r Sin[phis]]]/2 - phis/2]/Degree];
curves = Table[dTheta[w, r], {w, {0., 1.}}, {r, rs}];      (* curves[[w index, r index]] *)
Do[Print["w = ", {0, 1}[[i]], ", r = ", rs[[k]], " um:  max |theta - phi/2| = ",
   Max@Abs@curves[[i, k]], " deg at phi = ",
   phis[[First@Ordering[Abs@curves[[i, k]], -1]]]/Degree, " deg"], {i, 2}, {k, 3}];

(* ---------------------------------------------------------------- figure *)
<< MaTeX`
Get[FileNameJoin[{repoRoot, "plotgen", "frame_labels.wl"}]];
(* ticks and axis labels print at 9 pt, legends at 8 pt: the figure is 242 pt wide and printed
   at \linewidth = 246 pt *)
tex[s_, fs_ : 9*242/246.] := MaTeX[s, FontSize -> fs, "Preamble" -> {"\\usepackage{amsmath,upgreek}"}];
legFS = 8*242/246.;
cols = ColorData[97] /@ {1, 2, 3};
solid[c_] := Directive[c, AbsoluteThickness[1.3]];
dashed[c_] := Directive[c, AbsoluteThickness[1.1], AbsoluteDashing[{3.5, 2.}]];
yMax = 10.5;

pad = {{32, 8}, {23.5, 5}};            (* {{left, right}, {bottom, top}}, pt *)
{wFr, hFr} = {202, 131.5};             (* frame size, pt: one column with the padding *)
sample[sty_] := Graphics[{sty, Line[{{0, 0}, {1, 0}}]}, ImageSize -> {17, 6},
   PlotRange -> {{0, 1}, {-1, 1}}, ImagePadding -> 0, AspectRatio -> Full];
legend[rows_] := Grid[{sample[#1], tex[#2, legFS]} & @@@ rows, Alignment -> {Left, Center},
   Spacings -> {0.4, 0.15}];
(* the w = 1 curves are positive on (0, pi) and negative on (pi, 2 pi), the w = 0 curves stay
   within 0.9 deg: the lower left and the upper right of the frame are free *)
radiusLegend = legend[Table[{solid[cols[[k]]],
     "r=" <> {"0.1", "1", "2"}[[k]] <> "~\\upmu\\mathrm m"}, {k, 3}]];
styleLegend = legend[{{solid[Black], "\\text{with } T~(w=1)"},
    {dashed[Black], "\\text{without } T~(w=0)"}}];

(* MaTeX tick labels trimmed by 1 pt at the sides, as in the other figures (frame_labels.wl) *)
degTicks = {#, trimSides[tex[ToString[#] <> "^\\circ"], 1.]} & /@ Range[-8, 8, 4];
phiTicks = Join[{{0, tex["0"]}, {Pi, tex["\\pi"]}, {2 Pi, tex["2\\pi"]}},
   {#, ""} & /@ {Pi/2, 3 Pi/2}];

fig = Show[
   ListLinePlot[Join[Transpose[{phis, #}] & /@ curves[[2]], Transpose[{phis, #}] & /@ curves[[1]]],
    PlotStyle -> Join[solid /@ cols, dashed /@ cols],
    PlotRange -> {{0, 2 Pi}, {-yMax, yMax}}],   (* else Automatic cuts the w = 1 curves *)
   Graphics[{
     Inset[radiusLegend, Scaled[{0.03, 0.04}], {Left, Bottom}],
     Inset[styleLegend, Scaled[{0.985, 0.97}], {Right, Top}],
     frameLabels[{tex["\\text{Polar angle}~\\phi"], tex["\\theta-\\phi/2"]}, {13.4, 20.2}]}],
   Frame -> True, FrameStyle -> frameStyle[242/246.],
   Axes -> False,   (* ListLinePlot's axis at 0, black 0.5 pt, lay under the gray line *)
   GridLines -> {None, {0}}, GridLinesStyle -> gridStyle[242/246.],   (* under the curves *)
   FrameTicks -> {{degTicks, None}, {phiTicks, None}},   (* ticks on the left and bottom only *)
   PlotRange -> {{0, 2 Pi}, {-yMax, yMax}}, PlotRangePadding -> None,
   PlotRangeClipping -> False,   (* ListLinePlot clips by default, which hides the axis labels,
     insets in the margins; the curves stay within 8.8 deg *)
   AspectRatio -> Full, ImagePadding -> pad,
   ImageSize -> {wFr + Total[pad[[1]]], hFr + Total[pad[[2]]]}];

Export[FileNameJoin[{outDir, "drift_deformation.pdf"}], fig, "AllowRasterization" -> False];
Export[FileNameJoin[{outDir, "drift_deformation.png"}], fig, ImageResolution -> 600];
Print["exported to ", outDir];
