#!/usr/bin/env wolframscript
(* circle_average.wls -- the circle average of eq:p-layer-radial, drawn on the density
   variance around a +1/2 defect. Written 2026-10-12 for the prompt above
   eq:p-layer-radial (sec:nucleation-2).

   Two density plots of the quantity of fig:variance3d (varrho3D_high.png), the equal-point
   density variance over its pedestal, Var[rho]/Var[rho]_pedestal, from the production
   Lyapunov run (box 20, N_int = 180, B = 9e4, Dirichlet wall, delta Q != 0), in the same
   colour map (TemperatureMap), with the +1/2 director lines on top. Each panel carries a
   circle of radius r, centred
     (a) on the maximum of the variance, r_c = r_max (blue, as in fig:defects-layering(c)),
     (b) on the core, r_c = 0 (orange, as in fig:defects-layering(c)).
   P(r) of eq:p-layer-radial is the layering probability averaged over that circle.

   The variance comes from the run's companion file sigma_rho_<tag>.m (the 180 x 180
   diagonal of the covariance, first index y, second x), not from the 76 GB .bin. With the
   argument "analytic" the map is instead the analytic variance at delta Q = 0,
   eq:total-gauss, computed as in negative-defect/negative_defect.wls.

   The density map is one Raster (as in plotgen/DensityPlotRasterize/densityplot.wls);
   frame, director lines, circles and all text stay vector. All text is MaTeX: the
   headless front end garbles plain text in pdf exports.

   Usage:  plotgen/circle_average/run.sh [draft] [analytic] [key=value ...]
     r=<radius of the circle, um>, angle=<direction of the drawn radius, degrees>,
     out=<output directory>
   Outputs (default: next to this script): circle_average.pdf and circle_average.png
   (600 dpi); with "draft", circle_average_draft.png (150 dpi) only. *)

SetEnvironment["QT_QPA_PLATFORM" -> "offscreen"];

here = DirectoryName[$InputFileName];
repoRoot = ParentDirectory[ParentDirectory[here]];
draft = MemberQ[$ScriptCommandLine, "draft"];
analytic = MemberQ[$ScriptCommandLine, "analytic"];
opts = Association[Rule @@@ (StringSplit[#, "=", 2] & /@
      Select[Rest@$ScriptCommandLine, StringContainsQ[#, "="] &])];
num[key_, default_] := If[KeyExistsQ[opts, key], N@ToExpression@opts[key], default];
rCircle = num["r", 1.5];
angR = num["angle", -90] Degree;
outDir = Lookup[opts, "out", here];

(* ---------------------------------------------------------------- variance *)
tag = "box20_N180_Dir_lN0.5_z12.5_ld1_B90000_ss20";
sig = Import[FileNameJoin[{repoRoot, "data", "fluctuations-set2_2_comove_optimized_2",
     "sigma_rho_" <> tag <> ".m"}]];
modelParams = sig["modelParams"];
{lN, Bv, zetav, rho0v} = {lNoise, B, \[Zeta], \[Rho]0} /. modelParams;
pedestal = zetav/(2 Pi Bv rho0v lN^2);                 (* eq:pedestal *)

If[! analytic,
  grid = N@sig["gridInt"];
  varNum = ListInterpolation[Transpose[sig["sigmaRho"]]/pedestal, {grid, grid},
    InterpolationOrder -> 2];                           (* Transpose: x first *)
  var[x_?NumericQ, y_?NumericQ] := varNum[x, y],

  (* analytic variance at delta Q = 0, eq:total-gauss, as in negative_defect.wls:
     Var/pedestal - 1 = (lN^2/4) H_1(r) Cos[phi] for the accurate Pade profile *)
  fPade[z_] := Sqrt[(0.34010999 + 0.04516965 z + 0.00233656 z^2)/(1 +
       0.38387160 z + 0.04750621 z^2 + 0.00233656 z^3)];
  S0PadeOverR[r2_] := Sqrt[-a1/(2 b)] Sqrt[-a1/L1] fPade[-r2 a1/L1];
  Sfn[r_] = r S0PadeOverR[r^2] /. modelParams;
  g1[r_] = D[D[Sfn[r], r] + Sfn[r]/r, r] // Simplify;
  phiVec[t_] := Module[{s = UnitStep[1.*^-6 - t], td}, td = t + s;
     (1 - s) ((1 - Exp[-td])/td) + s (1 - t/2 + t^2/6)];
  filterH1[rs_List, nr_ : 4000, nA_ : 720, rmax_ : 400.] :=
   Module[{rp, rp2, gp, ca, wgt, dw, il2},
    rp = Exp[Subdivide[N@Log[0.0001], N@Log[rmax], nr - 1]]; rp2 = rp^2;
    gp = g1 /@ rp;
    ca = Cos[(Range[0, nA - 1] + 0.5) 2. Pi/nA];
    wgt = ca 2. Pi/nA;
    dw = Differences[rp]; il2 = 1./(Pi lN^2);
    Table[Module[{w2, k1, integ},
      w2 = (r^2 + rp2) - 2. r Outer[Times, rp, ca];
      k1 = il2 phiVec[2. w2/lN^2] . wgt;
      integ = rp gp k1;
      Total[(Most[integ] + Rest[integ])/2 dw]], {r, rs}]];
  hRadii = Join[{0.}, Exp[Subdivide[N@Log[0.01], N@Log[48.], 119]]];
  H1 = Interpolation[Transpose[{hRadii, Join[{0.}, filterH1[Rest[hRadii]]]}],
    InterpolationOrder -> 3];
  var[x_?NumericQ, y_?NumericQ] := With[{r = Sqrt[x^2 + y^2]},
    If[r == 0., 1., 1 + lN^2/4 H1[r] x/r]]];

{varMax, rMax} = Module[{s = FindMaximum[var[x, y], {{x, -1.3, -4., 0.}, {y, 0., -2., 2.}}]},
   {s[[1]], {x, y} /. s[[2]]}];
rMax = {rMax[[1]], 0.};   (* the map is mirror symmetric in y; the solver leaves ~1e-7 *)
Print[If[analytic, "analytic", "numerical (" <> tag <> ")"], " variance"];
Print["  peak ", 100 (varMax - 1), "% at r_max = ", rMax,
  " um;  mirror check var(-1, 0.7) - var(-1, -0.7) = ", var[-1., 0.7] - var[-1., -0.7]];
Print["  circle radius r = ", rCircle, " um"];

(* ---------------------------------------------------------------- graphics *)
<< MaTeX`
(* Computer Modern, as in varrho3D_high.png and the revtex body (plotgen/pick_font) *)
tex[s_, fs_ : 8] := MaTeX[s, FontSize -> fs, "Preamble" -> {"\\usepackage{amsmath,upgreek}"}];

L = 4;                      (* half-width of the plotted square, um *)
zRange = {0.967, 1.033};    (* symmetric about 1, so the pedestal is white; covers -3.25%..+2.4% *)
cmap = ColorData["TemperatureMap"];
nPix = If[draft, 200, 800];
s = 100;                    (* side of each plot square, pt; two panels fill one column *)

(* density map: pixel centres, rows bottom-up as Raster wants them *)
lut = List @@@ (ColorConvert[cmap[#], "RGB"] & /@ Subdivide[0., 1., 511]);
pix = -L + (Range[nPix] - 0.5) 2. L/nPix;
vals = Outer[var[#2, #1] &, pix, pix];                  (* vals[[iy, ix]] *)
Print["  map range ", MinMax[vals], " (colour range ", zRange, ")"];
idx = Round[511 Clip[Rescale[vals, zRange], {0., 1.}]] + 1;
densityRaster = Raster[Map[lut[[#]] &, idx, {2}], {{-L, -L}, {L, L}}];

(* +1/2 director lines: the level sets of Sqrt[r] Sin[phi/2], i.e. the parabolas
   x = (y^2 - 4c^4)/(4c^2) that bend around the core on the -x side, and the ray y = 0, x > 0.
   Cut to the runs inside the square by hand: PlotRangeClipping would also clip the axis
   labels, which are insets in the margins *)
inSquare[{x_, y_}] := Abs[x] <= L && Abs[y] <= L;
parabola[c_] := Line[Select[SplitBy[Table[{(y^2 - 4 c^4)/(4 c^2), y},
      {y, Subdivide[-L, L, 800]}], inSquare], inSquare@First@# && Length[#] > 1 &]];
directorLines = {Directive[Black, Opacity[0.6], AbsoluteThickness[0.6]],
   Line[{{0, 0}, {L, 0}}], parabola /@ Range[0.35, 2.1, 0.35]};

(* the circle of eq:p-layer-radial about c, its centre, and one radius labelled r; a white
   casing under the curves keeps them visible on the colour map (all casings first, so
   the radius does not notch the circle) *)
annot[c_, col_] := Module[{u = {Cos[angR], Sin[angR]}},
   {{White, AbsoluteThickness[2.8], Circle[c, rCircle], AbsoluteThickness[2.2],
     Line[{c, c + rCircle u}], Disk[c, 0.19]},
    {col, AbsoluteThickness[1.5], Circle[c, rCircle], AbsoluteThickness[1.],
     Line[{c, c + rCircle u}], Disk[c, 0.13]},
    Inset[tex["r", 9], c + 0.5 rCircle u + 0.3 {-u[[2]], u[[1]]}]}];

ticks = {#, tex[ToString[#]]} & /@ Range[-L, L, 2];
noLabels = {#, ""} & /@ Range[-L, L, 2];
(* axis labels as insets a fixed distance (pt) outside a frame of height h: FrameLabel
   leaves a gap sized for the default 12 pt tick font *)
xLabel[lab_, pt_, h_] := Inset[tex[lab, 9], Scaled[{0.5, -pt/h}]];
yLabel[lab_, pt_, w_] := Inset[tex[lab, 9], Scaled[{-pt/w, 0.5}], Center, Automatic, {0, 1}];
padA = {{25, 3}, {25, 4}}; padB = {{8, 3}, {25, 4}};   (* {{left, right}, {bottom, top}} *)
panel[c_, col_, label_, first_] := Graphics[{densityRaster, directorLines, annot[c, col],
    xLabel["x~(\\upmu\\mathrm{m})", 19, s],
    If[first, yLabel["y~(\\upmu\\mathrm{m})", 19, s], {}]},
   Frame -> True, FrameStyle -> Directive[Black, AbsoluteThickness[0.5]],
   FrameTicks -> {{If[first, ticks, noLabels], None}, {ticks, None}},   (* left and bottom only *)
   PlotLabel -> tex[label, 9],
   PlotRange -> {{-L, L}, {-L, L}}, PlotRangePadding -> None, AspectRatio -> 1,
   ImagePadding -> If[first, padA, padB],
   ImageSize -> {s + Total[If[first, padA, padB][[1]]], s + Total[If[first, padA, padB][[2]]]}];

(* shared horizontal colour bar under both panels, its frame aligned with theirs *)
barTicks = {#, tex[ToString[#]]} & /@ {0.97, 0.98, 0.99, 1, 1.01, 1.02, 1.03};
barH = 7;
padBar = {{padA[[1, 1]], padB[[1, 2]]}, {26, 2}};
barW = 2 s + padA[[1, 2]] + padB[[1, 1]];
colorBar = Graphics[{Raster[{lut}, {{zRange[[1]], 0}, {zRange[[2]], 1}}],
    xLabel["\\mathrm{Var}[\\rho]\\,/\\,\\mathrm{Var}[\\rho]_\\mathrm{pedestal}", 19, barH]},
   Frame -> True, FrameStyle -> Directive[Black, AbsoluteThickness[0.5]],
   FrameTicks -> {{None, None}, {barTicks, None}},
   PlotRange -> {zRange, {0, 1}}, PlotRangePadding -> None, AspectRatio -> Full,
   ImagePadding -> padBar, ImageSize -> {barW + Total[padBar[[1]]], barH + Total[padBar[[2]]]}];

fig = Grid[{
    {panel[rMax, ColorData[97][1], "\\textbf{(a)}\\enspace\\mathbf r_\\mathrm c=\\mathbf r_\\mathrm{max}", True],
     panel[{0., 0.}, ColorData[97][2], "\\textbf{(b)}\\enspace\\mathbf r_\\mathrm c=0", False]},
    {colorBar, SpanFromLeft}}, Spacings -> {0, 0.2}, Alignment -> {Left, Top}];

If[draft,
  Export[FileNameJoin[{outDir, "circle_average_draft.png"}], fig, ImageResolution -> 150],
  Export[FileNameJoin[{outDir, "circle_average.pdf"}], fig, "AllowRasterization" -> False];
  Export[FileNameJoin[{outDir, "circle_average.png"}], fig, ImageResolution -> 600]];
Print["exported to ", outDir];
