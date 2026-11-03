#!/usr/bin/env wolframscript
(* variance_map.wls -- the equal-point density variance around a +1/2 defect over its
   pedestal, Var[rho]/Var[rho]_pedestal, as one figure* (fig:variance3d). Written 2026-10-12
   for the prompt that merges Fig. 8 (varrho3D_high.png) with the circle-average maps of
   plotgen/circle_average/circle_average.wls.

     (a) the surface over a disc of radius 6 um about the core, in the style of
         varrho3D_high.png;
     (b) the same quantity seen from above, with a circle of radius r about the maximum of
         the variance, r_c = r_max (blue, as in fig:defects-layering(c));
     (c) the same circle about the core, r_c = 0 (orange, as in fig:defects-layering(c)).
   P(r) of eq:p-layer-radial is the layering probability averaged over that circle.

   All three panels share ONE absolute colour scale (TemperatureMap over zRange, from the
   minimum to the maximum of the variance, his choice of 2026-10-12c; the pedestal is
   therefore not at the white centre of the map), shown by the colour bar under (b, c). The
   director lines are the same analytic +1/2 lines in all panels.

   The variance comes from the production Lyapunov run (box 20, N_int = 180, B = 9e4,
   Dirichlet wall, delta Q != 0), through its companion file sigma_rho_<tag>.m (the
   180 x 180 diagonal of the covariance, first index y, second x), as in circle_average.wls.
   With the argument "analytic" all panels show instead the analytic variance at
   delta Q = 0, eq:total-gauss, computed as in negative-defect/negative_defect.wls.

   The 3D panel is rasterized (600 dpi, 150 in draft) and placed as an image; the density
   maps are Rasters (as in plotgen/DensityPlotRasterize/densityplot.wls); frames, lines,
   circles and the text of (b, c) stay vector. All text is MaTeX (Computer Modern): the
   headless front end garbles plain text in pdf exports.
   Sizes are design points. The composite is 510 wide; the pdf export scales by 0.72, so
   \includegraphics[width=\textwidth] prints every size 1:1.

   Usage:  plotgen/variance_map/run.sh [draft] [analytic] [key=value ...]
     r=<radius of the circle, um>, angle=<direction of the drawn radius, degrees>,
     out=<output directory>
   Outputs (default: next to this script): variance_map.pdf and variance_map.png (600 dpi);
   with "draft", variance_map_draft.png (150 dpi) only. *)

SetEnvironment["QT_QPA_PLATFORM" -> "offscreen"];
SetEnvironment["GALLIUM_DRIVER" -> "softpipe"];   (* llvmpipe crashes the front end on 3D rendering *)

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
t0 = AbsoluteTime[];

(* ---------------------------------------------------------------- variance *)
tag = "box20_N200_Dir_lN0.5_z12.5_ld1_B90000_ss20";
sig = Import[FileNameJoin[{repoRoot, "data", "fluctuations-set2_2_comove_optimized_2_ssDir",
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
{varMin, rMin} = Module[{s = FindMinimum[var[x, y], {{x, 1.7, 0., 4.}, {y, 0., -2., 2.}}]},
   {s[[1]], {x, y} /. s[[2]]}];
Print[If[analytic, "analytic", "numerical (" <> tag <> ")"], " variance"];
Print["  peak ", 100 (varMax - 1), "% at r_max = ", rMax,
  " um;  mirror check var(-1, 0.7) - var(-1, -0.7) = ", var[-1., 0.7] - var[-1., -0.7]];
Print["  trough ", 100 (varMin - 1), "% at ", rMin, " um"];
Print["  circle radius r = ", rCircle, " um"];

(* ---------------------------------------------------------------- graphics *)
<< MaTeX`
(* Computer Modern, as in varrho3D_high.png and the revtex body (plotgen/pick_font) *)
tex[s_, fs_ : 8] := MaTeX[s, FontSize -> fs, "Preamble" -> {"\\usepackage{amsmath,upgreek}"}];
varLabel = "\\mathrm{Var}[\\rho]\\,/\\,\\mathrm{Var}[\\rho]_\\mathrm{pedestal}";

R3 = 6;                     (* radius of the disc of panel (a), um *)
L = 4;                      (* half-width of the squares of (b, c), um *)
zRange = {varMin, varMax};  (* the colour range is the range of the variance (his choice) *)
cmap = ColorData["TemperatureMap"];
colour[z_] := cmap[Clip[Rescale[z, zRange], {0., 1.}]];   (* the one colour scale of all panels *)
nPix = If[draft, 200, 800];
res3D = If[draft, 150, 600];                               (* dpi of the rasterized 3D panel *)
{plotPts, maxRec} = If[draft, {40, 0}, {90, 2}];

(* layout, in design points; the composite is totalW wide *)
totalW = 510;
s = 105;                    (* side of each map frame *)
titleY = 8;                 (* height of the panel titles' centres above the frames *)
padTop = 16;                (* room for the titles *)
padA = {{25, 3}, {25, padTop}}; padB = {{8, 3}, {25, padTop}};   (* {{left, right}, {bottom, top}} *)
barH = 7; padBar = {{padA[[1, 1]], padB[[1, 2]]}, {26, 6}};
rightW = 2 s + Total[padA[[1]]] + Total[padB[[1]]];
rightH = s + Total[padA[[2]]] + barH + Total[padBar[[2]]];
leftW = totalW - rightW;
Print["  layout: (a) ", leftW, " x ", rightH, ",  (b, c) + bar ", rightW, " x ", rightH];

(* range of the variance on the plotted regions; anything outside zRange would be clipped *)
diskVals = Select[Flatten[Table[{x, y}, {x, -R3, R3, 0.05}, {y, -R3, R3, 0.05}], 1],
   Norm[#] <= R3 &];
diskRange = MinMax[var @@@ diskVals];
Print["  disc of radius ", R3, " um: min ", diskRange[[1]], ", max ", diskRange[[2]],
  If[diskRange[[1]] < zRange[[1]] || diskRange[[2]] > zRange[[2]],
   "  -- OUTSIDE the colour range " <> ToString[zRange], "  (inside the colour range " <>
    ToString[zRange] <> ")"]];

(* +1/2 director lines: the level sets of Sqrt[r] Sin[phi/2] = c, i.e. the parabolas
   x = (y^2 - 4c^4)/(4c^2), r = c^2 + y^2/(4c^2), that bend around the core on the -x side,
   and the ray y = 0, x > 0. The same family in all panels. *)
cs = Range[0.35, 2.1, 0.35];

(* -- (a): the surface, its colours on the absolute scale, the lines lifted onto it -- *)
xyTicks = {#, tex[ToString[#]]} & /@ Range[-R3, R3, 2];
zTicks = {#, tex[ToString[#]]} & /@ {0.98, 1, 1.02};
dz = 0.0004;                (* lift of the lines above the surface, 0.7% of zRange *)
lift[pts_] := ({#1, #2, var[#1, #2] + dz} &) @@@ pts;
diskParabola[c_] := With[{ym = 2 c Sqrt[R3 - c^2]},   (* the part with r <= R3 *)
   Table[{(y^2 - 4 c^4)/(4 c^2), y}, {y, Subdivide[-ym, ym, 400]}]];
lines3D = {Directive[Black, Opacity[0.6], AbsoluteThickness[0.6]],
   Line[lift[Table[{x, 0.}, {x, Subdivide[0., R3, 200]}]]],
   Line[lift[diskParabola[#]]] & /@ Select[cs, #^2 < R3 &]};
surface = Plot3D[var[x, y], {x, y} \[Element] Disk[{0, 0}, R3],
   PlotRange -> {{-R3, R3}, {-R3, R3}, zRange + {-1, 1} 0.001},   (* margin: the extremes are not clipped *)
   PlotPoints -> plotPts, MaxRecursion -> maxRec,
   ColorFunction -> (colour[#3] &), ColorFunctionScaling -> False,
   Lighting -> "Neutral", MeshStyle -> Directive[Dashed, Opacity[0.15]],
   BoxRatios -> {1, 1, 0.45}, ViewPoint -> {1.2, -2.0, 2.4},
   SphericalRegion -> False];   (* else the panel is padded to a bounding sphere *)
g3D = Show[surface, Graphics3D[lines3D],
   Ticks -> {xyTicks, xyTicks, zTicks},
   AxesLabel -> {tex["x~(\\upmu\\mathrm{m})", 9], tex["y~(\\upmu\\mathrm{m})", 9],
     Rotate[tex[varLabel, 9], 95 Degree]},
   LabelStyle -> Directive[Black], BoxStyle -> Directive[GrayLevel[0.45], AbsoluteThickness[0.5]],
   AxesStyle -> Directive[Black, AbsoluteThickness[0.5]],
   ImageSize -> {leftW, rightH - padTop}];
img3D = Rasterize[g3D, "Image", ImageResolution -> res3D, Background -> White];
Print["  3D panel rasterized: ", ImageDimensions[img3D], " px  (", Round[AbsoluteTime[] - t0], " s)"];
panelA = Graphics[{
    Inset[img3D, {leftW/2, (rightH - padTop)/2}, Center, leftW],
    Inset[tex["\\textbf{(a)}", 9], {leftW/2, rightH - padTop + titleY}]},
   PlotRange -> {{0, leftW}, {0, rightH}}, PlotRangePadding -> None, AspectRatio -> Full,
   ImagePadding -> None, ImageSize -> {leftW, rightH}];

(* -- (b, c): density maps; pixel centres, rows bottom-up as Raster wants them -- *)
lut = List @@@ (ColorConvert[cmap[#], "RGB"] & /@ Subdivide[0., 1., 511]);
pix = -L + (Range[nPix] - 0.5) 2. L/nPix;
vals = Outer[var[#2, #1] &, pix, pix];                  (* vals[[iy, ix]] *)
Print["  map range ", MinMax[vals], " (colour range ", zRange, ")"];
idx = Round[511 Clip[Rescale[vals, zRange], {0., 1.}]] + 1;
densityRaster = Raster[Map[lut[[#]] &, idx, {2}], {{-L, -L}, {L, L}}];

(* director lines cut to the runs inside the square by hand: PlotRangeClipping would also
   clip the axis labels, which are insets in the margins *)
inSquare[{x_, y_}] := Abs[x] <= L && Abs[y] <= L;
parabola[c_] := Line[Select[SplitBy[Table[{(y^2 - 4 c^4)/(4 c^2), y},
      {y, Subdivide[-L, L, 800]}], inSquare], inSquare@First@# && Length[#] > 1 &]];
directorLines = {Directive[Black, Opacity[0.6], AbsoluteThickness[0.6]],
   Line[{{0, 0}, {L, 0}}], parabola /@ cs};

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
(* labels and titles as insets a fixed distance (pt) from a frame of side h: FrameLabel
   and PlotLabel leave gaps sized for the default fonts *)
xLabel[lab_, pt_, h_] := Inset[tex[lab, 9], Scaled[{0.5, -pt/h}]];
yLabel[lab_, pt_, w_] := Inset[tex[lab, 9], Scaled[{-pt/w, 0.5}], Center, Automatic, {0, 1}];
title[lab_, pt_, h_] := Inset[tex[lab, 9], Scaled[{0.5, 1 + pt/h}]];
panel[c_, col_, label_, first_] := Graphics[{densityRaster, directorLines, annot[c, col],
    title[label, titleY, s], xLabel["x~(\\upmu\\mathrm{m})", 19, s],
    If[first, yLabel["y~(\\upmu\\mathrm{m})", 19, s], {}]},
   Frame -> True, FrameStyle -> Directive[Black, AbsoluteThickness[0.5]],
   FrameTicks -> {{If[first, ticks, noLabels], None}, {ticks, None}},   (* left and bottom only *)
   PlotRange -> {{-L, L}, {-L, L}}, PlotRangePadding -> None, AspectRatio -> 1,
   ImagePadding -> If[first, padA, padB],
   ImageSize -> {s + Total[If[first, padA, padB][[1]]], s + Total[If[first, padA, padB][[2]]]}];

(* shared horizontal colour bar under (b, c), its frame aligned with theirs *)
barTicks = {#, tex[ToString[#]]} & /@
   Select[{0.97, 0.98, 0.99, 1, 1.01, 1.02, 1.03}, zRange[[1]] <= # <= zRange[[2]] &];
barW = 2 s + padA[[1, 2]] + padB[[1, 1]];
colorBar = Graphics[{Raster[{lut}, {{zRange[[1]], 0}, {zRange[[2]], 1}}],
    xLabel[varLabel, 19, barH]},
   Frame -> True, FrameStyle -> Directive[Black, AbsoluteThickness[0.5]],
   FrameTicks -> {{None, None}, {barTicks, None}},
   PlotRange -> {zRange, {0, 1}}, PlotRangePadding -> None, AspectRatio -> Full,
   ImagePadding -> padBar, ImageSize -> {barW + Total[padBar[[1]]], barH + Total[padBar[[2]]]}];

maps = Grid[{
    {panel[rMax, ColorData[97][1], "\\textbf{(b)}\\enspace\\mathbf r_\\mathrm c=\\mathbf r_\\mathrm{max}", True],
     panel[{0., 0.}, ColorData[97][2], "\\textbf{(c)}\\enspace\\mathbf r_\\mathrm c=0", False]},
    {colorBar, SpanFromLeft}}, Spacings -> {0, 0}, Alignment -> {Left, Top}];
fig = Grid[{{panelA, maps}}, Spacings -> {0, 0}, Alignment -> {Left, Top}];

If[draft,
  Export[FileNameJoin[{outDir, "variance_map_draft.png"}], fig, ImageResolution -> 150],
  Export[FileNameJoin[{outDir, "variance_map.pdf"}], fig, "AllowRasterization" -> False];
  Export[FileNameJoin[{outDir, "variance_map.png"}], fig, ImageResolution -> 600]];
Print["exported to ", outDir, "  (", Round[AbsoluteTime[] - t0], " s)"];
