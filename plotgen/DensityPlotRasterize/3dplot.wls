#!/usr/bin/env wolframscript
(* Generates an example Plot3D with contour lines on the surface and exports it as PDF and SVG.

   Only the surface and the color-bar gradient are bitmaps; contour lines, box edges, ticks,
   tick labels, axis labels and the title stay vector.

   Mathematica rasterizes every Graphics3D on PDF/SVG export, so the figure is rebuilt as 2D
   Graphics: the surface alone is rasterized (transparent background) with a fixed camera
   (ViewMatrix), and every vector element is projected to 2D with that same camera. Box edges
   between two back-facing faces are drawn under the raster, all other edges over it.
   Contour lines get hidden-line removal by ray-marching from each point to the camera against
   the height field. *)

outFiles = If[Length[$ScriptCommandLine] > 1, Rest[$ScriptCommandLine], {"3dplot.pdf", "3dplot.svg"}];

f[x_, y_] := Sin[x] Cos[y] + 0.3 Sin[x y / 2];
colors = ColorData["SunsetColors"];
{{xmin, xmax}, {ymin, ymax}} = N[{{-2 Pi, 2 Pi}, {-2 Pi, 2 Pi}}];
{zmin, zmax} = MinMax[f @@ Transpose[Tuples[{Subdivide[xmin, xmax, 400], Subdivide[ymin, ymax, 400]}]]];
pr = {{xmin, xmax}, {ymin, ymax}, {zmin, zmax}};
levels = Subdivide[zmin, zmax, 12][[2 ;; -2]];
rasterPixels = 1200; (* surface render resolution *)
n = rasterPixels;
opts3D = {PlotRange -> pr, BoxRatios -> {1, 1, 0.4}, Boxed -> False, Axes -> False, ImageSize -> {n, n}};

(* ---- Surface (raster) and camera ---- *)
surface = Plot3D[
  f[x, y],
  {x, xmin, xmax}, {y, ymin, ymax},
  PlotPoints -> 80,
  Mesh -> None,
  ColorFunction -> (colors[Rescale[#3, {zmin, zmax}]] &), ColorFunctionScaling -> False,
  Evaluate[Sequence @@ opts3D]
];
vm = ViewMatrix /. AbsoluteOptions[surface, ViewMatrix];
{tm, pm} = vm;
camera = With[{c = Inverse[tm] . {0, 0, 0, 1}}, c[[1 ;; 3]] / c[[4]]];
If[camera[[3]] <= zmax, Print["Camera must be above the surface for hidden-line removal."]; Exit[1]];

img = Rasterize[Show[surface, ViewMatrix -> vm], "Image",
  RasterSize -> {n, n}, ImageSize -> {n, n}, Background -> None];

(* ---- Calibration: camera projection -> pixel coordinates of the rendered image ----
   Mathematica fits the projected plot box into the image; rather than rely on that rule,
   render marker points with identical settings and fit pixel = a * proj + b per axis. *)
proj[pts_?MatrixQ] := With[{h = Transpose[(pm . tm) . Transpose[Join[pts, ConstantArray[1., {Length[pts], 1}], 2]]]},
  h[[All, 1 ;; 2]] / h[[All, 4]]];

markers = With[{mid = Mean /@ pr, q = (#2 - #1) / 4 & @@@ pr},
  Prepend[mid + q # & /@ Tuples[{-1, 1}, 3], mid]];
calib = Rasterize[
  Graphics3D[{Red, AbsolutePointSize[8], Point[markers]}, ViewMatrix -> vm, Background -> White, Sequence @@ opts3D],
  "Image", RasterSize -> {n, n}, ImageSize -> {n, n}];
mask = With[{d = ImageData[calib]}, UnitStep[d[[All, All, 1]] - 0.7] UnitStep[0.4 - d[[All, All, 2]]]];
centroids = Values[ComponentMeasurements[MorphologicalComponents[Image[mask]], "Centroid"]];
If[Length[centroids] =!= Length[markers], Print["Calibration failed: found ", Length[centroids], " markers."]; Exit[1]];

predicted = proj[markers];
standardize[m_] := Transpose[Standardize /@ Transpose[m]];
measured = centroids[[First[Nearest[standardize[centroids] -> "Index", #]] & /@ standardize[predicted]]];
fits = Table[LeastSquares[Transpose[{predicted[[All, k]], ConstantArray[1., Length[markers]]}], measured[[All, k]]], {k, 2}];
toPx[pts_?MatrixQ] := With[{q = proj[pts]}, Transpose[Table[fits[[k, 1]] q[[All, k]] + fits[[k, 2]], {k, 2}]]];
toPx[v_?VectorQ] := First[toPx[{v}]];
Print["Calibration max residual: ", Max[Abs[toPx[markers] - measured]], " px"];

(* ---- Box edges ---- *)
mid = Mean /@ pr;
frontFace[a_, s_] := If[s == 2, camera[[a]] > pr[[a, 2]], camera[[a]] < pr[[a, 1]]];
edgePts[a_, {b_, sb_}, {c_, sc_}] := Table[
  ReplacePart[ConstantArray[0., 3], {a -> pr[[a, s]], b -> pr[[b, sb]], c -> pr[[c, sc]]}], {s, 2}];
edges = Flatten[Table[
  With[{bc = Complement[{1, 2, 3}, {a}]}, {a, {bc[[1]], sb}, {bc[[2]], sc}}],
  {a, 3}, {sb, 2}, {sc, 2}], 2];
backEdge[{a_, {b_, sb_}, {c_, sc_}}] := ! frontFace[b, sb] && ! frontFace[c, sc];
edgeLine[e_] := Line[toPx[edgePts @@ e]];

(* ---- Axes: ticks, numbers and labels, built in 3D and projected ---- *)
unit[k_, s_] := ReplacePart[ConstantArray[0., 3], k -> If[s == 2, 1., -1.]];
range[k_] := pr[[k, 2]] - pr[[k, 1]];
nearSide[k_] := If[camera[[k]] < mid[[k]], 1, 2];
(* x and y axes on the bottom edges nearest the camera, ticks pointing horizontally outward;
   z axis on the leftmost vertical silhouette edge, ticks along its back-facing face normal. *)
zEdge = First[SortBy[
  Select[Tuples[{{1, 2}, {1, 2}}], Xor[frontFace[1, #[[1]]], frontFace[2, #[[2]]]] &],
  First[toPx[Mean[edgePts[3, {1, #[[1]]}, {2, #[[2]]}]]]] &]];
zOut = If[frontFace[1, zEdge[[1]]], {2, zEdge[[2]]}, {1, zEdge[[1]]}];
axisSpecs = {
  {1, ReplacePart[ConstantArray[0., 3], {2 -> pr[[2, nearSide[2]]], 3 -> zmin}], {2, nearSide[2]}, "x"},
  {2, ReplacePart[ConstantArray[0., 3], {1 -> pr[[1, nearSide[1]]], 3 -> zmin}], {1, nearSide[1]}, "y"},
  {3, ReplacePart[ConstantArray[0., 3], {1 -> pr[[1, zEdge[[1]]]], 2 -> pr[[2, zEdge[[2]]]]}], zOut, "z"}
};
axisGraphics[{k_, base_, {o_, os_}, label_}] := Module[
  {dir = unit[o, os] range[o], ticks, at},
  ticks = Select[FindDivisions[pr[[k]], 5], pr[[k, 1]] <= # <= pr[[k, 2]] &];
  at[v_, d_] := toPx[ReplacePart[base, k -> v] + d dir];
  {
    Line[{at[#, 0], at[#, 0.025]}] & /@ ticks,
    Text[If[IntegerQ[#], #, N[#]], at[#, 0.09]] & /@ ticks,
    Text[Style[label, Italic], at[Mean[pr[[k]]], 0.2]]
  }
];

(* ---- Contour lines on the surface with hidden-line removal ---- *)
contourPolylines = Flatten[
  Cases[Normal[ContourPlot[f[x, y], {x, xmin, xmax}, {y, ymin, ymax},
      PlotPoints -> 80, Contours -> levels, ContourShading -> None]],
    Line[l_] :> If[MatrixQ[l, NumericQ], {l}, l], Infinity], 1];
resample[pts_, h_] := Append[
  Join @@ MapThread[Function[{a, b}, With[{k = Max[1, Ceiling[Norm[b - a] / h]]}, Table[a + (b - a) i / k, {i, 0, k - 1}]]],
    {Most[pts], Rest[pts]}],
  Last[pts]];
lift[pts_] := Join[pts, Transpose[{f @@ Transpose[pts]}], 2];

(* A point P is hidden if the segment P -> camera crosses the surface before leaving the plot
   box, i.e. it has samples both above and below the surface. A ray that stays below the
   surface is fine: it sees the underside, which Plot3D draws too. Vectorized over all
   points; samples run from just off P to the box exit. *)
visibleQ[pts3_] := Module[{d, sExit, s0, above, below, diff, tol = 1.*^-3 (zmax - zmin)},
  d = ConstantArray[camera, Length[pts3]] - pts3;
  sExit = Min[1, ##] & @@@ Transpose[Join[
    Table[With[{dk = d[[All, k]]}, (If[camera[[k]] > mid[[k]], pr[[k, 2]], pr[[k, 1]]] - pts3[[All, k]]) / dk], {k, 2}],
    {(zmax - pts3[[All, 3]]) / d[[All, 3]]}]];
  s0 = 0.01 / Sqrt[Total[d^2, {2}]];
  above = below = ConstantArray[0, Length[pts3]];
  Do[
    With[{s = s0 + (sExit - s0) u},
      diff = (pts3[[All, 3]] + s d[[All, 3]]) - f[pts3[[All, 1]] + s d[[All, 1]], pts3[[All, 2]] + s d[[All, 2]]];
      above += UnitStep[diff - tol];
      below += UnitStep[-diff - tol]],
    {u, Subdivide[0., 1., 200]}];
  Thread[above below == 0]
];
visibleRuns[pts2_] := Module[{pts3 = lift[resample[pts2, 0.03]], vis},
  vis = visibleQ[pts3];
  Cases[Split[Transpose[{pts3, vis}], #1[[2]] === #2[[2]] &],
    run : {{_, True}, __} :> run[[All, 1]]]
];
contourLines = Line[toPx[#]] & /@ Flatten[visibleRuns /@ contourPolylines, 1];

(* ---- Compose in 2D: back edges, surface raster, contours, front edges, axes ---- *)
boxBounds = MinMax /@ Transpose[toPx[Tuples[pr]]];
{{px0, px1}, {py0, py1}} = {Max[0, Floor[#1]], Min[n, Ceiling[#2]]} & @@@ boxBounds;
surfaceRaster = Raster[Reverse[ImageData[img]][[py0 + 1 ;; py1, px0 + 1 ;; px1]], {{px0, py0}, {px1, py1}}];

edgeStyle = Directive[GrayLevel[0.4], AbsoluteThickness[0.6]];
plot = Graphics[{
    {edgeStyle, edgeLine /@ Select[edges, backEdge]},
    surfaceRaster,
    {Black, Opacity[0.6], AbsoluteThickness[0.5], contourLines},
    {edgeStyle, edgeLine /@ Select[edges, ! backEdge[#] &]},
    {GrayLevel[0.2], AbsoluteThickness[0.6], axisGraphics /@ axisSpecs}
  },
  PlotLabel -> "Sin[x] Cos[y] + 0.3 Sin[x y/2]",
  PlotRange -> All, PlotRangePadding -> Scaled[0.06],
  ImageSize -> 320
];

(* Color bar: a 1x256 raster gradient with vector value ticks on the right. *)
colorBar = Graphics[
  Raster[Table[{List @@ ColorConvert[colors[t], "RGB"]}, {t, 0, 1, 1/255}], {{0, zmin}, {1, zmax}}],
  Frame -> {{False, True}, {False, False}},
  FrameTicks -> {{None, FindDivisions[{zmin, zmax}, 5]}, {None, None}},
  FrameStyle -> GrayLevel[0.4], AspectRatio -> 12, PlotRangePadding -> None,
  ImageSize -> {Automatic, 150}
];

Do[
  Export[file, Legended[plot, Placed[colorBar, After]], "AllowRasterization" -> False];
  Print["Wrote ", AbsoluteFileName[file]],
  {file, outFiles}
];
