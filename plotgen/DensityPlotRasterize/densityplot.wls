#!/usr/bin/env wolframscript
(* Generates an example DensityPlot with contour lines on top and exports it as PDF and SVG.

   Only the density map and the color-bar gradient are bitmaps; contours, frame, ticks,
   tick labels, axis labels and the title stay vector. The density map is rasterized on
   its own (no frame, no padding) and placed back as a Raster primitive spanning the exact
   plot range: Raster is the one primitive that exports as an embedded image, so
   everything drawn around and over it remains vector. *)

outFiles = If[Length[$ScriptCommandLine] > 1, Rest[$ScriptCommandLine], {"densityplot.pdf", "densityplot.svg"}];

f[x_, y_] := Sin[x] Cos[y] + 0.3 Sin[x y / 2];
colors = ColorData["SunsetColors"];
{{xmin, xmax}, {ymin, ymax}} = {{-2 Pi, 2 Pi}, {-2 Pi, 2 Pi}};
rasterPixels = 1000; (* density map resolution, ~360 dpi at the 200 pt figure width *)

(* Bare density map: no frame, axes or padding, so the plot range fills the whole image. *)
densityMap = DensityPlot[
  f[x, y],
  {x, xmin, xmax}, {y, ymin, ymax},
  PlotPoints -> 80,
  ColorFunction -> colors,
  Frame -> False, Axes -> False,
  PlotRange -> {{xmin, xmax}, {ymin, ymax}, All},
  PlotRangePadding -> None, ImagePadding -> None, AspectRatio -> Full
];
{zmin, zmax} = MinMax[Cases[densityMap, GraphicsComplex[pts_, ___] :> f @@@ pts, Infinity]];

(* Image rows run top-down, Raster rows bottom-up, hence the Reverse. *)
img = Rasterize[densityMap, "Image",
  RasterSize -> {rasterPixels, rasterPixels}, ImageSize -> {rasterPixels, rasterPixels}];
densityRaster = Raster[Reverse @ ImageData[img], {{xmin, ymin}, {xmax, ymax}}];

(* Contour lines drawn over the density map (lines only, no shading). *)
contours = ContourPlot[
  f[x, y],
  {x, xmin, xmax}, {y, ymin, ymax},
  PlotPoints -> 80,
  Contours -> Subdivide[zmin, zmax, 12][[2 ;; -2]],
  ContourShading -> None,
  ContourStyle -> Directive[Black, Opacity[0.6], AbsoluteThickness[0.5]]
];

(* Vector frame, ticks and labels around the raster. *)
plot = Show[
  Graphics[densityRaster,
    Frame -> True,
    FrameLabel -> {"x", "y"},
    PlotLabel -> "Sin[x] Cos[y] + 0.3 Sin[x y/2]",
    PlotRange -> {{xmin, xmax}, {ymin, ymax}}, PlotRangePadding -> None,
    AspectRatio -> 1, ImageSize -> 200
  ],
  contours
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
