#!/usr/bin/env wolframscript
(* defect_illustration.wls -- two small icons of a +1/2 and a -1/2 defect: director lines
   over a background that fades from white at the core to a saturated colour (the order
   parameter S rises from 0 at the core to S0 far away). Written 2026-10-06.

   Director lines are exact: for theta = phi/2 the integral curves are the level sets of
   Im Sqrt[z] = Sqrt[r] Sin[phi/2] (parabolas with focus at the core), and for
   theta = -phi/2 those of Im z^(3/2) = r^(3/2) Sin[3 phi/2] (three hyperbola-like
   families), with phi in [0, 2 Pi). The radial rays where the director points away from
   the core (one for +1/2, three for -1/2) are drawn separately.

   Only the background is a bitmap (a Raster spanning the plot range, as in
   plotgen/DensityPlotRasterize/densityplot.wls); lines, frame and title stay vector.

   Usage:  module load mathematica/14.10; unset DISPLAY
           wolframscript -file plotgen/defect_illustration/defect_illustration.wls
   Outputs: defect_plus.{pdf,png}, defect_minus.{pdf,png} next to this script. *)

SetEnvironment["QT_QPA_PLATFORM" -> "offscreen"];
SetEnvironment["GALLIUM_DRIVER" -> "softpipe"];   (* llvmpipe crashes the front end in Rasterize *)
outDir = DirectoryName[$InputFileName];

L = 1.;              (* half-width of the square *)
rc = 0.45;           (* core size of the colour fade *)
rasterPixels = 600;
size = 100;          (* icon width in pt *)

phi[x_, y_] := Mod[ArcTan[x, y], 2 Pi];
fade[x_, y_] := With[{r2 = x^2 + y^2}, Sqrt[r2/(r2 + rc^2)]];

background[col_] := Module[{map, img},
  map = DensityPlot[fade[x, y], {x, -L, L}, {y, -L, L}, PlotPoints -> 100,
    ColorFunction -> (Blend[{White, Lighter[col, 0.15], Darker[col, 0.2]}, #] &),
    ColorFunctionScaling -> False, Frame -> False, Axes -> False,
    PlotRange -> {{-L, L}, {-L, L}, All}, PlotRangePadding -> None,
    ImagePadding -> None, AspectRatio -> Full];
  img = Rasterize[map, "Image", RasterSize -> {rasterPixels, rasterPixels},
    ImageSize -> {rasterPixels, rasterPixels}];
  Raster[Reverse@ImageData[img], {{-L, -L}, {L, L}}]];

(* titles via MaTeX: the headless front end writes Helvetica text with a shifted glyph
   mapping ("+1/2" comes out as "+/ /0" in pdf viewers) *)
<< MaTeX`
title[s_] := MaTeX["\\text{" <> s <> "}", FontSize -> 14,
   "Preamble" -> {"\\usepackage[scaled]{helvet}", "\\usepackage{sansmath}", "\\sansmath",
     "\\renewcommand{\\familydefault}{\\sfdefault}"}];

title[s_] := MaTeX["\\text{" <> s <> "}", FontSize -> 14];

lineStyle = Directive[Black, AbsoluteThickness[0.6]];
director[g_, levels_] := ContourPlot[g[x, y], {x, -L, L}, {y, -L, L},
   Contours -> levels, ContourShading -> None, ContourStyle -> lineStyle,
   PlotPoints -> 60, MaxRecursion -> 3,
   Exclusions -> None];   (* the ArcTan jump on the -x axis is not a jump of g *)
rays[angles_] := Graphics[{lineStyle, Line[{{0, 0}, 2 L {Cos[#], Sin[#]}}] & /@ angles}];

icon[col_, label_, g_, levels_, rayAngles_, marks_ : {}] := Show[
   Graphics[background[col], PlotRange -> {{-L, L}, {-L, L}}, PlotRangePadding -> None,
    Frame -> False, AspectRatio -> 1, ImageSize -> size, ImagePadding -> {{2, 2}, {2, 2}},
    PlotLabel -> title[label]],
   director[g, levels], rays[rayAngles], Graphics[marks],
   PlotRange -> {{-L, L}, {-L, L}}, PlotRangeClipping -> True];

(* +1/2: vertex of each parabola on the negative x axis at x = -c^2 *)
gPlus[x_, y_] := Sqrt[Sqrt[x^2 + y^2]] Sin[phi[x, y]/2];
(* head (rounded side, x < 0) and tail (x > 0). "tail" sits on the x axis, with the ray
   interrupted behind it; "head" runs vertically between the parabolas through x = -0.28
   and -0.62 *)
mark[s_] := MaTeX["\\textcolor{white}{\\text{" <> s <> "}}", FontSize -> 14,
   "Preamble" -> {"\\usepackage{xcolor}"}];
unitsPerPt = 2 L/(size - 4);                       (* plot width over its pt width *)
tailX = 0.6;
tailHalf = unitsPerPt (First@ImageDimensions@Rasterize[mark["tail"], ImageResolution -> 72])/2 + 0.06;
plusMarks = {
   Inset[mark["tail"], {tailX, 0}],
   Inset[Rotate[mark["head"], 90 Degree], {-0.42, 0}],
   lineStyle, Line[{{0, 0}, {tailX - tailHalf, 0}}], Line[{{tailX + tailHalf, 0}, {2 L, 0}}]};
plus = icon[RGBColor[0.75, 0.08, 0.08], "$+$1/2 defects", gPlus,
   Sqrt[{0.08, 0.28, 0.62, 0.92}], {}, plusMarks];

(* -1/2: on the negative x axis |Im z^(3/2)| = |x|^(3/2), crossings evenly spaced *)
gMinus[x_, y_] := (x^2 + y^2)^(3/4) Abs@Sin[3 phi[x, y]/2];
minus = icon[RGBColor[0.1, 0.35, 0.68], "$-$1/2 defects", gMinus,
   Range[0.3, 1.8, 0.25]^(3/2), {0, 2 Pi/3, 4 Pi/3}];

Scan[Function[{nf},
   Export[FileNameJoin[{outDir, nf[[1]] <> ".pdf"}], nf[[2]], "AllowRasterization" -> False];
   Export[FileNameJoin[{outDir, nf[[1]] <> ".png"}], nf[[2]], ImageResolution -> 600];
   Print["Wrote ", nf[[1]]]],
  {{"defect_plus", plus}, {"defect_minus", minus}}];
