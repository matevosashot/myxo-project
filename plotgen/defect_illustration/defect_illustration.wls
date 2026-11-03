#!/usr/bin/env wolframscript
(* defect_illustration.wls -- two small icons of a +1/2 and a -1/2 defect: director lines
   on a white background, with the defect symbol of fig:defects-layering(b) at the core.
   Written 2026-10-06; background replaced by the symbols 2026-10-26.

   Director lines are exact: for theta = phi/2 the integral curves are the level sets of
   Im Sqrt[z] = Sqrt[r] Sin[phi/2] (parabolas with focus at the core), and for
   theta = -phi/2 those of Im z^(3/2) = r^(3/2) Sin[3 phi/2] (three hyperbola-like
   families), with phi in [0, 2 Pi). The radial rays where the director points away from
   the core (one for +1/2, three for -1/2) are drawn separately.

   Symbols as in Copenhagen et al., Nat. Phys. 17, 211 (2021), Fig. 2d,e: a disk at the
   core and segments along the radial rays (one for +1/2, three for -1/2), outlined in
   white. Colours are read off the markers in figures/handmade/main.pdf p.2, panel (b).

   Usage:  module load mathematica/14.10; unset DISPLAY
           wolframscript -file plotgen/defect_illustration/defect_illustration.wls
   Outputs: defect_plus.{pdf,png}, defect_minus.{pdf,png} next to this script. *)

SetEnvironment["QT_QPA_PLATFORM" -> "offscreen"];
SetEnvironment["GALLIUM_DRIVER" -> "softpipe"];   (* llvmpipe crashes the front end in Rasterize *)
outDir = DirectoryName[$InputFileName];

L = 1.;              (* half-width of the square *)
size = 100;          (* icon width in pt *)

phi[x_, y_] := Mod[ArcTan[x, y], 2 Pi];

(* titles via MaTeX: the headless front end writes Helvetica text with a shifted glyph
   mapping ("+1/2" comes out as "+/ /0" in pdf viewers) *)
<< MaTeX`
title[s_] := MaTeX["\\text{" <> s <> "}", FontSize -> 14];

lineStyle = Directive[Black, AbsoluteThickness[0.6]];
director[g_, levels_] := ContourPlot[g[x, y], {x, -L, L}, {y, -L, L},
   Contours -> levels, ContourShading -> None, ContourStyle -> lineStyle,
   PlotPoints -> 60, MaxRecursion -> 3,
   Exclusions -> None];   (* the ArcTan jump on the -x axis is not a jump of g *)
rays[angles_] := Graphics[{lineStyle, Line[{{0, 0}, 2 L {Cos[#], Sin[#]}}] & /@ angles}];

(* defect symbol: disk of radius rDisk at the core and arms of length armLen and width
   armW along the given angles; the white outline (width rim) separates it from the lines *)
rDisk = 0.075; armW = 0.065; rim = 0.018;
arm[len_, a_, d_] := Rotate[Rectangle[{0, -armW/2 - d}, {len + d, armW/2 + d}], a, {0, 0}];
symbol[col_, len_, angles_] := {
   White, Disk[{0, 0}, rDisk + rim], arm[len, #, rim] & /@ angles,
   col, Disk[{0, 0}, rDisk], arm[len, #, 0] & /@ angles};
red = RGBColor[250/255, 2/255, 4/255];
blue = RGBColor[2/255, 106/255, 250/255];

icon[label_, g_, levels_, rayAngles_, marks_] := Show[
   Graphics[{}, PlotRange -> {{-L, L}, {-L, L}}, PlotRangePadding -> None,
    Frame -> False, AspectRatio -> 1, ImageSize -> size, ImagePadding -> {{2, 2}, {2, 2}},
    PlotLabel -> title[label]],
   director[g, levels], rays[rayAngles], Graphics[marks],
   PlotRange -> {{-L, L}, {-L, L}}, PlotRangeClipping -> True];

(* +1/2: vertex of each parabola on the negative x axis at x = -c^2 *)
gPlus[x_, y_] := Sqrt[Sqrt[x^2 + y^2]] Sin[phi[x, y]/2];
(* head (rounded side, x < 0) and tail (x > 0). "tail" sits on the x axis beyond the end
   of the symbol's bar, which replaces the ray; "head" runs vertically between the
   parabolas through x = -0.28 and -0.62 *)
mark[s_] := MaTeX["\\text{" <> s <> "}", FontSize -> 14];
unitsPerPt = 2 L/(size - 4);                       (* plot width over its pt width *)
tailHalf = unitsPerPt (First@ImageDimensions@Rasterize[mark["tail"], ImageResolution -> 72])/2;
tailX = 0.92 L - tailHalf;                         (* label ends 0.08 inside the edge *)
barLen = tailX - tailHalf - 0.06;
plusMarks = {
   symbol[red, barLen, {0}],
   Inset[mark["tail"], {tailX, 0}],
   Inset[Rotate[mark["head"], 90 Degree], {-0.42, 0}]};
plus = icon["$+$1/2 defects", gPlus, Sqrt[{0.08, 0.28, 0.62, 0.92}], {}, plusMarks];
Print["+1/2 bar length = ", barLen, " (half-width L = ", L, ")"];

(* -1/2: on the negative x axis |Im z^(3/2)| = |x|^(3/2), crossings evenly spaced; the
   arms have the length of the +1/2 bar *)
gMinus[x_, y_] := (x^2 + y^2)^(3/4) Abs@Sin[3 phi[x, y]/2];
minus = icon["$-$1/2 defects", gMinus, Range[0.3, 1.8, 0.25]^(3/2), {0, 2 Pi/3, 4 Pi/3},
   symbol[blue, barLen, {0, 2 Pi/3, 4 Pi/3}]];

Scan[Function[{nf},
   Export[FileNameJoin[{outDir, nf[[1]] <> ".pdf"}], nf[[2]], "AllowRasterization" -> False];
   Export[FileNameJoin[{outDir, nf[[1]] <> ".png"}], nf[[2]], ImageResolution -> 600];
   Print["Wrote ", nf[[1]]]],
  {{"defect_plus", plus}, {"defect_minus", minus}}];
