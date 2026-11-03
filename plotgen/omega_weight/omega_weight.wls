#!/usr/bin/env wolframscript
(* omega_weight.wls -- the weight Omega of eq:variance-neighbors (app:variance-neighbors), as a
   small single-column figure. Written 2026-10-23 for his prompt "Let's have a small
   visualisation of the weight Omega, write the code in plotgen, I prefer mathematica script".

   eq:variance-neighbors writes the modulation of the density variance at a point r as an
   average over its neighbours r + w:
     Var/Var_pedestal - 1 = Int d^2d Omega(d) S cos 2(theta_nb - theta_d),
     Omega(d) = (l_n^2/(2 Pi d^4)) [2 - e^-X (X^2 + 2X + 2)],   X = 2 d^2/l_n^2,  l_n = 0.5 um,
   theta_nb the director angle at the neighbour, omega the angle of w.
   The figure shows the radial weight 2 Pi d Omega(d) (um^-1), whose area is 1, against d.
   2026-10-23: offset renamed w -> d (his choice; bold w_ss is a velocity in the same section).

   With g(X) = 2 - e^-X (X^2 + 2X + 2), g'(X) = X^2 e^-X, and d^2w = 2 Pi w dw = (Pi l_n^2/2) dX,
   the weight within w is in closed form
     F(w) = Int_0^X g(t)/t^2 dt = 1 - (1 + X) e^-X - g(X)/X,
   so Int Omega d^2w = F(infinity) = 1. Near w = 0, g = X^3/3 + ..., Omega = 4 w^2/(3 Pi l_n^4)
   (the series avoids the cancellation in g at small X); for w >> l_n, Omega -> l_n^2/(Pi w^4).
   The script checks these numbers and the ones the text quotes (half of the weight within
   0.68 um, 94% within 2 um, a quarter beyond 1 um) and prints them.

   Shaded: the half of the weight within the median radius. Dashed: the median radius (50%)
   and w = 2 um (94%, the dashed circle of fig:aimed-rods). The tick l_n marks the cutoff.

   All text is MaTeX (Computer Modern): the headless front end garbles plain text in pdf.
   Sizes are design points. The figure is 246 wide (the revtex \linewidth); the pdf export
   scales by 0.72, so \includegraphics[width=\linewidth] prints every size 1:1.

   Usage:  plotgen/omega_weight/run.sh [out=<dir>]
   Outputs: omega_weight.pdf and omega_weight.png (600 dpi) next to this script;
   run.sh copies the pdf to manuscript/figures/. *)

SetEnvironment["QT_QPA_PLATFORM" -> "offscreen"];
$HistoryLength = 0;

here = DirectoryName[$InputFileName];
opts = Association[Rule @@@ (StringSplit[#, "=", 2] & /@
      Select[Rest@$ScriptCommandLine, StringContainsQ[#, "="] &])];
outDir = Lookup[opts, "out", here];
t0 = AbsoluteTime[];

(* ---------------------------------------------------------------- the weight *)
ln = 0.5;                                              (* l_n, um (tab:field-params) *)
xOf[w_] := 2 w^2/ln^2;
ex[x_] := If[x > 600, 0., Exp[-x]];                  (* no underflow messages far out *)
gX[x_] := If[x < 10^-2, x^3/3 - x^4/4 + x^5/10, 2 - ex[x] (x^2 + 2 x + 2)];
omega[w_?NumericQ] := With[{x = xOf[w]},
   If[x < 10^-2, 4 w^2/(3 Pi ln^4) (1 - 3 x/4 + 3 x^2/10), ln^2/(2 Pi w^4) gX[x]]];
radial[w_?NumericQ] := 2 Pi w omega[w];               (* um^-1, area 1 *)
frac[w_?NumericQ] := With[{x = xOf[w]},                (* weight within w, closed form *)
   If[x < 10^-2, x^2/6 - x^3/12, 1 - (1 + x) ex[x] - gX[x]/x]];

(* checks *)
fmt[v_, n_] := ToString[NumberForm[v, n]];
exactX = Simplify[2 - Exp[-x] (x^2 + 2 x + 2)];
Print["  closed form: Int_0^Infinity g(t)/t^2 dt = ",
  Integrate[exactX/x^2, {x, 0, Infinity}], ";  Omega(0) = ",
  Limit[ln^2/(2 Pi w^4) (exactX /. x -> 2 w^2/ln^2), w -> 0]];
norm = NIntegrate[radial[w], {w, 0, 0.5, 1, 2, 5, 20, Infinity}, PrecisionGoal -> 10];
Print["  normalization Int 2 Pi w Omega dw = ", fmt[norm, 10], "  (minus 1: ", norm - 1, ")"];
wMed = w /. FindRoot[frac[w] == 1/2, {w, 0.7}];
Print["  median radius (50% within): ", fmt[wMed, 5], " um;  NIntegrate up to it: ",
  fmt[NIntegrate[radial[w], {w, 0, wMed}], 6]];
Print["  within 1 um: ", fmt[100 frac[1.], 4], "%,  within 2 um: ",
  fmt[100 frac[2.], 4], "%  (NIntegrate: ", fmt[100 NIntegrate[radial[w], {w, 0, 2}], 4],
  "%),  within 3 um: ", fmt[100 frac[3.], 4], "%"];
With[{wS = ln Sqrt[0.01/2]},                          (* the switch to the series, X = 0.01 *)
  Print["  min of Omega on (0, 10] um: ", Min[omega /@ Subdivide[0.001, 10., 5000]],
   " (positive);  series/exact - 1 at X = 0.01: ",
   N[(4 w^2/(3 Pi ln^4) (1 - 3 x/4 + 3 x^2/10))/(ln^2/(2 Pi w^4) (2 - Exp[-x] (x^2 + 2 x + 2))) - 1 /.
      {x -> 1/100, w -> 1/2 Sqrt[1/200], ln -> 1/2}, 5]]];
Print["  Omega / (l_n^2/(Pi w^4)) at w = 1, 2, 3, 5 um: ",
  fmt[omega[#]/(ln^2/(Pi #^4)), 5] & /@ {1., 2., 3., 5.}];
{peakVal, wPeak} = Module[{s = FindMaximum[radial[w], {w, 0.5}]}, {s[[1]], w /. s[[2]]}];
Print["  peak of 2 Pi w Omega: ", fmt[peakVal, 4], " um^-1 at w = ", fmt[wPeak, 4],
  " um"];

(* ---------------------------------------------------------------- graphics *)
<< MaTeX`
tex[s_, fs_ : 8] := MaTeX[s, FontSize -> fs, "Preamble" -> {"\\usepackage{amsmath,upgreek}"}];

wMax = 3.; yMax = 1.5;
curveColor = RGBColor[0.20, 0.40, 0.70];
dashStyle = Directive[GrayLevel[0.45], AbsoluteThickness[1.], AbsoluteDashing[{3, 2}]];

(* layout, in design points; the figure is totalW wide *)
totalW = 246; totalH = 160;
pad = {{38, 8}, {29, 7}};                             (* {{left, right}, {bottom, top}} *)
{fw, fh} = {totalW - Total[pad[[1]]], totalH - Total[pad[[2]]]};
Print["  frame ", fw, " x ", fh, " pt"];

ws = Join[Subdivide[0., 1.5, 600], Rest@Subdivide[1.5, wMax, 300]];
curve = {#, radial[#]} & /@ ws;
shade = Join[{{0., 0.}}, Select[curve, #[[1]] <= wMed &], {{wMed, radial[wMed]}, {wMed, 0.}}];

(* labels as insets a fixed distance (pt) from the frame, as in aimed_rods_exact.wls *)
xLabel[lab_, pt_] := Inset[tex[lab, 9], Scaled[{0.5, -pt/fh}]];
yLabel[lab_, pt_] := Inset[tex[lab, 9], Scaled[{-pt/fw, 0.5}], Center, Automatic, {0, 1}];
mark[w_, lab_] := {{dashStyle, Line[{{w, 0}, {w, yMax}}]},
   Inset[tex[lab, 8], {w, 0.93 yMax}, {Left, Center}, Automatic, {1, 0}]};

xTicks = Join[{#, tex[ToString[#]]} & /@ Range[0, 3],
   {{ln, tex["\\ell_\\mathrm{n}"]}}];
yTicks = {#, tex[ToString[NumberForm[#, {2, 1}]]]} & /@ {0.5, 1., 1.5};
yTicks = Prepend[yTicks, {0, tex["0"]}];

fig = Graphics[{
    {curveColor, Opacity[0.22], EdgeForm[], Polygon[shade]},
    {curveColor, AbsoluteThickness[1.5], CapForm["Round"], JoinForm["Round"], Line[curve]},
    mark[wMed, "\\;50\\%"], mark[2., "\\;94\\%"],
    xLabel["d~(\\upmu\\mathrm{m})", 21],
    yLabel["2\\pi d\\,\\Omega(d)~(\\upmu\\mathrm{m}^{-1})", 27]},
   Frame -> True, FrameStyle -> Directive[Black, AbsoluteThickness[0.5]],
   FrameTicks -> {{yTicks, None}, {xTicks, None}},
   PlotRange -> {{0, wMax}, {0, yMax}}, PlotRangePadding -> None, PlotRangeClipping -> False,
   AspectRatio -> Full, ImagePadding -> pad, ImageSize -> {totalW, totalH}];

Export[FileNameJoin[{outDir, "omega_weight.pdf"}], fig, "AllowRasterization" -> False];
Export[FileNameJoin[{outDir, "omega_weight.png"}], fig, ImageResolution -> 600];
Print["exported to ", outDir, "  (", Round[AbsoluteTime[] - t0], " s)"];
