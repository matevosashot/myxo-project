#!/usr/bin/env wolframscript
(* aimed_rods.wls -- the physical picture behind the head/tail asymmetry of the density
   variance at a +1/2 defect (written 2026-10-22 from the figure of the "Four questions on
   the model" artifact, Question 3).

   eq:total-gauss, integrated by parts twice, is an average over the neighbours of a spot r:
     Var/Var_pedestal - 1 = Int d^2w W(w) w^.Q_ss(r+w).w^ = < S cos 2(theta_nb - omega) >_W,
   W = (l_n^2/4)(kappa'' - kappa'/w) > 0, Int W d^2w = 1 (half of W within 0.68 um, 94% within
   2 um). A rod whose axis line passes through the spot ("aimed", cos > 0) raises the variance
   there, a rod lying across the line to the spot ("side-on", cos < 0) lowers it.

     (a) uniform nematic, director along x: aimed and side-on rods balance, modulation 0;
     (b) a spot on the head side, 1.5 um from the core: toward the core the rods turn to lie
         along the axis and become aimed at the spot, so the variance rises;
     (c) a spot on the tail side, 1.5 um from the core: toward the core the rods turn to lie
         across the axis and become side-on, so the variance drops.

   Rods: the +1/2 profile theta = phi/2 with the simple Pade S(r) of eq:pade-approximation
   (S0 = 0.9, l_d = 1 um), on a square grid of spacing 0.34 um about the spot; length
   proportional to S, colour by the sign of cos 2(theta - omega) (omega the direction from
   the spot to the rod), opacity by its size. Rods outside the dashed circle of radius 2 um
   (94% of the weight W) are drawn faint. Colours are the two ends of TemperatureMap, as in
   variance_map.wls (fig:variance3d): red raises, blue lowers the variance.
   The modulation in the titles is computed here from the formula above on the same profile
   (expected +-3.06%, scratchpad feedback/checks.py; the paper's 2.94% is for the exact
   profile of eq:pade-accurate).

   All text is MaTeX (Computer Modern): the headless front end garbles plain text in pdf.
   Sizes are design points. The figure is 510 wide; the pdf export scales by 0.72, so
   \includegraphics[width=\textwidth] prints every size 1:1, as for variance_map.pdf.
   aimed_rods_margin.wls draws the same three panels stacked, for the margin of a book.

   Usage:  plotgen/aimed_rods/run.sh [out=<dir>]
   Outputs: aimed_rods.pdf and aimed_rods.png (600 dpi) next to this script. *)

SetEnvironment["QT_QPA_PLATFORM" -> "offscreen"];
$HistoryLength = 0;

here = DirectoryName[$InputFileName];
opts = Association[Rule @@@ (StringSplit[#, "=", 2] & /@
      Select[Rest@$ScriptCommandLine, StringContainsQ[#, "="] &])];
outDir = Lookup[opts, "out", here];
t0 = AbsoluteTime[];

(* ---------------------------------------------------------------- model *)
s0 = 0.9; ld = 1.; ln = 0.5;                    (* tab:field-params *)
sPade[r_] := With[{x = r/ld}, s0 x Sqrt[(0.34 + 0.07 x^2)/(1 + 0.41 x^2 + 0.07 x^4)]];
(* W(w) = (l_n^2/4) K(w), K = kappa'' - kappa'/w = (2/Pi)(2 - e^-x (x^2 + 2x + 2))/w^4,
   x = 2 w^2/l_n^2; the series K = 16 w^2/(3 Pi l_n^6) avoids the cancellation at small w *)
kW[w_?NumericQ] := With[{x = 2 w^2/ln^2},
   If[x < 10^-3, 16 w^2/(3 Pi ln^6), 2/Pi (2 - Exp[-x] (x^2 + 2 x + 2))/w^4]];
(* modulation at a spot (x0, 0) of the +1/2 defect: S cos(phi - 2 omega), phi the polar angle
   of the neighbour about the core, theta = phi/2 *)
nu[x0_] := ln^2/4 NIntegrate[
    w kW[w] sPade[Sqrt[(x0 + w Cos[om])^2 + (w Sin[om])^2]]
      Cos[ArcTan[x0 + w Cos[om], w Sin[om]] - 2 om],
    {w, 0, 0.5, 1, 2, 4, 10, 40, 200}, {om, 0, Pi/2, Pi, 3 Pi/2, 2 Pi},
    PrecisionGoal -> 6, MaxRecursion -> 20];
xSpot = 1.5;
nuHead = nu[-xSpot]; nuTail = nu[xSpot];
Print["  modulation at x = -", xSpot, " um: ", 100 nuHead, "%,  at x = +", xSpot, " um: ",
  100 nuTail, "%   (", Round[AbsoluteTime[] - t0], " s)"];

(* ---------------------------------------------------------------- graphics *)
<< MaTeX`
tex[s_, fs_ : 8] := MaTeX[s, FontSize -> fs, "Preamble" -> {"\\usepackage{amsmath,upgreek}"}];
cmap = ColorData["TemperatureMap"];
{cAim, cSide} = {cmap[0.97], cmap[0.03]};       (* red raises, blue lowers the variance *)

L = 3.2;                    (* half-width of each window about its spot, um *)
step = 0.34; lRod = 0.24;   (* rod grid spacing and full rod length at S = S0, um *)
rWeight = 2.;               (* the dashed circle: 94% of W lies inside *)

(* layout, in design points; the figure is totalW wide *)
totalW = 510;
titleY = 8; padTop = 16;
padA = {{25, 3}, {25, padTop}}; padB = {{8, 3}, {25, padTop}};   (* {{left, right}, {bottom, top}} *)
s = (totalW - Total[padA[[1]]] - 2 Total[padB[[1]]])/3;           (* side of each frame *)
Print["  frame side ", s, " pt"];

(* each rod sits on a white casing slightly wider than itself, so the streamlines pass behind
   the rods even where a faint rod is mostly transparent *)
rods[x0_, uniform_] := Module[{g = Range[-L + step/2, L, step], out = {}, casing = {}},
   Do[Module[{w = Norm[{wx, wy}], X = x0 + wx, Y = wy, th, sv, c, op, d, seg},
      If[w > 0.12,
       {th, sv} = If[uniform, {0., s0}, {ArcTan[X, Y]/2, sPade[Norm[{X, Y}]]}];
       c = Cos[2 (th - ArcTan[wx, wy])];
       op = If[w <= rWeight, 0.22 + 0.78 Abs[c], 0.10 + 0.22 Abs[c]];
       d = 0.5 lRod sv/s0 {Cos[th], Sin[th]};
       seg = Line[{{X, Y} - d, {X, Y} + d}];
       AppendTo[casing, seg];
       AppendTo[out, {If[c > 0, cAim, cSide], Opacity[op], seg}]]],
    {wx, g}, {wy, g}];
   {CapForm["Round"], {White, AbsoluteThickness[2.6], casing}, AbsoluteThickness[1.3], out}];

spot[x0_] := {{White, Disk[{x0, 0}, 0.17]}, {Black, AbsoluteThickness[1.1], Circle[{x0, 0}, 0.13]}};
(* label on a small translucent white plate with rounded corners (it sits on the
   streamlines and rods), drawn before the marker so the plate does not cover it; plate
   sizes in um, fitted to the 8 pt "core" *)
coreY = 0.40;
core = {{White, Opacity[0.75], Rectangle[{-0.38, coreY - 0.12}, {0.38, coreY + 0.12},
      RoundingRadius -> 0.07]},
   Inset[tex["\\text{core}"], {0, coreY}],
   Black, Polygon[{{0.11, 0}, {0, 0.11}, {-0.11, 0}, {0, -0.11}}]};
weightCircle[x0_] := {GrayLevel[0.45], AbsoluteThickness[0.6], AbsoluteDashing[{3, 2}],
   Circle[{x0, 0}, rWeight]};

(* gray director streamlines behind the rods. +1/2: the level sets of Sqrt[r] Sin[phi/2] = c,
   the parabolas x = (y^2 - 4c^4)/(4c^2) that bend around the core on the head side (-x),
   plus the ray y = 0, x > 0, as in variance_map.wls; uniform: lines y = const. Cut to the
   window by hand, since PlotRangeClipping would also clip the axis labels. *)
streamC = Range[0.3, 2.4, 0.3];  streamDy = 0.68;
streamStyle = {GrayLevel[0.72], AbsoluteThickness[0.5]};   (* light: background only *)
streams[x0_, uniform_] := Module[{inW, para},
   inW[{x_, y_}] := Abs[x - x0] <= L && Abs[y] <= L;
   para[c_] := Line[Select[SplitBy[Table[{(y^2 - 4 c^4)/(4 c^2), y}, {y, Subdivide[-L, L, 800]}],
        inW], inW@First@# && Length[#] > 1 &]];
   (* the style must sit in the same list as the lines: directives do not reach out of a
      nested list *)
   Join[streamStyle, If[uniform,
     Line[{{x0 - L, #}, {x0 + L, #}}] & /@ N@Select[streamDy Range[-6, 6], Abs[#] <= L &],
     {para /@ streamC, Line[{{Max[0., x0 - L], 0}, {x0 + L, 0}}]}]]];

xLabel[lab_, pt_, h_] := Inset[tex[lab, 9], Scaled[{0.5, -pt/h}]];
yLabel[lab_, pt_, w_] := Inset[tex[lab, 9], Scaled[{-pt/w, 0.5}], Center, Automatic, {0, 1}];
title[lab_, pt_, h_] := Inset[tex[lab, 9], Scaled[{0.5, 1 + pt/h}]];
ticksIn[{a_, b_}] := {#, tex[ToString[#]]} & /@ Select[Range[-6, 6, 2], a <= # <= b &];
noLabels[{a_, b_}] := {#, ""} & /@ Select[Range[-6, 6, 2], a <= # <= b &];
pct[v_] := If[v > 0, "{+}", "{-}"] <> ToString[NumberForm[Abs[100 v], {3, 1}]] <> "\\%";

(* showY: y tick labels and y label; showX: x label (x tick labels are always drawn, since
   the x range differs from panel to panel) *)
panel[x0_, kind_, label_, showY_, showX_, pad_] := Module[{xr = {x0 - L, x0 + L}, yr = {-L, L}},
   Graphics[{streams[x0, kind === "uniform"], rods[x0, kind === "uniform"], weightCircle[x0],
     If[kind === "uniform", {}, core], spot[x0], title[label, titleY, s],
     If[showX, xLabel["x~(\\upmu\\mathrm{m})", 19, s], {}],
     If[showY, yLabel["y~(\\upmu\\mathrm{m})", 19, s], {}]},
    Frame -> True, FrameStyle -> Directive[Black, AbsoluteThickness[0.5]],
    FrameTicks -> {{If[showY, ticksIn[yr], noLabels[yr]], None}, {ticksIn[xr], None}},
    PlotRange -> {xr, yr}, PlotRangePadding -> None, PlotRangeClipping -> False,
    AspectRatio -> 1, ImagePadding -> pad,
    ImageSize -> {s + Total[pad[[1]]], s + Total[pad[[2]]]}]];

titles = {"\\textbf{(a)}\\enspace\\text{uniform nematic: } 0",
   "\\textbf{(b)}\\enspace\\text{head side: } " <> pct[nuHead],
   "\\textbf{(c)}\\enspace\\text{tail side: } " <> pct[nuTail]};
panels = {panel[0., "uniform", titles[[1]], True, True, padA],
   panel[-xSpot, "head", titles[[2]], False, True, padB],
   panel[xSpot, "tail", titles[[3]], False, True, padB]};

sample[g_] := Graphics[g, ImageSize -> {18, 8}, PlotRange -> {{-1, 1}, {-1, 1}},
   ImagePadding -> 0, AspectRatio -> Full];
legendItems = {
   {sample[{CapForm["Round"], AbsoluteThickness[1.3], cAim, Line[{{-0.8, 0}, {0.8, 0}}]}],
    "\\text{aimed at the spot: its axis line passes through it}"},
   {sample[{CapForm["Round"], AbsoluteThickness[1.3], cSide, Line[{{-0.8, 0}, {0.8, 0}}]}],
    "\\text{side-on: its axis line passes by}"},
   {sample[{GrayLevel[0.45], AbsoluteThickness[0.6], AbsoluteDashing[{3, 2}],
      Line[{{-1, 0}, {1, 0}}]}],
    "2~\\upmu\\mathrm{m}\\text{, 94\\% of the weight}"}};
legend = Grid[{With[{it = {#1, tex[#2]} & @@@ legendItems},
     Join[it[[1]], {Spacer[14]}, it[[2]], {Spacer[14]}, it[[3]]]]},
   Alignment -> {Left, Center}, Spacings -> {0.4, 0}];

fig = Grid[{panels, {legend, SpanFromLeft, SpanFromLeft}},
   Spacings -> {0, {0, 0.3}}, Alignment -> {Center, Top}];

Export[FileNameJoin[{outDir, "aimed_rods.pdf"}], fig, "AllowRasterization" -> False];
Export[FileNameJoin[{outDir, "aimed_rods.png"}], fig, ImageResolution -> 600];
Print["exported to ", outDir, "  (", Round[AbsoluteTime[] - t0], " s)"];
