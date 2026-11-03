#!/usr/bin/env wolframscript
(* aimed_rods_exact.wls -- the manuscript version of aimed_rods.wls (fig:aimed-rods, App F,
   app:variance-neighbors; written 2026-10-22 at his request "let this be another script").
   Identical to aimed_rods.wls except for the profile and the wording: S(r) is the accurate
   Pade approximant eq:pade-accurate (within 1e-4 S0 of the exact solution of H_a = 0), so the
   modulation in the titles is the paper's +-2.94% at 1.5 um, and the text says "point" where
   aimed_rods.wls says "spot". aimed_rods.wls (simple Pade, +-3.06%) is left unchanged.

   eq:variance-neighbors (eq:total-gauss integrated by parts twice) is an average over the
   neighbours of a point r:
     Var/Var_pedestal - 1 = Int d^2w Omega(w) w^.Q_ss(r+w).w^ = < S cos 2(theta_nb - omega) >_Omega,
   Omega = (l_n^2/4)(kappa'' - kappa'/w) > 0, Int Omega d^2w = 1 (half of it within 0.68 um,
   94% within 2 um). A rod whose axis line passes through the point ("aimed", cos > 0) raises
   the variance there, a rod lying across the line to the point ("side-on", cos < 0) lowers it.
   The legend states the rule exactly: red = within 45 deg of the line to the point (cos 2 alpha > 0),
   blue = more than 45 deg from it (changed 2026-10-23 at his request; the old legend said "axis line
   passes through it", which holds only at alpha = 0).

     (a) uniform nematic, director along x: aimed and side-on rods balance, modulation 0;
     (b) a point on the head side, 1.5 um from the core: toward and beyond the core the rods
         turn to lie along the axis and become aimed at the point, so the variance rises;
     (c) a point on the tail side, 1.5 um from the core: beyond the core the rods lie across
         the axis and become side-on, so the variance drops.

   Rods: the +1/2 profile theta = phi/2 with S(r) of eq:pade-accurate (S0 = 0.9, l_d = 1 um),
   on a square grid of spacing 0.34 um about the point; length proportional to S, colour by
   the sign of cos 2(theta - omega), opacity by its size; faint outside the dashed 2 um circle.
   Colours are the two ends of TemperatureMap, as in variance_map.wls (fig:variance3d).
   The modulation in the titles is computed here (expected +-2.942%, scratchpad
   round1022/checks_round.py).

   All text is MaTeX (Computer Modern): the headless front end garbles plain text in pdf.
   The figure is 510 wide; \includegraphics[width=\textwidth] prints every size 1:1.
   2026-10-23, his prompt "Increase the font size in this figure. Note that i used 0.8
   textwidth": the paper prints the figure at 0.8\textwidth, so every size prints at 0.8 of its
   design value. The fonts were raised from 8-9 to 10-11 design points (titles, axis labels and
   legend 11 -> 8.8 pt printed, tick labels and "core" 10 -> 8 pt), the paddings grown to
   match, and the legend broken into two rows to keep the design width 510.

   Usage:  plotgen/aimed_rods/run_exact.sh [out=<dir>]
   Outputs: aimed_rods_exact.pdf and aimed_rods_exact.png (600 dpi) next to this script;
   run_exact.sh copies the pdf to manuscript/figures/. *)

SetEnvironment["QT_QPA_PLATFORM" -> "offscreen"];
$HistoryLength = 0;

here = DirectoryName[$InputFileName];
opts = Association[Rule @@@ (StringSplit[#, "=", 2] & /@
      Select[Rest@$ScriptCommandLine, StringContainsQ[#, "="] &])];
outDir = Lookup[opts, "out", here];
t0 = AbsoluteTime[];

(* ---------------------------------------------------------------- model *)
s0 = 0.9; ld = 1.; ln = 0.5;                    (* tab:field-params *)
(* eq:pade-accurate *)
sPade[r_] := With[{x = r/ld}, s0 x Sqrt[(0.34011 + 0.04517 x^2 + 0.00234 x^4)/
    (1 + 0.38387 x^2 + 0.04751 x^4 + 0.00234 x^6)]];
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
(* font sizes in design points; printed at 0.8\textwidth they shrink by 0.8 *)
fsLabel = 11;               (* titles with the panel letters, axis labels, legend: 8.8 pt printed *)
fsTick = 10;                (* tick labels and "core": 8 pt printed *)
tex[s_, fs_ : fsTick] := MaTeX[s, FontSize -> fs, "Preamble" -> {"\\usepackage{amsmath,upgreek}"}];
cmap = ColorData["TemperatureMap"];
{cAim, cSide} = {cmap[0.97], cmap[0.03]};       (* red raises, blue lowers the variance *)

L = 3.2;                    (* half-width of each window about its spot, um *)
step = 0.34; lRod = 0.24;   (* rod grid spacing and full rod length at S = S0, um *)
rWeight = 2.;               (* the dashed circle: 94% of W lies inside *)

(* layout, in design points; the figure is totalW wide *)
totalW = 510;
titleY = 10; padTop = 20;
labelOff = 24;              (* distance of the axis labels' centres from the frame *)
padA = {{33, 3}, {32, padTop}}; padB = {{8, 3}, {32, padTop}};   (* {{left, right}, {bottom, top}} *)
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
   sizes in um, fitted to the 10 pt "core" (scaled from the 8 pt fit 2026-10-23). The label
   is shifted by dx toward the point, since at 10 pt it would otherwise touch the dashed
   circle, which passes 0.45 um beside the core at the label's height *)
coreY = 0.45;
core[dx_] := {{White, Opacity[0.75], Rectangle[{dx - 0.47, coreY - 0.15}, {dx + 0.47, coreY + 0.15},
      RoundingRadius -> 0.08]},
   Inset[tex["\\text{core}"], {dx, coreY}],
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

xLabel[lab_, pt_, h_] := Inset[tex[lab, fsLabel], Scaled[{0.5, -pt/h}]];
yLabel[lab_, pt_, w_] := Inset[tex[lab, fsLabel], Scaled[{-pt/w, 0.5}], Center, Automatic, {0, 1}];
title[lab_, pt_, h_] := Inset[tex[lab, fsLabel], Scaled[{0.5, 1 + pt/h}]];
ticksIn[{a_, b_}] := {#, tex[ToString[#]]} & /@ Select[Range[-6, 6, 2], a <= # <= b &];
noLabels[{a_, b_}] := {#, ""} & /@ Select[Range[-6, 6, 2], a <= # <= b &];
pct[v_] := If[v > 0, "{+}", "{-}"] <> ToString[NumberForm[Abs[100 v], {3, 1}]] <> "\\%";

(* showY: y tick labels and y label; showX: x label (x tick labels are always drawn, since
   the x range differs from panel to panel) *)
panel[x0_, kind_, label_, showY_, showX_, pad_] := Module[{xr = {x0 - L, x0 + L}, yr = {-L, L}},
   Graphics[{streams[x0, kind === "uniform"], rods[x0, kind === "uniform"], weightCircle[x0],
     If[kind === "uniform", {}, core[0.1 Sign[x0]]], spot[x0], title[label, titleY, s],
     If[showX, xLabel["x~(\\upmu\\mathrm{m})", labelOff, s], {}],
     If[showY, yLabel["y~(\\upmu\\mathrm{m})", labelOff, s], {}]},
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

sample[g_] := Graphics[g, ImageSize -> {22, 10}, PlotRange -> {{-1, 1}, {-1, 1}},
   ImagePadding -> 0, AspectRatio -> Full];
legendItems = {
   {sample[{CapForm["Round"], AbsoluteThickness[1.3], cAim, Line[{{-0.8, 0}, {0.8, 0}}]}],
    "\\text{within }45^\\circ\\text{ of the line to the point: raises the variance}"},
   {sample[{CapForm["Round"], AbsoluteThickness[1.3], cSide, Line[{{-0.8, 0}, {0.8, 0}}]}],
    "\\text{beyond }45^\\circ\\text{: lowers it}"},
   {sample[{GrayLevel[0.45], AbsoluteThickness[0.6], AbsoluteDashing[{3, 2}],
      Line[{{-1, 0}, {1, 0}}]}],
    "2~\\upmu\\mathrm{m}\\text{, 94\\% of the weight}"}};
(* two rows (2026-10-23): at 11 pt the three items no longer fit in one row of 510; the
   colour rule (red, blue) on the first, the dashed circle under the red item *)
legend = With[{it = {#1, tex[#2, fsLabel]} & @@@ legendItems},
   Grid[{Join[it[[1]], {Spacer[18]}, it[[2]]], Join[it[[3]], {"", "", ""}]},
    Alignment -> {Left, Center}, Spacings -> {0.4, 0.15}]];

fig = Grid[{panels, {legend, SpanFromLeft, SpanFromLeft}},
   Spacings -> {0, {0, 0.3}}, Alignment -> {Center, Top}];

Export[FileNameJoin[{outDir, "aimed_rods_exact.pdf"}], fig, "AllowRasterization" -> False];
Export[FileNameJoin[{outDir, "aimed_rods_exact.png"}], fig, ImageResolution -> 600];
Print["exported to ", outDir, "  (", Round[AbsoluteTime[] - t0], " s)"];
