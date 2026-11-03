#!/usr/bin/env wolframscript
(* aimed_rods_qss.wls -- fig:aimed-rods of the manuscript since 2026-10-23 (his request:
   "regenerate it with Q_ss profile, the dotted line in Fig 8. Move to the main text").
   A copy of aimed_rods_exact.wls (left unchanged) with the nematic field of panels (b,c)
   replaced: the profile theta = phi/2, S(r) of eq:pade-accurate becomes Q_ss, the numerical
   steady state of the moving defect, the same field and the same quadrature as the dotted
   curves of fig:variance (plotgen/varrho/varrho.wls):
     steady_<tag>.m of data/fluctuations-set2_2_comove_optimized_2 (RhoNeumann density wall,
     u = 0.362 um/min, as in the published varrho.pdf; "ssDir" takes the Dirichlet one),
     101 x 101 graded grid on [-20, 20]^2, ListInterpolation order 3;
     eq:variance-neighbors with Q_ss inside |w| < Rc = 12.5 um and the static profile outside,
     Gauss-Legendre radial segments x 256 angles (makeRule/nuRule copied from varrho.wls).

   Panels: (a) uniform nematic, director along y as at the head-side point (since 2026-11-01;
   before: along x); (b) the point of largest variance on the head side
   and (c) of smallest variance on the tail side, found here with FindMaximum/FindMinimum on
   the axis (varrho.wls: +2.546% at x = -1.26, -3.394% at x = +1.66). The titles give the
   modulation computed here. Rods: director angle arg(Q1 + i Q2)/2 and length proportional
   to |Q| = S (full length at S0 = 0.9). Gray director lines: integrated through Q_ss (RK4,
   the sign of the director kept continuous along each line) from the points (-c^2, 0) on the
   head-side axis, c = 0.3, 0.6, ..., 2.4, where the parabolas of the static profile cross it,
   upward and downward, plus the ray y = 0, x > 0 (Q2 = 0 on the axis by mirror symmetry).
   The "core" label is moved toward the point just enough to clear the dashed circle.

   Core panel (2026-10-23, his request "add the forth one ... the center at the core. This
   should picture the asymmetry and possibly explain why on fig 8 we get non-zero nu at the
   origin", then "swap c and d", "include in the manuscript"): panel (c) has the point at the
   core, with a dotted line x = 0 and, in its upper corners, the parts of the core value from
   the neighbors on the head side (x < 0, blue) and on the tail side (x > 0, red); they cancel
   exactly on the static profile (+-11.98%), on Q_ss they are -11.64% and +12.26% (net
   +0.628%), mainly because S_ss is weaker on the head side (static angle with S_ss: +0.98%,
   angle_ss with static S: -0.36%; claude_experiments/core_asym.wls). The tail side is (d).
   The four panels are designed for \textwidth (design width 510 printed 1:1: titles, axis
   labels and legend 9 pt, ticks and "core" 8 pt, rods 1.1 pt). Option "nocore" draws the
   earlier three panels for 0.8\textwidth (aimed_rods_qss_3panel.pdf/png, not copied).

   Usage:  plotgen/aimed_rods/run_qss.sh [nocore] [ssDir] [Rc=<um>] [out=<dir>]
   Outputs: aimed_rods_qss.pdf and aimed_rods_qss.png (600 dpi) next to this script;
   run_qss.sh copies the pdf to manuscript/figures/. Everything else (layout, fonts, colours,
   legend, design width 510 for 0.8\textwidth) as in aimed_rods_exact.wls. *)

SetEnvironment["QT_QPA_PLATFORM" -> "offscreen"];
$HistoryLength = 0;

here = DirectoryName[$InputFileName];
args = Rest[$ScriptCommandLine];
ssDir = MemberQ[args, "ssDir"];
four = ! MemberQ[args, "nocore"];   (* the core panel (c); the tail side moves to (d) *)
opts = Association[Rule @@@ (StringSplit[#, "=", 2] & /@ Select[args, StringContainsQ[#, "="] &])];
rCut = If[KeyExistsQ[opts, "Rc"], N@ToExpression@opts["Rc"], 12.5];
repoRoot = ParentDirectory[ParentDirectory[here]];
outDir = Lookup[opts, "out", here];
t0 = AbsoluteTime[];

(* ---------------------------------------------------------------- model *)
s0 = 0.9; ln = 0.5;                             (* tab:field-params *)

(* Q_ss: the steady state behind the dotted curves of fig:variance (varrho.wls) *)
tag = "box20_N200_Dir_lN0.5_z12.5_ld1_B90000_ss20";
dataDir = FileNameJoin[{repoRoot, "data",
    "fluctuations-set2_2_comove_optimized_2" <> If[ssDir, "_ssDir", ""]}];
st = Import[FileNameJoin[{dataDir, "steady_" <> tag <> ".m"}]];
{lN, a1v, L1v, bv} = {lNoise, a1, L1, b} /. st["modelParams"];
Print["steady state ", FileNameTake[dataDir], ": density wall ", st["bc"], ", Q wall ",
  st["qWall"], ", u = ", st["u"][[1]], " um/min"];
{gx, gy} = N /@ {st["gx"], st["gy"]};
q1ss = ListInterpolation[st["Q1Grid"], {gx, gy}, InterpolationOrder -> 3];   (* [[i,j]] = f(gx_i, gy_j) *)
q2ss = ListInterpolation[st["Q2Grid"], {gx, gy}, InterpolationOrder -> 3];

(* the static profile eq:pade-accurate, used beyond Rc as in varrho.wls *)
fPade[z_] := Sqrt[(0.34010999 + 0.04516965 z + 0.00233656 z^2)/(1 + 0.38387160 z +
      0.04750621 z^2 + 0.00233656 z^3)];
SoR[r2_] := Sqrt[-a1v/(2 bv)] Sqrt[-a1v/L1v] fPade[-r2 a1v/L1v];   (* S(r)/r *)
q1Pade[X_, Y_] := SoR[X^2 + Y^2] X;
q2Pade[X_, Y_] := SoR[X^2 + Y^2] Y;

(* eq:variance-neighbors, Var/Var_ped - 1 = Int d^2w Omega(w) w^.Q(r + w).w^, as in varrho.wls *)
fX[X_] := Module[{s = UnitStep[0.01 - X], Xc = Clip[X, {0., 600.}]},
   (1 - s) (2 - Exp[-Xc] (Xc^2 + 2 Xc + 2)) + s (X^3/3 - X^4/4 + X^5/10 - X^6/36)];
omega[w_] := lN^2/(2 Pi w^4) fX[2 w^2/lN^2];
glNodes[n_] := glNodes[n] = Module[{x, ws},
   x = Sort[x /. NSolve[LegendreP[n, x] == 0, x, WorkingPrecision -> 30]];
   ws = Map[Function[x0, 2/((1 - x0^2) (D[LegendreP[n, xx], xx] /. xx -> x0)^2)], x];
   N[{x, ws}]];
breaks = {0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 10, 12.5, 15, 20, 30, 50, 100, 200, 400};
makeRule[r1_, r2_, nr_ : 24, nth_ : 256] := Module[{edges, x, ws, rr, rw, th},
   edges = Union[{r1, r2}, Select[breaks, r1 < # < r2 &]];
   {x, ws} = glNodes[nr];
   {rr, rw} = Transpose@Flatten[Table[With[{a = edges[[k]], c = edges[[k + 1]]},
        Transpose[{(a + c)/2 + (c - a)/2 x, (c - a)/2 ws}]], {k, Length[edges] - 1}], 1];
   th = 2. Pi (Range[nth] - 0.5)/nth;
   <|"r" -> rr, "wr" -> rw rr omega[rr] (2. Pi/nth), "c" -> Cos[th], "s" -> Sin[th],
     "c2" -> Cos[2 th], "s2" -> Sin[2 th]|>];
nuRule[q1_, q2_, rule_, {x0_, y0_}] := Module[{X, Y, n = Length[rule["c"]]},
   X = Flatten[x0 + Outer[Times, rule["r"], rule["c"]]];
   Y = Flatten[y0 + Outer[Times, rule["r"], rule["s"]]];
   rule["wr"] . (Partition[q1[X, Y], n] . rule["c2"] + Partition[q2[X, Y], n] . rule["s2"])];
ruleIn = makeRule[0., rCut]; ruleOut = makeRule[rCut, 400.];
nu[x0_?NumericQ] := nuRule[q1ss, q2ss, ruleIn, {x0, 0.}] + nuRule[q1Pade, q2Pade, ruleOut, {x0, 0.}];

(* the points: the extremes of the dotted curve along the axis *)
{nuHead, xHead} = Module[{s = FindMaximum[nu[x], {x, -1.3, -4., 0.}]}, {s[[1]], x /. s[[2]]}];
{nuTail, xTail} = Module[{s = FindMinimum[nu[x], {x, 1.7, 0., 4.}]}, {s[[1]], x /. s[[2]]}];
(* the core: the whole average, and its parts from the neighbors on the head side (cos t < 0)
   and on the tail side (cos t > 0); on the static profile they cancel exactly *)
nuCore = nu[0.];
halves[q1_, q2_, rule_] := Module[{n = Length[rule["c"]], X, Y, terms, head},
   X = Flatten[Outer[Times, rule["r"], rule["c"]]]; Y = Flatten[Outer[Times, rule["r"], rule["s"]]];
   terms = rule["wr"] . (Partition[q1[X, Y], n] (ConstantArray[rule["c2"], Length[rule["r"]]]) +
       Partition[q2[X, Y], n] (ConstantArray[rule["s2"], Length[rule["r"]]]));
   head = UnitStep[-rule["c"]];
   {terms . head, terms . (1 - head)}];
{hIn, tIn} = halves[q1ss, q2ss, ruleIn]; {hOut, tOut} = halves[q1Pade, q2Pade, ruleOut];
{hSt, tSt} = halves[q1Pade, q2Pade, ruleIn] + halves[q1Pade, q2Pade, ruleOut];
Print["  core: ", 100 nuCore, "% = head-side neighbors ", 100 (hIn + hOut), "% + tail-side neighbors ",
  100 (tIn + tOut), "%;  static profile: ", 100 hSt, "% + ", 100 tSt, "% = ", 100 (hSt + tSt), "%"];
Print["  head side: max ", 100 nuHead, "% at x = ", xHead, ";  tail side: min ", 100 nuTail,
  "% at x = ", xTail, "   (", Round[AbsoluteTime[] - t0], " s)"];

(* ---------------------------------------------------------------- graphics *)
<< MaTeX`
(* font sizes in design points; printed at 0.8\textwidth they shrink by 0.8 *)
fsLabel = If[four, 9, 11];  (* titles with the panel letters, axis labels, legend: 8.8 pt printed *)
fsTick = If[four, 8, 10];   (* tick labels and "core": 8 pt printed *)
tex[s_, fs_ : fsTick] := MaTeX[s, FontSize -> fs, "Preamble" -> {"\\usepackage{amsmath,upgreek}"}];
cmap = ColorData["TemperatureMap"];
{cAim, cSide} = {cmap[0.97], cmap[0.03]};       (* red raises, blue lowers the variance *)

L = 3.2;                    (* half-width of each window about its spot, um *)
step = 0.34; lRod = 0.24;   (* rod grid spacing and full rod length at S = S0, um *)
rWeight = 2.;               (* the dashed circle: 94% of W lies inside *)

(* layout, in design points; the figure is totalW wide *)
totalW = 510;
{titleY, padTop, labelOff, padBot, padLeft} = If[four, {8, 17, 20, 27, 28}, {10, 20, 24, 32, 33}];
(* labelOff: distance of the axis labels' centres from the frame *)
padA = {{padLeft, 3}, {padBot, padTop}}; padB = {{If[four, 7, 8], 3}, {padBot, padTop}};   (* {{left, right}, {bottom, top}} *)
nPanels = If[four, 4, 3];
s = (totalW - Total[padA[[1]]] - (nPanels - 1) Total[padB[[1]]])/nPanels;   (* side of each frame *)
{rodW, casingW} = If[four, {1.1, 2.2}, {1.3, 2.6}];
Print["  frame side ", s, " pt"];

(* each rod sits on a white casing slightly wider than itself, so the streamlines pass behind
   the rods even where a faint rod is mostly transparent *)
rods[x0_, uniform_] := Module[{g = Range[-L + step/2, L, step], out = {}, casing = {}},
   Do[Module[{w = Norm[{wx, wy}], X = x0 + wx, Y = wy, th, sv, c, op, d, seg},
      If[w > 0.12,
       {th, sv} = If[uniform, {N[Pi/2], s0}, With[{q = {q1ss[X, Y], q2ss[X, Y]}},
          If[Norm[q] < 10^-9, {0., 0.}, {ArcTan @@ q/2, Norm[q]}]]];
       c = Cos[2 (th - ArcTan[wx, wy])];
       op = If[w <= rWeight, 0.22 + 0.78 Abs[c], 0.10 + 0.22 Abs[c]];
       d = 0.5 lRod sv/s0 {Cos[th], Sin[th]};
       seg = Line[{{X, Y} - d, {X, Y} + d}];
       AppendTo[casing, seg];
       AppendTo[out, {If[c > 0, cAim, cSide], Opacity[op], seg}]]],
    {wx, g}, {wy, g}];
   {CapForm["Round"], {White, AbsoluteThickness[casingW], casing}, AbsoluteThickness[rodW], out}];

spot[x0_] := {{White, Disk[{x0, 0}, 0.17]}, {Black, AbsoluteThickness[1.1], Circle[{x0, 0}, 0.13]}};
(* label on a small translucent white plate with rounded corners (it sits on the
   streamlines and rods), drawn before the marker so the plate does not cover it; plate
   sizes in um, fitted to the 10 pt "core" (scaled from the 8 pt fit 2026-10-23). The label
   is shifted by dx toward the point, since at 10 pt it would otherwise touch the dashed
   circle, which passes 0.45 um beside the core at the label's height *)
coreY = 0.45;
(* shift toward the point so that the plate (half-width 0.47) clears the circle by 0.06 *)
coreDx[x0_] := Sign[x0] Max[0.1, Abs[x0] - Sqrt[rWeight^2 - coreY^2] + 0.53];
core[dx_] := {{White, Opacity[0.75], Rectangle[{dx - 0.47, coreY - 0.15}, {dx + 0.47, coreY + 0.15},
      RoundingRadius -> 0.08]},
   Inset[tex["\\text{core}"], {dx, coreY}],
   Black, Polygon[{{0.11, 0}, {0, 0.11}, {-0.11, 0}, {0, -0.11}}]};
weightCircle[x0_] := {GrayLevel[0.45], AbsoluteThickness[0.6], AbsoluteDashing[{3, 2}],
   Circle[{x0, 0}, rWeight]};

(* gray director streamlines behind the rods. Uniform: lines x = const (director along y). Defect: director lines
   of Q_ss, integrated by RK4 (step 0.01 um) from (-c^2, 0) upward and downward, the sign of
   the director chosen at every stage to continue the previous direction, until they leave
   [-7, 7] x [-4.5, 4.5]; plus the ray y = 0, x > 0. Cut to the window by hand, since
   PlotRangeClipping would also clip the axis labels. *)
streamC = Range[0.3, 2.4, 0.3];  streamDy = 0.68;
streamStyle = {GrayLevel[0.72], AbsoluteThickness[0.5]};   (* light: background only *)
dirQ[{x_, y_}, prev_] := Module[{q1 = q1ss[x, y], q2 = q2ss[x, y], th, n},
   th = If[q1 == 0. && q2 == 0., 0., ArcTan[q1, q2]/2]; n = {Cos[th], Sin[th]};
   If[n . prev < 0, -n, n]];
traceLine[p0_, d0_, h_ : 0.01, nMax_ : 3000] := Module[{p = p0, d = d0, k1, k2, k3, k4, out = {p0}, i = 0},
   While[i++ < nMax && Abs[p[[1]]] < 7 && Abs[p[[2]]] < 4.5,
    k1 = dirQ[p, d]; k2 = dirQ[p + h/2 k1, k1]; k3 = dirQ[p + h/2 k2, k2]; k4 = dirQ[p + h k3, k3];
    d = (k1 + 2 k2 + 2 k3 + k4)/6; p = p + h d; AppendTo[out, p]];
   out];
streamPaths = Join @@ Table[{traceLine[{-c^2, 0.}, {0., 1.}], traceLine[{-c^2, 0.}, {0., -1.}]}, {c, streamC}];
Print["  director lines traced  (", Round[AbsoluteTime[] - t0], " s)"];
streams[x0_, uniform_] := Module[{inW, cut},
   inW[{x_, y_}] := Abs[x - x0] <= L && Abs[y] <= L;
   cut[path_] := Line[Select[SplitBy[path, inW], inW@First@# && Length[#] > 1 &]];
   (* the style must sit in the same list as the lines: directives do not reach out of a
      nested list *)
   Join[streamStyle, If[uniform,
     Line[{{x0 + #, -L}, {x0 + #, L}}] & /@ N@Select[streamDy Range[-6, 6], Abs[#] <= L &],
     {cut /@ streamPaths, Line[{{Max[0., x0 - L], 0}, {x0 + L, 0}}]}]]];

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
     If[kind === "head" || kind === "tail", core[coreDx[x0]], {}], spot[x0], title[label, titleY, s],
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
If[four, titles = {titles[[1]], titles[[2]], "\\textbf{(c)}\\enspace\\text{core: } " <> pct[nuCore],
    StringReplace[titles[[3]], "(c)" -> "(d)"]}];
panels = {panel[0., "uniform", titles[[1]], True, True, padA],
   panel[xHead, "head", titles[[2]], False, True, padB],
   panel[xTail, "tail", titles[[3]], False, True, padB]};
(* (d): the parts of the core value from the neighbors on the head side (x < 0, blue: they
   lie across the line to the core) and on the tail side (x > 0, red: they point at it), on
   plates in the upper corners, and a dotted line x = 0 between the halves *)
halfLabel[v_, col_, {x_, y_}, al_] := {{White, Opacity[0.8], Rectangle[{x - If[al === Left, 0.08, 1.62], y - 0.3},
      {x + If[al === Left, 1.62, 0.08], y + 0.3}, RoundingRadius -> 0.08]},
   Inset[MaTeX["\\color[rgb]{" <> StringRiffle[ToString[NumberForm[#, {4, 3}]] & /@ List @@ ColorConvert[col, "RGB"], ","] <>
      "}" <> pct[v], FontSize -> fsTick, "Preamble" -> {"\\usepackage{amsmath,upgreek,xcolor}"}], {x, y}, {al, Center}]};
halfLabels = If[four, {{GrayLevel[0.45], AbsoluteThickness[0.6], AbsoluteDashing[{1, 2}],
     Line[{{0, -L}, {0, -0.25}}], Line[{{0, 0.25}, {0, L}}]},
    halfLabel[hIn + hOut, cSide, {-L + 0.15, L - 0.42}, Left],
    halfLabel[tIn + tOut, cAim, {L - 0.15, L - 0.42}, Right]}, {}];
If[four, panels = Insert[panels, Show[panel[0., "core", titles[[3]], False, True, padB],
    Graphics[halfLabels]], 3]; panels[[4]] = panel[xTail, "tail", titles[[4]], False, True, padB]];

sample[g_] := Graphics[g, ImageSize -> {22, 10}, PlotRange -> {{-1, 1}, {-1, 1}},
   ImagePadding -> 0, AspectRatio -> Full];
legendItems = {
   {sample[{CapForm["Round"], AbsoluteThickness[rodW], cAim, Line[{{-0.8, 0}, {0.8, 0}}]}],
    "\\text{within }45^\\circ\\text{ of the line to the point: raises the variance}"},
   {sample[{CapForm["Round"], AbsoluteThickness[rodW], cSide, Line[{{-0.8, 0}, {0.8, 0}}]}],
    "\\text{beyond }45^\\circ\\text{: lowers it}"},
   {sample[{GrayLevel[0.45], AbsoluteThickness[0.6], AbsoluteDashing[{3, 2}],
      Line[{{-1, 0}, {1, 0}}]}],
    "2~\\upmu\\mathrm{m}\\text{, 94\\% of the weight}"}};
(* two rows (2026-10-23): at 11 pt the three items no longer fit in one row of 510; the
   colour rule (red, blue) on the first, the dashed circle under the red item *)
legend = With[{it = {#1, tex[#2, fsLabel]} & @@@ legendItems},
   Grid[{Join[it[[1]], {Spacer[18]}, it[[2]]], Join[it[[3]], {"", "", ""}]},
    Alignment -> {Left, Center}, Spacings -> {0.4, 0.15}]];

fig = Grid[{panels, Join[{legend}, ConstantArray[SpanFromLeft, nPanels - 1]]},
   Spacings -> {0, {0, 0.3}}, Alignment -> {Center, Top}];

outName = If[four, "aimed_rods_qss", "aimed_rods_qss_3panel"];
Export[FileNameJoin[{outDir, outName <> ".pdf"}], fig, "AllowRasterization" -> False];
Export[FileNameJoin[{outDir, outName <> ".png"}], fig, ImageResolution -> 600];
Print["exported to ", outDir, "  (", Round[AbsoluteTime[] - t0], " s)"];
