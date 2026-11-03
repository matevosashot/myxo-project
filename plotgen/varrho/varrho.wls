#!/usr/bin/env wolframscript
(* varrho.wls -- fig:variance (manuscript/figures/varrho.pdf): the equal-point density
   variance over its pedestal, Var[rho]/Var[rho]_pedestal, along the defect axis
   ((x, y) = (xi, 0), blue) and across it ((x, y) = (0, xi), orange). Written 2026-10-23 for
   the prompt above fig:variance.

   It reproduces the figure of the notebook negative-defect/plotting_small_core+ln_2.nb
   (Export cell "varrho", near line 2732) -- same data, curves, colours, line styles, axis
   ranges, ticks and width -- and adds a third pair of curves:

     dashed  eq:total-gauss at delta Q = 0 on the static profile of H^a = 0, eq:pade-accurate
             (rhoVarianceTheory of the notebook: the H_1 convolution, out to r' = 400 um);
     dotted  NEW: the same formula, delta Q = 0, but with Q_ss the numerical steady state of
             the moving defect, the one the Lyapunov run linearises about (steady_<tag>.m,
             written by jobscripts/set2_2_comove_optimized_2/fluctuations.wls on the solver's
             101 x 101 Pade-graded grid);
     solid   the numerical Lyapunov solution, delta Q != 0 (sigma_rho_<tag>.m, the 180 x 180
             diagonal of the covariance .bin, first index y, second x; it equals the .bin
             diagonal the notebook read).

   The dotted curve uses the neighbour form of eq:total-gauss (app:variance-neighbors),
     Var/Var_ped - 1 = Int d^2w Omega(w) w^.Q_ss(r + w).w^,
     Omega(w) = (lN^2/(2 Pi w^4)) [2 - e^-X (X^2 + 2X + 2)],  X = 2 w^2/lN^2,
   so the grid data are interpolated, never differentiated: w^.Q.w^ = Q1 cos2t + Q2 sin2t,
   t the angle of w. Polar quadrature about r: Gauss-Legendre on radial segments, uniform
   angles. Inside the disc |w| < Rc (default 12.5 um, inside the box for |xi| <= 7) Q_ss is
   the numerical field (InterpolationOrder 3); outside it, out to 400 um, the static
   profile, as for the dashed curve, so the two differ only through Q_ss within Rc.
   "check" prints the validation: the neighbour form on the static profile against the
   dashed curve, the quadrature, the interpolation order, the grid, and the cutoff Rc.
   Results of 2026-10-23 (default run): the static profile reproduces the dashed curve to
   1e-6 (+-2.9423% at x = -+1.5); the dotted curve moves by < 4e-5 with the interpolation
   order (2 vs 3), < 5e-6 with the grid (101 graded vs 180 uniform), and by <= 2.5e-4 with
   the treatment of the far field (no tail, whole box only, Rc = 9). Peaks: dotted +2.546%
   at x = -1.26 and -3.394% at +1.66; solid +2.398% at -1.26 and -3.261% at +1.66.

   Drawing. The four original curves keep the notebook's styles (Thickness 0.006 Darker solid,
   0.008 Dashed, ColorData[97] 1 and 2). The new curves are round dots in the base colours,
   drawn on top as Disks at equal arc length (zero-length dashes with round caps do not
   survive the headless pdf export). PlotRange -> All now also covers the dotted trough, so
   the y range reaches 0.966 instead of 0.967. Two legends replace the notebook's four-row
   one: line style = method (grey, lower left), colour = direction (upper right), in MaTeX at
   9 pt (prints at about 8.5-9 pt at \linewidth; the old legend was 7 pt).
   2026-10-23, his comments: no y minor ticks, dashed Thickness 0.005 (was 0.008), method
   legend in the order dashed "analytic, Q_Pade", dotted "analytic, Q_ss", solid "numerical".
   All text is MaTeX in Computer Modern, as in the exported varrho.pdf (whose glyphs are CM,
   although the notebook's tex[] now loads newtx); the headless front end garbles plain text
   in pdf exports, so the tick labels, plain text in the notebook, are MaTeX CM too, trimmed
   so they sit where the old ones sat. Design width 8.6 cm; the pdf export scales by 0.72,
   as the notebook's did, so the pdf is 176 x 118 pt like the original, with the frame at the
   same place (ImagePadding measured on the original).

   2026-11-01 (the author's request: the axis labels sit too far from the axes): FrameLabel replaced by
   frame_labels.wl, which puts the labels about 4 pt from the tick labels; the padding shrinks
   to fit, so the frame grows at the same width and the page gets lower.
   2026-11-02 (the author's rule: no line thinner than about 1 px; a line meant to look
   thinner is gray): frame and ticks print at 0.75 pt (0.20 pt before), the grid
   lines gray at 0.75 pt (frame_labels.wl); the legend swatches take Thickness relative to the
   frame width, as the curves do (relative to the figure width before, so thicker).
   2026-11-02, later requests of the author: grid lines 0.6 pt in a darker gray (gridStyle of
   frame_labels.wl); the method legend moved right, from 2% to 5% of the frame width.

   Usage:  plotgen/varrho/run.sh [check] [ssDir] [key=value ...]
     ssDir     take the Lyapunov curve AND the steady state from the ssDir run (Dirichlet
               density wall in the steady state, as the manuscript states); default: the run
               of the current figure, whose steady state has the RhoNeumann density wall
     nonew     omit the new curves and draw the original legend layout (reproduction check)
     Rc=<um>   radius of the numerical part of the neighbour integral (default 12.5)
     out=<dir> output directory (default: next to this script)
   Outputs: varrho.pdf, varrho.png (600 dpi); with nonew, varrho_orig.pdf/.png. *)

SetEnvironment["QT_QPA_PLATFORM" -> "offscreen"];

here = DirectoryName[$InputFileName];
repoRoot = ParentDirectory[ParentDirectory[here]];
args = Rest[$ScriptCommandLine];
check = MemberQ[args, "check"];
ssDir = MemberQ[args, "ssDir"];
noNew = MemberQ[args, "nonew"];
opts = Association[Rule @@@ (StringSplit[#, "=", 2] & /@ Select[args, StringContainsQ[#, "="] &])];
num[key_, default_] := If[KeyExistsQ[opts, key], N@ToExpression@opts[key], default];
rCut = num["Rc", 12.5];
outDir = Lookup[opts, "out", here];
t0 = AbsoluteTime[];

(* ---------------------------------------------------------------- data *)
tag = "box20_N200_Dir_lN0.5_z12.5_ld1_B90000_ss20";
dataDir = FileNameJoin[{repoRoot, "data",
    "fluctuations-set2_2_comove_optimized_2" <> If[ssDir, "_ssDir", ""]}];
sig = Import[FileNameJoin[{dataDir, "sigma_rho_" <> tag <> ".m"}]];
st = Import[FileNameJoin[{dataDir, "steady_" <> tag <> ".m"}]];
modelParams = sig["modelParams"];
If[Sort[N@modelParams] =!= Sort[N@st["modelParams"]], Print["*** WARNING: the steady state's parameters differ"]];
{lN, Bv, zetav, rho0v, a1v, L1v, bv} = {lNoise, B, \[Zeta], \[Rho]0, a1, L1, b} /. modelParams;
pedestal = zetav/(2 Pi Bv rho0v lN^2);                 (* eq:pedestal *)
Print["run ", FileNameTake[dataDir], "/", tag, ":  steady density wall ", st["bc"],
  ", Q wall ", st["qWall"], ", u = ", st["u"][[1]], " um/min; steady grid ", st["nx"], " x ",
  st["ny"], " (graded, h ", st["hMin"], "..", st["hMax"], "), box ", st["box"]];

(* solid: the Lyapunov variance, as readCovRho builds it (x first, InterpolationOrder 2) *)
gridInt = N@sig["gridInt"];
lyap = ListInterpolation[Transpose[sig["sigmaRho"]]/pedestal, {gridInt, gridInt},
   InterpolationOrder -> 2];

(* dashed: rhoVarianceTheory of the notebook (getRhoVariance, rInterpolation 40) *)
fPade[z_] := Sqrt[(0.34010999 + 0.04516965 z + 0.00233656 z^2)/(1 + 0.38387160 z +
      0.04750621 z^2 + 0.00233656 z^3)];
SoR[r2_] := Sqrt[-a1v/(2 bv)] Sqrt[-a1v/L1v] fPade[-r2 a1v/L1v];   (* S(r)/r, eq:pade-accurate *)
Sfn[r_] = r SoR[r^2];
gFn[r_] = D[D[Sfn[r], r] + Sfn[r]/r, r] // Simplify;           (* d_i d_j Q_ij = gFn(r) cos(phi) *)
phiVec[t_] := Module[{s = UnitStep[1.*^-6 - t], td}, td = t + s;
   (1 - s) ((1 - Exp[-td])/td) + s (1 - t/2 + t^2/6)];
filterH1[rs_List, nr_ : 4000, nA_ : 720, rmax_ : 400.] := Module[{rp, rp2, gp, ca, wgt, dw, il2},
   rp = Exp[Subdivide[N@Log[1.*^-4], N@Log[rmax], nr - 1]]; rp2 = rp^2; gp = gFn /@ rp;
   ca = Cos[(Range[0, nA - 1] + 0.5) (2. Pi/nA)]; wgt = ca (2. Pi/nA);
   dw = Differences[rp]; il2 = 1./(Pi lN^2);
   Table[Module[{w2, k1, integ},
     w2 = (r^2 + rp2) - 2. r Outer[Times, rp, ca];
     k1 = il2 (phiVec[2. w2/lN^2] . wgt); integ = rp gp k1;
     Total[(Most[integ] + Rest[integ])/2 dw]], {r, rs}]];
hRadii = Join[{0.}, Exp[Subdivide[N@Log[1.*^-2], N@Log[1.2 40], 89]]];
H1 = Interpolation[Transpose[{hRadii, Join[{0.}, filterH1[Rest[hRadii]]]}], InterpolationOrder -> 3];
theory[x_?NumericQ, y_?NumericQ] := With[{r = Sqrt[x^2 + y^2]},
   If[r == 0., 1., 1 + lN^2/4 H1[r] x/r]];

(* ------------------------------------------- dotted: the neighbour form on Q_ss *)
(* f(X) = 2 - e^-X (X^2 + 2X + 2) = Int_0^X t^2 e^-t dt; its series below X = 0.01 *)
fX[X_] := Module[{s = UnitStep[0.01 - X], Xc = Clip[X, {0., 600.}]},   (* Clip: no underflow *)
   (1 - s) (2 - Exp[-Xc] (Xc^2 + 2 Xc + 2)) + s (X^3/3 - X^4/4 + X^5/10 - X^6/36)];
omega[w_] := lN^2/(2 Pi w^4) fX[2 w^2/lN^2];
glNodes[n_] := glNodes[n] = Module[{x, ws},
   x = Sort[x /. NSolve[LegendreP[n, x] == 0, x, WorkingPrecision -> 30]];
   ws = Map[Function[x0, 2/((1 - x0^2) (D[LegendreP[n, xx], xx] /. xx -> x0)^2)], x];
   N[{x, ws}]];
breaks = {0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 10, 12.5, 15, 20, 30, 50, 100, 200, 400};
(* polar rule on the annulus r1 < w < r2: radial GL nodes (nr per segment) with the weight
   w Omega(w) dw dt folded in, nth uniform angles *)
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
q1Pade[X_, Y_] := SoR[X^2 + Y^2] X;     (* the static profile, eq:pade-accurate *)
q2Pade[X_, Y_] := SoR[X^2 + Y^2] Y;

{gx, gy} = N /@ {st["gx"], st["gy"]};
q1ss = ListInterpolation[st["Q1Grid"], {gx, gy}, InterpolationOrder -> 3];   (* [[i,j]] = f(gx_i, gy_j) *)
q2ss = ListInterpolation[st["Q2Grid"], {gx, gy}, InterpolationOrder -> 3];
ruleIn = makeRule[0., rCut]; ruleOut = makeRule[rCut, 400.];
nuSS[q1_, q2_, ruleI_, ruleO_][p_] := nuRule[q1, q2, ruleI, p] + nuRule[q1Pade, q2Pade, ruleO, p];
nuSteady = nuSS[q1ss, q2ss, ruleIn, ruleOut];

(* sampled every 0.1 um and interpolated (cubic); "check" tests the midpoints *)
xiS = Range[-7.5, 7.5, 0.1];
If[! noNew || check,
  ssAlongI = Interpolation[Transpose[{xiS, 1 + (nuSteady[{#, 0.}] & /@ xiS)}], InterpolationOrder -> 3];
  ssAcrossI = Interpolation[Transpose[{xiS, 1 + (nuSteady[{0., #}] & /@ xiS)}], InterpolationOrder -> 3];
  Print["neighbour form on Q_ss: ", 2 Length[xiS], " points, Rc = ", rCut, " um  (",
    Round[AbsoluteTime[] - t0], " s)"]];

(* ---------------------------------------------------------------- numbers *)
pct[v_] := Round[100 (v - 1), 0.001];
ext[f_, x0_, lo_, hi_, mx_] := Module[{g, s, x},
   g[x_?NumericQ] := f[x];          (* keeps FindMaximum from evaluating f symbolically *)
   s = If[mx, FindMaximum, FindMinimum][g[x], {x, x0, lo, hi}];
   {s[[1]], x /. s[[2]]}];
report[name_, along_, across_] := Module[{mx, mn, my},
   mx = ext[along, -1.3, -4., 0., True]; mn = ext[along, 1.7, 0., 4., False];
   my = ext[across, 2.7, 1., 5., False];
   Print[StringPadRight[name, 38], " max ", pct[mx[[1]]], "% at x = ", Round[mx[[2]], 0.01],
     ";  min ", pct[mn[[1]]], "% at x = ", Round[mn[[2]], 0.01],
     ";  core ", pct[along[0.]], "%;  across: min ", pct[my[[1]]], "% at |y| = ",
     Round[my[[2]], 0.01], ", y-mirror ", across[-2.] - across[2.]]];
Quiet@report["dashed  (static Q, delta Q = 0)", theory[#, 0.] &, theory[0., #] &];
If[! noNew || check,
  Quiet@report["dotted  (numerical Q_ss, delta Q = 0)", nuSteady[{#, 0.}] + 1 &, nuSteady[{0., #}] + 1 &]];
Quiet@report["solid   (Lyapunov, delta Q != 0)", lyap[#, 0.] &, lyap[0., #] &];
If[! noNew || check, Module[{xiC = Range[-7., 7., 0.05]},
  Print["max |dotted - solid|: along ", 100 Max[Abs[ssAlongI[xiC] - (lyap[#, 0.] & /@ xiC)]],
    "%, across ", 100 Max[Abs[ssAcrossI[xiC] - (lyap[0., #] & /@ xiC)]], "%;  max |dotted - dashed|: along ",
    100 Max[Abs[ssAlongI[xiC] - (theory[#, 0.] & /@ xiC)]], "%, across ",
    100 Max[Abs[ssAcrossI[xiC] - (theory[0., #] & /@ xiC)]], "%"]]];

(* ---------------------------------------------------------------- validation *)
If[check, Module[{xs = Range[-7., 7., 0.25], pts, ref, nP, nPin, d, q1g, q2g, res, st2},
   pts = {#, 0.} & /@ xs;
   Print["--- check (along the axis, xi = -7..7 step 0.25; differences in Var/Var_ped) ---"];
   Print["dotted curve, cubic interpolant at the midpoints vs direct: max |diff| along ",
     Max[Abs[ssAlongI[# + 0.05] - 1 - nuSteady[{# + 0.05, 0.}] & /@ Range[-7.05, 6.95, 0.1]]], ", across ",
     Max[Abs[ssAcrossI[# + 0.05] - 1 - nuSteady[{0., # + 0.05}] & /@ Range[-7.05, 6.95, 0.1]]]];
   (* 1. the method on the static profile against the dashed curve *)
   ref = theory @@@ pts - 1;
   nP = nuSS[q1Pade, q2Pade, ruleIn, ruleOut] /@ pts;
   Print["static profile, neighbour form (to 400 um) vs H_1 route (dashed): max |diff| ",
     Max[Abs[nP - ref]]];
   Print["  at x = -1.5, +1.5: neighbour ", 100 nP[[{23, 35}]], " %,  dashed ", 100 ref[[{23, 35}]], " %"];
   Print["  across the axis: max |nu| ", Max[Abs[nuSS[q1Pade, q2Pade, ruleIn, ruleOut][{0., #}] & /@ xs]]];
   nPin = nuRule[q1Pade, q2Pade, ruleIn, #] & /@ pts;
   Print["  without the tail beyond Rc = ", rCut, ": max |diff| ", Max[Abs[nPin - ref]],
     " (at |x| = 7: ", Abs[nPin[[1]] - ref[[1]]], ", at |x| = 1.5: ", Abs[nPin[[23]] - ref[[23]]], ")"];
   Print["  quadrature 24 x 256 vs 48 x 512: max |diff| ",
     Max[Abs[nP - (nuSS[q1Pade, q2Pade, makeRule[0., rCut, 48, 512], makeRule[rCut, 400., 48, 512]] /@ pts)]]];
   (* 2. the static profile sampled on the steady grid: the interpolation error *)
   Do[q1g = ListInterpolation[Outer[q1Pade, gx, gy], {gx, gy}, InterpolationOrder -> o];
     q2g = ListInterpolation[Outer[q2Pade, gx, gy], {gx, gy}, InterpolationOrder -> o];
     Print["static profile on the 101 x 101 steady grid, InterpolationOrder ", o, ": max |diff| ",
       Max[Abs[nP - (nuSS[q1g, q2g, ruleIn, ruleOut] /@ pts)]]], {o, {1, 2, 3}}];
   (* 3. the numerical Q_ss: interpolation order, grid, cutoff *)
   res = nuSteady /@ pts;
   Do[q1g = ListInterpolation[st["Q1Grid"], {gx, gy}, InterpolationOrder -> o];
     q2g = ListInterpolation[st["Q2Grid"], {gx, gy}, InterpolationOrder -> o];
     Print["Q_ss, InterpolationOrder ", o, " vs 3: max |diff| ",
       Max[Abs[res - (nuSS[q1g, q2g, ruleIn, ruleOut] /@ pts)]]], {o, {1, 2}}];
   q1g = ListInterpolation[st["Q1Int"], {gridInt, gridInt}, InterpolationOrder -> 3];
   q2g = ListInterpolation[st["Q2Int"], {gridInt, gridInt}, InterpolationOrder -> 3];
   Print["Q_ss resampled on the uniform 180 x 180 Lyapunov grid vs the steady grid: max |diff| ",
     Max[Abs[res - (nuSS[q1g, q2g, ruleIn, ruleOut] /@ pts)]]];
   Do[d = nuSS[q1ss, q2ss, makeRule[0., rc], makeRule[rc, 400.]] /@ pts;
     Print["Q_ss, Rc = ", rc, " vs ", rCut, ": max |diff| ", Max[Abs[d - res]],
       " (|x| <= 2: ", Max[Abs[(d - res)[[21 ;; 37]]]], ")"], {rc, {4., 6., 9.}}];
   d = nuRule[q1ss, q2ss, ruleIn, #] & /@ pts;
   Print["Q_ss without the static tail beyond Rc = ", rCut, ": max |diff| ", Max[Abs[d - res]],
     " (|x| <= 2: ", Max[Abs[(d - res)[[21 ;; 37]]]], ")"];
   (* the whole box, nothing outside: a masked disc reaching the farthest corner *)
   d = Module[{rule = makeRule[0., 38., 24, 512], inBox},
      inBox[X_, Y_] := UnitStep[19.999 - Abs[X]] UnitStep[19.999 - Abs[Y]];
      nuRule[Function[{X, Y}, inBox[X, Y] q1ss[Clip[X, {-20., 20.}], Clip[Y, {-20., 20.}]]],
        Function[{X, Y}, inBox[X, Y] q2ss[Clip[X, {-20., 20.}], Clip[Y, {-20., 20.}]]], rule, #] & /@ pts];
   Print["Q_ss over the whole box, nothing outside: max |diff| ", Max[Abs[d - res]],
     " (|x| <= 2: ", Max[Abs[(d - res)[[21 ;; 37]]]], ")"];
   st2 = Import[FileNameJoin[{repoRoot, "data", "fluctuations-set2_2_comove_optimized_2" <>
        If[ssDir, "", "_ssDir"], "steady_" <> tag <> ".m"}]];
   q1g = ListInterpolation[st2["Q1Grid"], {gx, gy}, InterpolationOrder -> 3];
   q2g = ListInterpolation[st2["Q2Grid"], {gx, gy}, InterpolationOrder -> 3];
   Print["Q_ss of the other steady density wall (", st2["bc"], ", u = ", st2["u"][[1]],
     "): max |diff| ", Max[Abs[res - (nuSS[q1g, q2g, ruleIn, ruleOut] /@ pts)]]];
   Print["--- check done (", Round[AbsoluteTime[] - t0], " s)"]]];

(* ---------------------------------------------------------------- figure *)
<< MaTeX`
Get[FileNameJoin[{repoRoot, "plotgen", "frame_labels.wl"}]];
w = 8.6*72/2.54;                    (* design width, as in the notebook *)
tex[s_, fs_ : 11] := MaTeX[s, FontSize -> fs,
   "Preamble" -> {"\\usepackage{amsmath,siunitx,upgreek}"}];
colors = ColorData[97, "ColorList"];
(* 2026-11-01, the author's request for one font and one set of sizes in all figures: tick
   and axis labels print at 9 pt, legends at 8 pt (at \linewidth = 246 pt; until then ticks 8,
   axis labels 11, legends 9) *)
{fsT, legFS} = {9, 8} w/246.;

(* ticks: the frame ticks the notebook asked for (x: Range[-20, 20, 2]; y: Automatic), with
   MaTeX labels at the 8 pt of FrameStyle; the automatic y ticks are read off a first
   rendering with AbsoluteOptions *)
(* MaTeX boxes carry about 1 pt of margin on each side, which pushes the labels away from the
   frame; trimmed by 1 pt they sit where the notebook's plain-text tick labels sat *)
trim[g_, d_] := Module[{pr = PlotRange /. Options[g], is = First@Cases[Options[g], (ImageSize -> v_) :> v]},
   Graphics[First[g], PlotRange -> {pr[[1]] + {d, -d}, pr[[2]]}, ImageSize -> is - {2 d, 0},
    AspectRatio -> Automatic, BaselinePosition -> (BaselinePosition /. Options[g])]];
xLabel[v_] := trim[tex[ToString[Round[v]], fsT], 1.];
yLabel[v_] := trim[tex[ToString[NumberForm[v, {3, 2}]], fsT], 1.];

styleSolid[c_] := Directive[Darker[c], Thickness[0.006]];        (* as in the notebook *)
styleDashed[c_] := Directive[Dashed, c, Thickness[0.005]];    (* 0.008 in the notebook; thinner, his request 2026-10-23 *)
(* the new curves: round dots in the base colour, drawn on top as Disks at equal arc length
   (zero-length dashes with round caps do not survive the headless pdf export) *)
dotD = 2.4; dotStep = 4.8;          (* dot diameter and spacing, design pt *)

funcs = {lyap[\[Xi], 0], theory[\[Xi], 0], lyap[0, \[Xi]], theory[0, \[Xi]]};
styles = {styleSolid[colors[[1]]], styleDashed[colors[[1]]], styleSolid[colors[[2]]],
   styleDashed[colors[[2]]]};
funcsAll = If[noNew, funcs, Join[funcs, {ssAlongI[\[Xi]], ssAcrossI[\[Xi]]}]];

(* until 2026-11-01 the frame sat where the original's automatic padding put it,
   {{40.5, 1.78}, {39.1, 0.25}}; now the padding only leaves room for the labels of
   frameLabels (4 pt from the tick labels) *)
pad = {{32.9, 1.78}, {24.1, 0.25}};
labelDist = {13.8, 21.1};           (* frame to the label boxes, x and y, design pt *)
plotW = w - Total[pad[[1]]]; plotH = plotW/GoldenRatio;
plotOpts = {GridLines -> {{0}, {1}}, GridLinesStyle -> gridStyle[w/246.], Axes -> False, Frame -> True,
   ImagePadding -> pad, ImageSize -> {w, plotH + Total[pad[[2]]]}};
(* a first rendering with PlotRange -> All over all curves gives the range (the new curves dip
   below the old ones at the trough) and the automatic y ticks, relabelled in MaTeX *)
prelim = Plot[Evaluate[funcsAll], {\[Xi], -7, 7}, PlotRange -> All, Evaluate[Sequence @@ plotOpts],
   FrameTicks -> {{Automatic, None}, {Range[-20, 20, 2], None}}];
pRange = PlotRange /. AbsoluteOptions[prelim, PlotRange];
pPad = PlotRangePadding /. AbsoluteOptions[prelim, PlotRangePadding];
yTicks = First@First[FrameTicks /. AbsoluteOptions[prelim, FrameTicks]];
yTicks = Cases[yTicks, {v_, lab_, rest___} /; !MatchQ[lab, "" | Null] :> {v, yLabel[v], rest}];
   (* labelled ticks only: no y minor ticks (his request 2026-10-23) *)
xTicks = {#, xLabel[#]} & /@ Range[-20, 20, 2];
Print["y range ", pRange[[2]], If[noNew, "", " (the four original curves alone: " <>
     ToString[MinMax[Table[funcs, {\[Xi], -7., 7., 0.005}]]] <> ")"]];

(* dots every dotStep pt of arc length, in the plot's own scale (pt per unit) *)
{sx, sy} = {plotW, plotH}/(Subtract @@@ Reverse /@ pRange + Total /@ pPad);
dotted[f_] := Module[{xs = Subdivide[-7., 7., 6000], ys, s, at},
   ys = f[xs];
   s = Prepend[Accumulate[Sqrt[(sx Differences[xs])^2 + (sy Differences[ys])^2]], 0.];
   at = Interpolation[Transpose[{s, xs}], InterpolationOrder -> 1] /@
     Range[Mod[Last[s], dotStep]/2, Last[s], dotStep];      (* centred: same margin at both ends *)
   Transpose[{at, f[at]}]];
(* Disks with radii in display units (round on the page); Points come out square in the png *)
dots[pts_] := Disk[#, (dotD/2) {1/sx, 1/sy}] & /@ pts;
dotPrims = If[noNew, {}, {colors[[1]], dots[dotted[ssAlongI]], colors[[2]], dots[dotted[ssAcrossI]]}];

(* legends: a swatch and a MaTeX label per row, as insets in the plot *)
swOpts[len_] := Sequence[PlotRange -> {{0, len}, {-4, 4}}, ImageSize -> {len, 8},
   AspectRatio -> Full, ImagePadding -> 0];
swatch[style_, len_ : 13] := Graphics[{style /. Thickness[t_] :> AbsoluteThickness[t plotW],
    Line[{{0, 0}, {len, 0}}]}, swOpts[len]];
swatchDots[c_, len_ : 13] := Module[{n = Floor[(len - dotD)/dotStep] + 1},   (* whole dots, centred *)
   Graphics[{c, Disk[{(len - (n - 1) dotStep)/2 + # dotStep, 0}, dotD/2] & /@ Range[0, n - 1]},
    swOpts[len]]];
legendBox[rows_] := Grid[{#[[1]], tex[#[[2]], legFS]} & /@ rows,
   Alignment -> {Left, Center}, Spacings -> {0.3, 0.15}];
grey = GrayLevel[0.5];               (* the method legend in grey, Darker for solid as in the plot *)
(* the longest label in the bottom row, where the curves lie furthest to the right *)
legendMethod = legendBox[{
    {swatch[styleDashed[grey]], "\\text{analytic, }\\mathbf Q_\\text{Pad\\'e}"},
    {swatchDots[grey], "\\text{analytic, }\\mathbf Q_\\mathrm{ss}"},
    {swatch[styleSolid[grey]], "\\text{numerical}"}}];   (* order and labels: his request 2026-10-23 *)
legendDir = legendBox[{
    {swatch[Directive[colors[[1]], Thickness[0.008]]], "\\text{along: }(x,\\!y)\\!=\\!(\\xi,\\!0)"},
    {swatch[Directive[colors[[2]], Thickness[0.008]]], "\\text{across: }(x,\\!y)\\!=\\!(0,\\!\\xi)"}}];
(* the original 4-row legend, for the reproduction check only *)
legendOrig = LineLegend[{styleDashed[colors[[1]]], styleSolid[colors[[1]]], styleDashed[colors[[2]]],
    styleSolid[colors[[2]]]},
   {tex["(x,\\!y)\\!=\\!(\\xi,\\!0)~~\\text{Theory}", 7], tex["~~~~~\\text{Lyapunov solver}", 7],
    tex["(x,\\!y)\\!=\\!(0,\\!\\xi)~~\\text{Theory}", 7], tex["~~~~~\\text{Lyapunov solver}", 7]},
   LegendMarkerSize -> 15, Spacings -> {0.4, -0.6}];

legendOpt = If[noNew, {PlotLegends -> Placed[legendOrig, {0.25, 0.22}]}, {}];
fig = Plot[Evaluate[funcs], {\[Xi], -7, 7},
   PlotStyle -> styles,
   Evaluate[Sequence @@ legendOpt],
   Epilog -> {If[noNew, {}, {dotPrims,
      Inset[legendMethod, Scaled[{0.05, 0.015}], {Left, Bottom}],   (* clear of the y ticks *)
      Inset[legendDir, Scaled[{0.995, 0.985}], {Right, Top}]}],
     (* the axis labels, close to the frame (frame_labels.wl) *)
     frameLabels[{tex["\\xi~(\\upmu \\rm{m})", fsT],
       tex["\\text{Var}[\\rho]~/~\\text{Var}[\\rho]_\\text{pedestal}", fsT]}, labelDist]},
   PlotRange -> pRange, PlotRangePadding -> pPad,   (* explicit range: the padding must be too *)
   Evaluate[Sequence @@ plotOpts],
   FrameTicks -> {{yTicks, None}, {xTicks, None}},
   FrameStyle -> Directive[frameStyle[w/246.], FontFamily -> "Times", FontSize -> 8],
   PlotRangeClipping -> False];     (* the labels lie outside the frame; the curves stay inside *)

name = If[noNew, "varrho_orig", "varrho"];
Export[FileNameJoin[{outDir, name <> ".pdf"}], fig];
Export[FileNameJoin[{outDir, name <> ".png"}], fig, ImageResolution -> 600];
Print["exported ", FileNameJoin[{outDir, name <> ".pdf"}], " and .png  (", Round[AbsoluteTime[] - t0], " s)"];
