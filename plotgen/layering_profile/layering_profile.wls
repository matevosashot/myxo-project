#!/usr/bin/env wolframscript
(* layering_profile.wls -- panel (c) of fig:defects-layering (manuscript/figures/final.pdf,
   placed by hand into the Illustrator composite figures/handmade/main.pdf, page 2): the
   layering probability averaged over a circle of radius r, eq:p-layer-radial, relative to
   its far-field value, about the point of maximum probability r_max (blue) and about the
   defect core (orange), with the measured points of Fig. 3f of Copenhagen et al. 2021 (red,
   error bars as given there). Written 2026-10-26 to regenerate the panel from the N_int = 200
   Lyapunov run.

   Source: remote_kernel/plotting_small_core+ln_2.nb (read only), the cells
     dataExperiment = ... top_plot.csv                                 (In[248], line ~3457)
     Qapprox[z] = Exp[-z^2/2]/(z Sqrt[2 Pi])                            (In[356], line ~3502)
     p = Compile[{x, y, pcrit}, 2 Qapprox[Sqrt[pcrit^2/(B^2 sigRhoFnDirichlet[x, y])]]],
     x00 = -x /. NMaximize[sigRhoFnDirichlet[-x, 0], {x, 0, 5}]         (In[360], line ~3590)
     pRad[r, x0, y0, pcrit, nPoints] = Sum[p[x0 + r Cos[phi], y0 + r Sin[phi], pcrit]/nPoints,
         {phi, 0., 2 Pi, 2 Pi/nPoints}]                                (In[362], line ~3627)
     the figure Block (x0 = x00, pcrit = 35258, anglePoints = 100, legFS = 9) (In[364], ~3790)
     export["final.pdf", fig]                                           (In[365], line ~4002)
   and the initialization cells (w = 8.6 cm, tex[] = MaTeX with amsmath,siunitx,upgreek,
   colors = ColorData[97]). The notebook reads the variance as readCovRho[dataPath]["rho"]
   (binaryReadTools.m: the diagonal of the covariance .bin, ListInterpolation with
   InterpolationOrder 2, x first) and modelParams = readSteady[dataPath]["modelParams"].

   What it computes, as the notebook does:
     P(x, y)  = 2 Qapprox(z),  z = pcrit/(B sqrt(Var[rho](x, y)))   (eq:p-layer, no fit)
     r_max    = (x00, 0), x00 the maximum of Var[rho] on the negative x axis (NMaximize;
                a 2D maximisation, printed, finds the same point, y = 0 to 1e-7)
     curve(r) = pRad[r, r_c] / pRad[12, r_max]: P (not the variance) averaged over 100 angles
                on the circle of radius r about r_c = r_max (blue) or r_c = 0 (orange),
                r = 0.1 .. 12 um, log scale.
   So the far-field value P(inf) of the panel is NOT the pedestal value of eq:sigma-p: it is
   the circle average of P at r = 12 um about r_max (the blue curve is 1 there by
   construction); it is 0.867 of the pedestal value for N = 180 and 200. norm=pedestal
   divides by the pedestal value instead (all curves move down by that factor).
   The angular sum runs over phi = 0, 2 Pi/100, ..., 2 Pi: 101 terms over 100, phi = 0 counted
   twice. Kept as in the notebook; a proper 100-angle mean changes the curves by < 1%
   (printed).

   pcrit. The figure cell sets pcrit = 35258 (35.3 kPa um, eq:pcrit with the parameters of
   app:parameters-extrusion, z_inf = 20.831 as in the manuscript), but the old final.pdf was
   drawn with 35000: the notebook compiles p with "// Evaluate" while the parameter cell
   (In[56]) has set the GLOBAL pcrit = 35000, so 35000 is substituted into the body before
   Compile binds its argument pcrit, and the 35258 passed by the figure cell is ignored.
   Checked: N=180 pcrit=35000 reproduces the curves of the old final.pdf to the pdf's
   coordinate precision (blue at r = 0.1: 174.06 vs 174.02, orange maximum 34.22 vs 34.21,
   orange at 12: 1.0583 vs 1.058); with 35258 the same run gives 187.8 and 36.65. This script
   compiles p with pcrit unset, so the pcrit it is given is the one used; default 35258, the
   value the cell and the manuscript intend.

   Changes from the notebook: the data source, the pcrit fix above, and the text rendering.
   - Data: the companion file sigma_rho_<tag>.m of the ssDir run (the Nint x Nint diagonal
     of the density covariance, first index y, hence the Transpose; it equals the .bin
     diagonal readCovRho reads), never the .bin; modelParams from steady_<tag>.m as
     readSteady does (checked equal to the companion file's). Default tag
     box20_N200_Dir_lN0.5_z12.5_ld1_B90000_ss20; the notebook used N180.
   - Tick labels: MaTeX (Computer Modern) instead of plain "Times" 8 pt text. The headless
     front end sets plain text in a different serif (and garbles it in some exports); the
     old final.pdf has CM-like digits 6 pt high on the page. The MaTeX labels are 8.4 pt,
     trimmed by 1 pt and moved down slightly, so that size and position match the old digits
     to about 0.1 pt (measured on the rasterized pdfs). Tick positions are the automatic ones,
     read from a first rendering with AbsoluteOptions; the labelled x ticks are lengthened to
     the 2.0 pt they have in the old pdf.
   - ImagePadding is set to the frame position measured on the old final.pdf (frame at
     x = 27.74..174.34, y = 12.15..102.75 pt on its 176 x 130 pt page; the pdf export scales
     design points by 0.72), so the new panel drops into main.pdf at the same place. The
     frame, ticks and grid line of N=180 pcrit=35000 coincide with the old pdf's.
   Everything else -- curves, colours, thicknesses, PointSize, error-bar style, legend text
   and position, frame labels, grid line at 1, ranges, width -- is the notebook's.

   Usage:  plotgen/layering_profile/run.sh [N=200|180|150] [pcrit=<value>] [norm=pedestal] [out=<dir>]
     N=<n>          the ssDir run with Nint = n (default 200)
     pcrit=<value>  critical pressure in Pa um (default 35258; 35000 reproduces the old pdf)
     norm=pedestal  normalise by P at the pedestal variance instead of pRad[12, r_max]
     out=<dir>      output directory (default: next to this script)
   Outputs: layering_profile.pdf and layering_profile.png (600 dpi); with non-default options
   the name gets _N<n>, _pc<pcrit>, _ped. run.sh (no arguments) copies layering_profile.pdf
   to manuscript/figures/final.pdf. Numbers are printed: r_max, peak enhancement, maxima of
   the curves, monotonicity of the blue curve. About 30 s.
   Results 2026-10-26, pcrit 35258 (z_inf = 20.831), N200 [N180]: r_max = (-1.250, 0)
   [(-1.268, 0)] um, nu0 = 2.424% [2.400%]; P(r_max)/pRad[12, r_max] = 198.1 [188.6],
   P(r_max)/P(pedestal) = 171.8 [163.5]; blue largest at the plot start r = 0.1, 197.0 [187.8];
   orange maximum 38.28 [36.65] at r = 1.225 [1.223]; blue decreases to 0.967 at r = 7.7 and
   rises back to 1 at r = 12 (not monotonic beyond 7.7 um, by 3%). Log of the runs:
   claude_experiments/round_2026_10_26/layering_profile_run.log.
   2026-11-02, the author's request: the x label is "Distance from the center r" (was "from
   the core": r is measured from the centre of the circle, which for the blue curve is
   r_max), and the legend reads "Center at P_max", "Center at the core" (was "Core at
   P_max", "Core at origin"), as the caption words it; the legend is anchored in the upper
   right corner of the frame.
   2026-11-02, the author's request: the axis labels sit closer to the axes, about 3.8 design
   pt (0.45 of the tick font size, as in the restyled figures) from the tick labels, placed by
   plotgen/frame_labels.wl instead of FrameLabel. The page and the frame stay where they were,
   so the panel still drops into main.pdf at the same place; the space the labels free stays
   white at the left and bottom of the page. *)

SetEnvironment["QT_QPA_PLATFORM" -> "offscreen"];

here = DirectoryName[$InputFileName];
repoRoot = ParentDirectory[ParentDirectory[here]];
args = Rest[$ScriptCommandLine];
opts = Association[Rule @@@ (StringSplit[#, "=", 2] & /@ Select[args, StringContainsQ[#, "="] &])];
nInt = Lookup[opts, "N", "200"];
pcritV = ToExpression@Lookup[opts, "pcrit", "35258"];
normPed = Lookup[opts, "norm", "r12"] === "pedestal";
outDir = Lookup[opts, "out", here];
name = "layering_profile" <> If[nInt === "200", "", "_N" <> nInt] <>
   If[pcritV === 35258, "", "_pc" <> ToString[pcritV]] <> If[normPed, "_ped", ""];
t0 = AbsoluteTime[];

(* ---------------------------------------------------------------- data *)
tag = "box20_N" <> nInt <> "_Dir_lN0.5_z12.5_ld1_B90000_ss20";
dataDir = FileNameJoin[{repoRoot, "data", "fluctuations-set2_2_comove_optimized_2_ssDir"}];
sig = Import[FileNameJoin[{dataDir, "sigma_rho_" <> tag <> ".m"}]];
modelParams = Import[FileNameJoin[{dataDir, "steady_" <> tag <> ".m"}]]["modelParams"];
If[Sort[N@modelParams] =!= Sort[N@sig["modelParams"]],
  Print["*** WARNING: the steady state's parameters differ from the companion file's"]];
grid = N@sig["gridInt"];
(* readCovRho[...]["rho"]: covFn = ListInterpolation[mat, {grid, grid}, InterpolationOrder -> 2],
   mat x first; the companion file stores it y first *)
sigRhoFnDirichlet = ListInterpolation[Transpose[sig["sigmaRho"]], {grid, grid},
   InterpolationOrder -> 2];
Print["run ", FileNameTake[dataDir], "/", tag, ": grid ", Length[grid], " x ", Length[grid],
  " on ", MinMax[grid], ", B = ", B /. modelParams];

dataExperiment = Apply[Function[{r, low, mean, high},
     {r, Around[mean, {mean - low, high - mean}]}],
    Rest@Import[FileNameJoin[{repoRoot, "top_plot.csv"}]], 1];

(* ---------------------------------------------------------------- the notebook's functions *)
(* pcrit has no value here, so the compiled p uses its argument (see the header) *)
Qapprox[z_] = Exp[-z^2/2]/(z Sqrt[2 \[Pi]]);
p = Compile[{{x, _Real}, {y, _Real}, {pcrit, _Real}},
   2 Qapprox[Sqrt[(pcrit^2)/(B^2 sigRhoFnDirichlet[x, y])]] /. modelParams // Simplify // Evaluate];
x00 = -x /. First@Rest@NMaximize[sigRhoFnDirichlet[-x, 0] /. modelParams // Evaluate, {x, 0, 5}];
pRad[r_, x0_, y0_, pcrit_, nPoints_ : 1000] := Sum[
   With[{x = x0 + r Cos[\[Phi]], y = y0 + r Sin[\[Phi]]}, p[x, y, pcrit]/nPoints],
   {\[Phi], 0., 2 \[Pi], 2 \[Pi]/nPoints}];

pcrit = pcritV;         (* the figure cell's 35258 by default; set only after p is compiled *)
anglePoints = 100;
x0 = x00;
pedestal = \[Zeta]/(2 \[Pi] B \[Rho]0 lNoise^2) /. modelParams;   (* eq:pedestal *)
pPed = 2 Qapprox[N[pcrit/((B /. modelParams) Sqrt[pedestal])]];  (* P at the pedestal variance *)
norm = If[normPed,
   Sum[pPed/anglePoints, {\[Phi], 0., 2 \[Pi], 2 \[Pi]/anglePoints}],   (* pRad of a uniform field *)
   pRad[12, x0, 0, pcrit, anglePoints]];                                (* the notebook's *)
blue[r_?NumericQ] := pRad[r, x0, 0, pcrit, anglePoints]/norm;
orange[r_?NumericQ] := pRad[r, 0, 0, pcrit, anglePoints]/norm;

(* ---------------------------------------------------------------- numbers *)
Module[{Bv, ped, zInf, nu0, pMax2D, rs, bs, os, db, iMin, ext, bMax, oMax, mean100, norm100, dev},
  Bv = B /. modelParams;
  ped = pedestal;
  zInf = N[pcrit/(Bv Sqrt[ped])];
  nu0 = sigRhoFnDirichlet[x0, 0]/ped - 1;
  Print["pcrit = ", pcrit, " (kPa um x 1e3), anglePoints = ", anglePoints, ", z_inf = ", zInf,
    ", sigma_p,inf = ", Bv Sqrt[ped]];
  Print["r_max (notebook: max of Var on the -x axis): (x, y) = (", x0, ", 0) um;  nu0 = Var/ped - 1 = ",
    100 nu0, " %"];
  pMax2D = NMaximize[{sigRhoFnDirichlet[x, y], -5 <= x <= 5, -5 <= y <= 5}, {x, y}];
  Print["2D check: max of Var over |x|,|y| <= 5 at (x, y) = ", {x, y} /. Last[pMax2D],
    ",  nu0 = ", 100 (First[pMax2D]/ped - 1), " %"];
  Print["far-field normalisation: ", If[normPed, "P(pedestal) (norm=pedestal)", "pRad[12, r_max] (notebook)"],
    ";  pRad[12, r_max] / P(pedestal) = ", pRad[12, x0, 0, pcrit, anglePoints]/pPed];
  Print["peak enhancement P(r_max)/pRad[12, r_max] = ", p[x0, 0., pcrit]/pRad[12, x0, 0, pcrit, anglePoints],
    ";  P(r_max)/P(pedestal) = ", p[x0, 0., pcrit]/pPed,
    ";  eq:R-value exp[z^2/2 nu0/(1+nu0)] = ", Exp[zInf^2/2 nu0/(1 + nu0)]];
  rs = Range[0.1, 12., 0.01];
  bs = blue /@ rs; os = orange /@ rs;
  ext[f_, r1_] := Module[{g, s, rr}, g[rr_?NumericQ] := f[rr];
     s = Quiet@FindMaximum[{g[rr], 0.1 <= rr <= 12}, {rr, r1}]; {s[[1]], rr /. s[[2]]}];
  bMax = {Max[bs], rs[[First@Ordering[bs, -1]]]};
  oMax = ext[orange, rs[[First@Ordering[os, -1]]]];
  Print["blue (about r_max):  max ", bMax[[1]], " at r = ", bMax[[2]], " (plot start r = 0.1);  r = 0: ",
    blue[0.]];
  Print["orange (about core): max ", oMax[[1]], " at r = ", oMax[[2]], ";  r = 0.1: ", os[[1]]];
  db = Differences[bs];
  If[AllTrue[db, # < 0 &],
   Print["blue decreases monotonically on r = 0.1..12 (step 0.01)"],
   iMin = First@Ordering[bs, 1];
   Print["blue is NOT monotonic: decreases to ", bs[[iMin]], " at r = ", rs[[iMin]],
     ", then rises to ", Last[bs], " at r = 12 (first increase at r = ",
     rs[[First@FirstPosition[db, _?(# >= 0 &)]]], ")"]];
  Print["values at r = 1, 2, 3, 6, 11:  blue ", blue /@ {1., 2., 3., 6., 11.}, ",  orange ",
    orange /@ {1., 2., 3., 6., 11.}];
  (* the angular sum counts phi = 0 twice; a proper 100-angle mean, for the record *)
  mean100[r_, xc_] := Mean[p[xc + r Cos[#], r Sin[#], N@pcrit] & /@ (2 \[Pi] Range[0, 99]/100.)];
  norm100 = If[normPed, pPed, mean100[12, x0]];
  dev = Max[Abs[(mean100[#, x0]/norm100)/blue[#] - 1] & /@ Range[0.1, 12, 0.1]];
  Print["proper 100-angle mean instead of the notebook's 101-term sum: max relative change of blue ",
    dev, ", of orange ", Max[Abs[(mean100[#, 0.]/norm100)/orange[#] - 1] & /@ Range[0.1, 12, 0.1]]];
  Print["numbers done (", Round[AbsoluteTime[] - t0], " s)"]];

(* ---------------------------------------------------------------- figure *)
<< MaTeX`
Get[FileNameJoin[{repoRoot, "plotgen", "frame_labels.wl"}]];
w = 8.6*72/2.54;
tex[s_, fs_ : 11] := MaTeX[s, FontSize -> fs,
   "Preamble" -> {"\\usepackage{amsmath,siunitx,upgreek}"}];
colors = ColorData[97, "ColorList"];
legFS = 9;
darkred = Darker[Red, 0.15];

(* MaTeX boxes carry about 1 pt of margin on each side; trimmed by d so the labels sit where
   plain text would (as in plotgen/varrho/varrho.wls), the content moved down by s and the
   box cut by c at the bottom, in design pt. Size and offsets match the digits of the old
   final.pdf to about 0.1 pt (its plain-text labels print 6 pt high and are centred on the
   y ticks): 8.4 pt instead of 8, s = 0.67 on y, 0.21 on x; c = 1 on x puts the x frame
   label, which sits below the tick-label boxes, where it was. *)
trim[g_, d_, s_ : 0, c_ : 0] := Module[{pr = PlotRange /. Options[g], is = First@Cases[Options[g], (ImageSize -> v_) :> v]},
   Graphics[First[g], PlotRange -> {pr[[1]] + {d, -d}, pr[[2]] + s + {c, 0}}, ImageSize -> is - {2 d, c},
    AspectRatio -> Automatic, BaselinePosition -> (BaselinePosition /. Options[g])]];
tickText[v_] := If[Abs[v - Round[v]] < 1.*^-9, ToString[Round[v]], ToString[v]];
tickFS = 8.4;
tickLabel[v_, s_, c_] := trim[tex[tickText[v], tickFS], 1., s, c];

(* the frame of the old final.pdf, in design points (pdf points / 0.72): left 27.738, right
   175.52 - 174.34, top 12.145; the bottom keeps the page at 130 pt *)
pad = {{27.738, 175.52 - 174.34}, {26.88, 12.145}}/0.72;

curves = LogPlot[{blue[r], orange[r]}, {r, 0.1, 12}, PlotRange -> All,
   GridLines -> {None, {1}},
   Frame -> True,
   FrameStyle -> Directive[Black, FontFamily -> "Times", FontSize -> 8],
   PlotStyle -> Thickness[0.008],
   FrameTicks -> {{Automatic, None}, {Automatic, None}},
   PlotLegends -> Placed[LineLegend[
      {Directive[colors[[1]], Thickness[0.004]],
       Directive[colors[[2]], Thickness[0.004]],
       Directive[darkred, Dashed]},
      {tex["\\text{Center at }{\\cal P}_\\text{max}", legFS], tex["\\text{Center at the core}", legFS],
       tex["\\text{Experiment}", legFS]},
      LegendMarkers -> {None, None,
        Graphics[{darkred, EdgeForm[None], Disk[{0, 0}, 0.3]}, PlotRange -> {{-1, 1}, {-1, 1}}]},
      LabelStyle -> {FontSize -> 8, FontFamily -> "Times New Roman"},
      LegendMarkerSize -> 15,
      Spacings -> {0.4, -0.5}], {{0.98, 0.97}, {1, 1}}]];   (* centred at Scaled[{0.80, 0.84}]
     until 2026-11-02; the longer labels ran past the frame, so the corner is anchored *)
points = ListLogPlot[dataExperiment, PlotRange -> All,
   PlotStyle -> Directive[PointSize[0.02], darkred],
   IntervalMarkersStyle -> Directive[Dashing[None], Thickness[0.004]]];
Print["curves sampled (", Round[AbsoluteTime[] - t0], " s)"];

(* the notebook's automatic ticks, read off a first rendering, relabelled in MaTeX *)
(* (same ImagePadding as the final figure: the tick lengths come back scaled to the frame) *)
prelim = Show[curves /. Legended[g_, ___] :> g, points, ImageSize -> w, ImagePadding -> pad];
{yAuto, xAuto} = First /@ (FrameTicks /. AbsoluteOptions[prelim, FrameTicks]);
(* unlabelled (minor) ticks come back as "" or as a Spacer2 box *)
unlabelled[lab_] := MatchQ[lab, "" | Null | {}] || ! FreeQ[lab, "Spacer2" | Spacer];
relabel[ticks_, toValue_, s_, c_] := Replace[ticks,
   {pos_, lab_, rest___} :> {pos, If[unlabelled[lab], "", tickLabel[toValue[pos], s, c]], rest}, {1}];
xTicks = relabel[xAuto, Identity, 0.21, 1.];
yTicks = relabel[yAuto, Exp, 0.67, 0.];
(* the labelled x ticks print 2.0 pt long in the old pdf (automatic ticks); the explicit list
   AbsoluteOptions returns draws them 1.73 pt long, so they are lengthened by 2.0/1.727 *)
xTicks = Replace[xTicks, {pos_, lab_, {l1_, l2_}, rest___} /; lab =!= "" :>
     {pos, lab, {l1 2.0/1.727, l2}, rest}, {1}];
Print["labelled ticks: x ", Cases[xAuto, {v_, l_, ___} /; ! unlabelled[l] :> v],
  ",  y ", Cases[yAuto, {v_, l_, ___} /; ! unlabelled[l] :> Round[Exp[v], 0.01]]];

(* the axis labels as insets, labelDist from the frame to the label boxes (design pt) *)
labelDist = {13.3, 17.5};
fig = Show[curves, points, ImageSize -> w,
   FrameTicks -> {{yTicks, None}, {xTicks, None}},
   Epilog -> frameLabels[{tex["\\text{Distance from the center}~~r~(\\upmu \\rm{m})"],
       tex["\\text{Layering probability}~~{\\cal P}(r) / {\\cal P}(\\infty)"]}, labelDist],
   PlotRangeClipping -> False,      (* the labels lie outside the frame; the data stay inside *)
   ImagePadding -> pad];

Export[FileNameJoin[{outDir, name <> ".pdf"}], fig];
Export[FileNameJoin[{outDir, name <> ".png"}], fig, ImageResolution -> 600];
Print["exported ", FileNameJoin[{outDir, name <> ".pdf"}], " and .png  (", Round[AbsoluteTime[] - t0], " s)"];
