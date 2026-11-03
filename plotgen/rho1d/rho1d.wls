#!/usr/bin/env wolframscript
(* rho1d.wls -- three figures of the steady comoving +1/2 defect, written 2026-11-01:
     fig:rho1d-dirichlet          rho1d_Dirichlet.pdf          (rho - 1)/(zeta/B) along y = 0,
                                  reservoir wall: numerical and closed form eq:rho-analytic
     fig:rho1d-dirichlet-neumann  rho1d_Dirichlet_Neumann.pdf  the same with the no-flux wall
                                  and eq:rho-noflux
     fig:nematicaxis              nematicaxis.pdf              director lines, numerical and Pade

   Source: the notebook remote_kernel/plotting_small_core+ln_2.nb (the copy of 2026-10-26; the
   older copy in negative-defect/ has newtx labels and other sizes), whose cells made the
   figures in manuscript/figures; this script runs them, with the style changes listed below:
     tex, w, colors, ArcTan2, export   cells 4-8 (MaTeX in Computer Modern; export: pdf through
                                       RemoteExport, png at 900 dpi as <name>_high.png and at
                                       600 dpi as <name>.png)
     modelParams, dataPath             readSteady of the ssDir N_int = 180 run (B = 9e4)
     Pade profile                      S0Pade, Q1Pade, Q2Pade
     closed forms                      rhoSteady, bcDirichlet, bc0J (the no-flux wall)
     reservoir wall                    the saved steady state of that run (solSteadyDirichletSaved:
                                       Dirichlet density wall, free Q), for rho and the director
     no-flux wall                      comoving_steady_Qfree_solver.m, "Neumann" (= no flux), with
                                       delta mesh 0.09, hMax 1, box 20 (the notebook's
                                       solSteadyNeumann; its solve cell is no longer in the notebook)
   Check of 2026-11-01: run with the notebook's own options, the script reproduced the pdfs of
   Figs. 5 and 6 in manuscript/figures pixel for pixel. Fig. 11 there was an older export, with
   newtx labels and the reservoir curve from a fresh solve (2e-4 off the saved state at the
   minimum).

   Style, changed on the author's requests of 2026-11-01:
   - the axis labels sat too far from the axes: they are placed by plotgen/frame_labels.wl,
     about 4 pt (printed) from the tick labels, with an explicit ImagePadding;
   - one font and one set of sizes in all figures: all text in MaTeX Computer Modern, ticks
     included, at 9 pt for tick and axis labels and 8 pt for legends as printed (design sizes
     scaled by the width in manuscript.tex: 0.85 \linewidth, 0.7 for nematicaxis);
   - nematicaxis: no title ("Nematic axis" repeated the caption);
   - rho1d_Dirichlet_Neumann: the legend names the walls as the caption does, reservoir and no
     flux (it said "Dirichlet b.c.", "Neumann b.c.").
   2026-11-02, the author's requests:
   - no line prints thinner than 0.75 pt (frame_labels.wl): frame and ticks 0.75 pt (0.17 pt
     before), the axes x = 0 and y = 0 gray at 0.75 pt (black 0.43 pt before);
   - legend samples as thick as the curves (Thickness in a legend is relative to the sample);
   - rho1d_Dirichlet_Neumann: the legend stays clear of the axis x = 0 (shorter samples,
     smaller margins);
   - nematicaxis: director lines 0.75 pt (0.35 pt before); the legend box tight around its
     rows, in the upper right corner.
   2026-11-02, later requests of the author:
   - the lines x = 0 and y = 0 are grid lines (gridStyle of frame_labels.wl, 0.6 pt, gray 0.5),
     drawn under the curves; as axes they lay on top;
   - the y label on one line, (rho - 1)/(zeta/B), without the \tfrac; it is lower, so the
     left padding shrinks by 2.8 pt and the label moves 1 pt closer, to keep the 4 pt gap;
   - legend rows 9.5 pt apart as printed, the line spacing of 8 pt text in LaTeX (12.1 pt
     before).

   The no-flux solve takes about a minute; its result is cached in rho1d_neumann.mx next to
   this script and reused ("resolve" redoes it).
   Usage:  plotgen/rho1d/run.sh [resolve] [out=<dir>]
   Outputs: rho1d_Dirichlet, rho1d_Dirichlet_Neumann, nematicaxis (.pdf, .png, _high.png) in
   out=<dir>, default next to this script. *)

SetEnvironment["QT_QPA_PLATFORM" -> "offscreen"];
$HistoryLength = 0;

here = DirectoryName[$InputFileName];
repoRoot = ParentDirectory[ParentDirectory[here]];
args = Rest[$ScriptCommandLine];
resolve = MemberQ[args, "resolve"];
opts = Association[Rule @@@ (StringSplit[#, "=", 2] & /@ Select[args, StringContainsQ[#, "="] &])];
figDir = Lookup[opts, "out", here];
If[! DirectoryQ[figDir], CreateDirectory[figDir]];
t0 = AbsoluteTime[];
stamp[] := ToString[Round[AbsoluteTime[] - t0]] <> " s";

(* ---------------------------------------------------------------- notebook set-up *)
<< MaTeX`
w = 8.6*72/2.54;
tex[s_, fs_ : 11] := MaTeX[s, FontSize -> fs,
   "Preamble" -> {"\\usepackage{amsmath,siunitx,upgreek}"}];
colors = ColorData[97, "ColorList"];
ArcTan2[x_, y_] := Pi - ArcTan[-x, y];

Import[FileNameJoin[{repoRoot, "RemoteExport.m"}]];
export[name_, fig_] := Block[{base = FileNameJoin[{figDir, name}]},
   RemoteExport[base <> ".pdf", fig];
   Export[base <> "_high.png", fig, ImageResolution -> 900];
   Export[base <> ".png", fig, ImageResolution -> 600];
   Print["exported ", base, ".pdf  (", stamp[], ")"]];
Get[FileNameJoin[{repoRoot, "plotgen", "frame_labels.wl"}]];

Get[FileNameJoin[{repoRoot, "binaryReadTools.m"}]];
dataPath = FileNameJoin[{repoRoot,
    "data/fluctuations-set2_2_comove_optimized_2_ssDir/covariance_box20_N180_Dir_lN0.5_z12.5_ld1_B90000_ss20.bin"}];
solSteadyDirichletSaved = readSteady[dataPath];
modelParams = solSteadyDirichletSaved["modelParams"];
Print["modelParams ", modelParams];

(* the Pade profile *)
fPade[z_] := Sqrt[(0.34010999 + 0.04516965*z + 0.00233656*z^2)/(1 + 0.3838716*z +
      0.04750621*z^2 + 0.00233656*z^3)];
S0PadeOverR[r2_] := Sqrt[-a1/(2*b)]*Sqrt[-a1/L1]*fPade[(-r2)*(a1/L1)];
S0Pade[r_] := r*S0PadeOverR[r^2];
Q0Pade[xv_, yv_] := S0PadeOverR[xv^2 + yv^2]*{{xv, yv}, {yv, -xv}};
Q1Pade[xx_, yy_] = Q0Pade[xx, yy][[1, 1]];
Q2Pade[xx_, yy_] = Q0Pade[xx, yy][[1, 2]];

(* the closed forms of the steady density on a disc of radius R *)
P1[r_] = r^2*Derivative[1][S][r] - r*S[r];
Qu1[r_] = Derivative[1][S][R] - Derivative[1][S][r] + S[R]/R - S[r]/r;
u[r_] := A1*r + B1/r - P1[r]/(2*r) - (r*Qu1[r])/2;
rhoSteady[x_, y_] := With[{r = Sqrt[x^2 + y^2], \[Phi] = ArcTan[x, y]},
   (-(\[Zeta]/B))*(u[r]*Cos[\[Phi]] + h)];
bcDirichlet = {B1 -> 0, A1 -> P1[R]/(2*R^2), h -> -(B/\[Zeta])};
bc0J = Simplify[{B1 -> 0, A1 -> p1 - P1[R]/(2*R^2), h -> -(B/\[Zeta])} /.
    {p1 -> Derivative[1][S][R] + S[R]/R}];

(* ---------------------------------------------------------------- the no-flux solve *)
solverParams = {\[Delta]mesh -> 0.09, hMax -> 1, box -> 20, zetaSteps -> 2, maxSteps -> 25,
   relativeDiff -> 1.*^-8, "ReuseFactorization" -> False};
cacheFile = FileNameJoin[{here, "rho1d_neumann.mx"}];
If[FileExistsQ[cacheFile] && ! resolve,
  solSteadyNeumann = Import[cacheFile];
  Print["no-flux solution from the cache ", cacheFile],
  Import[FileNameJoin[{repoRoot, "comoving_defect", "comoving_steady_Qfree_solver.m"}]];
  solSteadyNeumann = SolveActiveNematicComovingSteadyQfree[modelParams, solverParams, "Neumann"];
  Print["Neumann solved (", stamp[], ")"];
  Export[cacheFile, KeyTake[solSteadyNeumann, {"rho", "Q1", "Q2", "u"}], "MX"]];
Print["u: Dirichlet (saved) ", solSteadyDirichletSaved["u"], ", Neumann ", solSteadyNeumann["u"]];
zB = \[Zeta]/B /. modelParams;
Print["min (rho - 1)/(zeta/B) on y = 0: Dirichlet ",
  FindMinimum[(solSteadyDirichletSaved["rho"][x, 0] - 1)/zB, {x, 2.5, 0.5, 10}],
  ";  at x = -20, 20: Neumann ", (solSteadyNeumann["rho"][#, 0] - 1)/zB & /@ {-20, 20}];

(* ---------------------------------------------------------------- figures *)
(* sizes: design pt that print at pt when the figure (width w) is included at frac \linewidth *)
printed[pt_, frac_] := pt w/(frac 246.);
{fsT, fsL} = printed[{9, 8}, 0.85];      (* tick and axis labels 9 pt, legends 8 pt *)
texT[str_] := tex[str, fsT];
fR = printed[1, 0.85];                   (* design pt per printed pt *)
curveT = AbsoluteThickness[1.2];         (* 1.03 pt printed, as Thickness[0.006] of the frame *)
rhoLabels = texT /@ {"x~(\\upmu \\rm{m})", "(\\rho-1)/(\\zeta/B)"};
rhoOpts[yTicks_] := Sequence[
   FrameTicks -> {{texTicks[yTicks, texT], None}, {texTicks[Range[-20, 20, 10], texT], None}},
   FrameStyle -> frameStyle[fR], Axes -> False,
   GridLines -> {{0}, {0}}, GridLinesStyle -> gridStyle[fR],   (* under the curves *)
   Epilog -> frameLabels[rhoLabels, {15.6, 17.9}],
   PlotRangeClipping -> False, ImagePadding -> {{31.5, 8}, {29, 4}}];
(* legends: samples as thick as the curves, no margin beyond the rows, rows 9.5 pt apart *)
rhoLegend[styles_, labels_, size_] := LineLegend[Directive[curveT, #] & /@ styles,
   tex[#, fsL] & /@ labels, LegendMarkerSize -> size, Spacings -> {0.3, -0.2},
   LegendMargins -> 0];

fig = With[{box = box /. solverParams},
   Plot[Evaluate[({solSteadyDirichletSaved["rho"][x, 0],
         rhoSteady[x, 0] /. bcDirichlet /. {S -> S0Pade} /. {R -> box} /. modelParams} - 1)/
      (\[Zeta]/B) /. modelParams], {x, -box, box},
    PlotStyle -> (Directive[curveT, #1] &) /@ {colors[[1]], Directive[Black, Dashed]},
    PlotLegends -> Placed[rhoLegend[{colors[[1]], Directive[Black, Dashed]},
       {"\\text{Numerical}", "\\text{Analytical}"}, {16, 9}], {{0.98, 0.97}, {1, 1}}],
    ImageSize -> w, PlotRange -> {{-21, 21}, {-1.31, 1.1}}, Frame -> True,
    Evaluate[rhoOpts[{-1, 0, 1}]]]];
export["rho1d_Dirichlet", fig];

fig = With[{box = box /. solverParams},
   Plot[Evaluate[({solSteadyDirichletSaved["rho"][x, 0],
         rhoSteady[x, 0] /. bcDirichlet /. {S -> S0Pade} /. {R -> box} /. modelParams,
         solSteadyNeumann["rho"][x, 0],
         rhoSteady[x, 0] /. bc0J /. {S -> S0Pade} /. {R -> box} /. modelParams} - 1)/
      (\[Zeta]/B) /. modelParams], {x, -box, box},
    PlotStyle -> (Directive[curveT, #1] &) /@ {colors[[1]],
       Directive[colors[[1]], Dashed], colors[[2]], Directive[colors[[2]], Dashed]},
    PlotLegends -> Placed[rhoLegend[{colors[[1]], colors[[2]], Directive[Black, Dashed]},
       {"\\text{Reservoir, }\\rho=1", "\\text{No flux}", "\\text{Analytical}"}, {16, 9}],
      {{0.995, 0.98}, {1, 1}}],
    ImageSize -> w, PlotRange -> {{-21, 21}, {-2.1, 2.1}}, Frame -> True,
    Evaluate[rhoOpts[Range[-2, 2, 1]]]]];
export["rho1d_Dirichlet_Neumann", fig];

(* the director lines run on to the edge of the StreamPlot region, |x|, |y| <= 12, beyond the
   plot range, |x|, |y| <= 10, so without clipping (needed for the labels) they are cut here *)
clipSeg[{p1_, p2_}, W_] := Module[{t0 = 0., t1 = 1., d = p2 - p1, ok = True},
   Do[With[{pp = c[[1]], qq = c[[2]]},
     If[pp == 0., If[qq < 0., ok = False],
      If[pp < 0., t0 = Max[t0, qq/pp], t1 = Min[t1, qq/pp]]]],
    {c, {{-d[[1]], p1[[1]] + W}, {d[[1]], W - p1[[1]]}, {-d[[2]], p1[[2]] + W},
      {d[[2]], W - p1[[2]]}}}];
   If[ok && t0 <= t1, {p1 + t0 d, p1 + t1 d, t0 > 0., t1 < 1.}, None]];
clipLine[pts_, W_] := Module[{runs = {}, cur = {}, s},
   Do[s = clipSeg[seg, W];
    If[s === None, If[Length[cur] > 1, AppendTo[runs, cur]]; cur = {},
     If[cur === {} || s[[3]], If[Length[cur] > 1, AppendTo[runs, cur]]; cur = {s[[1]]}];
     AppendTo[cur, s[[2]]];
     If[s[[4]], AppendTo[runs, cur]; cur = {}]],
    {seg, Partition[pts, 2, 1]}];
   If[Length[cur] > 1, AppendTo[runs, cur]];
   Line /@ runs];
clipGraphics[g_Graphics, W_] := MapAt[ReplaceAll[#, {
      (Line | Arrow)[pts : {{_?NumericQ, _?NumericQ} ..}] :> clipLine[pts, W],
      (Line | Arrow)[l : {{{_?NumericQ, _?NumericQ} ..} ..}] :> (clipLine[#, W] & /@ l)}] &,
   Normal[g], 1];
clipGraphics[Legended[g_, l___], W_] := Legended[clipGraphics[g, W], l];

{fsT6, fsL6} = printed[{9, 8}, 0.7];
f6 = printed[1, 0.7];
dirStyles = {Directive[Black, Dashed, thinLine[f6]], Directive[Red, thinLine[f6]]};
fig = With[{W = 12, W2 = 10, Q1 = solSteadyDirichletSaved["Q1"], Q2 = solSteadyDirichletSaved["Q2"]},
   Module[{seeds, g},
    seeds = Join[({#1, 0.001} &) /@ Range[-11, -3, 2], {{-1.5, 0}, {10, 0.0001}, {9, 1.5},
       {9, -1.5}, {1, 1.5}, {1, -1.5}}];
    g = StreamPlot[{With[{\[Theta] = ArcTan2[Q1Pade[x, y], Q2Pade[x, y]]/2 /. modelParams},
        {Cos[\[Theta]], Sin[\[Theta]]}], With[{\[Theta] = ArcTan2[Q1[x, y], Q2[x, y]]/2},
        {Cos[\[Theta]], Sin[\[Theta]]}]}, {x, -W, W}, {y, -W, W},
      StreamScale -> None, StreamColorFunction -> None,
      StreamStyle -> dirStyles, StreamPoints -> {seeds, 0.2},
      PlotRange -> {{-W2, W2}, {-W2, W2}}, ImageSize -> w, FrameStyle -> frameStyle[f6],
      FrameTicks -> {{texTicks[Range[-10, 10, 5], tex[#, fsT6] &], None},
        {texTicks[Range[-10, 10, 5], tex[#, fsT6] &], None}},
      (* the box tight around the rows: no legend margin, 3 pt of frame margin *)
      PlotLegends -> Placed[LineLegend[dirStyles,
         tex[#, fsL6] & /@ {"\\text{Pad\\'e approximation}", "\\text{Numerical solution}"},
         LegendMarkerSize -> {20, 10}, LegendMargins -> 0,
         LegendFunction -> (Framed[#1, RoundingRadius -> 3, Background -> RGBColor[1, 1, 1, 0.8],
             FrameStyle -> refStyle[f6], FrameMargins -> 3] &), Spacings -> {0.3, 0.1}],
        {{0.98, 0.98}, {1, 1}}],
      Epilog -> frameLabels[tex[#, fsT6] & /@ {"x~(\\upmu \\rm{m})", "y~(\\upmu \\rm{m})"},
        {18.5, 27.5}],
      PlotRangeClipping -> False, ImagePadding -> {{44, 8}, {33.8, 6}}];
    clipGraphics[g, W2]]];
export["nematicaxis", fig];

Print["done (", stamp[], ")"];
