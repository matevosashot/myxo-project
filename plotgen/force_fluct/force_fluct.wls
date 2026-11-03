#!/usr/bin/env wolframscript
(* force_fluct.wls -- the two figures on the fluctuating active force of the Lyapunov appendix
   (app:lyapunov-results), regenerated from the N_int = 200 run. Written 2026-10-26.

     fig:force-axis  fActiveMeanAndStd.pdf  mean <f_x> and Std[f_x], Std[f_y] along y = 0, um^-1
     fig:force-map   f2.png, f2_high.png    <|delta f|^2> = <(d_b dQ_xb)^2> + <(d_b dQ_yb)^2>
                                            around the core, um^-2
   with f_a = d_b Q_ab (the active force is -zeta f).

   Source: the subsection "Active Force Fluctuations" of the notebook
   remote_kernel/plotting_small_core+ln_2.nb, which this script runs cell for cell:
     data   res = readActiveStress[dataPath]          (In[346], near line 3146)
            modelParams = readSteady[dataPath]["modelParams"]   (In[43], near line 563)
     mean   fActive[x, y] from Q1Pade, Q2Pade          (In[350], near line 3270; the Pade
            profile of the static core, In[17] near line 330, with the run's modelParams)
     figs   the DensityPlot loop over fxfx, fyfy, f2   (In[349], near line 3160; export near
            line 3232) and the Plot of fig              (In[351]; export near line 3438)
     tex, w, export                                    (In[1], In[41]: pdf through
            RemoteExport, png at 900 dpi as <name>_high.png and at 600 dpi as <name>.png)
   Only dataPath changes: the notebook reads the ssDir N_int = 180 run, this script by default
   the ssDir N_int = 200 run, data/fluctuations-set2_2_comove_optimized_2_ssDir/
   covariance_box20_N200_Dir_lN0.5_z12.5_ld1_B90000_ss20.bin (h = 0.199 um against 0.221 um).

   How the derivative variances are read. <(d_al dQ_a)(d_be dQ_b)> at one point needs the
   covariance of dQ between neighbouring nodes, not only its diagonal. readActiveStress.m
   (repo root, the notebook's reader) applies the solver's 2-point central difference to both
   sides, so it touches only the nine shifted diagonals Delta in {0, +-2, +-2N, N+-1, -N+-1}
   of the Q1Q1, Q2Q2 and Q1Q2 blocks: one contiguous read of 4N+3 reals per matrix row, about
   0.8 GB in all and 0.4 GB of memory at N = 200, instead of three full N^4 blocks (12.8 GB
   each) for readCorrQ. The fields live on gridF = gridInt[[2 ;; -2]]. Reading the 115 GB
   file takes about 7 min from cold (/data); "cache=<file.mx>" keeps the arrays (1 MB).
   Check of 2026-10-26: run on the N = 180 file (N=180), it reproduces the old figures, e.g.
   Std[f_x](x, 0) = 2.676 at x = -0.8 and 2.761 at x = 1.3, max <|delta f|^2> = 13.57.

   Export. Every figure is exported in a parallel subkernel, as the notebook's RemoteExport
   did for the pdf: exported from the script's own kernel, the headless front end drops the
   PlotLegends (LineLegend) and crashes on the BarLegend. Ticks stay plain Times text (it
   survives this route, as in the old pdf). The pdf comes out as the old one did (176 x 132 pt,
   the same fonts). The pngs have the old layout at 0.8 of its pixel size (f2.png 2462 x 2037
   against 3084 x 2550): the old ones came from the author's desktop front end, which
   rasterised at magnification 1.25; included at a fixed width, they print the same.

   The script also prints the numbers behind the text of app:lyapunov-results (std/mean along
   the axis, the maxima of Std[f_x] and Std[f_y], the ring of <|delta f|^2> and the share of
   <f_x^2> on it, f_x against f_y at the core and far away) and compares the mean of the
   notebook (Pade profile of the static core) with the mean from the numerical steady state
   the Lyapunov run is linearised about (steady_<tag>.m, derivatives of its order-3
   interpolants); "meanSS" draws the latter instead.

   2026-11-01 (the author's request: the axis labels sit too far from the axes): in fActiveMeanAndStd
   the x label is placed by frame_labels.wl, about 4 pt from the tick labels, with an explicit
   ImagePadding that leaves room for it; the frame grows at the same width.
   2026-11-02 (the author's rule: no line thinner than about 1 px; a line meant to look
   thinner is gray): in fActiveMeanAndStd frame
   and ticks print at 0.75 pt (0.19 pt before), the grid line gray at 0.75 pt (frame_labels.wl);
   the curves keep their width (Thickness 0.004 of the frame) as an AbsoluteThickness, so that
   the legend samples match them (1.14 pt for 0.88 pt before).

   Usage:  plotgen/force_fluct/run.sh [N=<Nint>] [meanSS] [cache=<file.mx>] [out=<dir>]
     N=<Nint>        the ssDir run box20_N<Nint>_..._B90000_ss20 (default 200; 180 = old figures)
     meanSS          the mean <f_x> from the numerical steady state instead of the Pade profile
     cache=<file>    read the force arrays from <file> if it exists, else write them there
     out=<dir>       output directory (default: next to this script)
   Outputs: fActiveMeanAndStd.pdf; f2.png (600 dpi), f2_high.png (900 dpi); likewise fxfx and
   fyfy (the notebook's loop exports all three; only f2 is in the manuscript).
   run.sh copies fActiveMeanAndStd.pdf, f2.png and f2_high.png to manuscript/figures/ unless
   out= or N= is given. *)

SetEnvironment["QT_QPA_PLATFORM" -> "offscreen"];
SetEnvironment["GALLIUM_DRIVER" -> "softpipe"];   (* llvmpipe crashes the front end on raster export *)

here = DirectoryName[$InputFileName];
repoRoot = ParentDirectory[ParentDirectory[here]];
args = Rest[$ScriptCommandLine];
opts = Association[Rule @@@ (StringSplit[#, "=", 2] & /@ Select[args, StringContainsQ[#, "="] &])];
meanSS = MemberQ[args, "meanSS"];
nInt = If[KeyExistsQ[opts, "N"], ToExpression[opts["N"]], 200];
outDir = Lookup[opts, "out", here];
cacheFile = Lookup[opts, "cache", None];
If[! DirectoryQ[outDir], CreateDirectory[outDir]];
t0 = AbsoluteTime[];
stamp[] := ToString[Round[AbsoluteTime[] - t0]] <> " s";

Get[FileNameJoin[{repoRoot, "binaryReadTools.m"}]];
Get[FileNameJoin[{repoRoot, "readActiveStress.m"}]];

(* ---------------------------------------------------------------- data *)
tag = "box20_N" <> ToString[nInt] <> "_Dir_lN0.5_z12.5_ld1_B90000_ss20";
dataPath = FileNameJoin[{repoRoot, "data", "fluctuations-set2_2_comove_optimized_2_ssDir",
    "covariance_" <> tag <> ".bin"}];
steady = readSteady[dataPath];
modelParams = steady["modelParams"];
Print["run ", tag, ": steady density wall ", steady["bc"], ", u = ", steady["u"][[1]],
  " um/min;  modelParams ", modelParams];

res = If[StringQ[cacheFile] && FileExistsQ[cacheFile],
   Print["force arrays from the cache ", cacheFile]; Import[cacheFile],
   Module[{r = readActiveStress[dataPath], keep},
     keep = KeyTake[r, {"gridF", "h", "Nint", "fxfxValues", "fyfyValues", "fxfyValues", "f2Values"}];
     If[StringQ[cacheFile], Export[cacheFile, keep]];
     keep]];
Print["readActiveStress: Nint = ", res["Nint"], ", h = ", N@res["h"], " um  (", stamp[],
  ", MaxMemoryUsed ", Round[MaxMemoryUsed[]/2.^30, 0.01], " GiB)"];
(* the interpolants exactly as readActiveStress returns them: covFn = ListInterpolation,
   InterpolationOrder 2, x first, on gridF *)
Do[res[k] = covFn[res[k <> "Values"], res["gridF"]], {k, {"fxfx", "fyfy", "fxfy", "f2"}}];

(* ---------------------------------------------------------------- mean force *)
(* the Pade profile of the static core, In[17] of the notebook; a1, b, L1 stay symbols until
   /. modelParams, as there *)
fPade[z_] := Sqrt[(0.34010999 + 0.04516965 z + 0.00233656 z^2)/(1 + 0.38387160 z +
      0.04750621 z^2 + 0.00233656 z^3)];
S0PadeOverR[r2_] := Sqrt[-a1/(2 b)] Sqrt[-a1/L1] fPade[-r2 a1/L1];
Q0Pade[xv_, yv_] := S0PadeOverR[xv^2 + yv^2] {{xv, yv}, {yv, -xv}};
Q1Pade[xx_, yy_] = Q0Pade[xx, yy][[1, 1]];
Q2Pade[xx_, yy_] = Q0Pade[xx, yy][[1, 2]];
fActive[x_, y_] := {Derivative[1, 0][Q1Pade][x, y] + Derivative[0, 1][Q2Pade][x, y],
   Derivative[1, 0][Q2Pade][x, y] - Derivative[0, 1][Q1Pade][x, y]};    (* In[350] *)

meanPade[x_?NumericQ] := fActive[x, 0.][[1]] /. modelParams;
(* the numerical steady state (comoving, u > 0) the Lyapunov run is linearised about *)
meanNum[x_?NumericQ] := Derivative[1, 0][steady["Q1"]][x, 0.] + Derivative[0, 1][steady["Q2"]][x, 0.];
meanFx = If[meanSS, meanNum, meanPade];

(* ---------------------------------------------------------------- numbers *)
fmt[v_Integer] := ToString[v];
fmt[v_Real] := ToString[NumberForm[Chop[v, 1.*^-6] + 0., {8, 3}]];
fmt[v_List] := "{" <> StringRiffle[fmt /@ v, ", "] <> "}";
fmt[v_] := ToString[v];
r3[v_] := fmt[N[v]];      (* for printing: 3 decimals *)
printTable[head_List, rows_List] := (Print[StringRiffle[head, "  |  "]];
   Scan[Print[StringRiffle[fmt /@ #, "  |  "]] &, rows]);
sx[x_?NumericQ] := Sqrt[res["fxfx"][x, 0.]];
sy[x_?NumericQ] := Sqrt[res["fyfy"][x, 0.]];
xMax = Floor[Max[res["gridF"]], 0.5];
xsPlot = Range[-7., 7., 0.005];
xsBox = Range[-xMax, xMax, 0.01];
mP = meanPade /@ xsPlot; mN = meanNum /@ xsPlot;
Print["\n--- mean <f_x>(x, 0): Pade (notebook) vs numerical steady state, |x| <= 7 ---"];
Print["max |num - Pade| = ", r3@Max[Abs[mN - mP]], " at x = ", r3@xsPlot[[Ordering[Abs[mN - mP], -1][[1]]]],
  ";  at x = 0: Pade ", r3@meanPade[0.], ", num ", r3@meanNum[0.],
  ";  max: Pade ", r3@Max[mP], " at ", r3@xsPlot[[Ordering[mP, -1][[1]]]], ", num ", r3@Max[mN],
  " at ", r3@xsPlot[[Ordering[mN, -1][[1]]]]];
Print["drawn: ", If[meanSS, "numerical steady state", "Pade profile (as the notebook)"]];

Print["\n--- along the axis y = 0 (um^-1) ---"];
Module[{vs, m, s1, s2},
  vs = {-19., -15., -10., -7., -5., -3., -2., -1., -0.5, 0., 0.5, 1., 2., 3., 5., 7., 10., 15., 19.};
  printTable[{"x", "<fx> Pade", "<fx> num", "Std fx", "Std fy", "Std fx/<fx>", "Std fy/<fx>",
     "Std fx/Std fy"},
   Table[{x, meanPade[x], meanNum[x], sx[x], sy[x], sx[x]/meanFx[x], sy[x]/meanFx[x], sx[x]/sy[x]},
    {x, vs}]];
  Do[m = meanFx /@ xs; s1 = sx /@ xs; s2 = sy /@ xs;
    Print[If[Length[xs] == Length[xsPlot], "|x| <= 7:  ", "|x| <= " <> ToString[xMax] <> ": "],
      "min Std fx/<fx> = ", r3@Min[s1/m], " at x = ", r3@xs[[Ordering[s1/m, 1][[1]]]],
      ";  min Std fy/<fx> = ", r3@Min[s2/m], " at x = ", r3@xs[[Ordering[s2/m, 1][[1]]]],
      ";  min (Std fx - <fx>) = ", r3@Min[s1 - m], ", min (Std fy - <fx>) = ", r3@Min[s2 - m]],
    {xs, {xsPlot, xsBox}}]];
localExtrema[f_, xs_] := Module[{v = f /@ xs, d},
   d = Sign[Differences[v]];
   <|"max" -> ({r3@xs[[# + 1]], r3@v[[# + 1]]} & /@ Flatten@Position[Partition[d, 2, 1], {1, -1}]),
     "min" -> ({r3@xs[[# + 1]], r3@v[[# + 1]]} & /@ Flatten@Position[Partition[d, 2, 1], {-1, 1}])|>];
Print["Std fx on |x| <= 7, local extrema {x, value}: ", localExtrema[sx, xsPlot]];
Print["Std fy on |x| <= 7, local extrema {x, value}: ", localExtrema[sy, xsPlot]];
Print["mirror check: max |fxfx(x, 0.6) - fxfx(x, -0.6)| / fxfx = ",
  Max[Abs[res["fxfx"][#, 0.6] - res["fxfx"][#, -0.6]]/res["fxfx"][#, 0.6] & /@ Range[-3., 3., 0.1]]];

Print["\n--- <|delta f|^2> in the plane (um^-2) ---"];
Module[{f2 = res["f2"], c, f2c, th, rs, prof, k, ring, sq, n, gF = res["gridF"]},
  sq = Select[gF, Abs[#] <= 3.9 &];
  Print["range on the plotted square |x|, |y| <= 3.9 (grid nodes): ", r3[MinMax[Outer[f2, sq, sq]]]];
  c = Module[{s = FindMinimum[f2[x, y], {{x, 0.3, -1., 1.}, {y, 0., -1., 1.}}]}, {x, y} /. s[[2]]];
  f2c = f2 @@ c;
  Print["at the core (0, 0): ", r3@f2[0., 0.], " (fxfx ", r3@res["fxfx"][0., 0.], ", fyfy ",
    r3@res["fyfy"][0., 0.], ", fxfx/fyfy ", r3[res["fxfx"][0., 0.]/res["fyfy"][0., 0.]],
    ");  local minimum near the core at ", r3[c], ": ", r3@f2c, " (fxfx/fyfy ",
    r3[(res["fxfx"] @@ c)/(res["fyfy"] @@ c)], ")"];
  (* rays from that minimum: the radius and value of the largest f2 along each, and the share
     of <f_x^2> there *)
  rs = Range[0., 3.6, 0.005];
  ring = Table[
     prof = f2 @@@ (Function[r, c + r {Cos[th], Sin[th]}] /@ rs);
     k = Ordering[prof, -1][[1]];
     With[{p = c + rs[[k]] {Cos[th], Sin[th]}},
      {Round[th/Degree], rs[[k]], p, prof[[k]], (res["fxfx"] @@ p)/prof[[k]],
       (res["fxfx"] @@ p)/(res["fyfy"] @@ p), k < Length[rs]}],
     {th, Range[0., 345., 15.] Degree}];
  printTable[{"angle (deg)", "r_max", "point", "f2 max", "fxfx/f2", "fxfx/fyfy", "interior max"}, ring];
  Print["ring closed (an interior maximum above the centre on every ray): ",
    And @@ ring[[All, -1]] && Min[ring[[All, 4]]] > f2c, ";  f2 max on the ring: ",
    fmt@Max[ring[[All, 4]]], ", min over the rays: ", fmt@Min[ring[[All, 4]]], ", centre ", fmt@f2c];
  Print["fxfx/f2 on the ring: mean over the rays ", fmt@Mean[ring[[All, 5]]], ", range ",
    fmt@MinMax[ring[[All, 5]]], ";  fxfx/fyfy range ", fmt@MinMax[ring[[All, 6]]]];
  (* the global maximum *)
  Module[{s = FindMaximum[f2[x, y], {{x, c[[1]], -3., 3.}, {y, 1., 0.1, 3.}}], p},
   p = {x, y} /. s[[2]];
   Print["global maximum ", fmt@s[[1]], " at ", fmt@p, " (and its mirror image): fxfx ",
     fmt[res["fxfx"] @@ p], ", fyfy ", fmt[res["fyfy"] @@ p], ", fxfx/f2 ", fmt[(res["fxfx"] @@ p)/s[[1]]]]];
  (* far from the core *)
  Print["fxfx/fyfy at radius r (rows) and angle (columns, deg, from +x):"];
  printTable[Prepend[ToString /@ Range[0, 315, 45], "r"],
   Table[Prepend[Table[res["fxfx"][r Cos[a], r Sin[a]]/res["fyfy"][r Cos[a], r Sin[a]],
      {a, Range[0., 315., 45.] Degree}], r], {r, {3., 5., 7., 10., 15.}}]];
  Print["f2 at radius r (rows) and angle (columns):"];
  printTable[Prepend[ToString /@ Range[0, 315, 45], "r"],
   Table[Prepend[Table[f2[r Cos[a], r Sin[a]], {a, Range[0., 315., 45.] Degree}], r],
    {r, {3., 5., 7., 10., 15.}}]]];
Print["(", stamp[], ")"];

(* ---------------------------------------------------------------- figures *)
<< MaTeX`
Get[FileNameJoin[{repoRoot, "plotgen", "frame_labels.wl"}]];
w = 8.6*72/2.54;                                             (* In[1] *)
tex[s_, fs_ : 11] := MaTeX[s, FontSize -> fs,
   "Preamble" -> {"\\usepackage{amsmath,siunitx,upgreek}"}];

(* export in a subkernel (the notebook's RemoteExport for the pdf; here for the pngs too) *)
exportKernel = First@LaunchKernels[1];
subExport[file_String, expr_, o___] := With[{f = file, e = expr, a = {o}},
   ParallelEvaluate[Export[f, e, Sequence @@ a], exportKernel]];
export[file_String, fig_] := Module[{ext = ToLowerCase@FileExtension[file], base, r},  (* In[41] *)
   base = StringDrop[file, -(StringLength[ext] + 1)];
   r = Switch[ext,
     "pdf", {subExport[FileNameJoin[{outDir, file}], fig]},
     "png", {subExport[FileNameJoin[{outDir, base <> "_high.png"}], fig, ImageResolution -> 900],
       subExport[FileNameJoin[{outDir, file}], fig, ImageResolution -> 600]}];
   Print["exported ", r, "  (", stamp[], ")"];
   If[! AllTrue[r, StringQ], Print["*** EXPORT FAILED for ", file]];
   r];

(* In[349]: the three maps; only f2 is in the manuscript *)
With[{W = 3.9},
  Apply[(fig = DensityPlot[res[#1][x, y], {x, -W, W}, {y, -W, W}, PlotPoints -> 50,
       ImageSize -> w,
       PlotRangePadding -> None,
       PlotRange -> All,
       FrameStyle -> Directive[Black, FontFamily -> "Times", FontSize -> 8],
       FrameTicks -> {{Range[-5, 5, 1], None}, {Range[-5, 5, 1], None}},
       FrameLabel -> {tex["x~(\\upmu \\rm{m})"], tex["y~(\\upmu \\rm{m})"]},
       PlotLegends ->
        BarLegend[Automatic, LegendLabel -> tex[StringTemplate["`1`"][#2]],
         LegendMarkerSize -> 200],
       LabelStyle -> {FontSize -> 7, FontFamily -> "Times"}
       ]; export[#1 <> ".png", fig]; fig) &, {
    {"fxfx", "\\text{Var}\\left[ f_x \\right]"},
    {"fyfy", "\\text{Var}\\left[  f_y\\right]"},
    {"f2", "\\,\\langle |\\delta f|^2\\rangle \\,"}
    }, 2]];

(* In[351]: the mean and the standard deviations along the axis.
   2026-11-01, the author's request for one font and one set of sizes in all figures: all text
   in MaTeX Computer Modern (plain Times ticks until then); tick and axis labels print at 9 pt,
   the legend and the formula at 8 pt, at 0.95 \linewidth (246 pt; until then 7.7 pt ticks and
   legend, 10.5 pt label, 8.6 pt formula) *)
{fsT, fsL} = {9, 8} w/(0.95*246.);
fig = With[{box = 7, legFS = fsL}, Plot[{
     meanFx[x],
     Sqrt[res["fxfx"][x, 0]],
     Sqrt[res["fyfy"][x, 0]]}, {x, -box, box},
    PlotStyle -> AbsoluteThickness[0.92],   (* Thickness[0.004] of the 231 pt frame *)
    ImageSize -> w,
    GridLines -> {{0}, None}, GridLinesStyle -> refStyle[w/(0.95*246.)],
    Axes -> False,
    PlotRange -> {All, {0, 3.2}},
    FrameTicks -> {{texTicks[Range[0, 3, 1], tex[#, fsT] &], None},
      {texTicks[Range[-6, 6, 2], tex[#, fsT] &], None}},
    Frame -> True,
    FrameStyle -> frameStyle[w/(0.95*246.)],
    ImagePadding -> {{9, 4}, {25.5, 3}},  (* room for the x label of frameLabels *)
    PlotRangeClipping -> False,           (* the label lies outside; the curves stay inside *)
    PlotLegends -> Placed[LineLegend[Automatic,
       {tex["\\big\\langle f_x(x,0) \\big\\rangle", legFS], tex["\\text{Std}\\left[\\, f_x(x,0) \\,\\right]", legFS],
        tex["\\text{Std}\\left[\\, f_y(x,0) \\,\\right]", legFS]},
       LegendMarkerSize -> {15, 9},
       Spacings -> {0.4, 0.1}], Scaled[{0.2, 0.84}]],
    Epilog -> {Text[tex["f_\\alpha=\\partial_\\beta Q_{\\alpha\\beta}", fsL], Scaled[{0.85, 0.9}]],
      frameLabels[{tex["x~(\\upmu \\rm{m})", fsT], None}, {14.45, 0}]}
    ]];
Print["Std f_x, Std f_y on |x| <= 7 reach ", r3@Max[sx /@ xsPlot], ", ", r3@Max[sy /@ xsPlot],
  " (y range of the figure 0..3.2)"];
export["fActiveMeanAndStd.pdf", fig];

CloseKernels[];
Print["done (", stamp[], ")"];
