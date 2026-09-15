#!/usr/bin/env wolframscript
(* Plots from evolve_run.wls output.  No solving.

     evolve_plots.wls Dirichlet [plotBox] [nFrames]
     evolve_plots.wls both

   Writes out_<bc>/ :
     u_vs_t.png        u_x(t), u_y(t), log t
     rho_slice.png     rho(x,0) at log-spaced times
     rho_spacetime.png rho(x,0) as x vs log t
     Q_slice.png       Q1(x,0), Q2(x,0) at the same times
     frames/           defect shape, strength + director
     history.csv                                                        *)

(* Rasterising plots needs the front end, which crashes on a headless box
   ("No matching fbConfigs", "glx: failed to create drisw screen") unless
   DISPLAY is cleared -- the same reason run.sh and the jobscripts all do
   `unset DISPLAY`.  Clear it here so the script works either way. *)
SetEnvironment["DISPLAY" -> ""];

Get[FileNameJoin[{DirectoryName[$InputFileName], "comoving_core.m"}]];
dir = DirectoryName[$InputFileName];
args = Drop[$ScriptCommandLine, 1];
argN[k_, d_] := If[Length[args] >= k, ToExpression[args[[k]]], d];
argS[k_, d_] := If[Length[args] >= k, args[[k]], d];
BCARG  = argS[1, "both"];
PLOTBOX = argN[2, 3.];
NFRAMES = argN[3, 24];
bcList = If[ToLowerCase[BCARG] === "both", {"Dirichlet", "Neumann"}, {BCARG}];

(* the user's own definitions *)
strength[q1_, q2_] := q1^2 + q2^2;                  (* = ddot[W,W]/2 = S^2 *)
director[q1_, q2_] := Module[{th = ArcTan[q1, q2]/2}, {Cos[th], Sin[th]}];

Do[
 Module[{mxf, d, gg, op, n, nx, ny, gx, hist, snaps, S0, ts, pick, cols,
         outdir, framedir, slice, rhoC, q1C, q2C, tmin, tmax, cf, uS, uR,
         legend, k, fr, X, q1f, q2f, dens, vec, expo, rhoPanel, q1Panel},

  mxf = FileNameJoin[{dir, "evolve_" <> bc <> ".mx"}];
  If[!FileExistsQ[mxf],
    Print["missing ", FileNameTake[mxf], " -- run evolve_run.wls first"];
    Continue[]];
  d = Import[mxf];
  gg = d["grid"]; n = d["n"];
  op = cdOps[gg, gg, 4];
  nx = op["nx"]; ny = op["ny"]; gx = op["gx"];
  hist = d["history"]; snaps = d["snapshots"];
  S0 = d["dp"]["S0"];
  uS = Lookup[d, "uSteady", Last[hist]["u"][[1]]];
  uR = Lookup[d, "uRef", None];

  outdir = FileNameJoin[{dir, "out_" <> bc}];
  framedir = FileNameJoin[{outdir, "frames"}];
  Quiet@CreateDirectory[outdir]; Quiet@CreateDirectory[framedir];

  (* Line plots go out as PDF as well as PNG: dean.tex is a5paper and every
     figure it includes is vector.  Density plots and frames stay PNG. *)
  expo[name_, g_] := (Export[FileNameJoin[{outdir, name <> ".png"}], g,
                             ImageResolution -> 120];
                      Export[FileNameJoin[{outdir, name <> ".pdf"}], g]);

  Print["=== ", bc, " ===  n = ", n, " (", nx, "^2),  ",
        Length[hist], " steps,  ", Length[snaps], " snapshots"];
  Print["    t range ", cdFmt@hist[[1]]["t"], " .. ", cdFmt@hist[[-1]]["t"],
        "   u_steady = ", NumberForm[uS, 10],
        If[uR =!= None,
           "   (reference " <> ToString@NumberForm[uR, 8] <> ")", ""]];

  (* ---------------- history csv ---------------- *)
  Export[FileNameJoin[{outdir, "history.csv"}],
    Prepend[{#["t"], #["dt"], #["u"][[1]], #["u"][[2]], #["rate"],
             #["detG"], #["newton"]} & /@ hist,
            {"t", "dt", "ux", "uy", "rate", "detG", "newton"}]];

  (* ---------------- u over time ---------------- *)
  expo["u_vs_t",
   ListLogLinearPlot[
     {{#["t"], #["u"][[1]]} & /@ hist,
      {#["t"], Abs[#["u"][[2]]]} & /@ hist},
     Joined -> True, PlotRange -> All, ImageSize -> 700,
     PlotStyle -> {Directive[Thick, ColorData[97][1]],
                   Directive[Thick, Dashed, ColorData[97][4]]},
     Frame -> True, GridLines -> Automatic,
     FrameLabel -> {"t  [min]", "u  [\[Mu]m/min]"},
     PlotLabel -> bc <> ":  comoving frame velocity, u_steady = " <>
                  ToString@NumberForm[uS, 8],
     Epilog -> {Gray, Dashed, Line[{{Log[hist[[1]]["t"]], uS},
                                    {Log[hist[[-1]]["t"]], uS}}]},
     PlotLegends -> Placed[{"u_x", "|u_y|  (symmetry check)"}, {0.25, 0.82}]]];

  (* ---------------- log-spaced snapshot selection ---------------- *)
  ts = #["t"] & /@ snaps;
  tmin = ts[[1]]; tmax = ts[[-1]];
  pick = DeleteDuplicates@Table[
     First@Ordering[Abs[ts - tt], 1],
     {tt, Exp@Subdivide[Log[tmin], Log[tmax], NFRAMES - 1]}];
  Print["    frames: ", Length[pick], " at t = ",
        Row[Riffle[cdFmt /@ ts[[pick]], ", "]]];

  (* "Rainbow", not "TemperatureMap": the latter is white in the middle, so
     the mid-time curves vanish against the plot background. *)
  tmin = ts[[pick[[1]]]]; tmax = ts[[pick[[-1]]]];
  cf[t_] := ColorData["Rainbow"][
     (Log[t] - Log[tmin])/(Log[tmax] - Log[tmin])];
  cols = Directive[Thick, cf[#]] & /@ ts[[pick]];
  legend = BarLegend[{"Rainbow", {Log10[tmin], Log10[tmax]}},
              LegendLabel -> "log10 t"];

  slice[X_, blk_] := Transpose@{gx,
     Partition[X[[(blk - 1) n + 1 ;; blk n]], nx][[(ny + 1)/2]]};

  (* ---------------- rho(x,0) ---------------- *)
  rhoC = slice[#["X"], 1] & /@ snaps[[pick]];
  rhoPanel = ListLinePlot[rhoC, PlotStyle -> cols,
     PlotRange -> {{-PLOTBOX, PLOTBOX}, All},
     Frame -> True, GridLines -> Automatic, ImageSize -> 700,
     FrameLabel -> {"x  [\[Mu]m]", "\[Rho](x, 0)"},
     PlotLabel -> bc <> ":  density slice, early (blue) to late (red)"];
  expo["rho_slice", Legended[rhoPanel, legend]];

  (* ---------------- rho space-time ---------------- *)
  Export[FileNameJoin[{outdir, "rho_spacetime.png"}],
   ListDensityPlot[
     Flatten[Table[{gx[[i]], Log10[ts[[pick[[j]]]]], rhoC[[j, i, 2]]},
        {j, Length[pick]}, {i, nx}], 1],
     PlotRange -> {{-PLOTBOX, PLOTBOX}, All, All}, ImageSize -> 700,
     ColorFunction -> "TemperatureMap", PlotLegends -> Automatic,
     Frame -> True, FrameLabel -> {"x  [\[Mu]m]", "log10 t"},
     PlotLabel -> bc <> ":  \[Rho](x, 0) relaxation", InterpolationOrder -> 1],
   ImageResolution -> 120];

  (* ---------------- Q1, Q2 slices ---------------- *)
  q1C = slice[#["X"], 2] & /@ snaps[[pick]];
  q2C = slice[#["X"], 3] & /@ snaps[[pick]];
  q1Panel = ListLinePlot[q1C, PlotStyle -> cols,
     PlotRange -> {{-PLOTBOX, PLOTBOX}, All}, Frame -> True,
     GridLines -> Automatic, ImageSize -> 700,
     FrameLabel -> {"x  [\[Mu]m]", "Q1(x, 0)"},
     PlotLabel -> "(b)  nematic slice"];
  (* combined two-panel figure for the paper: dean.tex bakes multi-panel
     figures into one file and labels them (a)/(b) in the caption *)
  expo["slices", Legended[Column[{
        Show[rhoPanel, PlotLabel -> "(a)  density slice"], q1Panel}], legend]];

  Export[FileNameJoin[{outdir, "Q_slice.png"}],
   Legended[Column[{
     ListLinePlot[q1C, PlotStyle -> cols,
       PlotRange -> {{-PLOTBOX, PLOTBOX}, All}, Frame -> True,
       GridLines -> Automatic, ImageSize -> 700,
       FrameLabel -> {None, "Q1(x, 0)"},
       PlotLabel -> bc <> ":  nematic slices, early (blue) to late (red)"],
     (* Q2 is ODD in y, so Q2(x,0) is identically zero in exact arithmetic.
        This panel therefore shows the numerical symmetry breaking, not
        physics -- it tracks u_y and shrinks with refinement. *)
     ListLinePlot[q2C, PlotStyle -> cols,
       PlotRange -> {{-PLOTBOX, PLOTBOX}, All}, Frame -> True,
       GridLines -> Automatic, ImageSize -> 700,
       FrameLabel -> {"x  [\[Mu]m]", "Q2(x, 0)"},
       PlotLabel -> Style[
         "Q2(x,0) = 0 by symmetry: this panel is the numerical error, "
         <> "max |Q2(x,0)| = " <> cdFmt[Max@Abs[q2C[[All, All, 2]]]],
         Italic, Gray]]}], legend],
   ImageResolution -> 120];

  (* ---------------- defect-shape frames ---------------- *)
  (* PlotRange on strength is FIXED at {0, S0^2} across all frames, so the
     colour scale does not rescale frame to frame and the sequence does not
     appear to pulse. *)
  Do[
   k = pick[[fr]]; X = snaps[[k]]["X"];
   q1f = cdInterp[X, "Q1", op]; q2f = cdInterp[X, "Q2", op];
   dens = DensityPlot[strength[q1f[x, y], q2f[x, y]],
      {x, -PLOTBOX, PLOTBOX}, {y, -PLOTBOX, PLOTBOX},
      PlotRange -> {Full, Full, {0, S0^2}}, ColorFunction -> "TemperatureMap",
      PlotPoints -> 40, MaxRecursion -> 1, PlotLegends -> Automatic,
      Frame -> True, FrameLabel -> {"x  [\[Mu]m]", "y  [\[Mu]m]"},
      PlotLabel -> Row[{bc, "   t = ", NumberForm[snaps[[k]]["t"], 4],
                        " min   u = ", NumberForm[snaps[[k]]["u"][[1]], 6]}]];
   vec = VectorPlot[director[q1f[x, y], q2f[x, y]],
      {x, -PLOTBOX, PLOTBOX}, {y, -PLOTBOX, PLOTBOX},
      (* Arrowheads[0]: the nematic director is headless, it has no sign.
         VectorColorFunction -> None keeps them black rather than letting
         VectorPlot colour by magnitude (which is constant here anyway). *)
      VectorStyle -> {Black, Arrowheads[0]}, VectorPoints -> 16,
      VectorColorFunction -> None,
      VectorScaling -> None, VectorSizes -> {0.35, 0.35}];
   Export[FileNameJoin[{framedir,
      "frame_" <> IntegerString[fr, 10, 3] <> ".png"}],
     Show[dens, vec, ImageSize -> 600], ImageResolution -> 110];
   (* the last frame is the converged state -- also save it under a stable
      name, since that is the one the paper figure uses *)
   If[fr == Length[pick],
     Export[FileNameJoin[{outdir, "defect_steady.png"}],
       Show[dens, vec, ImageSize -> 600], ImageResolution -> 150]],
   {fr, Length[pick]}];

  Print["    wrote ", Length[pick], " frames to out_", bc, "/frames/"];

  (* ---------------- checks ---------------- *)
  Print["    core (Q1,Q2)(0,0) over snapshots, max |.| = ",
        cdFmt@Max@Abs@Flatten[
          {#["X"][[n + op["k0"]]], #["X"][[2 n + op["k0"]]]} & /@ snaps]];
  Print["    max |u_y| over run = ",
        cdFmt@Max@Abs[#["u"][[2]] & /@ hist]];
  If[bc === "Neumann",
    Print["    mass (after cdShiftRhoMean) spread over snapshots = ",
          cdFmt@Max@Abs[(op["w"] . cdShiftRhoMean[#["X"], op][[1 ;; n]]/
                         op["area"]) - 1. & /@ snaps]]];
 ], {bc, bcList}];

Print["done"];
