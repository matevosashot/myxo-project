#!/usr/bin/env wolframscript
(* ::Package:: *)

(* export_fields.wls -- derived fields of the full-model comoving solutions for
   the report.

   Reads ../comoving_defect_full/out/<tag>.mx (cfRun of run_full.wls), rebuilds
   the solver's operators from the saved grid (cdOps + cdMirrorOps, order 4) and
   writes, for every milestone of every run, data/<tag>_m<k>.json with
     - the milestone record: weights, u, diagnostics;
     - profiles along y = 0 and x = 0 on the solver nodes: rho, (rho-1)B/zeta,
       Q1, Q2, S, H1, H2 and the x-components of the three force densities of
       the density flux, xi0 J = -B grad rho + Q:grad H + f_a (eq:velocity):
         pressure  -B d_x rho,  elastic  2 (Q1 d_x H1 + Q2 d_x H2),
         active    -zeta (d_x Q1 + d_y Q2);
     - the director deviation theta - phi/2 on circles r = 0.5, 1, 2, 5 um,
       theta = Arg(Q1 + i Q2)/2 from the cubic interpolants;
     - maps on a uniform 161^2 grid of [-4, 4]^2: S, |H| = Sqrt[H1^2 + H2^2],
       (rho-1)B/zeta, theta;
     - the director field on a uniform 81^2 grid of [-10, 10]^2 (for director
       lines).
   plus data/<tag>_runs.json with the u trace of every converged continuation
   point and the run metadata.

   Usage:  wolframscript -file export_fields.wls <tag> [<tag> ...] *)

$HistoryLength = 0;
thisDir = DirectoryName[$InputFileName];
repoRoot = ParentDirectory[thisDir];
Get[FileNameJoin[{repoRoot, "comoving_defect_full", "comoving_full_core.m"}]];

tags = Rest[$ScriptCommandLine];
If[tags === {}, Print["usage: export_fields.wls tag [tag ...]"]; Exit[1]];
dataDir = FileNameJoin[{thisDir, "data"}];
Quiet@CreateDirectory[dataDir];

num[x_] := N[x];
wrapHalfPi[a_] := Mod[a + Pi/2, Pi] - Pi/2;

Do[
  Module[{file, run, mpR, gg, op, n, nx, ny, i0, j0, B0, zeta, rsc, ms, rowY0, colX0, sig0,
          axisData, circles, maps, dirField, k, X, r, q1, q2, h1, h2, S, fx, el, pr,
          q1f, q2f, h1f, h2f, rf, uni, uni2, out},
    file = FileNameJoin[{repoRoot, "comoving_defect_full", "out", tag <> ".mx"}];
    If[!FileExistsQ[file], Print["missing ", file]; Continue[]];
    Clear[cfRun]; Get[file]; run = cfRun;
    mpR = Association[run["modelParams"]];
    B0 = mpR[B]; zeta = mpR[\[Zeta]];
    (* rhoScaled = (rho - 1) B/zeta; at zeta = 0 it is the pressure B(rho - 1) instead *)
    rsc = If[zeta != 0, B0/zeta, B0];
    gg = run["mesh"]["grid"];
    op = cfOps[cdMirrorOps[cdOps[gg, gg, 4]]];
    n = op["n"]; nx = op["nx"]; ny = op["ny"];
    i0 = (nx + 1)/2; j0 = (ny + 1)/2;
    rowY0 = Table[i + (j0 - 1) nx, {i, nx}];       (* nodes of y = 0 *)
    colX0 = Table[i0 + (j - 1) nx, {j, ny}];       (* nodes of x = 0 *)
    uni  = Subdivide[-4., 4., 160];
    uni2 = Subdivide[-10., 10., 80];
    ms = run["milestones"];
    Print[tag, ": n = ", n, ", ", Length[ms], " milestones"];

    (* state layout: Join[rho - 1, ...] since 2026-10-12; older runs carry rho *)
    sig0 = If[KeyExistsQ[run, "stateLayout"], 0., -1.];
    Do[
      With[{m = ms[[k]]},
        X = m["X"];
        r  = sig0 + X[[1 ;; n]];  (* r = rho - 1 *)
        q1 = X[[n + 1 ;; 2 n]];  q2 = X[[2 n + 1 ;; 3 n]];
        h1 = X[[3 n + 1 ;; 4 n]]; h2 = X[[4 n + 1 ;; 5 n]];
        S  = Sqrt[q1^2 + q2^2];
        fx = -zeta ((op["Dx"] . q1) + (op["Dy"] . q2));
        el = 2 (q1 (op["Dx"] . h1) + q2 (op["Dx"] . h2));
        pr = -B0 (op["Dx"] . r);
        axisData[v_, idx_] := num[v[[idx]]];
        {rf, q1f, q2f, h1f, h2f} = Table[cfInterp[X, f, op, If[f == 1, sig0, 0.]], {f, 5}];
        circles = Table[
          Module[{phis = Subdivide[-Pi, Pi, 360], xs, ys, th},
            xs = rad Cos[phis]; ys = rad Sin[phis];
            th = Arg[q1f[xs, ys] + I q2f[xs, ys]]/2;
            <|"r" -> rad, "phi" -> num[phis],
              "dev" -> num[wrapHalfPi[th - phis/2]]|>],
          {rad, {0.5, 1., 2., 5.}}];
        maps = Module[{gx = Flatten@Table[x, {y, uni}, {x, uni}],
                       gy = Flatten@Table[y, {y, uni}, {x, uni}], a1, a2, b1, b2, rr},
          a1 = q1f[gx, gy]; a2 = q2f[gx, gy]; b1 = h1f[gx, gy]; b2 = h2f[gx, gy];
          rr = rf[gx, gy];
          <|"x" -> num[uni], "y" -> num[uni],
            "S" -> Partition[num[Sqrt[a1^2 + a2^2]], 161],
            "Habs" -> Partition[num[Sqrt[b1^2 + b2^2]], 161],
            "rhoScaled" -> Partition[num[rr rsc], 161],
            "theta" -> Partition[num[Arg[a1 + I a2]/2], 161]|>];
        dirField = Module[{gx = Flatten@Table[x, {y, uni2}, {x, uni2}],
                           gy = Flatten@Table[y, {y, uni2}, {x, uni2}], a1, a2},
          a1 = q1f[gx, gy]; a2 = q2f[gx, gy];
          <|"x" -> num[uni2], "y" -> num[uni2],
            "Q1" -> Partition[num[a1], 81], "Q2" -> Partition[num[a2], 81]|>];
        out = <|
          "tag" -> tag, "milestone" -> k, "w" -> num[m["w"]], "ok" -> m["ok"],
          "u" -> num[m["u"]],
          "diag" -> <|"residual" -> num[m["diag"]["residual"]],
                      "lineForce" -> num[m["diag"]["lineForce"]],
                      "wallS" -> num /@ m["diag"]["wallS"],
                      "SInfUniform" -> num[m["diag"]["SInfUniform"]],
                      "rhoRange" -> num[m["diag"]["rhoRange"]],
                      "symmetry" -> num /@ m["diag"]["symmetry"]|>,
          "B" -> B0, "zeta" -> zeta, "box" -> num[run["mesh"]["box"]],
          "nx" -> nx, "hMin" -> num[run["mesh"]["hMin"]],
          "axisX" -> <|"x" -> num[op["gx"]],
                       "rho" -> num[1 + r[[rowY0]]],
                       "rhoScaled" -> num[r[[rowY0]] rsc],
                       "Q1" -> axisData[q1, rowY0], "Q2" -> axisData[q2, rowY0],
                       "S" -> axisData[S, rowY0],
                       "H1" -> axisData[h1, rowY0], "H2" -> axisData[h2, rowY0],
                       "fActive" -> axisData[fx, rowY0],
                       "fElastic" -> axisData[el, rowY0],
                       "fPressure" -> axisData[pr, rowY0]|>,
          "axisY" -> <|"y" -> num[op["gy"]],
                       "S" -> axisData[S, colX0],
                       "Q1" -> axisData[q1, colX0], "Q2" -> axisData[q2, colX0],
                       "rhoScaled" -> num[r[[colX0]] rsc]|>,
          "circles" -> circles, "maps" -> maps, "director" -> dirField|>;
        Export[FileNameJoin[{dataDir, tag <> "_m" <> ToString[k] <> ".json"}], out, "RawJSON"];
        Print["   milestone ", k, " w = ", m["w"], " u = ", m["u"][[1]], " -> ",
              tag <> "_m" <> ToString[k] <> ".json"]],
      {k, Length[ms]}];

    Export[FileNameJoin[{dataDir, tag <> "_runs.json"}],
      <|"tag" -> tag, "args" -> run["args"], "time" -> num[run["time"]],
        "host" -> run["host"], "date" -> run["date"], "g2Mode" -> run["g2Mode"],
        "nx" -> nx, "hMin" -> num[run["mesh"]["hMin"]], "box" -> num[run["mesh"]["box"]],
        "zeta" -> zeta, "B" -> B0,
        "milestones" -> (<|"w" -> num[#["w"]], "u" -> num[#["u"]], "ok" -> #["ok"],
                           "residual" -> num[#["diag"]["residual"]],
                           "wallS" -> num /@ #["diag"]["wallS"],
                           "SInfUniform" -> num[#["diag"]["SInfUniform"]],
                           "rhoRange" -> num[#["diag"]["rhoRange"]],
                           "lineForce" -> num[#["diag"]["lineForce"]],
                           "newton" -> Lookup[#, "newton", 0],
                           "halvings" -> Lookup[#, "halvings", 0],
                           "wall" -> num[#["wall"]]|> & /@ ms),
        "trace" -> (<|"p" -> num[#["p"]], "u" -> num[#["u"]],
                      "iterations" -> #["iterations"]|> & /@ run["trace"])|>,
      "RawJSON"];
    Print["   ", tag, "_runs.json"]],
  {tag, tags}];
