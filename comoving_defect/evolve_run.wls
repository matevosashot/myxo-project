#!/usr/bin/env wolframscript
(* Comoving-frame transient: evolve the rho-Q system with the defect core
   pinned at the origin, keeping field snapshots, then polish to the fixed
   point.  Writes evolve_<bc>.mx for evolve_plots.wls.

     evolve_run.wls                     -> L=8, refine=1, both BCs
     evolve_run.wls 8 1 Dirichlet       -> one BC
     evolve_run.wls 8 0.7 both          -> cheaper grid

   Long-running and meant to be watched: prints every accepted step and
   checkpoints every 20 steps, so the .mx is usable before it finishes. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "comoving_core.m"}]];
dir = DirectoryName[$InputFileName];
T0 = AbsoluteTime[];

args  = Drop[$ScriptCommandLine, 1];
argN[k_, d_] := If[Length[args] >= k, ToExpression[args[[k]]], d];
argS[k_, d_] := If[Length[args] >= k, args[[k]], d];
LL    = argN[1, 8.];
REF   = argN[2, 1];
BCARG = argS[3, "both"];
bcList = If[ToLowerCase[BCARG] === "both", {"Dirichlet", "Neumann"}, {BCARG}];

(* ---- parameters: the user's block, verbatim ---- *)
ClearAll[ell, S0v, xi0, zet, Bv, bv, gam, Lam, xir, av, Kv];
psol = First@Solve[{
    ell == Sqrt[(Kv + zet xir/(4 xi0))/(-av - Lam xir)],
    S0v == Sqrt[(-av - Lam xir)/(2 bv)],
    gam == 8 bv/xir} /.
   {ell -> 0.8, S0v -> 0.9, xi0 -> 5., zet -> 10, Bv -> 10000.,
    bv -> 0.1, gam -> 8., Lam -> 0.1}, {xir, av, Kv}];
mk[z_] := <|"B" -> 10000., "a" -> (av /. psol), "b" -> 0.1, "K" -> (Kv /. psol),
            "zeta" -> z, "xi0" -> 5., "xir" -> (xir /. psol),
            "rho0" -> 0.3, "Lambda" -> 0.1, "rhoGauge" -> "Point"|>;
dp = cdDerived[mk[10.]]; dp = Append[dp, "rhoGauge" -> "Point"];
pf = cdPadeFuns[dp];

(* already-converged reference values, for the cross-check at the end *)
uRef = <|"Dirichlet" -> 0.2842160764525491, "Neumann" -> 0.6980352321428196|>;

Print["==============================================================="];
Print["comoving transient    L = ", LL, "   refine = ", REF,
      "   bc = ", bcList];
Print["ld = ", dp["ld"], "   S0 = ", dp["S0"],
      "   B/xi0 = ", dp["Bx"], "   4K'/xir = ", cdFmt@dp["kap"]];
Print["alpha = ", cdFmt[10. (xir /. psol) dp["S0"]/(4 5. dp["Kp"])],
      "   relaxation time xir/(4|a'|) = ",
      cdFmt[(xir /. psol)/(4 Abs[dp["ap"]])], " min"];
Print["==============================================================="];

Do[
 Module[{bcType = bc, gg, op, n, bcs, ne, X0, ev, st, mxf, ckpt, uEnd, gap},
  gg = cdGridPade[LL, dp, 0.05, 1.0, REF];
  op = cdOps[gg, gg, 4]; n = op["n"];
  bcs = cdBC[op, dp, pf, bcType];
  ne  = If[bcType === "Neumann" && cdRhoGauge[dp] =!= "Replace", 3, 2];
  mxf = FileNameJoin[{dir, "evolve_" <> bcType <> ".mx"}];

  Print[""];
  Print["### ", bcType, " ###   n = ", n, " nodes (", op["nx"], "^2)",
        "   hMin = ", cdFmt@gg["hMin"], "   ld/hMin = ",
        cdFmt[dp["ld"]/gg["hMin"]]];

  (* Initial density.  Uniform rho = 1 as requested for Dirichlet, where it
     satisfies the wall.  For the sealed wall uniform rho does NOT satisfy
     n.(B grad rho + zeta div Q) = 0, and at B/xi0 = 2000 that inconsistency
     blows up the first implicit step, so Neumann starts from the
     quasi-static seed instead. *)
  X0 = cdInitial[op, pf];
  If[bcType === "Neumann",
    X0 = cdInitialRho[X0, op, dp, bcs, bcType];
    Print["    rho seed: quasi-static (uniform violates the sealed wall)"],
    Print["    rho seed: uniform 1"]];
  Print["    rho in ", cdFmt[{Min@X0[[1 ;; n]], Max@X0[[1 ;; n]]}]];
  Print[""];

  ckpt = Function[{kk, tt, XX, exx, hh, ss},
    Export[mxf, <|"bc" -> bcType, "L" -> LL, "refine" -> REF, "n" -> n,
                  "grid" -> gg, "modelParams" -> mk[10.], "dp" -> dp,
                  "history" -> hh, "snapshots" -> ss, "final" -> XX,
                  "ex" -> exx, "partial" -> True|>];
    Print["    [checkpoint at k = ", kk, " -> ", FileNameTake[mxf], "]"]];

  (* dtMax = 0.5, not 0.02: implicit Euler is L-stable, the relaxation time
     is only 0.154 min, and at 0.02 the run stalled at t = 5.2 after 400
     steps without ever reaching tolerance. *)
  ev = cdEvolve[X0, ConstantArray[0., ne], op, dp, bcs, bcType,
        "dt0" -> 1.*^-5, "dtMax" -> 0.5, "MaxSteps" -> 400,
        "Tolerance" -> 1.*^-7, "PrintEvery" -> 1,
        "SnapshotEvery" -> 1, "Checkpoint" -> ckpt, "CheckpointEvery" -> 25];

  Print[""];
  Print["    evolve: ", ev["steps"], " steps, t = ", cdFmt@ev["t"],
        ", converged = ", ev["converged"],
        ", ", Round[ev["wall"], 0.1], " s, ",
        Length@ev["snapshots"], " snapshots"];

  st = cdSteady[ev["X"], ev["ex"], op, dp, bcs, bcType,
                "MaxIterations" -> 30, "Verbose" -> False];
  uEnd = st["ex"][[1]];
  gap  = (uEnd - uRef[bcType])/uRef[bcType];
  Print["    steady polish: u = ", NumberForm[uEnd, 10],
        "   |R| = ", cdFmt@st["residual"], "   converged = ", st["converged"]];
  Print["    reference (refine 1, L 8): ", NumberForm[uRef[bcType], 10],
        "   relative gap = ", cdFmt[gap]];
  Print["    u_y = ", cdFmt@st["ex"][[2]],
        "   core = ", cdFmt[{st["X"][[n + op["k0"]]], st["X"][[2 n + op["k0"]]]}]];

  Export[mxf, <|"bc" -> bcType, "L" -> LL, "refine" -> REF, "n" -> n,
                "grid" -> gg, "modelParams" -> mk[10.], "dp" -> dp,
                "history" -> ev["history"], "snapshots" -> ev["snapshots"],
                "final" -> st["X"], "ex" -> st["ex"],
                "uSteady" -> uEnd, "uRef" -> uRef[bcType],
                "evolveConverged" -> ev["converged"],
                "steadyConverged" -> st["converged"], "partial" -> False|>];
  Print["    saved ", FileNameTake[mxf], "  (",
        Round[FileByteCount[mxf]/1.*^6, 0.1], " MB)"];
 ], {bc, bcList}];

Print[""];
Print["total wall ", Round[AbsoluteTime[] - T0, 0.1], " s"];
Print["next: wolframscript -file evolve_plots.wls <bcType>"];
