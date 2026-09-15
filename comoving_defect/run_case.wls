#!/usr/bin/env wolframscript
(* One case per invocation.  Appends a line to results.csv so a kill never
   loses progress.
     run_case.wls jac
     run_case.wls passive  <L> <refine>
     run_case.wls active   <L> <refine> <bcType>
     run_case.wls drift    <L> <refine> <bcType>
     run_case.wls sinh     <L> <n>      <bcType>    (grading cross-check)  *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "comoving_core.m"}]];
dir = DirectoryName[$InputFileName];
csv = FileNameJoin[{dir, "results.csv"}];
t00 = AbsoluteTime[];

args = Drop[$ScriptCommandLine, 1];
what = args[[1]];
argN[k_, d_] := If[Length[args] >= k, ToExpression[args[[k]]], d];
argS[k_, d_] := If[Length[args] >= k, args[[k]], d];

(* ---------------- parameters (user's block, verbatim) ---------------- *)
ClearAll[ell, S0v, xi0, zet, Bv, bv, gam, Lam, xir, av, Kv];
psol = First@Solve[{
    ell == Sqrt[(Kv + zet xir/(4 xi0))/(-av - Lam xir)],
    S0v == Sqrt[(-av - Lam xir)/(2 bv)],
    gam == 8 bv/xir} /.
   {ell -> 0.8, S0v -> 0.9, xi0 -> 5., zet -> 10, Bv -> 10000.,
    bv -> 0.1, gam -> 8., Lam -> 0.1}, {xir, av, Kv}];
mk[z_] := <|"B" -> 10000., "a" -> (av /. psol), "b" -> 0.1, "K" -> (Kv /. psol),
            "zeta" -> z, "xi0" -> 5., "xir" -> (xir /. psol),
            "rho0" -> 0.3, "Lambda" -> 0.1|>;

DMESH = 0.05; HMAX = 1.0;
(* 5th arg selects the sealed-wall density gauge: "Point" (default, sparse)
   or "Mean" (dense row, physical normalisation).  Both give the same u. *)
GAUGE = argS[5, "Point"];
dpOf[z_] := Append[cdDerived[mk[z]], "rhoGauge" -> GAUGE];

emit[assoc_] := Module[{line},
  line = StringRiffle[
    ToString[#, InputForm] & /@ {DateString["ISO8601"], what,
      Lookup[assoc, "L", ""], Lookup[assoc, "refine", ""],
      Lookup[assoc, "bc", ""], Lookup[assoc, "n", ""],
      Lookup[assoc, "hMin", ""], Lookup[assoc, "u", ""],
      Lookup[assoc, "uy", ""], Lookup[assoc, "ok", ""],
      Lookup[assoc, "R", ""], Lookup[assoc, "extra", ""],
      Round[AbsoluteTime[] - t00, 0.1]}, ","];
  If[!FileExistsQ[csv],
    Export[csv, "time,what,L,refine,bc,n,hMin,u,uy,ok,R,extra,wall\n", "Text"]];
  WriteString[csv, line <> "\n"];
  Print["CSV: ", line]];

(* ---------------- shared solve ---------------- *)
build[LL_, refine_, bcType_, dpp_, pff_] := Module[{gg, op, bc},
  gg = cdGridPade[LL, dpp, DMESH, HMAX, refine];
  op = cdOps[gg, gg, 4];
  bc = cdBC[op, dpp, pff, bcType];
  {gg, op, bc}];

(* Continuation in zeta.  At alpha = 0.43 a cold Newton from the Pade seed
   does NOT converge (measured: 40 wasted iterations + a 30-step evolve +
   more, ~200 linear solves, 57 s at n = 75).  Ramping zeta from 0 keeps
   every Newton a small perturbation of the previous solution: ~3
   iterations per step, ~20 solves total.  The grid is built once at the
   TARGET zeta and held fixed; only the equations and the Pade far field
   move with zeta. *)
solveCase[LL_, refine_, bcType_, dpTarget_, pfTarget_, gridSpec_: Automatic,
          zTarget_: 10., zSteps_: 5] :=
 Module[{gg, op, ne, X, ex, st = Null, zs, dpz, pfz, bcz, bcOut, k},
  gg = If[gridSpec === Automatic,
          cdGridPade[LL, dpTarget, DMESH, HMAX, refine], gridSpec];
  op = cdOps[gg, gg, 4];
  ne = If[bcType === "Neumann" && GAUGE =!= "Replace", 3, 2];
  zs = If[zTarget == 0., {0.}, Rest@Subdivide[0., zTarget, zSteps]];
  zs = Prepend[zs, 0.];

  (* seed at zeta = 0 *)
  dpz = dpOf[0.]; pfz = cdPadeFuns[dpz];
  bcz = cdBC[op, dpz, pfz, bcType];
  X  = cdInitialRho[cdInitial[op, pfz], op, dpz, bcz, bcType];
  ex = ConstantArray[0., ne];

  Do[
   dpz = dpOf[zs[[k]]]; pfz = cdPadeFuns[dpz];
   bcz = cdBC[op, dpz, pfz, bcType];
   st  = cdSteady[X, ex, op, dpz, bcz, bcType,
                  "MaxIterations" -> 25, "Verbose" -> False];
   If[!TrueQ@st["converged"],
     Print["   continuation stalled at zeta = ", zs[[k]],
           "  |R| = ", cdFmt@st["residual"]]];
   X = st["X"]; ex = st["ex"]; bcOut = bcz,
   {k, Length[zs]}];

  <|"gg" -> gg, "op" -> op, "bc" -> bcOut, "st" -> st, "ev" -> Null,
    "u" -> ex[[1]], "uy" -> ex[[2]], "ok" -> st["converged"],
    "R" -> st["residual"], "n" -> op["n"], "zs" -> zs|>];

(* ==================================================================== *)
Switch[what,

 "jac",
  Module[{dp, pf, gg, op, n, bc, X, ex, Xold, Rf, Jan, Jnum, v, eps = 1.*^-6,
          e, errs = {}},
   dp = dpOf[10.]; pf = cdPadeFuns[dp];
   Do[
    gg = cdGrid[2., 11, 1.0]; op = cdOps[gg, gg, 4]; n = op["n"];
    bc = cdBC[op, dp, pf, bcType]; e = If[bcType === "Neumann" && GAUGE =!= "Replace", 3, 2];
    SeedRandom[11]; X = cdInitial[op, pf] + 0.05 RandomReal[{-1, 1}, 3 n];
    ex = RandomReal[{-1, 1}, e]; Xold = X + 0.1 RandomReal[{-1, 1}, 3 n];
    (* Compare the UNSCALED Jacobian against finite differences of the
       UNSCALED residual.  cdRJ row-equilibrates, and its scale factors
       depend on the state, so differencing the scaled residual would be
       differencing a different function at every perturbation. *)
    Rf[xx_, ee_] := cdRes[xx, ee, Xold, 0.001, 1, op, dp, bc, bcType];
    Jan = Module[{q = cdRJ[X, ex, Xold, 0.001, 1, op, dp, bc, bcType]},
            q[[3]] q[[2]]];     (* undo the row scaling *)
    Jnum = Transpose@Join[
      Table[v = ConstantArray[0., 3 n]; v[[k]] = eps;
        (Rf[X + v, ex] - Rf[X - v, ex])/(2 eps), {k, 3 n}],
      Table[v = ConstantArray[0., e]; v[[k]] = eps;
        (Rf[X, ex + v] - Rf[X, ex - v])/(2 eps), {k, e}]];
    AppendTo[errs, Max@Abs[Normal[Jan] - Jnum]/Max@Abs@Normal[Jan]];
    Print["   ", bcType, "  rel err = ", cdFmt@Last[errs]],
    {bcType, {"Dirichlet", "Neumann"}}];
   emit[<|"extra" -> "jacErr=" <> cdFmt[Max[errs]],
          "ok" -> (Max[errs] < 1.*^-7)|>]],

 "passive",
  Module[{dp0, pf0, LL, rf, r, op, n, S, line, xp, err, ref},
   dp0 = dpOf[0.]; pf0 = cdPadeFuns[dp0];
   LL = argN[2, 6.]; rf = argN[3, 1];
   r = solveCase[LL, rf, "Dirichlet", dp0, pf0, Automatic, 0.];
   op = r["op"]; n = r["n"];
   ref = cdRadialBVP[dp0["ld"], dp0["S0"], 100., 20000];
   S = Sqrt[r["st"]["X"][[n + 1 ;; 2 n]]^2 + r["st"]["X"][[2 n + 1 ;; 3 n]]^2];
   line = Partition[S, op["nx"]][[(op["ny"] + 1)/2]];
   xp = Select[Range[op["nx"]], 0.08 < op["gx"][[#]] <= 2.5 &];
   err = Max@Abs[line[[xp]] - (ref["f"] /@ op["gx"][[xp]])]/dp0["S0"];
   Print["   n = ", n, "  hMin = ", cdFmt@r["gg"]["hMin"],
         "  |u| = ", cdFmt@Norm[{r["u"], r["uy"]}],
         "  profileErr = ", cdFmt@err,
         "  max|rho-1| = ", cdFmt@Max@Abs[r["st"]["X"][[1 ;; n]] - 1.]];
   emit[<|"L" -> LL, "refine" -> rf, "bc" -> "Dirichlet", "n" -> n,
          "hMin" -> r["gg"]["hMin"], "u" -> Norm[{r["u"], r["uy"]}],
          "uy" -> r["uy"], "ok" -> r["ok"], "R" -> r["R"],
          "extra" -> "err=" <> cdFmt[err]|>]],

 "active" | "sinh",
  Module[{dp, pf, LL, rf, bcType, r, op, n, est, sym, gspec, X},
   dp = dpOf[10.]; pf = cdPadeFuns[dp];
   LL = argN[2, 8.]; rf = argN[3, 1]; bcType = argS[4, "Dirichlet"];
   gspec = If[what === "sinh", cdGrid[LL, rf, 2.5], Automatic];
   r = solveCase[LL, rf, bcType, dp, pf, gspec];
   op = r["op"]; n = r["n"]; X = r["st"]["X"];
   est = cdCoreEstimates[X, op, dp]; sym = cdSymmetry[X, op];
   Print["   ", bcType, "  L = ", LL, "  n = ", n,
         "  hMin = ", cdFmt@r["gg"]["hMin"],
         "  ld/hMin = ", cdFmt[dp["ld"]/r["gg"]["hMin"]]];
   Print["   u  = ", NumberForm[r["u"], 10], "   uy = ", cdFmt@r["uy"],
         "   ok = ", r["ok"], "   |R| = ", cdFmt@r["R"]];
   Print["   naive = ", cdFmt@est["naive"][[1]],
         "   naive/u = ", cdFmt[est["naive"][[1]]/r["u"]],
         "   coreExact = ", cdFmt@est["exact"][[1]],
         "   detG = ", cdFmt@est["detG"]];
   Print["   symmetry = ", cdFmt /@ Values@sym,
         "   core = ", cdFmt[{X[[n + op["k0"]]], X[[2 n + op["k0"]]]}]];
   Print["   rho in ", cdFmt[{Min@X[[1 ;; n]], Max@X[[1 ;; n]]}],
         If[bcType === "Neumann",
            "   gauge = " <> GAUGE <>
            "   mu = " <> If[Length[r["st"]["ex"]] >= 3,
                             cdFmt[r["st"]["ex"][[3]]], "n/a"] <>
            "   mass = " <> cdFmt[op["w"] . X[[1 ;; n]]/op["area"]] <>
            "   mass after shift = " <>
            cdFmt[op["w"] . cdShiftRhoMean[X, op][[1 ;; n]]/op["area"]], ""]];
   If[r["ev"] =!= Null, Print["   (evolve fallback: ", r["ev"]["steps"], " steps)"]];
   emit[<|"L" -> LL, "refine" -> rf, "bc" -> bcType, "n" -> n,
          "hMin" -> r["gg"]["hMin"], "u" -> r["u"], "uy" -> r["uy"],
          "ok" -> r["ok"], "R" -> r["R"],
          "extra" -> "naiveRatio=" <> cdFmt[est["naive"][[1]]/r["u"]] <>
                     ";sym=" <> cdFmt[Max@Values@sym]|>]],

 "drift",
  Module[{dp, pf, LL, rf, bcType, r, op, n, ne, bc, X, dt, traj, corePos,
          f1, f2, slope, nStep = 8},
   dp = dpOf[10.]; pf = cdPadeFuns[dp];
   LL = argN[2, 8.]; rf = argN[3, 1]; bcType = argS[4, "Dirichlet"];
   r = solveCase[LL, rf, bcType, dp, pf];
   op = r["op"]; n = r["n"]; bc = r["bc"]; ne = If[bcType === "Neumann" && GAUGE =!= "Replace", 3, 2];
   corePos[XX_] := Module[{a, b, s},
     a = cdInterp[XX, "Q1", op]; b = cdInterp[XX, "Q2", op];
     s = Quiet@FindRoot[{a[xx, yy] == 0, b[xx, yy] == 0},
                        {{xx, 0.}, {yy, 0.}}, AccuracyGoal -> 12];
     {xx, yy} /. s];
   X = r["st"]["X"]; dt = 0.002; traj = {{0., corePos[X]}};
   Do[Module[{Xi = X, R, M, rsv, keep, d, j = 0},
      While[j < 8,
       {R, M, rsv} = cdRJ[Xi, ConstantArray[0., ne], X, dt, 1, op, dp, bc, bcType];
       keep = If[ne >= 3, Append[Range[3 n], 3 n + 3], Range[3 n]];
       d = LinearSolve[M[[keep, keep]], -R[[keep]]];
       Xi += d[[1 ;; 3 n]];
       If[Max@Abs@d < 1.*^-12 Max[1., Max@Abs@Xi], j = 99, j++]];
      X = Xi];
     AppendTo[traj, {Last[traj][[1]] + dt, corePos[X]}], {nStep}];
   slope = Coefficient[Fit[{#[[1]], #[[2, 1]]} & /@ traj, {1, tt}, tt], tt];
   Print["   lab drift dx/dt = ", NumberForm[slope, 10]];
   Print["   bordered u      = ", NumberForm[r["u"], 10]];
   Print["   rel diff        = ", cdFmt[(slope - r["u"])/r["u"]]];
   Print["   y drift         = ", cdFmt[Last[traj][[2, 2]]]];
   emit[<|"L" -> LL, "refine" -> rf, "bc" -> bcType, "n" -> n,
          "u" -> r["u"], "ok" -> r["ok"],
          "extra" -> "drift=" <> ToString[slope, InputForm] <>
                     ";rel=" <> cdFmt[(slope - r["u"])/r["u"]]|>]],

 _, Print["unknown: ", what]
];

Print["[", Round[AbsoluteTime[] - t00, 0.1], "s] done"];
