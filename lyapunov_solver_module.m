(* ::Package:: *)

(* ::Section:: *)
(*Gaussian fluctuations around the \[Rho]\[Dash]Q steady state: discrete Lyapunov solve.*)

(* psi = (d\[Rho], dQ1, dQ2)^T on a uniform interior grid with psi|_bdry = 0;
   psi_t = A psi + noise, Sigma solves A Sigma + Sigma A^T + Q = 0.
   Noise uses the EXACT max-entropy moments; the leading-order closure is not
   PSD above S = 0.618.

   Layout: the public functions get their usage messages here, in Global`, and
   everything else lives in the ActiveNematic`Private` context at the bottom --
   so the helpers (kappaOfS, localFluxCovariance, ...) need no name prefix and
   are invisible outside.  iterative_solver_module.m opens the same private
   context and owns the steady state: SolveActiveNematicSteady,
   visualizeSteadyState, and the shared prepareDir/savePlot. *)

Needs["NDSolve`FEM`"];
SetSystemOptions["ParallelOptions" -> "ParallelThreadNumber" -> $ProcessorCount];

If[Names["Global`SolveActiveNematicSteady"] === {},
   Get[FileNameJoin[{DirectoryName[$InputFileName], "iterative_solver_module.m"}]]];

Clear[SolveFluctuationsLyapunov, visualizeFluctuations,
      exportData, MergeParams, ScriptParamOverrides];


SolveFluctuationsLyapunov::usage =
  "SolveFluctuationsLyapunov[lyapunovParams, modelParams, steadyState] solves " <>
  "A \[CapitalSigma] + \[CapitalSigma] A^T + Q = 0 for the equal-time covariance of " <>
  "(\[Delta]\[Rho], \[Delta]Q1, \[Delta]Q2) linearized about steadyState, on a uniform " <>
  "Nint x Nint grid over [-fbox, fbox]^2 with Dirichlet (zero) boundary values.\n\n" <>
  "  lyapunovParams : rules {Nint->..., fbox->...}; both positive integers, " <>
  "h = 2 fbox/(Nint+1) is kept exact so odd Nint puts the defect core on a node\n" <>
  "  modelParams    : same rules as SolveActiveNematicSteady\n" <>
  "  steadyState    : <|\"rho\"->..., \"Q1\"->..., \"Q2\"->..., \"mesh\"->...|>\n\n" <>
  "Returns <|\"C\", \"sigmaRho\", \"sigmaRhoFn\", \"sigmaQ1\", \"sigmaQ2\", " <>
  "\"gridInt\", \"h\", \"Nint\", \"fbox\", \"eigenvalues\", \"maxReLambda\", " <>
  "\"residual\", \"residualRefined\", \"modelParams\", \"lyapunovParams\", " <>
  "\"timings\"|>.";

SolveFluctuationsLyapunov::badgrid =
  "Nint and fbox must be positive integers (h = 2 fbox/(Nint+1) is kept exact); got `1`.";
SolveFluctuationsLyapunov::badss =
  "steadyState must be an association with keys \"rho\", \"Q1\", \"Q2\"; got `1`.";
SolveFluctuationsLyapunov::unstable =
  "Linearization is not strictly stable: max Re(\[Lambda]) = `1`; results may be unphysical.";
SolveFluctuationsLyapunov::notpsd =
  "Condition (ii) Q >= 0 is violated even with the exact moments at `1` grid points; \
check the closure assembly.";
SolveFluctuationsLyapunov::norefine =
  "Refinement did not improve the residual (`1` -> `2`); consider LyapunovSolve.";
SolveFluctuationsLyapunov::kappa =
  "\[Kappa](S) root-find did not converge at some grid points (max |I1/I0 - S| = `1`).";

visualizeFluctuations::usage =
  "visualizeFluctuations[correlations, outputDir, plotBox] returns {plot2D, plot3D} " <>
  "of Var(\[Delta]\[Rho])(r) over [-plotBox, plotBox]^2 and writes them as PNGs into " <>
  "outputDir. outputDir = None plots without saving.";

exportData::usage =
  "exportData[correlations, outputDir, saveCovariance] writes sigma_rho_<tag>.m " <>
  "(grid + Var(\[Delta]\[Rho]) + spectrum) into outputDir and, when saveCovariance is " <>
  "True, the full covariance as row-major Real64 covariance_<tag>.bin plus a meta " <>
  "file carrying the layout and the timings (lyapunovTime, steadyStateTime, " <>
  "totalTime, and the per-stage breakdown; steadyStateTime is None unless the " <>
  "caller added it to correlations). " <>
  "Returns the list of written paths; outputDir = None writes nothing.";

MergeParams::usage =
  "MergeParams[defaults, overrides] returns defaults with the keys given in " <>
  "overrides replaced, keeping every key not mentioned. " <>
  "MergeParams[defaults, overrides, label] aborts on a key absent from defaults.";

ScriptParamOverrides::usage =
  "ScriptParamOverrides[str] parses \"name1={k->v,...}; name2={...}\" into " <>
  "<|\"name1\"->{k->v,...}, ...|> without evaluating the left-hand sides. " <>
  "ScriptParamOverrides[] reads $ScriptCommandLine[[2]], or returns <||> when no " <>
  "argument was given.";


Begin["ActiveNematic`Private`"];

(* ::Subsection:: *)
(*Parameter merging*)

MergeParams[defaults : {___Rule}, overrides : {___Rule}] :=
  Normal @ Association @ Join[defaults, overrides];

MergeParams[defaults : {___Rule}, overrides : {___Rule}, label_String] :=
  Module[{bad = Complement[First /@ overrides, First /@ defaults]},
    If[bad =!= {},
      Print["ERROR: unknown key(s) ", bad, " in ", label,
            "; known keys: ", First /@ defaults];
      Abort[]];
    MergeParams[defaults, overrides]];

ScriptParamOverrides[str_String] :=
  Module[{held, stmts},
    held = Quiet @ ToExpression[str, InputForm, Hold];
    If[Head[held] =!= Hold,
      Print["ERROR: cannot parse the parameter string:\n  ", str]; Abort[]];
    (* one held statement per assignment; a trailing ";" leaves a Null *)
    stmts = Replace[held, {
       Hold[CompoundExpression[s__]] :> List @@ Map[Hold, Hold[s]],
       Hold[s_] :> {Hold[s]}}];
    If[!AllTrue[stmts, MatchQ[#, Hold[Set[_Symbol, _]] | Hold[Null]] &],
      Print["ERROR: the parameter string must be assignments separated by \";\", ",
            "e.g. \"modelParams={B->2000.}; lyapunovSolverParams={Nint->51}\"\n  got: ",
            str];
      Abort[]];
    Association @ Cases[stmts,
      Hold[Set[s_Symbol, v_]] :> SymbolName[Unevaluated[s]] -> v]];

ScriptParamOverrides[] :=
  If[TrueQ[$Notebooks] || Length[$ScriptCommandLine] < 2, <||>,
     ScriptParamOverrides[$ScriptCommandLine[[2]]]];


(* ::Subsection:: *)
(*Closure: exact max-entropy moments*)

(* g2 = I0 I2/I1^2, g3 = I0^2 I3/I1^3, kappa fixed pointwise by I1/I0 = S.
   Below S = 0.02 both are numerically 0/0, so use the series. *)
kappaOfS[s_?NumericQ] := If[s < 1.*^-6, 0.,
   kk /. Quiet @ FindRoot[BesselI[1, kk]/BesselI[0, kk] == s,
                          {kk, 2 s/(1 - s^2)}, MaxIterations -> 200]];
g2Exact[s_?NumericQ] := If[s < 0.02, 1/2 + s^2/6 + 11 s^4/180,
   With[{k = kappaOfS[s]}, BesselI[0, k] BesselI[2, k]/BesselI[1, k]^2]];
g3Exact[s_?NumericQ] := If[s < 0.02, 1/6 + s^2/8 + 23 s^4/480,
   With[{k = kappaOfS[s]}, BesselI[0, k]^2 BesselI[3, k]/BesselI[1, k]^3]];

(* Local 6x6 flux covariance, ordered (rho x, rho y, Q1 x, Q1 y, Q2 x, Q2 y);
   its PSD implies PSD of the divergence part.  TWO EXACT ZEROS are expected:
   Qhat(theta).p = p forces two relations among the fluxes, so rank <= 4. *)
localFluxCovariance[r_, q1_, q2_, g2_, g3_] := Module[
  {a2 = g2/(2 r), a3 = g3/(4 r^2),
   k2c = q1^2 - q2^2, k2s = 2 q1 q2,
   k3c = q1^3 - 3 q1 q2^2, k3s = 3 q1^2 q2 - q2^3},
  ArrayFlatten[{
   {{{q1 + r, q2}, {q2, -q1 + r}},
    {{q1 + r/2 + a2 k2c, a2 k2s}, {a2 k2s, q1 - r/2 - a2 k2c}},
    {{q2 + a2 k2s, r/2 - a2 k2c}, {r/2 - a2 k2c, q2 - a2 k2s}}},
   {{{q1 + r/2 + a2 k2c, a2 k2s}, {a2 k2s, q1 - r/2 - a2 k2c}},
    {{3 q1/4 + r/2 + a2 k2c + a3 k3c, q2/4 + a3 k3s},
     {q2/4 + a3 k3s, -3 q1/4 + r/2 + a2 k2c - a3 k3c}},
    {{q2/4 + a2 k2s + a3 k3s, q1/4 - a3 k3c},
     {q1/4 - a3 k3c, -q2/4 + a2 k2s - a3 k3s}}},
   {{{q2 + a2 k2s, r/2 - a2 k2c}, {r/2 - a2 k2c, q2 - a2 k2s}},
    {{q2/4 + a2 k2s + a3 k3s, q1/4 - a3 k3c},
     {q1/4 - a3 k3c, -q2/4 + a2 k2s - a3 k3s}},
    {{q1/4 + r/2 - a2 k2c - a3 k3c, 3 q2/4 - a3 k3s},
     {3 q2/4 - a3 k3s, -q1/4 + r/2 - a2 k2c + a3 k3c}}}}]];


(* ::Subsection:: *)
(*Output helpers*)

runTag[corr_] := "box" <> ToString[corr["fbox"]] <> "_N" <> ToString[corr["Nint"]];

(* prepareDir and savePlot come from iterative_solver_module.m -- same private
   context, so no qualification is needed. *)

(* Plain text: Grid/NumberForm do not render under wolframscript. *)
fmtSeconds[x_] := StringPadLeft[ToString@NumberForm[N@x, {7, 2}], 8];

printTimings[timings_, nInt_] := Module[{tot, slow},
  tot = Total[timings[[All, 2]]];
  Print[];
  Print["=== timing breakdown (Nint = ", nInt, ", n = ", 3 nInt^2, ") ==="];
  Print[StringPadRight["stage", 46], fmtSeconds["  s"], "   %tracked"];
  Print[StringJoin @ ConstantArray["-", 68]];
  Do[Print[StringPadRight[r[[1]], 46], fmtSeconds[r[[2]]],
           StringPadLeft[ToString@NumberForm[N[100 r[[2]]/tot], {5, 1}], 8],
           If[r[[2]] > 4, "  <== >4 s", ""]],
     {r, SortBy[timings, -#[[2]] &]}];
  Print[StringJoin @ ConstantArray["-", 68]];
  Print[StringPadRight["TOTAL tracked", 46], fmtSeconds[tot], StringPadLeft["100.0", 8]];
  slow = Select[timings, #[[2]] > 4 &][[All, 1]];
  Print["Stages above 4 s: ", If[slow === {}, "none", StringRiffle[slow, ", "]]];
  tot];


(* ::Subsection:: *)
(*Main solver*)

SolveFluctuationsLyapunov[
    lyapunovParams : {___Rule},
    modelParams : {___Rule},
    steadyState_Association] :=
Module[
  {
    nInt, fboxV, timings = {}, rec,
    Bv, xi0, xir, zet, rho0v, Lam, av, bv, Lv, cL,
    rhoSs, Q1Ss, Q2Ss,
    h, gridFull, gridInt, iZero,
    d1Full, l1Full, D1, L1, Idn, Dx, Dy, Lx, Ly, Dxy, Lap,
    xflat, yflat, sample, diag,
    rhoV, Q1V, Q2V, rhoLap, Q1x, Q1y, Q2x, Q2y,
    Svals, g2V, g3V, kapResid, tSamp, tClo, tKchk, tCond,
    minEigLoc, scaleLoc,
    A11, A12, A13, A21, A22, A23, A31, A32, A33, Atilde, tAsmA,
    kap, pf, locPf, divBlock, h2c, h2s, h3c, h3s, e2, e3,
    Qrr, QQ1Q1, QQ2Q2, QQ1Q2, QQ2Q1, QQ1r, QQ2r, QrQ1, QrQ2, Qtilde, tAsmQ,
    Adense, Qdense, tDense, eigvals, eigvecs, tEs, maxReLam, Pev,
    eigLyapSolve, Cmat, tSol, resid, tRes1, tRef, residRef, tRes3,
    sigRhoVec, sigRhoMat, sigRhoFn, sigQ1, sigQ2, tExtr
  },

  rec[lbl_, t_] := (AppendTo[timings, {lbl, t}]; t);

  (* --- unpack + validate ---------------------------------------------
     The keys on the right MUST be Global`: the caller builds the rule lists
     in Global`, and a bare name here would resolve into this private
     context instead and match nothing. *)
  {nInt, fboxV} = {Global`Nint, Global`fbox} /. lyapunovParams;
  If[!IntegerQ[nInt] || !IntegerQ[fboxV] || nInt <= 0 || fboxV <= 0,
    Message[SolveFluctuationsLyapunov::badgrid, {fboxV, nInt}]; Return[$Failed]];
  If[!AllTrue[{"rho", "Q1", "Q2"}, KeyExistsQ[steadyState, #] &],
    Message[SolveFluctuationsLyapunov::badss, Keys[steadyState]]; Return[$Failed]];

  {Bv, xi0, xir, zet, rho0v, Lam, av, bv, Lv} =
    {Global`B, Global`\[Xi]0, Global`\[Xi]r, Global`\[Zeta], Global`\[Rho]0,
     Global`\[CapitalLambda], Global`a, Global`b, Global`L} /. modelParams;
  cL = Lv + zet xir/(4 xi0);   (* elastic constant in \[Delta]H *)

  {rhoSs, Q1Ss, Q2Ss} = steadyState /@ {"rho", "Q1", "Q2"};

  (* --- grid ---------------------------------------------------------- *)
  (* x_k = -fbox + k h needs k = (Nint+1)/2 to hit 0: ODD Nint puts the origin
     (defect core, S = 0) on the grid, EVEN Nint does not. *)
  h = (2 fboxV)/(nInt + 1);
  gridFull = Table[-fboxV + k h, {k, 0, nInt + 1}];
  gridInt = gridFull[[2 ;; -2]];
  Print["[2] Grid: fbox=", fboxV, ", Nint=", nInt,
        ", h=", ToString @ NumberForm[N[h], 5], ", state-dim=", 3 nInt^2];
  If[OddQ[nInt],
     iZero = (nInt + 1)/2;
     Print["    Nint odd  -> origin ON the grid at gridInt[[", iZero,
           "]] = ", gridInt[[iZero]], ", flat index k = ", iZero + (iZero - 1) nInt,
           " of ", nInt^2],
     iZero = None;
     Print["    Nint even -> origin NOT on the grid; nearest nodes at x = ",
           ToString @ NumberForm[N @ gridInt[[nInt/2]], 4], " and ",
           ToString @ NumberForm[N @ gridInt[[nInt/2 + 1]], 4]]];

  (* --- 1D FD matrices; Dirichlet = entries outside the interior dropped *)
  d1Full = NDSolve`FiniteDifferenceDerivative[Derivative[1], gridFull,
     "DifferenceOrder" -> 2]["DifferentiationMatrix"];
  l1Full = NDSolve`FiniteDifferenceDerivative[Derivative[2], gridFull,
     "DifferenceOrder" -> 2]["DifferentiationMatrix"];
  D1 = SparseArray @ d1Full[[2 ;; -2, 2 ;; -2]];
  L1 = SparseArray @ l1Full[[2 ;; -2, 2 ;; -2]];

  (* State ordering k = i + (j-1) N, x fastest.  Hence Dx = I (x) D1, Dy = D1 (x) I. *)
  Idn = IdentityMatrix[nInt, SparseArray];
  Dx = KroneckerProduct[Idn, D1];
  Dy = KroneckerProduct[D1, Idn];
  Lx = KroneckerProduct[Idn, L1];
  Ly = KroneckerProduct[L1, Idn];
  Dxy = KroneckerProduct[D1, D1];
  Lap = Lx + Ly;

  (* --- sample the steady state -------------------------------------- *)
  (* Batched sampling, ~200x faster than scalar calls.  N@ is load-bearing:
     exact rationals bypass the fast path and are 19x slower. *)
  xflat = Flatten @ ConstantArray[N @ gridInt, nInt];
  yflat = Flatten @ Transpose @ ConstantArray[N @ gridInt, nInt];
  sample[f_] := f[xflat, yflat];

  Print["[3] Sampling steady state on interior grid ..."];
  tSamp = rec["[3] sample steady state", First @ AbsoluteTiming[
     rhoV = sample[rhoSs];
     Q1V = sample[Q1Ss];
     Q2V = sample[Q2Ss];
     rhoLap = sample[Derivative[2, 0][rhoSs][##] + Derivative[0, 2][rhoSs][##] &];
     Q1x = sample[Derivative[1, 0][Q1Ss][##] &];
     Q1y = sample[Derivative[0, 1][Q1Ss][##] &];
     Q2x = sample[Derivative[1, 0][Q2Ss][##] &];
     Q2y = sample[Derivative[0, 1][Q2Ss][##] &];]];
  Print["    sampling (8 batched FEM calls): ", tSamp, " s"];
  Print["    \[Rho]_ss \[Element] ", MinMax@rhoV,
        "   Q1_ss \[Element] ", MinMax@Q1V,
        "   Q2_ss \[Element] ", MinMax@Q2V];

  Svals = Sqrt[Q1V^2 + Q2V^2]/rhoV;
  {tClo, {g2V, g3V}} = AbsoluteTiming[{g2Exact /@ Svals, g3Exact /@ Svals}];
  rec["[3] kappa(S) inversion + g2,g3", tClo];
  (* kappaOfS is not memoized, so this repeats all Nint^2 root-finds. *)
  {tKchk, kapResid} = AbsoluteTiming[
     Max @ Abs[(BesselI[1, #]/BesselI[0, #] & /@ (kappaOfS /@ Svals)) - Svals]];
  rec["[3] kappa(S) verification (diagnostic)", tKchk];
  Print["    S = |Q_ss|/\[Rho]_ss \[Element] ", MinMax@Svals, ",  mean ", Mean@Svals];
  Print["    g2 \[Element] ", MinMax@g2V, "   g3 \[Element] ", MinMax@g3V,
        "   (small-S limits 1/2, 1/6)"];
  Print["    \[Kappa](S) inversion: ", tClo, " s, max |I1/I0 - S| = ", kapResid];
  If[kapResid > 1.*^-8, Message[SolveFluctuationsLyapunov::kappa, kapResid]];

  (* --- drift matrix A ------------------------------------------------ *)
  (* Q = {{Q1,Q2},{Q2,-Q1}}.  The comoving term cancels the one from
     Div(dQ Grad rho_ss).  Sign a -> a + Lambda xir: rotational diffusion
     DESTROYS order. *)
  Print["[4] Assembling A ..."];
  diag[c_] := DiagonalMatrix[SparseArray @ c];
  tAsmA = rec["[4] assemble A", First @ AbsoluteTiming[
     A11 = (Bv/xi0) Lap;
     A12 = (zet/xi0) (Lx - Ly);
     A13 = (2 zet/xi0) Dxy;

     A21 = (Bv/xi0) (Q1V Lap + Q1x Dx + Q1y Dy);
     A22 = (Bv/xi0) diag[rhoLap] +
           (4/xir) (diag[-(av + Lam xir) - bv (6 Q1V^2 + 2 Q2V^2)] + cL Lap);
     A23 = -(16 bv/xir) diag[Q1V Q2V];

     A31 = (Bv/xi0) (Q2V Lap + Q2x Dx + Q2y Dy);
     A32 = A23;
     A33 = (Bv/xi0) diag[rhoLap] +
           (4/xir) (diag[-(av + Lam xir) - bv (2 Q1V^2 + 6 Q2V^2)] + cL Lap);

     Atilde = ArrayFlatten[{{A11, A12, A13},
                            {A21, A22, A23},
                            {A31, A32, A33}}];]];
  Print["    A dims = ", Dimensions[Atilde], ",  assembled in ", tAsmA, " s"];

  (* --- noise matrix Q ------------------------------------------------ *)
  (* f = d_i F_i, <F_i F_j> = C_ij delta(r-r') -> (1/h^2) sum_ij D_i.diag(C_ij).D_j^T.
     EXACT max-entropy moments:
       <QQ> = rho J + (g2/2rho) U,  <Q'Q'> = 4 rho J - (2 g2/rho) U,
       <QQQ> = (1/2) T + (g3/rho^2) W;   g2 = g3 = 0 gives the old closure. *)
  Print["[5] Assembling Q ..."];
  kap = 2 zet/(xi0 rho0v);
  pf = kap/h^2;
  locPf = 2 Lam/(rho0v h^2);
  divBlock[cxx_, cyy_, cxy_, cyx_] :=
    Dx . (cxx * Transpose[Dx]) + Dy . (cyy * Transpose[Dy]) +
    Dx . (cxy * Transpose[Dy]) + Dy . (cyx * Transpose[Dx]);

  tAsmQ = rec["[5] assemble Q (24 sparse products)", First @ AbsoluteTiming[
     h2c = Q1V^2 - Q2V^2;        h2s = 2 Q1V Q2V;
     h3c = Q1V^3 - 3 Q1V Q2V^2;  h3s = 3 Q1V^2 Q2V - Q2V^3;
     e2 = g2V/(2 rhoV);          (* U terms *)
     e3 = g3V/(4 rhoV^2);        (* W terms *)

     (* rho-rho: no closure enters *)
     Qrr = pf * divBlock[Q1V + rhoV, -Q1V + rhoV, Q2V, Q2V];

     QQ1Q1 = (pf * divBlock[3 Q1V/4 + rhoV/2 + e2 h2c + e3 h3c,
                            -3 Q1V/4 + rhoV/2 + e2 h2c - e3 h3c,
                            Q2V/4 + e3 h3s,
                            Q2V/4 + e3 h3s] +
              locPf * diag[2 rhoV - 4 e2 h2c]);

     QQ2Q2 = (pf * divBlock[Q1V/4 + rhoV/2 - e2 h2c - e3 h3c,
                            -Q1V/4 + rhoV/2 - e2 h2c + e3 h3c,
                            3 Q2V/4 - e3 h3s,
                            3 Q2V/4 - e3 h3s] +
              locPf * diag[2 rhoV + 4 e2 h2c]);

     (* NB Q1-Q2 carries a LOCAL part: <Q'Q'> is no longer proportional to J. *)
     QQ1Q2 = (pf * divBlock[Q2V/4 + e2 h2s + e3 h3s,
                            -Q2V/4 + e2 h2s - e3 h3s,
                            Q1V/4 - e3 h3c,
                            Q1V/4 - e3 h3c] +
              locPf * diag[-4 e2 h2s]);
     QQ2Q1 = Transpose[QQ1Q2];

     QQ1r = pf * divBlock[Q1V + rhoV/2 + e2 h2c,
                          Q1V - rhoV/2 - e2 h2c,
                          e2 h2s,
                          e2 h2s];
     QQ2r = pf * divBlock[Q2V + e2 h2s,
                          Q2V - e2 h2s,
                          rhoV/2 - e2 h2c,
                          rhoV/2 - e2 h2c];
     QrQ1 = Transpose[QQ1r];
     QrQ2 = Transpose[QQ2r];

     Qtilde = ArrayFlatten[{{Qrr,  QrQ1,  QrQ2},
                            {QQ1r, QQ1Q1, QQ1Q2},
                            {QQ2r, QQ2Q1, QQ2Q2}}];
     Qtilde = (Qtilde + Transpose[Qtilde])/2;    (* symmetrize *)
     ]];
  Print["    Q dims = ", Dimensions[Qtilde], ",  assembled in ", tAsmQ, " s"];

  (* --- condition (ii): local PSD of the flux covariance --------------- *)
  tCond = rec["[5] condition (ii): N^2 local 6x6 eigensolves", First @ AbsoluteTiming[
     minEigLoc = Min /@ Eigenvalues /@
        MapThread[localFluxCovariance, {rhoV, Q1V, Q2V, g2V, g3V}];]];
  scaleLoc = Max @ Abs @ Flatten @ {rhoV, Q1V, Q2V};
  Print["    condition (ii) check: ", tCond, " s"];
  Print["    local 6x6 min eig (exact closure)  \[Element] ", MinMax@minEigLoc,
        ";  negative at ", Count[minEigLoc, x_ /; x < -1.*^-10 scaleLoc],
        "/", Length@minEigLoc, " grid points"];
  If[Min@minEigLoc < -1.*^-10 scaleLoc,
     Message[SolveFluctuationsLyapunov::notpsd,
             Count[minEigLoc, x_ /; x < -1.*^-10 scaleLoc]]];

  (* --- Lyapunov solve ------------------------------------------------- *)
  (* LyapunovSolve is single-threaded; diagonalize instead: A = P D P^-1 =>
     C = P Y P^T, Y_ij = -(P^-1 Q P^-T)_ij/(lam_i + lam_j).  A is NOT self-adjoint
     (dropped Onsager partners), but the error is set by kappa(P)^2 eps ~ 1e-10;
     one refinement step then reaches ~1e-14.  The printed residual verifies it. *)
  Print["[6] Lyapunov solve via eigendecomposition + refinement ..."];
  tDense = rec["[6] densify A, Q (Normal[])", First @ AbsoluteTiming[
     Adense = Normal @ Atilde;
     Qdense = Normal @ Qtilde;]];
  Print["    [0] densify A,Q        ", tDense, " s   (",
        ToString @ NumberForm[N[2 * 8 * Length[Adense]^2/2^30], 3], " GiB)"];

  {tEs, {eigvals, eigvecs}} = AbsoluteTiming @ Eigensystem[Adense];
  rec["[6] Eigensystem (dense, nonsymmetric)", tEs];
  Print["    [a] Eigensystem        ", tEs, " s"];

  maxReLam = Max @ Re @ eigvals;
  Print["    max Re(\[Lambda]) = ", maxReLam,
        If[maxReLam < 0, "  (stable)", "  (UNSTABLE!)"]];
  If[maxReLam >= 0, Message[SolveFluctuationsLyapunov::unstable, maxReLam]];
  Print["    min |\[Lambda]_i+\[Lambda]_j| = ", 2 Min @ Abs @ Re @ eigvals,
        ",  complex fraction = ",
        N[100 Count[eigvals, z_ /; Abs[Im[z]] > 1.*^-8 Abs[z]]/Length[eigvals]], "%"];

  Pev = Transpose @ eigvecs;
  (* Solve A X + X A^T = -Rhs in the eigenbasis; reused by the refinement. *)
  eigLyapSolve[Rhs_] := Re[Pev . (
     -LinearSolve[Pev, Transpose @ LinearSolve[Pev, Transpose[Rhs]]] /
        Outer[Plus, eigvals, eigvals]) . Transpose[Pev]];

  {tSol, Cmat} = AbsoluteTiming @ eigLyapSolve[Qdense];
  rec["[6] back-substitute", tSol];
  Print["    [b] back-substitute    ", tSol, " s"];

  (* Each residual is two dense n x n matmuls, same O(n^3) as the Eigensystem. *)
  {tRes1, resid} = AbsoluteTiming[
     Max @ Abs[Adense . Cmat + Cmat . Transpose[Adense] + Qdense]];
  rec["[6] residual eval #1 (unrefined)", tRes1];
  Print["    residual (unrefined)   = ", resid, "   [", tRes1, " s]"];

  {tRef, Cmat} = AbsoluteTiming[
     Cmat + eigLyapSolve[Adense . Cmat + Cmat . Transpose[Adense] + Qdense]];
  rec["[6] refinement step (incl. residual eval #2)", tRef];
  {tRes3, residRef} = AbsoluteTiming[
     Max @ Abs[Adense . Cmat + Cmat . Transpose[Adense] + Qdense]];
  rec["[6] residual eval #3 (verification)", tRes3];
  Print["    [c] refinement step    ", tRef, " s"];
  Print["    residual (refined)     = ", residRef, "   [", tRes3, " s]",
        "   (improved ", ToString@NumberForm[N[resid/residRef], 4], "x)"];
  If[residRef > resid,
     Message[SolveFluctuationsLyapunov::norefine, resid, residRef]];
  Print["    total solve time = ", tEs + tSol + tRef, " s"];

  (* --- extract -------------------------------------------------------- *)
  Print["[7] Extracting diagonal blocks ..."];
  (* ListInterpolation, NOT {x,y,z} triples: the triple form makes
     ListDensityPlot triangulate.  Wants arr[[i,j]] = f(x_i,y_j). *)
  tExtr = rec["[7] extract diagonals + build interpolant", First @ AbsoluteTiming[
     sigRhoVec = Diagonal @ Cmat[[1 ;; nInt^2, 1 ;; nInt^2]];
     sigQ1 = Partition[Diagonal @ Cmat[[nInt^2 + 1 ;; 2 nInt^2,
                                        nInt^2 + 1 ;; 2 nInt^2]], nInt];
     sigQ2 = Partition[Diagonal @ Cmat[[2 nInt^2 + 1 ;; 3 nInt^2,
                                        2 nInt^2 + 1 ;; 3 nInt^2]], nInt];
     sigRhoMat = Partition[sigRhoVec, nInt];  (* sigRhoMat[[j,i]] = Sigma(x_i,y_j) *)
     sigRhoFn = ListInterpolation[Transpose[sigRhoMat], {gridInt, gridInt},
                                  InterpolationOrder -> 2];]];
  Print["    extract + interpolant: ", tExtr, " s"];
  Print["Sigma_rho stats: min=", Min@sigRhoVec, "  max=", Max@sigRhoVec,
        "  mean=", Mean@sigRhoVec];

  printTimings[timings, nInt];

  <| "C" -> Cmat,
     "sigmaRho" -> sigRhoMat,
     "sigmaRhoFn" -> sigRhoFn,
     "sigmaQ1" -> sigQ1,
     "sigmaQ2" -> sigQ2,
     "gridInt" -> gridInt,
     "h" -> N[h],
     "Nint" -> nInt,
     "fbox" -> fboxV,
     "eigenvalues" -> eigvals,
     "maxReLambda" -> maxReLam,
     "residual" -> resid,
     "residualRefined" -> residRef,
     "minEigLocal" -> MinMax@minEigLoc,
     "modelParams" -> modelParams,
     "lyapunovParams" -> lyapunovParams,
     "timings" -> timings |>
];


(* ::Subsection:: *)
(*Visualization and export*)

visualizeFluctuations[correlations_Association, outputDir_, plotBox_?NumericQ] :=
Module[{dir, gridInt, sigRhoFn, xLo, xHi, plotPts, lbl, plot2D, plot3D, tag, x, y},
  dir = prepareDir[outputDir];
  gridInt = correlations["gridInt"];
  sigRhoFn = correlations["sigmaRhoFn"];
  tag = runTag[correlations];
  (* plotBox crops the Dirichlet boundary layer out of the figures; the
     interpolant only spans +-(fbox - h), so clamp to that. *)
  {xLo, xHi} = {Max[-plotBox, First@gridInt], Min[plotBox, Last@gridInt]};
  plotPts = Min[correlations["Nint"], 200];
  lbl = ("Var(\[Delta]\[Rho])(r)  [Dirichlet BC, h=" <>
         ToString @ NumberForm[correlations["h"], 3] <> "]");

  plot2D = DensityPlot[sigRhoFn[x, y], {x, xLo, xHi}, {y, xLo, xHi},
     PlotPoints -> plotPts, MaxRecursion -> 0,
     PlotRange -> {{-plotBox, plotBox}, {-plotBox, plotBox}},
     PlotLegends -> Automatic, FrameLabel -> {"x", "y"},
     PlotLabel -> lbl, ColorFunction -> "SunsetColors"];

  plot3D = Plot3D[sigRhoFn[x, y], {x, xLo, xHi}, {y, xLo, xHi},
     PlotPoints -> plotPts, MaxRecursion -> 0, Mesh -> None,
     PlotRange -> {{-plotBox, plotBox}, {-plotBox, plotBox}, Automatic},
     AxesLabel -> {"x", "y", "Var(\[Delta]\[Rho])"},
     PlotLabel -> lbl, ColorFunction -> "SunsetColors"];

  savePlot[dir, "sigma_rho_dirichlet_" <> tag <> ".png", plot2D];
  savePlot[dir, "sigma_rho_surface_" <> tag <> ".png", plot3D];
  {plot2D, plot3D}];

exportData[correlations_Association, outputDir_, saveCovariance_ : False] :=
Module[{dir, tag, out, dumpPath, binPath, metaPath, Cmat, n, lyapTime, ssTime,
        paths = {}},
  dir = prepareDir[outputDir];
  If[dir === None, Return[{}]];
  tag = runTag[correlations];
  out[stem_, ext_] := FileNameJoin[{dir, stem <> tag <> ext}];

  dumpPath = out["sigma_rho_", ".m"];
  Export[dumpPath,
     <|"fbox" -> correlations["fbox"], "Nint" -> correlations["Nint"],
       "h" -> correlations["h"], "gridInt" -> correlations["gridInt"],
       "sigmaRho" -> correlations["sigmaRho"],
       "sigmaQ1" -> correlations["sigmaQ1"],
       "sigmaQ2" -> correlations["sigmaQ2"],
       "maxReLambda" -> correlations["maxReLambda"],
       "residualRefined" -> correlations["residualRefined"],
       "modelParams" -> correlations["modelParams"]|>];
  Print["Saved: ", dumpPath];
  AppendTo[paths, dumpPath];

  (* Full covariance, row-major Real64, n = 3 Nint^2 with blocks (drho, dQ1, dQ2).
     Read back with ArrayReshape[BinaryReadList[f, "Real64"], {n, n}].  8 n^2
     bytes -- 487 MB at Nint = 51. *)
  If[TrueQ[saveCovariance],
     Cmat = correlations["C"];
     n = Length[Cmat];
     binPath = out["covariance_", ".bin"];
     metaPath = out["covariance_", "_meta.m"];
     Export[binPath, Cmat, {"Binary", "Real64"}];
     lyapTime = Total[correlations["timings"][[All, 2]]];
     ssTime = Lookup[correlations, "steadyStateTime", None];
     Export[metaPath,
        <|"n" -> n, "Nint" -> correlations["Nint"], "fbox" -> correlations["fbox"],
          "h" -> correlations["h"], "gridInt" -> correlations["gridInt"],
          "blocks" -> {"rho", "Q1", "Q2"},
          "ordering" -> "row-major Real64; grid index k = i + (j-1) Nint, x fastest",
          "lyapunovTime" -> lyapTime,
          "steadyStateTime" -> ssTime,
          "totalTime" -> If[NumericQ[ssTime], ssTime + lyapTime, lyapTime],
          "timings" -> correlations["timings"]|>];
     Print["Saved: ", binPath, "  (",
           ToString@NumberForm[N[FileByteCount[binPath]/2^20], {6, 1}], " MiB)"];
     Print["Saved: ", metaPath];
     paths = Join[paths, {binPath, metaPath}]];
  paths];

End[];
