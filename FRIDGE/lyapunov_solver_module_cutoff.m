(* ::Package:: *)

(* ::Section:: *)
(*Gaussian fluctuations around the \[Rho]\[Dash]Q steady state: discrete Lyapunov solve.*)

(* psi = (d\[Rho], dQ1, dQ2)^T on a uniform interior grid with psi|_bdry = 0;
   psi_t = A psi + noise, Sigma solves A Sigma + Sigma A^T + Q = 0.
   Noise uses the EXACT max-entropy moments; the leading-order closure is not
   PSD above S = 0.618.

   ========================== CUTOFF VARIANT ==========================
   This is lyapunov_solver_module_R1.m (consistent stencils, R == 1) with the
   grid REPLACED as the UV regulator by a PHYSICAL one: the noise is given a
   finite spatial correlation length lNoise, a new model parameter.

       <F_i(r,t) F_j(r',t')> = M_ij(r) G_l(r-r') delta(t-t'),
       G_l(r) = exp(-r^2/(2 l^2))/(2 Pi l^2),   Ghat_l(k) = exp(-k^2 l^2/2)

   Only the delta in TIME is needed for the Lyapunov equation (it is what makes
   psi Markovian); the spatial kernel is unconstrained, so this costs nothing
   and changes no step of the derivation.  The k^2 of the divergence-form noise
   still cancels against the k^2 of the relaxation rate, so

       S(k) = (zeta/(B rho0)) exp(-k^2 l^2/2)

   which is INTEGRABLE, and the equal-point variance is finite with no cutoff
   imposed by hand:

       Sigma_iso(r) = zeta rho_ss(r)/(2 Pi B rho0 lNoise^2).

   That is dean.tex eq:iso-final term for term: l_UV there is the noise
   CORRELATION LENGTH, not a smoothing width and not the grid spacing.  Setting
   lNoise = 0 recovers the R1 result exactly.  See dean.tex sec-noise-cutoff.

   IMPLEMENTATION.  The kernel enters as a congruence on the noise matrix that
   eq:Q-assembly already builds:

       Q_l = W . Q . W,      W = MatrixExp[(lNoise^2/4) Lap]

   whose symbol exp(-k^2 l^2/4), applied twice, is Ghat_l.  Four reasons for
   this form rather than a dense kernel sandwiched inside the flux products:
     - W = W^T, so Q_l stays symmetric;
     - a congruence preserves PSD, so the realizability structure of
       dean.tex sec-realizability carries over verbatim;
     - W is a function of the DIRICHLET Laplacian, so the BC is exact and there
       is no kernel clipping at the boundary;
     - Lx.Ly = Ly.Lx = Kron(L1,L1), so the exponential splits EXACTLY as
       W = Kron(E,E) with E = MatrixExp[(lNoise^2/4) L1], one Nint x Nint
       matrix.  E is built in closed form from the sine basis that diagonalises
       L1 (eigenvalues -4 sin^2(Pi p/(2(Nint+1)))/h^2), so no numerical matrix
       exponential is ever needed and the mode weight underflows gracefully to 0
       when lNoise >> h.

   ACCURACY.  The lattice weight decays more slowly than the continuum Gaussian
   (sin^2(phi/2) < phi^2/4), so it keeps slightly too much short-wavelength
   noise.  The error depends on lNoise/h ALONE and is +(1/4)(h/l)^2 at leading
   order: 36% at l = h, 7.7% at l = 2h, 1.0% at l = 5h, 0.25% at l = 10h.
   Require lNoise >~ 5 h.  Below lNoise ~ h the Gaussian is invisible to the
   grid and the answer reverts to the grid-limited R1 value -- such a run
   measures the mesh, not the model.

   ---------------- inherited from the R1 variant ----------------
   Consistent stencil set, so that the discrete fluctuation-dissipation balance
   holds mode by mode and the (grid-regulated) UV pedestal is exactly

       Sigma_iso(r) = zeta rho_ss(r)/(B rho0 h^2),      i.e. R == 1.

   The original module divides the divergence-form noise with the CENTRAL
   difference D1 (symbol i sin(phi)/h) while the drift uses the COMPACT
   Laplacian (symbol -4 sin^2(phi/2)/h^2).  They are not an adjoint pair, so
   each mode gets C_k = (zeta/(B rho0 h^2)) R(p,q) with

       R = (sin^2 p + sin^2 q)/(4(sin^2(p/2) + sin^2(q/2))),   <R> = 0.363380,

   i.e. the pedestal comes out 1/<R> = 2.752x too small.  See dean.tex,
   appendix "Linear OU field theory" (sec-ou).

   WHAT CHANGED, and only this:

   (1) G1 = forward difference from the Nint interior NODES to the Nint+1 cell
       EDGES.  Then -(Gx^T Gx + Gy^T Gy) = Lap EXACTLY for Dirichlet -- the
       boundary edge supplies the stencil entry an Nint x Nint collocated
       backward difference would drop.  Asserted at run time in [2].

   (2) The DIAGONAL noise flux terms use Gx, Gy with their coefficients sampled
       on the corresponding cell edges; the CROSS terms keep the collocated
       central differences Dx, Dy with node-sampled coefficients.  This split is
       not a compromise -- it is forced.  A corner-interpolated cross term picks
       up a phase cos((q-p)/2) from the half-cell offset between an x-edge and a
       y-edge, which destroys the parity in p behind the continuum angular
       cancellation for traceless Q^ss, and injects a SPURIOUS anisotropic
       pedestal 0.087 Q2(r) zeta/(B rho0 h^2) -- an ~8% h^-2 modulation at
       S ~ 0.86.  The collocated symbol i sin(phi)/h is phase-free, so
       sin p sin q stays odd in p and the cancellation is exact.  Measured BZ
       averages (claude_experiments/theory_numeric_discrepancy/uv_regulator.wls):

           scheme                             c_iso    c_Q1   c_Q2
           central everywhere                 0.36338    0     0
           fully staggered, corner cross      1          0     0.0870
           staggered diagonal + central cross 1          0     0

   (3) A21, A31 use the divergence form -(B/xi0) sum_i Gi^T diag(Q|edge) Gi in
       place of the product-rule expansion (Q1V Lap + Q1x Dx + Q1y Dy).  The two
       differ by 2.8% at Nint = 41 and only the divergence form is exactly
       self-adjoint, as Div(c Grad .) is in the continuum.  A11, A12 and A13 are
       already exactly the staggered operators and are UNCHANGED -- verified in
       claude_experiments/theory_numeric_discrepancy/drift_consistency.wls:
       Lx = -Gx^T Gx, Ly = -Gy^T Gy, and Dxy is exactly the staggered composite
       div_x(d_y Q_xy|corner), so its Nyquist zero is intrinsic to d_x d_y and
       not a stencil defect.

   The locPf (rotational, Lambda) noise terms are not in divergence form, so the
   R1 stencil change left them alone -- but W DOES act on them, because the
   rotational noise is particle discreteness too and carries the same correlation
   length; Q_l = W.Q.W is applied to the whole assembled matrix.
   ---------------------------------------------------------------

   Entry point is SolveFluctuationsLyapunovCutoff and runTag carries BOTH the
   grid and an "_lN<value>" suffix, so no two (Nint, lNoise) runs can collide in
   one directory.  Do NOT load this module together with
   lyapunov_solver_module.m or lyapunov_solver_module_R1.m in the same kernel:
   the private helpers (runTag, exportData, ...) are shared and the last one
   loaded wins.
   ====================================================================

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

Clear[SolveFluctuationsLyapunovCutoff, visualizeFluctuations,
      exportData, MergeParams, ScriptParamOverrides];


SolveFluctuationsLyapunovCutoff::usage =
  "SolveFluctuationsLyapunovCutoff[lyapunovParams, modelParams, steadyState] solves " <>
  "A \[CapitalSigma] + \[CapitalSigma] A^T + Q = 0 for the equal-time covariance of " <>
  "(\[Delta]\[Rho], \[Delta]Q1, \[Delta]Q2) linearized about steadyState, on a uniform " <>
  "Nint x Nint grid over [-fbox, fbox]^2 with Dirichlet (zero) boundary values.\n\n" <>
  "The noise carries a Gaussian spatial correlation of length lNoise (a " <>
  "modelParams key), so the equal-point variance is UV-finite: " <>
  "\[CapitalSigma]_iso = \[Zeta] \[Rho]_ss/(2 \[Pi] B \[Rho]0 lNoise^2). " <>
  "lNoise = 0 reproduces the grid-regulated R1 result.\n\n" <>
  "  lyapunovParams : rules {Nint->..., fbox->...}; both positive integers, " <>
  "h = 2 fbox/(Nint+1) is kept exact so odd Nint puts the defect core on a node\n" <>
  "  modelParams    : same rules as SolveActiveNematicSteady, PLUS lNoise >= 0 " <>
  "(the noise correlation length, in the same length units as fbox; lNoise >= 5 h " <>
  "is required for the continuum value to within 1%)\n" <>
  "  steadyState    : <|\"rho\"->..., \"Q1\"->..., \"Q2\"->..., \"mesh\"->...|>\n\n" <>
  "Returns <|\"C\", \"sigmaRho\", \"sigmaRhoFn\", \"sigmaQ1\", \"sigmaQ2\", " <>
  "\"gridInt\", \"h\", \"Nint\", \"fbox\", \"lNoise\", \"Wsymbol\", " <>
  "\"eigenvalues\", \"maxReLambda\", \"residual\", \"residualRefined\", " <>
  "\"minEigLocal\", \"minEigQ\", \"lapIdentity\", \"scheme\", \"modelParams\", " <>
  "\"lyapunovParams\", \"timings\"|>.";

SolveFluctuationsLyapunovCutoff::badgrid =
  "Nint and fbox must be positive integers (h = 2 fbox/(Nint+1) is kept exact); got `1`.";
SolveFluctuationsLyapunovCutoff::badss =
  "steadyState must be an association with keys \"rho\", \"Q1\", \"Q2\"; got `1`.";
SolveFluctuationsLyapunovCutoff::badlnoise =
  "modelParams must carry a numeric lNoise >= 0 (the noise correlation length); got `1`.";
SolveFluctuationsLyapunovCutoff::subgrid =
  "lNoise = `1` is only `2` grid spacings: the Gaussian is barely resolved and the \
result is grid-limited, not the continuum value.  Use lNoise >= 5 h.";
SolveFluctuationsLyapunovCutoff::unstable =
  "Linearization is not strictly stable: max Re(\[Lambda]) = `1`; results may be unphysical.";
SolveFluctuationsLyapunovCutoff::notpsd =
  "Condition (ii) Q >= 0 is violated even with the exact moments at `1` grid points; \
check the closure assembly.";
SolveFluctuationsLyapunovCutoff::norefine =
  "Refinement did not improve the residual (`1` -> `2`); consider LyapunovSolve.";
SolveFluctuationsLyapunovCutoff::kappa =
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

(* The tag must separate BOTH the grid and the noise correlation length: a sweep
   writes 18 x 5 runs into one directory and nothing may be overwritten.  Four
   decimals, never scientific, so 0.1000 / 0.5000 / 1.0000 / 1.7321 / 5.0000 are
   all distinct and sort sensibly.  The "." is harmless in a filename and to the
   ".bin" -> "_meta.m" replacement readers use. *)
lTag[l_] := "_lN" <> ToString@NumberForm[N[l], {8, 4},
   ExponentFunction -> (Null &), NumberPadding -> {"", "0"}];

runTag[corr_] := "box" <> ToString[corr["fbox"]] <> "_N" <> ToString[corr["Nint"]] <>
                 lTag[corr["lNoise"]];

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

SolveFluctuationsLyapunovCutoff[
    lyapunovParams : {___Rule},
    modelParams : {___Rule},
    steadyState_Association] :=
Module[
  {
    nInt, fboxV, timings = {}, rec,
    Bv, xi0, xir, zet, rho0v, Lam, av, bv, Lv, cL,
    rhoSs, Q1Ss, Q2Ss,
    h, gridFull, gridInt, iZero, xFace,
    d1Full, l1Full, D1, L1, Idn, Dx, Dy, Lx, Ly, Dxy, Lap,
    G1, Gx, Gy, lapIdentity,
    xflat, yflat, sample, diag,
    xeX, xeY, yeX, yeY, sampleXE, sampleYE,
    rhoV, Q1V, Q2V, rhoLap, Q1x, Q1y, Q2x, Q2y,
    rhoXE, Q1XE, Q2XE, rhoYE, Q1YE, Q2YE,
    SXE, SYE, g2XE, g3XE, g2YE, g3YE, e2XE, e3XE, e2YE, e3YE,
    h2cXE, h2sXE, h3cXE, h3sXE, h2cYE, h2sYE, h3cYE, h3sYE,
    Svals, g2V, g3V, kapResid, tSamp, tClo, tKchk, tCond,
    minEigLoc, scaleLoc, minEigQ, tQpsd,
    A11, A12, A13, A21, A22, A23, A31, A32, A33, Atilde, tAsmA,
    kap, pf, locPf, divBlockR1, h2c, h2s, h3c, h3s, e2, e3,
    Qrr, QQ1Q1, QQ2Q2, QQ1Q2, QQ2Q1, QQ1r, QQ2r, QrQ1, QrQ2, Qtilde, tAsmQ,
    Adense, Qdense, tDense, eigvals, eigvecs, tEs, maxReLam, Pev,
    eigLyapSolve, Cmat, tSol, resid, tRes1, tRef, residRef, tRes3,
    sigRhoVec, sigRhoMat, sigRhoFn, sigQ1, sigQ2, tExtr,
    lNoiseV, Wsym, Emat, Esym, smoothBlock, tSmooth, blk, Qasym
  },

  rec[lbl_, t_] := (AppendTo[timings, {lbl, t}]; t);

  (* --- unpack + validate ---------------------------------------------
     The keys on the right MUST be Global`: the caller builds the rule lists
     in Global`, and a bare name here would resolve into this private
     context instead and match nothing. *)
  {nInt, fboxV} = {Global`Nint, Global`fbox} /. lyapunovParams;
  If[!IntegerQ[nInt] || !IntegerQ[fboxV] || nInt <= 0 || fboxV <= 0,
    Message[SolveFluctuationsLyapunovCutoff::badgrid, {fboxV, nInt}]; Return[$Failed]];
  If[!AllTrue[{"rho", "Q1", "Q2"}, KeyExistsQ[steadyState, #] &],
    Message[SolveFluctuationsLyapunovCutoff::badss, Keys[steadyState]]; Return[$Failed]];

  {Bv, xi0, xir, zet, rho0v, Lam, av, bv, Lv} =
    {Global`B, Global`\[Xi]0, Global`\[Xi]r, Global`\[Zeta], Global`\[Rho]0,
     Global`\[CapitalLambda], Global`a, Global`b, Global`L} /. modelParams;
  cL = Lv + zet xir/(4 xi0);   (* elastic constant in \[Delta]H *)

  (* the cutoff parameter.  A missing key would silently stay the symbol lNoise
     and poison MatrixExp, so check it here rather than failing deep in [5b]. *)
  lNoiseV = Global`lNoise /. modelParams;
  If[!NumericQ[lNoiseV] || N[lNoiseV] < 0,
    Message[SolveFluctuationsLyapunovCutoff::badlnoise, lNoiseV]; Return[$Failed]];
  lNoiseV = N[lNoiseV];

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

  (* --- R1: forward difference, interior NODES -> the nInt+1 cell EDGES.
     Edge m lies between nodes m-1 and m, at x = -fbox + (m - 1/2) h, with
     u_0 = u_{nInt+1} = 0 (Dirichlet), so (G1 u)_m = (u_m - u_{m-1})/h.
     The (nInt+1)-th row is what makes -G^T G = L1 hold EXACTLY: it carries the
     stencil entry that an nInt x nInt collocated backward difference drops. *)
  G1 = SparseArray[Join[Table[{m, m} -> 1/h, {m, nInt}],
                        Table[{m + 1, m} -> -1/h, {m, nInt}]], {nInt + 1, nInt}];
  Gx = KroneckerProduct[Idn, G1];
  Gy = KroneckerProduct[G1, Idn];
  xFace = Table[-fboxV + (m - 1/2) h, {m, nInt + 1}];

  (* The identity the whole scheme rests on.  Exact in rationals, so test == 0. *)
  lapIdentity = Max @ Abs @ Normal[
     Lap + (Transpose[Gx] . Gx + Transpose[Gy] . Gy)];
  Print["    R1: |Lap + (Gx^T Gx + Gy^T Gy)| = ", lapIdentity,
        If[lapIdentity == 0, "   (exact -- adjoint pair)",
           "   *** NOT AN ADJOINT PAIR, R != 1 ***"]];
  If[lapIdentity != 0,
     Print["    ABORT: the staggered/compact identity failed; the R1 scheme is void."];
     Abort[]];

  (* --- the noise-correlation operator W = Kron[Emat, Emat] -------------
     W = MatrixExp[(lNoise^2/4) Lap] and Lap = Lx + Ly with Lx.Ly = Ly.Lx, so the
     exponential splits exactly into one nInt x nInt factor.  L1 is diagonalised
     ANALYTICALLY by the Dirichlet sine basis
         v_p(i) = Sqrt[2/(nInt+1)] Sin[Pi p i/(nInt+1)],
         mu_p   = -4 Sin[Pi p/(2(nInt+1))]^2/h^2,
     so Emat = V . diag(Wsym) . V^T with Wsym_p = Exp[-(lNoise/h)^2 Sin[...]^2].
     Closed form beats MatrixExp here: it is one symmetric matmul, it is exact,
     and when lNoise >> h the weights underflow to 0 -- which is the correct
     behaviour (those modes carry no noise) rather than an overflow. *)
  If[lNoiseV == 0.,
     Wsym = ConstantArray[1., nInt];
     Emat = None;   (* short-circuit: W = I, so this run reproduces R1 exactly *)
     Print["    lNoise = 0 -> W = I; this run is the grid-regulated R1 scheme"],
     (* else *)
     Wsym = Exp[-(lNoiseV/N[h])^2 Sin[N[Pi] Range[nInt]/(2 (nInt + 1))]^2];
     Emat = With[{V = Sqrt[2./(nInt + 1)] Table[
                       Sin[N[Pi] p i/(nInt + 1)], {i, nInt}, {p, nInt}]},
                 V . (Wsym * Transpose[V])];
     Esym = Max @ Abs[Emat - Transpose[Emat]];
     Print["    lNoise = ", lNoiseV, " = ", ToString@NumberForm[lNoiseV/N[h], 4],
           " h;  W = Kron[E,E], E symmetry residual = ", Esym,
           ", max |E| = ", Max@Abs@Emat];
     Print["    W mode weights \[Element] ", MinMax@Wsym,
           "  (Exp[-(lNoise/h)^2 Sin[Pi p/(2(Nint+1))]^2])"];
     If[Esym > 1.*^-12, Print["    WARNING: E is not symmetric to 1e-12."]];
     If[lNoiseV < 5 N[h],
        Message[SolveFluctuationsLyapunovCutoff::subgrid, lNoiseV,
                ToString@NumberForm[lNoiseV/N[h], 3]]]];

  (* --- sample the steady state -------------------------------------- *)
  (* Batched sampling, ~200x faster than scalar calls.  N@ is load-bearing:
     exact rationals bypass the fast path and are 19x slower. *)
  xflat = Flatten @ ConstantArray[N @ gridInt, nInt];
  yflat = Flatten @ Transpose @ ConstantArray[N @ gridInt, nInt];
  sample[f_] := f[xflat, yflat];

  (* R1: the diagonal flux coefficients live on the cell edges, so sample the FEM
     steady state THERE rather than averaging node values -- the interpolant is
     defined everywhere, so this costs nothing in accuracy.
       x-edges: (face_x, node_y), (nInt+1) x nInt, face index fastest
       y-edges: (node_x, face_y), nInt x (nInt+1), node index fastest
     The layout matches Gx = Idn (x) G1 and Gy = G1 (x) Idn respectively. *)
  xeX = Flatten @ ConstantArray[N @ xFace, nInt];
  xeY = Flatten @ Transpose @ ConstantArray[N @ gridInt, nInt + 1];
  yeX = Flatten @ ConstantArray[N @ gridInt, nInt + 1];
  yeY = Flatten @ Transpose @ ConstantArray[N @ xFace, nInt];
  sampleXE[f_] := f[xeX, xeY];
  sampleYE[f_] := f[yeX, yeY];

  Print["[3] Sampling steady state on interior grid ..."];
  tSamp = rec["[3] sample steady state", First @ AbsoluteTiming[
     rhoV = sample[rhoSs];
     Q1V = sample[Q1Ss];
     Q2V = sample[Q2Ss];
     rhoLap = sample[Derivative[2, 0][rhoSs][##] + Derivative[0, 2][rhoSs][##] &];
     Q1x = sample[Derivative[1, 0][Q1Ss][##] &];
     Q1y = sample[Derivative[0, 1][Q1Ss][##] &];
     Q2x = sample[Derivative[1, 0][Q2Ss][##] &];
     Q2y = sample[Derivative[0, 1][Q2Ss][##] &];
     (* R1: the same three fields on the two edge sets *)
     rhoXE = sampleXE[rhoSs];  Q1XE = sampleXE[Q1Ss];  Q2XE = sampleXE[Q2Ss];
     rhoYE = sampleYE[rhoSs];  Q1YE = sampleYE[Q1Ss];  Q2YE = sampleYE[Q2Ss];]];
  Print["    sampling (14 batched FEM calls: 8 nodal + 6 on the edges): ", tSamp, " s"];
  Print["    \[Rho]_ss \[Element] ", MinMax@rhoV,
        "   Q1_ss \[Element] ", MinMax@Q1V,
        "   Q2_ss \[Element] ", MinMax@Q2V];

  Svals = Sqrt[Q1V^2 + Q2V^2]/rhoV;
  SXE = Sqrt[Q1XE^2 + Q2XE^2]/rhoXE;
  SYE = Sqrt[Q1YE^2 + Q2YE^2]/rhoYE;
  (* R1: the closure is needed wherever a coefficient is needed, so three times
     over.  ~3 Nint^2 extra root-finds; at Nint = 131 this stage goes from ~11 s
     to ~33 s, against a multi-hour Eigensystem. *)
  {tClo, {g2V, g3V, g2XE, g3XE, g2YE, g3YE}} = AbsoluteTiming[
     {g2Exact /@ Svals, g3Exact /@ Svals,
      g2Exact /@ SXE,   g3Exact /@ SXE,
      g2Exact /@ SYE,   g3Exact /@ SYE}];
  rec["[3] kappa(S) inversion + g2,g3 (nodes + both edge sets)", tClo];
  (* kappaOfS is not memoized, so this repeats all Nint^2 root-finds. *)
  {tKchk, kapResid} = AbsoluteTiming[
     Max @ Abs[(BesselI[1, #]/BesselI[0, #] & /@ (kappaOfS /@ Svals)) - Svals]];
  rec["[3] kappa(S) verification (diagnostic)", tKchk];
  Print["    S = |Q_ss|/\[Rho]_ss \[Element] ", MinMax@Svals, ",  mean ", Mean@Svals];
  Print["    S on x-edges \[Element] ", MinMax@SXE,
        ",  on y-edges \[Element] ", MinMax@SYE];
  Print["    g2 \[Element] ", MinMax@g2V, "   g3 \[Element] ", MinMax@g3V,
        "   (small-S limits 1/2, 1/6)"];
  Print["    \[Kappa](S) inversion: ", tClo, " s, max |I1/I0 - S| = ", kapResid];
  If[kapResid > 1.*^-8, Message[SolveFluctuationsLyapunovCutoff::kappa, kapResid]];

  (* --- drift matrix A ------------------------------------------------ *)
  (* Q = {{Q1,Q2},{Q2,-Q1}}.  The comoving term cancels the one from
     Div(dQ Grad rho_ss).  Sign a -> a + Lambda xir: rotational diffusion
     DESTROYS order.

     R1 note.  A11, A12, A13 are ALREADY exactly the staggered operators and are
     left alone:  Lx = -Gx^T Gx and Ly = -Gy^T Gy separately, so
        A11 = (B/xi0) Lap        = -(B/xi0)(Gx^T Gx + Gy^T Gy)
        A12 = (zeta/xi0)(Lx-Ly)  = (zeta/xi0)(-Gx^T Gx + Gy^T Gy)
     and Dxy is exactly div_x(d_y Q_xy|corner) + div_y(d_x Q_yx|corner) with Q2
     averaged from nodes to corners, so A13 is the staggered double divergence
     too and its Nyquist zero is intrinsic to d_x d_y rather than a defect.
     Verified in drift_consistency.wls (all three identities exact to 0).

     Only A21, A31 change: Div(Q_ss Grad drho) was expanded by the product rule
     at the nodes, which is NOT self-adjoint (asymmetry 1.6e-2 at Nint = 41)
     although Div(c Grad .) is in the continuum.  The divergence form below is
     exactly self-adjoint and reuses the edge samples the noise already needs.
     The two forms differ by 2.8% at Nint = 41; the effect on Var(drho) is
     bounded by the whole delta-Q sector, measured at ~0.34%. *)
  Print["[4] Assembling A ..."];
  diag[c_] := DiagonalMatrix[SparseArray @ c];
  tAsmA = rec["[4] assemble A", First @ AbsoluteTiming[
     A11 = (Bv/xi0) Lap;
     A12 = (zet/xi0) (Lx - Ly);
     A13 = (2 zet/xi0) Dxy;

     A21 = -(Bv/xi0) (Transpose[Gx] . (Q1XE * Gx) + Transpose[Gy] . (Q1YE * Gy));
     A22 = (Bv/xi0) diag[rhoLap] +
           (4/xir) (diag[-(av + Lam xir) - bv (6 Q1V^2 + 2 Q2V^2)] + cL Lap);
     A23 = -(16 bv/xir) diag[Q1V Q2V];

     A31 = -(Bv/xi0) (Transpose[Gx] . (Q2XE * Gx) + Transpose[Gy] . (Q2YE * Gy));
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
       <QQQ> = (1/2) T + (g3/rho^2) W;   g2 = g3 = 0 gives the old closure.

     R1: the ii terms are the ADJOINT PARTNERS of the drift Laplacian and use
     Gx, Gy with edge-sampled coefficients; the i != j terms keep the collocated
     Dx, Dy with node-sampled coefficients.  Setting every coefficient to 1 makes
     the first two terms exactly -Lap, which is what forces R == 1.  See the
     header for why the cross terms must NOT be staggered. *)
  Print["[5] Assembling Q ..."];
  kap = 2 zet/(xi0 rho0v);
  pf = kap/h^2;
  locPf = 2 Lam/(rho0v h^2);
  divBlockR1[cxxE_, cyyE_, cxy_, cyx_] :=
    Transpose[Gx] . (cxxE * Gx) + Transpose[Gy] . (cyyE * Gy) +
    Dx . (cxy * Transpose[Dy]) + Dy . (cyx * Transpose[Dx]);

  tAsmQ = rec["[5] assemble Q (24 sparse products)", First @ AbsoluteTiming[
     h2c = Q1V^2 - Q2V^2;        h2s = 2 Q1V Q2V;
     h3c = Q1V^3 - 3 Q1V Q2V^2;  h3s = 3 Q1V^2 Q2V - Q2V^3;
     e2 = g2V/(2 rhoV);          (* U terms *)
     e3 = g3V/(4 rhoV^2);        (* W terms *)
     (* the same combinations on the two edge sets, for the ii flux terms *)
     h2cXE = Q1XE^2 - Q2XE^2;        h2sXE = 2 Q1XE Q2XE;
     h3cXE = Q1XE^3 - 3 Q1XE Q2XE^2; h3sXE = 3 Q1XE^2 Q2XE - Q2XE^3;
     e2XE = g2XE/(2 rhoXE);          e3XE = g3XE/(4 rhoXE^2);
     h2cYE = Q1YE^2 - Q2YE^2;        h2sYE = 2 Q1YE Q2YE;
     h3cYE = Q1YE^3 - 3 Q1YE Q2YE^2; h3sYE = 3 Q1YE^2 Q2YE - Q2YE^3;
     e2YE = g2YE/(2 rhoYE);          e3YE = g3YE/(4 rhoYE^2);

     (* rho-rho: no closure enters *)
     Qrr = pf * divBlockR1[Q1XE + rhoXE, -Q1YE + rhoYE, Q2V, Q2V];

     QQ1Q1 = (pf * divBlockR1[3 Q1XE/4 + rhoXE/2 + e2XE h2cXE + e3XE h3cXE,
                              -3 Q1YE/4 + rhoYE/2 + e2YE h2cYE - e3YE h3cYE,
                              Q2V/4 + e3 h3s,
                              Q2V/4 + e3 h3s] +
              locPf * diag[2 rhoV - 4 e2 h2c]);

     QQ2Q2 = (pf * divBlockR1[Q1XE/4 + rhoXE/2 - e2XE h2cXE - e3XE h3cXE,
                              -Q1YE/4 + rhoYE/2 - e2YE h2cYE + e3YE h3cYE,
                              3 Q2V/4 - e3 h3s,
                              3 Q2V/4 - e3 h3s] +
              locPf * diag[2 rhoV + 4 e2 h2c]);

     (* NB Q1-Q2 carries a LOCAL part: <Q'Q'> is no longer proportional to J. *)
     QQ1Q2 = (pf * divBlockR1[Q2XE/4 + e2XE h2sXE + e3XE h3sXE,
                              -Q2YE/4 + e2YE h2sYE - e3YE h3sYE,
                              Q1V/4 - e3 h3c,
                              Q1V/4 - e3 h3c] +
              locPf * diag[-4 e2 h2s]);
     QQ2Q1 = Transpose[QQ1Q2];

     QQ1r = pf * divBlockR1[Q1XE + rhoXE/2 + e2XE h2cXE,
                            Q1YE - rhoYE/2 - e2YE h2cYE,
                            e2 h2s,
                            e2 h2s];
     QQ2r = pf * divBlockR1[Q2XE + e2XE h2sXE,
                            Q2YE - e2YE h2sYE,
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

  (* --- [5b] impose the physical cutoff: Q -> W . Q . W ------------------
     W = Kron[Emat, Emat] acting on the grid index of each of the 9 blocks.  Never
     form W: with the flat index k = i + (j-1) nInt (x fastest), an nInt^2 x nInt^2
     block reshapes to the 4-tensor X[[j,i,j',i']] and W.X.W is Emat contracted
     onto each of the four slots.  That is 4 nInt^5 flops per block, 36 nInt^5 in
     all -- 0.04% of the 27 nInt^6 Eigensystem below, i.e. free.

     IN PLACE, ONE BLOCK PAIR AT A TIME.  Qtilde is densified once here (stage [6]
     needs it dense anyway, so this moves that cost rather than adding it) and each
     block is then overwritten via Part assignment.  Peak extra memory is ~3 blocks,
     2.4 GiB each at nInt = 131.  Building all nine blocks and ArrayFlatten-ing them
     would instead hold 9 + 9 = 42 GiB of temporaries, and a whole-matrix
     (Q + Q^T)/2 would add two more 21 GiB copies -- hence the pairwise
     symmetrisation below, which is block-sized and enforces exact global symmetry.
     Folding W into Gx/Dx instead would make them dense and cost O(nInt^6). *)
  If[Emat =!= None,
     smoothBlock[X_] := Module[{T = ArrayReshape[X, {nInt, nInt, nInt, nInt}]},
        T = Emat . T;
        T = Transpose[Emat . Transpose[T, {2, 1, 3, 4}], {2, 1, 3, 4}];
        T = Transpose[Emat . Transpose[T, {3, 2, 1, 4}], {3, 2, 1, 4}];
        T = Transpose[Emat . Transpose[T, {4, 2, 3, 1}], {4, 2, 3, 1}];
        ArrayReshape[T, {nInt^2, nInt^2}]];
     rng[a_] := (a - 1) nInt^2 + 1 ;; a nInt^2;
     Print["[5b] Correlating the noise: Q -> W.Q.W, W = Kron[E,E] ..."];
     tSmooth = rec["[5b] noise correlation W.Q.W (9 blocks, in place)",
       First @ AbsoluteTiming[
        Qasym = 0.;
        Qtilde = Normal @ Qtilde;
        Do[Module[{Sab, Sba, M},
           Sab = smoothBlock[Qtilde[[rng[a], rng[b]]]];
           If[a === b,
              Qasym = Max[Qasym, Max @ Abs[Sab - Transpose[Sab]]];
              Qtilde[[rng[a], rng[a]]] = (Sab + Transpose[Sab])/2,
              (* else: pair (a,b) with (b,a) so the result is exactly symmetric *)
              Sba = smoothBlock[Qtilde[[rng[b], rng[a]]]];
              Qasym = Max[Qasym, Max @ Abs[Sab - Transpose[Sba]]];
              M = (Sab + Transpose[Sba])/2;
              Qtilde[[rng[a], rng[b]]] = M;
              Qtilde[[rng[b], rng[a]]] = Transpose[M]]],
          {a, 3}, {b, a, 3}];]];
     Print["    smoothed in ", tSmooth, " s;  asymmetry before re-symmetrising = ",
           Qasym],
     Print["[5b] lNoise = 0: skipping the noise correlation (W = I)"];
     Qasym = 0.];

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
     Message[SolveFluctuationsLyapunovCutoff::notpsd,
             Count[minEigLoc, x_ /; x < -1.*^-10 scaleLoc]]];

  (* R1: the node-wise 6x6 above is still the RIGHT physical condition -- at symbol
     level the ii terms give 4(cxx sp^2 + cyy sq^2) and the cross term
     8 cxy sp cp sq cq, and 4(cxx sp^2 + cyy sq^2) >= 8 Sqrt[cxx cyy] sp sq by
     AM-GM while cp cq <= 1, so PSD follows from Sqrt[cxx cyy] >= |cxy|, which is
     exactly pointwise PSD of the local flux covariance.  But the mixed scheme is
     no longer a single Gram matrix D M D^T, so check the assembled Qtilde too.
     Arnoldi on the sparse matrix; skipped above n = 20000 where it gets slow.

     NOTE this now runs AFTER [5b], so it tests the SMOOTHED Q_l = W.Q.W.  That is
     the matrix the Lyapunov solve actually uses, and since a congruence preserves
     PSD the check should pass whenever the unsmoothed one would -- W only shrinks
     the spectrum (all its mode weights are in (0,1]), so eigenvalues move toward
     zero and a marginal Q can read as marginally negative at round-off level. *)
  If[3 nInt^2 <= 20000,
     {tQpsd, minEigQ} = AbsoluteTiming[
        Quiet @ Check[
           Min @ Re @ Eigenvalues[N @ Qtilde, -4, Method -> "Arnoldi"],
           Indeterminate]];
     rec["[5] condition (ii): min eig of assembled Q (Arnoldi)", tQpsd];
     Print["    assembled Q min eig = ", minEigQ, "   [", tQpsd, " s]",
           If[NumericQ[minEigQ] && minEigQ < -1.*^-8 Max@Abs@Diagonal@Qtilde,
              "   *** Q IS NOT PSD ***", "   (PSD)"]],
     minEigQ = Missing["skipped, n > 20000"];
     Print["    assembled Q min eig: skipped (n = ", 3 nInt^2, " > 20000)"]];

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
  If[maxReLam >= 0, Message[SolveFluctuationsLyapunovCutoff::unstable, maxReLam]];
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
     Message[SolveFluctuationsLyapunovCutoff::norefine, resid, residRef]];
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
     "lNoise" -> lNoiseV,
     "lNoiseOverH" -> lNoiseV/N[h],
     "Wsymbol" -> Wsym,
     "eigenvalues" -> eigvals,
     "maxReLambda" -> maxReLam,
     "residual" -> resid,
     "residualRefined" -> residRef,
     "minEigLocal" -> MinMax@minEigLoc,
     "minEigQ" -> minEigQ,
     "lapIdentity" -> lapIdentity,
     "QasymBeforeSym" -> Qasym,
     "scheme" -> "cutoff: R1 stencils (adjoint-pair diagonal noise flux, collocated \
cross terms, divergence-form A21/A31) + Gaussian noise correlation length lNoise, \
imposed as Q -> W.Q.W with W = MatrixExp[(lNoise^2/4) Lap] = Kron[E,E]",
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
         ToString @ NumberForm[correlations["h"], 3] <>
         ", l_noise=" <> ToString @ NumberForm[correlations["lNoise"], 4] <>
         " = " <> ToString @ NumberForm[correlations["lNoiseOverH"], 3] <> " h]");

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
       (* lNoise is in modelParams too, but every reader of a sweep wants it
          without unpacking a rule list, and lNoiseOverH is the number that says
          whether the run is grid-limited or converged. *)
       "lNoise" -> correlations["lNoise"],
       "lNoiseOverH" -> correlations["lNoiseOverH"],
       "sigmaRho" -> correlations["sigmaRho"],
       "sigmaQ1" -> correlations["sigmaQ1"],
       "sigmaQ2" -> correlations["sigmaQ2"],
       "maxReLambda" -> correlations["maxReLambda"],
       "residualRefined" -> correlations["residualRefined"],
       "minEigQ" -> correlations["minEigQ"],
       "scheme" -> correlations["scheme"],
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
          "lNoise" -> correlations["lNoise"],
          "lNoiseOverH" -> correlations["lNoiseOverH"],
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
