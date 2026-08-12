#!/usr/bin/env wolframscript
(* ::Package:: *)




(* Fluctuations around the Dirichlet steady state of the rho-Q active nematic:
   discrete Lyapunov solve for psi = (drho, dQ1, dQ2)^T, psi|_bdry = 0.
   Theory: dean_corrected_2.tex, "Linear OU field theory".
   Noise uses the EXACT max-entropy moments; the leading-order closure is not
   PSD above S = 0.618 and this steady state reaches 0.89. *)

Needs["NDSolve`FEM`"];
SetSystemOptions["ParallelOptions" -> "ParallelThreadNumber" -> $ProcessorCount];
Print["ProcessorCount: ", $ProcessorCount];

(* Timing ledger; O(n^3) stages are Eigensystem, back-substitute, refinement,
   and EACH of the three residual evaluations.  n = 3 Nint^2. *)
timings = {};
rec[lbl_, t_] := (AppendTo[timings, {lbl, t}]; t);


(* NotebookDirectory[] fails under wolframscript; $ScriptCommandLine is {} in the
   front end, where Rest[{}] does not evaluate (ListQ[{}] is True, so ListQ misses it). *)
inNotebook = TrueQ[$Notebooks];

scriptDir = If[inNotebook, NotebookDirectory[], DirectoryName[$InputFileName]];
If[!StringQ[scriptDir],
   Print["ERROR: cannot locate the script directory (inNotebook = ",
         inNotebook, "). Save the notebook, or run via wolframscript -file."];
   Abort[]];
repoRoot = ParentDirectory[scriptDir];
repoRoot = scriptDir;

Import[FileNameJoin[{repoRoot, "iterative_solver_module.m"}]];


(* ell and S0 use the RENORMALIZED relations: a -> a + Lambda xir,
   L -> L + zeta xir/(4 xi0).  The bare ones would give ell = 0.80, not 2.5. *)

(* known *)
ellv = 2.5; S0v = 0.9; rho0in = 0.3; xi0in = 5.; zetin = 10.;
(* guessed *)
Bin = 1000.; bin = 1.; gamin = 1.; Lamin = 1.;
(* derived *)
xirin   = 8 bin/gamin;                     (* 8      *)
aRenorm = -2 bin S0v^2;                    (* a + \[CapitalLambda] \[Xi]r = -1.62  *)
LRenorm = ellv^2 (-aRenorm);               (* L + \[Zeta] \[Xi]r/(4 \[Xi]0) = 10.125 *)
ain     = aRenorm - Lamin xirin;           (* -9.62  *)
Lin     = LRenorm - zetin xirin/(4 xi0in); (* 6.125  *)

modelParams = {B -> Bin, a -> ain, b -> bin, \[Xi]0 -> xi0in, \[Xi]r -> xirin,
               \[Zeta] -> zetin, \[Rho]0 -> rho0in, L -> Lin,
               \[CapitalLambda] -> Lamin};

solverParams = {\[Delta]mesh -> 0.061, hMax -> 3, box -> 20,
                maxSteps -> 8, relativeDiff -> 1.*^-4};

Bv      = B                    /. modelParams;
xi0     = \[Xi]0               /. modelParams;
xir     = \[Xi]r               /. modelParams;
zet     = \[Zeta]              /. modelParams;
rho0v   = \[Rho]0              /. modelParams;
Lam     = \[CapitalLambda]     /. modelParams;
av      = a                    /. modelParams;
bv      = b                    /. modelParams;
Lv      = L                    /. modelParams;
cL      = Lv + zet xir/(4 xi0);   (* elastic constant in \[Delta]H *)


Print["[1] Steady state with Dirichlet BC ..."];
{tSS, res} = AbsoluteTiming[
   SolveActiveNematicSteady[modelParams, solverParams, "Dirichlet"]];
rec["[1] steady state (FEM, N-independent)", tSS];
Print["    done in ", tSS, " s"];

rhoSs = res["rho"];
Q1Ss  = res["Q1"];
Q2Ss  = res["Q2"];


(* Options as a single rule-list string, identical in notebook and script:
     wolframscript -file <this> "{fbox -> 20, Nint -> 21, saveData -> True}"
   ClearAll before ToExpression is essential -- these symbols are the working
   variables too, so on a re-run "fbox -> 20" would parse as "20 -> 20". *)
optsString = "{fbox -> 20, Nint -> 21, plotBox -> 8, saveData -> False}";
If[!inNotebook && Length[$ScriptCommandLine] > 1,
   optsString = $ScriptCommandLine[[2]]];

ClearAll[fbox, Nint, plotBox, outputPath, saveData];
optsGiven = ToExpression[optsString];
If[!MatchQ[optsGiven, {(_Rule | _RuleDelayed) ...}],
   Print["ERROR: options must be a list of rules, e.g.\n",
         "  \"{fbox -> 20, Nint -> 21, plotBox -> 8, saveData -> True}\"\n",
         "got: ", optsString];
   Abort[]];
optsKnown = {fbox, Nint, plotBox, outputPath, saveData};
optsBad = Complement[First /@ optsGiven, optsKnown];
If[optsBad =!= {},
   Print["ERROR: unknown option(s) ", optsBad, "; known: ", optsKnown]; Abort[]];

opts = Join[optsGiven,
   {fbox -> 20, Nint -> 21, plotBox -> 8, outputPath -> Directory[],
    saveData -> False}];
{fboxV, nIntV, plotBoxV, outPathV, saveDataV} = optsKnown /. opts;
fbox = fboxV; Nint = nIntV; plotBox = plotBoxV;
outputPath = outPathV; saveData = saveDataV;

If[!IntegerQ[fbox] || !IntegerQ[Nint] || fbox <= 0 || Nint <= 0,
   Print["ERROR: fbox and Nint must be positive integers (h = 2 fbox/(Nint+1) ",
         "is kept exact); got ", {fbox, Nint}]; Abort[]];
If[!NumericQ[plotBox] || plotBox <= 0,
   Print["ERROR: plotBox must be a positive number; got ", plotBox]; Abort[]];
If[!BooleanQ[saveData],
   Print["ERROR: saveData must be True or False; got ", saveData]; Abort[]];
If[!StringQ[outputPath] || !DirectoryQ[outputPath],
   Print["ERROR: outputPath is not an existing directory: ", outputPath]; Abort[]];

Print["    ", If[inNotebook, "NOTEBOOK", "SCRIPT"], " mode; fbox=", fbox,
      ", Nint=", Nint, ", plotBox=", plotBox, ", saveData=", saveData];
Print["    outputPath = ", outputPath];

(* x_k = -fbox + k h needs k = (Nint+1)/2 to hit 0: ODD Nint puts the origin
   (defect core, S = 0) on the grid, EVEN Nint does not.  h is kept exact so
   that node is identically 0. *)
h    = (2 fbox)/(Nint + 1);
gridFull = Table[-fbox + k h, {k, 0, Nint + 1}];
gridInt  = gridFull[[2 ;; -2]];
Print["[2] Grid: fbox=", fbox, ", Nint=", Nint,
      ", h=", ToString @ NumberForm[N[h], 5],
      ", state-dim=", 3 Nint^2];
If[OddQ[Nint],
   iZero = (Nint + 1)/2;
   If[gridInt[[iZero]] =!= 0,
      Print["WARNING: Nint odd but gridInt[[", iZero, "]] = ", gridInt[[iZero]],
            " is not exactly 0."],
      Print["    Nint odd  -> origin ON the grid at gridInt[[", iZero,
            "]] = 0, flat index k = ", iZero + (iZero - 1) Nint,
            " of ", Nint^2]],
   iZero = None;
   Print["    Nint even -> origin NOT on the grid; nearest nodes at x = ",
         ToString @ NumberForm[N @ gridInt[[Nint/2]], 4], " and ",
         ToString @ NumberForm[N @ gridInt[[Nint/2 + 1]], 4]]];


(* 1D FD matrices; Dirichlet = entries outside the interior are dropped. *)
d1Full = NDSolve`FiniteDifferenceDerivative[Derivative[1], gridFull,
   "DifferenceOrder" -> 2]["DifferentiationMatrix"];
l1Full = NDSolve`FiniteDifferenceDerivative[Derivative[2], gridFull,
   "DifferenceOrder" -> 2]["DifferentiationMatrix"];
D1 = SparseArray @ d1Full[[2 ;; -2, 2 ;; -2]];
L1 = SparseArray @ l1Full[[2 ;; -2, 2 ;; -2]];

(* State ordering k = i + (j-1) N, x fastest.  Hence Dx = I (x) D1, Dy = D1 (x) I. *)
Idn = IdentityMatrix[Nint, SparseArray];
Dx  = KroneckerProduct[Idn, D1];
Dy  = KroneckerProduct[D1, Idn];
Lx  = KroneckerProduct[Idn, L1];
Ly  = KroneckerProduct[L1, Idn];
Dxy = KroneckerProduct[D1, D1];
Lap = Lx + Ly;
Ig  = IdentityMatrix[Nint^2, SparseArray];


(* Batched sampling, ~200x faster than scalar calls.  N@ is load-bearing:
   exact rationals bypass the fast path and are 19x slower. *)
xflat = Flatten @ ConstantArray[N @ gridInt, Nint];
yflat = Flatten @ Transpose @ ConstantArray[N @ gridInt, Nint];
sample[f_] := f[xflat, yflat];

Print["[3] Sampling steady state on interior grid ..."];
tSamp = rec["[3] sample steady state", First @ AbsoluteTiming[
   rhoV   = sample[ rhoSs ];
   Q1V    = sample[ Q1Ss  ];
   Q2V    = sample[ Q2Ss  ];
   rhoLap = sample[ Derivative[2, 0][rhoSs][##] + Derivative[0, 2][rhoSs][##] & ];
   Q1x    = sample[ Derivative[1, 0][Q1Ss][##] & ];
   Q1y    = sample[ Derivative[0, 1][Q1Ss][##] & ];
   Q2x    = sample[ Derivative[1, 0][Q2Ss][##] & ];
   Q2y    = sample[ Derivative[0, 1][Q2Ss][##] & ];]];
Print["    sampling (8 batched FEM calls): ", tSamp, " s"];

Print["    \[Rho]_ss \[Element] ", MinMax@rhoV,
      "   Q1_ss \[Element] ", MinMax@Q1V,
      "   Q2_ss \[Element] ", MinMax@Q2V];


(* g2 = I0 I2/I1^2, g3 = I0^2 I3/I1^3, kappa fixed pointwise by I1/I0 = S.
   Below S = 0.02 both are numerically 0/0, so use the series.
   Local iterator is kk -- kap is the noise prefactor below. *)
kapOfS[s_?NumericQ] := If[s < 1.*^-6, 0.,
   kk /. Quiet @ FindRoot[BesselI[1, kk]/BesselI[0, kk] == s,
                          {kk, 2 s/(1 - s^2)}, MaxIterations -> 200]];
g2ex[s_?NumericQ] := If[s < 0.02, 1/2 + s^2/6 + 11 s^4/180,
   With[{k = kapOfS[s]}, BesselI[0, k] BesselI[2, k]/BesselI[1, k]^2]];
g3ex[s_?NumericQ] := If[s < 0.02, 1/6 + s^2/8 + 23 s^4/480,
   With[{k = kapOfS[s]}, BesselI[0, k]^2 BesselI[3, k]/BesselI[1, k]^3]];

Svals = Sqrt[Q1V^2 + Q2V^2]/rhoV;
{tClo, {g2V, g3V}} = AbsoluteTiming[{g2ex /@ Svals, g3ex /@ Svals}];
rec["[3] kappa(S) inversion + g2,g3", tClo];
(* kapOfS is not memoized, so this repeats all Nint^2 root-finds; diagnostic only. *)
{tKchk, kapResid} = AbsoluteTiming[
   Max @ Abs[(BesselI[1, #]/BesselI[0, #] & /@ (kapOfS /@ Svals)) - Svals]];
rec["[3] kappa(S) verification (diagnostic)", tKchk];
Print["    S = |Q_ss|/\[Rho]_ss \[Element] ", MinMax@Svals, ",  mean ", Mean@Svals];
Print["    g2 \[Element] ", MinMax@g2V, "   g3 \[Element] ", MinMax@g3V,
      "   (small-S limits 1/2, 1/6)"];
Print["    \[Kappa](S) inversion: ", tClo, " s, max |I1/I0 - S| = ", kapResid];
If[kapResid > 1.*^-8,
   Print["WARNING: kappa(S) root-find did not converge at some grid points."]];


(* Q = {{Q1,Q2},{Q2,-Q1}}.  The comoving term cancels the one from
   Div(dQ Grad rho_ss).  Sign a -> a + Lambda xir: rotational diffusion
   DESTROYS order; dean.tex had it the other way. *)
Print["[4] Assembling A ..."];
diag[c_] := DiagonalMatrix[SparseArray @ c];

tAsmA = rec["[4] assemble A", First @ AbsoluteTiming[
   A11 = (Bv/xi0) Lap;
   A12 = (zet/xi0) (Lx - Ly);
   A13 = (2 zet/xi0) Dxy;

   A21 = (Bv/xi0) (Q1V Lap + Q1x Dx + Q1y Dy);
   A22 = (Bv/xi0) diag[rhoLap] +
         (4/xir) ( diag[ -(av + Lam xir) - bv (6 Q1V^2 + 2 Q2V^2) ] + cL Lap );
   A23 = -(16 bv/xir) diag[ Q1V Q2V ];

   A31 = (Bv/xi0) (Q2V Lap + Q2x Dx + Q2y Dy);
   A32 = A23;
   A33 = (Bv/xi0) diag[rhoLap] +
         (4/xir) ( diag[ -(av + Lam xir) - bv (2 Q1V^2 + 6 Q2V^2) ] + cL Lap );

   Atilde = ArrayFlatten[{{A11, A12, A13},
                          {A21, A22, A23},
                          {A31, A32, A33}}];]];
Print["    A dims = ", Dimensions[Atilde], ",  assembled in ", tAsmA, " s"];


(* Noise f = d_i F_i, <F_i F_j> = C_ij delta(r-r') -> (1/h^2) sum_ij D_i.diag(C_ij).D_j^T.
   EXACT max-entropy moments \[Eqref]{exact-moments}:
     <QQ> = rho J + (g2/2rho) U,  <Q'Q'> = 4 rho J - (2 g2/rho) U,
     <QQQ> = (1/2) T + (g3/rho^2) W;   g2 = g3 = 0 gives the old closure.
   Verified vs von Mises quadrature to 1e-16..4e-13 in
   comparison_draft/{derive_exact_blocks,check_exact_blocks}.wls. *)
Print["[5] Assembling Q ..."];
kap = 2 zet / (xi0 rho0v);
pf  = kap / h^2;
locPf = 2 Lam / (rho0v h^2);
divBlock[cxx_, cyy_, cxy_, cyx_] :=
  Dx . (cxx * Transpose[Dx]) + Dy . (cyy * Transpose[Dy]) +
  Dx . (cxy * Transpose[Dy]) + Dy . (cyx * Transpose[Dx]);

tAsmQ = rec["[5] assemble Q (24 sparse products)", First @ AbsoluteTiming[
h2c = Q1V^2 - Q2V^2;          h2s = 2 Q1V Q2V;
h3c = Q1V^3 - 3 Q1V Q2V^2;    h3s = 3 Q1V^2 Q2V - Q2V^3;
e2  = g2V/(2 rhoV);           (* U terms *)
e3  = g3V/(4 rhoV^2);         (* W terms *)

(* rho-rho: no closure enters *)
Qrr = pf * divBlock[Q1V + rhoV, -Q1V + rhoV, Q2V, Q2V];

QQ1Q1 = (pf * divBlock[ 3 Q1V/4 + rhoV/2 + e2 h2c + e3 h3c,
                       -3 Q1V/4 + rhoV/2 + e2 h2c - e3 h3c,
                        Q2V/4 + e3 h3s,
                        Q2V/4 + e3 h3s]
         + locPf * diag[ 2 rhoV - 4 e2 h2c ]);

QQ2Q2 = (pf * divBlock[ Q1V/4 + rhoV/2 - e2 h2c - e3 h3c,
                       -Q1V/4 + rhoV/2 - e2 h2c + e3 h3c,
                        3 Q2V/4 - e3 h3s,
                        3 Q2V/4 - e3 h3s]
         + locPf * diag[ 2 rhoV + 4 e2 h2c ]);

(* NB Q1-Q2 now carries a LOCAL part: <Q'Q'> is no longer proportional to J. *)
QQ1Q2 = (pf * divBlock[ Q2V/4 + e2 h2s + e3 h3s,
                       -Q2V/4 + e2 h2s - e3 h3s,
                        Q1V/4 - e3 h3c,
                        Q1V/4 - e3 h3c]
         + locPf * diag[ -4 e2 h2s ]);
QQ2Q1 = Transpose[QQ1Q2];

QQ1r = pf * divBlock[ Q1V + rhoV/2 + e2 h2c,
                      Q1V - rhoV/2 - e2 h2c,
                      e2 h2s,
                      e2 h2s];
QQ2r = pf * divBlock[ Q2V + e2 h2s,
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


(* Local 6x6 flux covariance, ordered (rho x, rho y, Q1 x, Q1 y, Q2 x, Q2 y);
   its PSD implies PSD of the divergence part.  TWO EXACT ZEROS are expected:
   Qhat(theta).p = p forces two relations among the fluxes, so rank <= 4 -- the
   exact M sits ON the cone boundary, which is why the old closure fell off it. *)
localC[r_, q1_, q2_, g2_, g3_] := Module[
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

tCond = rec["[5] condition (ii): 2 N^2 local 6x6 eigensolves", First @ AbsoluteTiming[
   minEigLoc = Min /@ Eigenvalues /@ MapThread[localC, {rhoV, Q1V, Q2V, g2V, g3V}];
   minEigOld = Min /@ Eigenvalues /@
      MapThread[localC[#1, #2, #3, 0, 0] &, {rhoV, Q1V, Q2V}];]];
scaleLoc  = Max @ Abs @ Flatten @ {rhoV, Q1V, Q2V};
Print["    condition (ii) check: ", tCond, " s"];
Print["    local 6x6 min eig (exact closure)  \[Element] ", MinMax@minEigLoc,
      ";  negative at ", Count[minEigLoc, x_ /; x < -1.*^-10 scaleLoc],
      "/", Length@minEigLoc, " grid points"];
Print["    local 6x6 min eig (old O(S) one)   \[Element] ", MinMax@minEigOld,
      ";  negative at ", Count[minEigOld, x_ /; x < -1.*^-10 scaleLoc],
      "/", Length@minEigOld, " grid points   [S > ",
      ToString@NumberForm[N[(Sqrt[5] - 1)/2], 4], "]"];
If[Min@minEigLoc < -1.*^-10 scaleLoc,
   Print["    WARNING: condition (ii) Q >= 0 is violated even with the exact ",
         "moments -- this should not happen; check the closure assembly."]];


(* LyapunovSolve is single-threaded; diagonalize instead: A = P D P^-1 =>
   C = P Y P^T, Y_ij = -(P^-1 Q P^-T)_ij/(lam_i + lam_j).  A is NOT self-adjoint
   (dropped Onsager partners), but the error is set by kappa(P)^2 eps ~ 1e-10,
   not by the asymmetry; one refinement step then reaches ~1e-14.  The printed
   residual is the verification -- if it stops dropping, fall back. *)
Print["[6] Lyapunov solve via eigendecomposition + refinement ..."];
(* memory-bandwidth bound, 8 n^2 bytes each *)
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
If[maxReLam >= 0,
   Print["WARNING: linearization not strictly stable; results may be unphysical."]];
Print["    min |\[Lambda]_i+\[Lambda]_j| = ", 2 Min @ Abs @ Re @ eigvals,
      ",  complex fraction = ",
      N[100 Count[eigvals, z_ /; Abs[Im[z]] > 1.*^-8 Abs[z]]/Length[eigvals]],
      "%"];

Pev = Transpose @ eigvecs;


(* Solve A X + X A^T = -Rhs in the eigenbasis; reused by the refinement. *)
eigLyapSolve[Rhs_] := Re[Pev . (
   -LinearSolve[Pev, Transpose @ LinearSolve[Pev, Transpose[Rhs]]] /
      Outer[Plus, eigvals, eigvals]) . Transpose[Pev]];

{tSol, Cmat} = AbsoluteTiming @ eigLyapSolve[Qdense];
rec["[6] back-substitute", tSol];
Print["    [b] back-substitute    ", tSol, " s"];


(* Each residual is two dense n x n matmuls, same O(n^3) as the Eigensystem,
   evaluated three times: here, inside the refinement, and to verify. *)
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
Print["    residual (refined)     = ", residRef,
      "   [", tRes3, " s]",
      "   (improved ", ToString@ScientificForm[N[resid/residRef], 3], "x)"];
If[residRef > resid,
   Print["WARNING: refinement did not improve the residual; ",
         "consider falling back to LyapunovSolve."]];

Print["    total solve time = ", tEs + tSol + tRef, " s"];


Print["[7] Extracting Sigma_rho_rho and plotting ..."];
(* ListInterpolation, NOT {x,y,z} triples: the triple form makes ListDensityPlot
   triangulate (1.5 s vs 0.44 s at Nint=151).  Wants arr[[i,j]] = f(x_i,y_j). *)
tExtr = rec["[7] extract diagonal + build interpolant", First @ AbsoluteTiming[
   sigRhoVec = Diagonal @ Cmat[[1 ;; Nint^2, 1 ;; Nint^2]];
   sigRhoMat = Partition[sigRhoVec, Nint];  (* sigRhoMat[[j,i]] = Sigma(x_i, y_j) *)

   sigRhoFn = ListInterpolation[Transpose[sigRhoMat], {gridInt, gridInt},
                                InterpolationOrder -> 2];]];
Print["    extract + interpolant: ", tExtr, " s"];




(* plotBox crops the Dirichlet boundary layer out of the figures; the
   interpolant only spans +-(fbox - h), so clamp to that to avoid extrapolation. *)
{xLo, xHi} = {Max[-plotBox, First@gridInt], Min[plotBox, Last@gridInt]};
plotPts = Min[Nint, 200];
lbl = ("Var(\[Delta]\[Rho])(r)  [Dirichlet BC, h="
       <> ToString @ NumberForm[N[h], 3] <> "]");

   plot2D = DensityPlot[sigRhoFn[x, y], {x, xLo, xHi}, {y, xLo, xHi},
      PlotPoints -> plotPts, MaxRecursion -> 0, PlotRange -> {{-plotBox,plotBox},{-plotBox,plotBox}},
      PlotLegends -> Automatic, FrameLabel -> {"x", "y"},
      PlotLabel -> lbl, ColorFunction -> "SunsetColors"];

   plot3D = Plot3D[sigRhoFn[x, y], {x, xLo, xHi}, {y, xLo, xHi},
      PlotPoints -> plotPts, MaxRecursion -> 0, Mesh -> None, PlotRange -> {{-plotBox,plotBox},{-plotBox,plotBox},Automatic},
      AxesLabel -> {"x", "y", "Var(\[Delta]\[Rho])"},
      PlotLabel -> lbl, ColorFunction -> "SunsetColors"];


(* Displayed in the notebook only.  A bare "plot3D" would make wolframscript
   print the whole Graphics expression -- every mesh vertex -- to stdout. *)
If[inNotebook, plot3D]

If[inNotebook, plot2D]


tPlot = rec["[7] plot + export", First @ AbsoluteTiming[
   tag = "box" <> ToString[fbox] <> "_N" <> ToString[Nint];
   out[stem_, ext_] := FileNameJoin[{outputPath, stem <> tag <> ext}];
   outPath   = out["sigma_rho_dirichlet_", ".png"];
   outPath3D = out["sigma_rho_surface_",   ".png"];
   Export[outPath,   plot2D, ImageResolution -> 150];
   Export[outPath3D, plot3D, ImageResolution -> 150];

   dumpPath = out["sigma_rho_", ".m"];
   Export[dumpPath,
      <|"fbox" -> fbox, "Nint" -> Nint, "h" -> N[h],
        "gridInt" -> gridInt, "sigmaRho" -> sigRhoMat,
        "maxReLambda" -> maxReLam|>];]];
Print["Saved: ", outPath];
Print["Saved: ", outPath3D];
Print["Saved: ", dumpPath];
Print["    plot + export: ", tPlot, " s   (PlotPoints -> ", plotPts, ")"];

(* Full covariance, row-major Real64, n = 3 Nint^2 with blocks (drho, dQ1, dQ2).
   Read back with ArrayReshape[BinaryReadList[f, "Real64"], {n, n}].  8 n^2
   bytes -- 487 MB at Nint = 51, hence off by default. *)
If[TrueQ[saveData],
   binPath = out["covariance_", ".bin"];
   tBin = rec["[7] save covariance (binary)", First @ AbsoluteTiming[
      Export[binPath, Cmat, {"Binary", "Real64"}];
      Export[out["covariance_", "_meta.m"],
         <|"n" -> Length[Cmat], "Nint" -> Nint, "fbox" -> fbox, "h" -> N[h],
           "gridInt" -> gridInt, "blocks" -> {"rho", "Q1", "Q2"},
           "ordering" -> "row-major Real64; grid index k = i + (j-1) Nint, x fastest"|>];]];
   Print["Saved: ", binPath, "  (",
         ToString@NumberForm[N[FileByteCount[binPath]/2^20], {6, 1}], " MiB, ",
         tBin, " s)"]];

Print["Sigma_rho stats: min=", Min@sigRhoVec, "  max=", Max@sigRhoVec,
      "  mean=", Mean@sigRhoVec];


(* Plain text: Grid/NumberForm do not render under wolframscript. *)
fmtS[x_] := StringPadLeft[ToString@NumberForm[N@x, {7, 2}], 8];
Print[];
Print["=== timing breakdown (Nint = ", Nint, ", n = ", 3 Nint^2, ") ==="];
totT = Total[timings[[All, 2]]];
Print[StringPadRight["stage", 46], fmtS["  s"], "   %tracked"];
Print[StringJoin@ConstantArray["-", 68]];
Do[Print[StringPadRight[r[[1]], 46], fmtS[r[[2]]],
         StringPadLeft[ToString@NumberForm[N[100 r[[2]]/totT], {5, 1}], 8],
         If[r[[2]] > 4, "  <== >4 s", ""]],
   {r, SortBy[timings, -#[[2]] &]}];
Print[StringJoin@ConstantArray["-", 68]];
Print[StringPadRight["TOTAL tracked", 46], fmtS[totT], StringPadLeft["100.0", 8]];
slow = Select[timings, #[[2]] > 4 &][[All, 1]];
Print["Stages above 4 s: ",
      If[slow === {}, "none", StringRiffle[slow, ", "]]];
