#!/usr/bin/env wolframscript
(* -------------------------------------------------------------------
   Fluctuations around the Dirichlet steady state of the coupled
   \[Rho]-Q active-nematic system.

   Implements the discrete Lyapunov procedure from dean.tex
   "Linear OU field theory" (sec.~702ff), for state
       \[Psi] = (\[Delta]\[Rho], \[Delta]Q1, \[Delta]Q2)^T,
   with Dirichlet BC \[Psi]|_boundary = 0.

   Steady state is solved on the large box (via iterative_solver_module.m
   with bcType = "Dirichlet"). The fluctuation problem is set on a
   smaller inner box discretized on a uniform NxN interior grid.
   The steady-state InterpolatingFunction is evaluated (with its exact
   derivatives) on that grid.

   Output:  sigma_rho_rho.png  ---  ListDensityPlot of Var(\[Delta]\[Rho])(r).
   -------------------------------------------------------------------- *)

Needs["NDSolve`FEM`"];
SetSystemOptions["ParallelOptions" -> "ParallelThreadNumber" -> 10];

(* ---------- locate parent repo ----------------------------------------- *)
scriptDir = If[TrueQ[$InputFileName == "" || $InputFileName === $Failed],
                NotebookDirectory[],
                DirectoryName[$InputFileName]];
repoRoot = ParentDirectory[scriptDir];

Import[FileNameJoin[{repoRoot, "iterative_solver_module.m"}]];

(* ---------- model + solver parameters (from fluctuations.wls) ---------- *)
preParams   = {ell -> 2.5, S0 -> 0.9, A -> -1};
modelParams = {B -> 1000., a -> A, b -> -A/(2 S0^2), \[Xi]0 -> 5.,
               \[Xi]r -> 0.1, \[Zeta] -> 10., \[Rho]0 -> 1.,
               L -> ell^2 (-A), \[CapitalLambda] -> 0.1} /. preParams;

solverParams = {\[Delta]mesh -> 0.061, hMax -> 3, box -> 20,
                maxSteps -> 8, relativeDiff -> 1.*^-4};

(* Numeric parameter values that appear in the linearized equations. *)
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

(* ---------- steady state (Dirichlet) ----------------------------------- *)
Print["[1] Steady state with Dirichlet BC ..."];
{tSS, res} = AbsoluteTiming[
   SolveActiveNematicSteady[modelParams, solverParams, "Dirichlet"]];
Print["    done in ", tSS, " s"];

rhoSs = res["rho"];
Q1Ss  = res["Q1"];
Q2Ss  = res["Q2"];

(* ---------- fluctuation grid: uniform, interior only ------------------- *)
(* Override from the command line:  wolframscript -file <this> <fbox> <Nint> *)
cliArgs = If[ListQ[$ScriptCommandLine], Rest[$ScriptCommandLine], {}];
fbox = If[Length[cliArgs] >= 1, ToExpression[cliArgs[[1]]], 10];
Nint = If[Length[cliArgs] >= 2, ToExpression[cliArgs[[2]]], 30];
h    = (2 fbox)/(Nint + 1);
gridFull = Table[-fbox + k h, {k, 0, Nint + 1}];
gridInt  = gridFull[[2 ;; -2]];
Print["[2] Grid: fbox=", fbox, ", Nint=", Nint,
      ", h=", ToString @ NumberForm[N[h], 5],
      ", state-dim=", 3 Nint^2];

(* ---------- 1D FD matrices (Dirichlet: entries outside dropped) -------- *)
d1Full = NDSolve`FiniteDifferenceDerivative[Derivative[1], gridFull,
   "DifferenceOrder" -> 2]["DifferentiationMatrix"];
l1Full = NDSolve`FiniteDifferenceDerivative[Derivative[2], gridFull,
   "DifferenceOrder" -> 2]["DifferentiationMatrix"];
D1 = SparseArray @ d1Full[[2 ;; -2, 2 ;; -2]];
L1 = SparseArray @ l1Full[[2 ;; -2, 2 ;; -2]];

(* 2D Kronecker operators.  State ordering: k = i + (j-1) N (i-fastest, x-fastest).
   With this convention:  Dx = I \[CircleTimes] D1, Dy = D1 \[CircleTimes] I,
                          Lx = I \[CircleTimes] L1, Ly = L1 \[CircleTimes] I. *)
Idn = IdentityMatrix[Nint, SparseArray];
Dx  = KroneckerProduct[Idn, D1];
Dy  = KroneckerProduct[D1, Idn];
Lx  = KroneckerProduct[Idn, L1];
Ly  = KroneckerProduct[L1, Idn];
Dxy = KroneckerProduct[D1, D1];
Lap = Lx + Ly;
Ig  = IdentityMatrix[Nint^2, SparseArray];

(* ---------- sample steady state on interior grid ----------------------- *)
(* Vectorized sampling: one batched call per field instead of Nint^2 scalar
   calls.  FEM InterpolatingFunctions accept array arguments and do the mesh
   point-location in a single compiled pass -- bit-identical results, ~200x
   faster (measured 155x at Nint=30, 251x at Nint=60; see
   draft_fast_sampling.wls).  Two traps:
     - N@ is load-bearing: exact rationals bypass the fast path (19x slower
       for the same answer), and gridInt is exact since h = 2 fbox/(Nint+1).
     - ParallelTable is SLOWER here (0.63x): every kernel has to be sent the
       whole InterpolatingFunction, and per-call overhead dominates.
   Ordering matches k(i,j) = i + (j-1) Nint, i.e. x fastest. *)
xflat = Flatten @ ConstantArray[N @ gridInt, Nint];
yflat = Flatten @ Transpose @ ConstantArray[N @ gridInt, Nint];
sample[f_] := f[xflat, yflat];

Print["[3] Sampling steady state on interior grid ..."];
rhoV   = sample[ rhoSs ];
Q1V    = sample[ Q1Ss  ];
Q2V    = sample[ Q2Ss  ];
rhoLap = sample[ Derivative[2, 0][rhoSs][##] + Derivative[0, 2][rhoSs][##] & ];
Q1x    = sample[ Derivative[1, 0][Q1Ss][##] & ];
Q1y    = sample[ Derivative[0, 1][Q1Ss][##] & ];
Q2x    = sample[ Derivative[1, 0][Q2Ss][##] & ];
Q2y    = sample[ Derivative[0, 1][Q2Ss][##] & ];

Print["    \[Rho]_ss \[Element] ", MinMax@rhoV,
      "   Q1_ss \[Element] ", MinMax@Q1V,
      "   Q2_ss \[Element] ", MinMax@Q2V];

(* Convenience: multiply "diag(c) . operator" via broadcast; result is Sparse. *)

(* ---------- assemble A ------------------------------------------------- *)
(* Derivation notes (see conversation).
   Row 1 (\[Delta]\[Rho]):
       A11 = (B/xi0) Lap
       A12 = (zet/xi0)(Lx - Ly)          [from (\[Partial]x^2 - \[Partial]y^2)\[Delta]Q1]
       A13 = (2 zet/xi0) Dxy             [from 2 \[Partial]x\[Partial]y \[Delta]Q2]

   Row 2 (\[Delta]Q1):
       A21 = (B/xi0)( Q1_ss Lap + Q1x Dx + Q1y Dy )
       A22 = (B/xi0) diag(\[Del]^2 \[Rho]_ss)
             + (4/xir)( -(a - \[Lambda] xir) - b(6 Q1_ss^2 + 2 Q2_ss^2) + cL Lap )
       A23 = -(16 b/xir) diag(Q1_ss Q2_ss)

   Row 3 (\[Delta]Q2):
       A31 = (B/xi0)( Q2_ss Lap + Q2x Dx + Q2y Dy )
       A32 = A23
       A33 = (B/xi0) diag(\[Del]^2 \[Rho]_ss)
             + (4/xir)( -(a - \[Lambda] xir) - b(2 Q1_ss^2 + 6 Q2_ss^2) + cL Lap )
*)

Print["[4] Assembling A ..."];
diag[c_] := DiagonalMatrix[SparseArray @ c];

A11 = (Bv/xi0) Lap;
A12 = (zet/xi0) (Lx - Ly);
A13 = (2 zet/xi0) Dxy;

A21 = (Bv/xi0) (Q1V Lap + Q1x Dx + Q1y Dy);
A22 = (Bv/xi0) diag[rhoLap] +
      (4/xir) ( diag[ -(av - Lam xir) - bv (6 Q1V^2 + 2 Q2V^2) ] + cL Lap );
A23 = -(16 bv/xir) diag[ Q1V Q2V ];

A31 = (Bv/xi0) (Q2V Lap + Q2x Dx + Q2y Dy);
A32 = A23;
A33 = (Bv/xi0) diag[rhoLap] +
      (4/xir) ( diag[ -(av - Lam xir) - bv (2 Q1V^2 + 6 Q2V^2) ] + cL Lap );

Atilde = ArrayFlatten[{{A11, A12, A13},
                       {A21, A22, A23},
                       {A31, A32, A33}}];
Print["    A dims = ", Dimensions[Atilde]];

(* ---------- assemble Q (noise covariance) ------------------------------ *)
(* For a divergence-form noise  f = \[Partial]_i F_i  with  <F_i F_j> = C_ij \[Delta](r-r'),
   the discrete block is  (1/h^2) sum_{ij} D_i . diag(C_ij) . D_j^T.

   Local (\[Lambda]) part gives 2 B_loc \[Delta]_kk'/h^2 on the diagonal.

   Flux covariances C^{ab, cd}_ij for the three tensor pairs (11,11),
   (12,12), (11,12) come from  T_{ab mn ij}(Q_ss) + \[Rho]_ss J_{ab mn} \[Delta]_ij,
   and for the \[Rho]-Q cross pairs from  \[Rho]_ss J_{ab ij} + Q^ss_{ab} \[Delta]_ij.
   Overall prefactor  \[Kappa] = 2 \[Zeta] / (\[Xi]0 \[Rho]0). *)

Print["[5] Assembling Q ..."];
kap = 2 zet / (xi0 rho0v);    (* continuous coefficient *)
pf  = kap / h^2;              (* discrete divergence prefactor *)
locFactor = 4 Lam / (rho0v h^2);   (* 2 B_loc / h^2 with B_loc = 2\[Lambda]\[Rho]/\[Rho]_0 *)

(* helper: assemble D_i . diag(c_ij) . D_j^T.
   Uses broadcast identity  diag(c) . M  =  c * M  (row-wise), sparse-safe. *)
divBlock[cxx_, cyy_, cxy_, cyx_] :=
  Dx . (cxx * Transpose[Dx]) + Dy . (cyy * Transpose[Dy]) +
  Dx . (cxy * Transpose[Dy]) + Dy . (cyx * Transpose[Dx]);

(* rho-rho  (no local part) *)
Qrr = pf * divBlock[Q1V + rhoV, -Q1V + rhoV, Q2V, Q2V];

(* Q1-Q1 *)
QQ1Q1 = (pf * divBlock[3 Q1V/2 + rhoV/2, -3 Q1V/2 + rhoV/2, Q2V/2, Q2V/2]
         + locFactor * diag[ rhoV ]);

(* Q2-Q2 *)
QQ2Q2 = (pf * divBlock[Q1V/2 + rhoV/2, -Q1V/2 + rhoV/2, 3 Q2V/2, 3 Q2V/2]
         + locFactor * diag[ rhoV ]);

(* Q1-Q2 cross (no local part; J_{1112} = 0) *)
QQ1Q2 = pf * divBlock[Q2V/2, -Q2V/2, Q1V/2, Q1V/2];
QQ2Q1 = Transpose[QQ1Q2];

(* Q_ab - rho cross:
     C^{Q_ab, rho}_ij = \[Kappa] ( \[Rho]_ss J_{ab ij} + Q^ss_{ab} \[Delta]_ij ).
   For ab = 11:  diag(1/2,-1/2) J-part, plus Q1 * I.
   For ab = 12:  off-diag 1/2, plus Q2 * I.                                     *)
QQ1r = pf * divBlock[rhoV/2 + Q1V, -rhoV/2 + Q1V, 0 rhoV, 0 rhoV];
QQ2r = pf * divBlock[Q2V, Q2V, rhoV/2, rhoV/2];
QrQ1 = Transpose[QQ1r];
QrQ2 = Transpose[QQ2r];

Qtilde = ArrayFlatten[{{Qrr,  QrQ1,  QrQ2},
                       {QQ1r, QQ1Q1, QQ1Q2},
                       {QQ2r, QQ2Q1, QQ2Q2}}];
Qtilde = (Qtilde + Transpose[Qtilde])/2;    (* symmetrize *)
Print["    Q dims = ", Dimensions[Qtilde]];

(* ---------- Lyapunov solve  A C + C A^T = -Q --------------------------- *)
(* Mma's LyapunovSolve (Bartels-Stewart) is single-threaded and dominates the
   wall time.  We diagonalize instead:  A = P Diag P^-1  =>  C = P Y P^T with
   Y_{ij} = -(P^-1 Q P^-T)_{ij}/(\[Lambda]_i + \[Lambda]_j).  All heavy steps
   (Eigensystem, LinearSolve, Dot) go through multithreaded MKL.

   CAVEAT on accuracy.  A is NOT self-adjoint: the \[Rho]<->Q coupling blocks
   violate Onsager reciprocity because two reciprocal terms were dropped in
   dean.tex (the \[Zeta] Grad Grad \[Rho] term in the Q equation, by the B>>1
   approximation at L598; and Div(Q:Grad H) in the \[Rho] equation, to avoid
   4th-order derivatives, at L625).  See diagnose_symmetry.wls /
   verify_onsager.wls -- the skew part is ~35% of the symmetric part, and the
   one-sided \[Zeta] coupling is what makes ~18% of the spectrum complex.

   What controls the error is NOT the asymmetry but kappa(P), the eigenvector
   condition number, via  err ~ kappa(P)^2 eps.  Measured: kappa(P) ~ 3e2-8e2
   with no systematic growth in N, giving ~1e-10 -- matching the observed
   residual.  Two things keep this safe:
     (i)  min|\[Lambda]_i+\[Lambda]_j| >= 2 min|Re \[Lambda]| ~ 7, so no
          near-cancellation in the denominators;
     (ii) ONE step of iterative refinement (below) reuses the cached
          eigendecomposition and drives the residual to ~1e-14, i.e. ~30x
          BETTER than LyapunovSolve, at a fraction of its cost.
   The refinement residual printed below is the verification -- if it ever
   fails to drop, fall back to LyapunovSolve. *)
Print["[6] Lyapunov solve via eigendecomposition + refinement ..."];
Adense = Normal @ Atilde;
Qdense = Normal @ Qtilde;

{tEs, {eigvals, eigvecs}} = AbsoluteTiming @ Eigensystem[Adense];
Print["    [a] Eigensystem        ", tEs, " s"];

(* Stability check, free from the full spectrum. *)
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

(* Solve  A X + X A^T = -Rhs  in the eigenbasis (reused by refinement). *)
eigLyapSolve[Rhs_] := Re[Pev . (
   -LinearSolve[Pev, Transpose @ LinearSolve[Pev, Transpose[Rhs]]] /
      Outer[Plus, eigvals, eigvals]) . Transpose[Pev]];

{tSol, Cmat} = AbsoluteTiming @ eigLyapSolve[Qdense];
Print["    [b] back-substitute    ", tSol, " s"];

resid = Max @ Abs[Adense . Cmat + Cmat . Transpose[Adense] + Qdense];
Print["    residual (unrefined)   = ", resid];

(* Iterative refinement: solve A dC + dC A^T = -R for the residual R and
   correct.  Cheap -- the eigendecomposition is already cached. *)
{tRef, Cmat} = AbsoluteTiming[
   Cmat + eigLyapSolve[Adense . Cmat + Cmat . Transpose[Adense] + Qdense]];
residRef = Max @ Abs[Adense . Cmat + Cmat . Transpose[Adense] + Qdense];
Print["    [c] refinement step    ", tRef, " s"];
Print["    residual (refined)     = ", residRef,
      "   (improved ", ToString@ScientificForm[N[resid/residRef], 3], "x)"];
If[residRef > resid,
   Print["WARNING: refinement did not improve the residual; ",
         "consider falling back to LyapunovSolve."]];

Print["    total solve time = ", tEs + tSol + tRef, " s"];

(* ---------- variance of \[Delta]\[Rho] and plot ------------------------ *)
Print["[7] Extracting Sigma_rho_rho and plotting ..."];
sigRhoVec = Diagonal @ Cmat[[1 ;; Nint^2, 1 ;; Nint^2]];
sigRhoMat = Partition[sigRhoVec, Nint];  (* sigRhoMat[[j, i]] = Sigma(x_i, y_j) *)

data = Flatten[
   Table[{gridInt[[i]], gridInt[[j]], sigRhoMat[[j, i]]},
         {j, Nint}, {i, Nint}], 1];

plot = ListDensityPlot[data,
   PlotLegends -> Automatic,
   FrameLabel -> {"x", "y"},
   PlotLabel  -> "Var(\[Delta]\[Rho])(r)  [Dirichlet BC, h="
                 <> ToString @ NumberForm[N[h], 3] <> "]",
   ColorFunction -> "SunsetColors",
   InterpolationOrder -> 0];

tag = "box" <> ToString[fbox] <> "_N" <> ToString[Nint];
outPath = FileNameJoin[{scriptDir, "sigma_rho_dirichlet_" <> tag <> ".png"}];
Export[outPath, plot, ImageResolution -> 150];
Print["Saved: ", outPath];

(* Save the raw variance field so it can be re-plotted / compared without
   re-running the solve. *)
dumpPath = FileNameJoin[{scriptDir, "sigma_rho_" <> tag <> ".m"}];
Export[dumpPath,
   <|"fbox" -> fbox, "Nint" -> Nint, "h" -> N[h],
     "gridInt" -> gridInt, "sigmaRho" -> sigRhoMat,
     "maxReLambda" -> maxReLam|>];
Print["Saved: ", dumpPath];

(* also dump a small summary *)
Print["Sigma_rho stats: min=", Min@sigRhoVec, "  max=", Max@sigRhoVec,
      "  mean=", Mean@sigRhoVec];
