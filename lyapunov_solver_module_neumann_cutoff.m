(* ::Package:: *)

(* ::Section:: *)
(*Gaussian fluctuations around the \[Rho]\[Dash]Q steady state: discrete Lyapunov solve.*)

(* psi = (d\[Rho], dQ1, dQ2)^T on a uniform cell-centred grid with NO-FLUX walls;
   psi_t = A psi + noise, Sigma solves A Sigma + Sigma A^T + Q = 0.
   Noise uses the EXACT max-entropy moments; the leading-order closure is not
   PSD above S = 0.618.

   ====================== NEUMANN CUTOFF VARIANT ======================
   This is lyapunov_solver_module_neumann.m -- the SEALED, ANCHORING wall -- with
   the grid REPLACED as the UV regulator by a PHYSICAL one: the noise is given a
   finite spatial correlation length lNoise, a new model parameter.

       <F_i(r,t) F_j(r',t')> = M_ij(r) G_l(r-r') delta(t-t'),
       G_l(r) = exp(-r^2/(2 l^2))/(2 Pi l^2),   Ghat_l(k) = exp(-k^2 l^2/2)

   Only the delta in TIME is needed for the Lyapunov equation (it is what makes
   psi Markovian); the spatial kernel is unconstrained, so this costs nothing and
   changes no step of the derivation.  The k^2 of the divergence-form noise still
   cancels against the k^2 of the relaxation rate, so

       S(k) = (zeta/(B rho0)) exp(-k^2 l^2/2)

   which is INTEGRABLE, and the equal-point variance is finite with no cutoff
   imposed by hand:

       Sigma_iso(r) = zeta rho_ss(r)/(2 Pi B rho0 lNoise^2).

   That is dean.tex eq:iso-final term for term: l_UV there is the noise
   CORRELATION LENGTH, not a smoothing width and not the grid spacing.  Setting
   lNoise = 0 recovers the sealed-wall result exactly.  See sec-noise-cutoff.

   IMPLEMENTATION.  The kernel enters as a congruence on the noise matrix that
   eq:Q-assembly already builds:

       Q_l = W . Q . W,      W = MatrixExp[(lNoise^2/4) Lap]

   whose symbol exp(-k^2 l^2/4), applied twice, is Ghat_l.  Reasons for this form
   rather than a dense kernel sandwiched inside the flux products:
     - W = W^T, so Q_l stays symmetric;
     - a congruence preserves PSD, so the realizability structure of
       dean.tex sec-realizability carries over verbatim;
     - Lx.Ly = Ly.Lx = Kron(L1,L1), so the exponential splits EXACTLY as
       W = Kron(E,E) with E = MatrixExp[(lNoise^2/4) L1], one Nint x Nint matrix.

   Lap here is the NEUMANN Laplacian, for ALL sectors including dQ.  The
   correlation length is a property of the rods, and a rod cannot straddle a wall
   it cannot cross, so the kernel REFLECTS.  Two consequences:

     - E is built in closed form from the DCT-II basis that diagonalises the
       cell-centred Neumann L1,

           v_p(i) = c_p Cos[Pi p (i - 1/2)/Nint],  p = 0..Nint-1,
           c_0 = Sqrt[1/Nint],  c_p = Sqrt[2/Nint],
           mu_p = -4 Sin[Pi p/(2 Nint)]^2/h^2,
           Wsym_p = Exp[-(lNoise/h)^2 Sin[Pi p/(2 Nint)]^2],

       so no numerical matrix exponential is ever needed and the mode weight
       underflows gracefully to 0 when lNoise >> h.  Verified against L1 directly
       in claude_experiments/fluctuations_neumann/rho_only_probe.wls.

     - Wsym_0 = 1 EXACTLY, because p = 0 is the constant mode and mu_0 = 0.  Hence
       W.1 = 1, and 1^T (W Q W) = 1^T Q W = 0: the congruence preserves the
       conservation structure of the sealed wall to machine precision.  Asserted
       in [2] as Max|E.1 - 1|.  Contrast the RESERVOIR cutoff variant
       (lyapunov_solver_module_cutoff.m), where W is the DIRICHLET heat kernel and
       its diagonal VANISHES at the wall -- the rods shrink as they approach it.
       Reflecting off a sealed wall the image charge ADDS instead, so the variance
       is ENHANCED there, toward a factor 2 for lNoise >> h.  That factor, and its
       sign, is the sharpest single check that the boundary condition is real.

   ACCURACY.  The lattice weight decays more slowly than the continuum Gaussian
   (sin^2(phi/2) < phi^2/4), so it keeps slightly too much short-wavelength noise.
   The error depends on lNoise/h ALONE and is +(1/4)(h/l)^2 at leading order: 36%
   at l = h, 7.7% at l = 2h, 1.0% at l = 5h, 0.25% at l = 10h.  Require
   lNoise >~ 5 h.  Below lNoise ~ h the Gaussian is invisible to the grid and the
   answer reverts to the grid-limited value -- such a run measures the mesh, not
   the model.

   ---------------- inherited from the sealed-wall variant ----------------
   The boundary condition is changed from the RESERVOIR case to the SEALED case.  dean.tex appendix sec-ou, subsection
   "Boundary conditions" (sec-bc), states that the boundary condition is THREE
   independent choices -- on the steady state, on the fluctuating field, and on
   the noise flux -- and that only two of the eight combinations are consistent:

     reservoir : rho_ss pinned,      psi|_bdry = 0,        flux FREE
                 -> lyapunov_solver_module_R1.m            (set2_1)
     sealed    : n.Grad rho_ss = 0,  n.Grad psi = 0,       flux BLOCKED
                 -> THIS FILE                              (set3_1)

   The pairing is not a convention.  Whether psi is pinned at the wall and
   whether F.n vanishes there are the same question -- can a bacterium cross --
   asked of the mean and of the fluctuation; answering it differently for the
   two puts the wall out of detailed balance by fiat.

   WHAT CHANGED relative to R1, and only this:

   (1) GRID.  Cell-centred, so the wall is a FACE and not a node:

           h = 2 fbox/Nint,   x_i = -fbox + (i - 1/2) h,   i = 1..Nint
           faces at -fbox + m h,  m = 0..Nint;  m = 0 and m = Nint ARE the walls

       Only the Nint-1 INTERIOR faces exist, because the flux is blocked on the
       other two -- so no coefficient is ever sampled on the boundary itself.
       Odd Nint still puts the defect core on a node, at index (Nint+1)/2.

       Useful consequence: Neumann Nint = M+1 has EXACTLY the same h as
       Dirichlet Nint = M (both 2 fbox/(M+1)), so the set2_1 grid list already
       contains an h-matched partner for every odd run here.

   (2) OPERATORS.  G1 is the forward difference from the Nint cell centres to
       the Nint-1 interior faces, and the Laplacian is DEFINED as its Gram
       matrix, so the adjoint-pair identity holds by construction:

           L1 = -G1^T G1 = (1/h^2) tridiag(1, -2, 1) with -1/h^2 in the CORNERS

       which is dean.tex eq:G-neumann.  The collocated difference is the
       node<-face average of the same gradient, D1 = P^T G1, i.e. the reflective
       central difference (row 1 = (-1,1,0,...)/2h, row Nint = (0,...,-1,1)/2h).
       Both satisfy G1.1 = 0 and D1.1 = 0: the gradient of a constant vanishes
       everywhere INCLUDING at the wall.  That is what "reflective" buys, and
       everything below follows from it.

   (3) DIVERGENCES MUST BE TRANSPOSED GRADIENTS.  R1 writes the cross terms as
       Dx.(c * Dy^T), using Dx itself as the divergence.  Under Dirichlet that is
       harmless because D1^T = -D1.  Under Neumann it is fatal: the COLUMN sums
       of D1 are (-1,0,...,0,1)/h, so a flux divergence would create mass in the
       four corner cells at O(h^-2) -- a defect that does not go away with
       refinement.  The same flaw sits in A13 = (2 zeta/xi0) Kron(D1,D1), whose
       column sums are 1/h^2.

       The fix is one rule applied uniformly: every discrete divergence is MINUS
       THE TRANSPOSE OF A GRADIENT, B_i in {-Gi^T, -Di^T}, so that
       1^T B_i = -(Gi.1)^T = 0 identically.  Two expressions change:

           divBlockNeu   cross term   Dx^T.(c * Dy)      was  Dx.(c * Dy^T)
           Dxy           = -(1/2)(Kron[D1,D1^T] + Kron[D1^T,D1])
                                                        was  Kron[D1,D1]

       Both REDUCE TO THE R1 FORMS EXACTLY when D1^T = -D1, which is what the
       Dirichlet central difference is -- so this is a generalisation, not a
       competing scheme (verified in
       claude_experiments/fluctuations_neumann/ops_probe.wls).  Dxy stays a
       symmetric matrix, the two Kroneckers being each other's transpose.  No
       corner interpolation is introduced anywhere, so the phase cos((q-p)/2)
       that would inject a spurious 0.087 Q2 anisotropic pedestal -- the reason
       R1 keeps its cross terms collocated -- never appears here either.

   (4) THE CONSERVED MODE.  With (1)-(3), v = (1,0,0) is an exact RIGHT and LEFT
       null vector of Atilde, and Qtilde v = 0 exactly:

           Atilde . v  = 0    A11.1 = Lap.1 = 0;  A21.1 = A31.1 = 0 via Gx.1 = 0
           v^T . Atilde= 0    A11, A12 symmetric with 1 in the kernel;
                              A13 by the new Dxy
           v^T . Qtilde= 0    every rho-row block is B_i diag(c) B_j^T with
                              1^T B_i = 0

       All four are asserted at run time in [4a] and [5a] and abort on failure;
       they are the sealed-wall analogue of R1's lapIdentity gate.

       So 0 is an exact, semi-simple eigenvalue and the Lyapunov equation is
       SOLVABLE BUT NOT UNIQUE: C is fixed only up to c v v^T.  Physically this
       is dean.tex eq:canonical -- a sealed box conserves the total exactly, so
       the uniform mode has ZERO variance and the right solution is the one with
       the null component projected out.  eigLyapSolve zeroes row and column
       i_null of the transformed right-hand side before dividing, which makes the
       refinement step safe as well.  maxReLambda and min|lam_i+lam_j| are then
       reported EXCLUDING i_null, or the stability warning would fire on the null
       mode every run.

       Closed-form check (claude_experiments/fluctuations_neumann/rho_only_probe.wls):
       for Q_ss = 0 and rho_ss = rho constant the rho sector alone gives

           Sigma = (zeta rho/(B rho0 h^2)) [I - 11^T/Nint^2]

       EXACTLY -- flat pedestal (eq:fdt-split) and the canonical rank-one
       subtraction (eq:canonical) in one object.  Reproduced to 5e-16 relative.
       The Dirichlet counterpart is pedestal * I; the 1/Nint^2 is the whole
       difference between the two consistent boundary conditions.

   ---------------- inherited from the R1 variant ----------------
   Consistent stencil set, so that the discrete fluctuation-dissipation balance
   holds mode by mode and the (grid-regulated) UV pedestal is exactly

       Sigma_iso(r) = zeta rho_ss(r)/(B rho0 h^2),      i.e. R == 1.

   The DIAGONAL noise flux terms use Gx, Gy with their coefficients sampled on
   the corresponding cell faces; the CROSS terms keep the collocated differences
   with node-sampled coefficients.  Setting every coefficient to 1 makes the
   first two terms exactly -Lap, which is what forces R == 1.

   A21, A31 use the divergence form -(B/xi0) sum_i Gi^T diag(Q|face) Gi in place
   of the product-rule expansion; only the divergence form is exactly
   self-adjoint, as Div(c Grad .) is in the continuum.

   The locPf (rotational, Lambda) noise terms are not in divergence form and are
   untouched.  They appear only in the dQ blocks, never in a rho block, so they
   cannot drive the conserved mode.
   ---------------------------------------------------------------

   (5) THE WALL ANCHORS THE DIRECTOR.  Sealing the wall is a statement about
       MASS, and only rho is conserved -- so only the rho sector's boundary
       condition is fixed by the flux argument above.  dQ is not conserved and
       its wall condition is a separate physical question: does the surface exert
       a torque on the director?  It does, and it must.

       Measured, in claude_experiments/fluctuations_neumann/instability_probe.wls:
       giving dQ Neumann as well leaves Atilde with exactly THREE unstable
       eigenvalues at every grid, and they CONVERGE rather than refining away --

           Nint            21        31        41        61
           dQ Neumann   +0.1043   +0.0781   +0.0696   +0.0635   (2 modes)
                        +0.0044   +0.0028   +0.0022   +0.0019   (1 mode)
           dQ Dirichlet   none      none      none      none

       The two degenerate ones are the defect TRANSLATION modes and the third has
       0.994 overlap with the global director ROTATION (dQ1,dQ2) ~ (-2 Q2, 2 Q1).
       That is textbook: a free-anchoring wall images a +1/2 defect
       ATTRACTIVELY, so the centred defect is a saddle and drifts into the wall,
       while free anchoring also un-freezes the rotational Goldstone mode.  A
       wall that blocks mass but exerts no torque simply cannot confine the
       defect, there is no stable steady state, and no equal-time covariance
       exists to compute.

       So the dQ sector gets STRONG ANCHORING, dQ|_bdry = 0, which is exactly
       what the steady state already assumes -- SolveActiveNematicSteady applies
       its "Neumann" only to rho and keeps DirichletCondition[Q == Q_Pade, True]
       on the whole boundary (rho is pinned to 1 at the single boundary node
       x == 0, y = -box, which fixes the additive constant of an otherwise
       singular pure-Neumann problem).  Steady state and fluctuations therefore
       agree, which they did not in the all-Neumann version.

       On a cell-centred grid the wall face sits h/2 from the first centre, so
       its flux is 2 u_1/h and it acts over half a cell.  That gives a clean
       adjoint pair for the anchored sector too,

           L1q = -G1q^T diag(wq) G1q,   wq = (1/2, 1, ..., 1, 1/2),
           G1q : (Nint+1) x Nint, wall rows +-2/h, interior rows (u_m-u_{m-1})/h

       i.e. tridiag(1,-2,1)/h^2 with -3/h^2 in the CORNERS -- the Neumann L1 with
       its corners shifted by -2/h^2, negative definite with no null mode
       (dirichlet_pair_probe.wls).  L1q enters A22 and A33 and NOTHING ELSE.

       Three things deliberately stay Neumann, and none of them is an
       approximation:
         - the rho ROW of the drift (A11, A12, A13).  A12 and A13 are the flux of
           rho driven by grad Q, and THAT flux is blocked at the wall whatever dQ
           does, so the Neumann divergence is the correct statement, not a
           compromise -- and it is what keeps mass conservation exact;
         - every noise block, including the dQ ones.  The flux noise is carried
           by bacteria, which cannot cross the wall, so its boundary condition is
           blocked flux regardless of the anchoring;
         - hence the dQ sector pairs a Dirichlet drift with a blocked-flux noise.
           That is not a stencil error: strong anchoring is the infinite-stiffness
           limit of a surface free energy, so Var(dQ) -> 0 at the wall is the
           correct answer, and dean.tex sec-bc already says the dQ sector is the
           one place a boundary layer survives h -> 0.

   NOTE.  set3_1 therefore realises "sealed for mass, anchored for the director":
   dean.tex sec-bc's sealed row in the rho sector, which is what the flux argument
   constrains, with the director condition inherited from the steady state.

   Entry point is SolveFluctuationsLyapunovNeumannCutoff and runTag carries BOTH
   the grid and a "_Neu_lN<value>" suffix, so no two (Nint, lNoise) runs can
   collide in one directory.  Do NOT load this module together with
   lyapunov_solver_module.m, _R1.m, _cutoff.m or _neumann.m in the same kernel:
   the private helpers are shared and the last one loaded wins.
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

Clear[SolveFluctuationsLyapunovNeumannCutoff, visualizeFluctuations,
      exportData, MergeParams, ScriptParamOverrides];


SolveFluctuationsLyapunovNeumannCutoff::usage =
  "SolveFluctuationsLyapunovNeumannCutoff[lyapunovParams, modelParams, steadyState] solves " <>
  "A \[CapitalSigma] + \[CapitalSigma] A^T + Q = 0 for the equal-time covariance of " <>
  "(\[Delta]\[Rho], \[Delta]Q1, \[Delta]Q2) linearized about steadyState, on a uniform " <>
  "CELL-CENTRED Nint x Nint grid over [-fbox, fbox]^2 with a SEALED, ANCHORING " <>
  "wall: Neumann (zero normal derivative) on \[Delta]\[Rho] and zero noise flux " <>
  "through \[PartialD]\[CapitalOmega], with the director strongly anchored " <>
  "(\[Delta]Q|_bdry = 0) as the steady state already assumes.  Giving \[Delta]Q " <>
  "Neumann too leaves the +1/2 defect unconfined and the linearization unstable; " <>
  "see point (5) of the file header.\n\n" <>
  "The total \[Delta]\[Rho] is then exactly conserved, so \[CapitalSigma] is the " <>
  "CANONICAL covariance: the uniform mode carries zero variance and is projected out " <>
  "(dean.tex eq:canonical).\n\n" <>
  "The noise carries a Gaussian spatial correlation of length lNoise (a " <>
  "modelParams key), so the equal-point variance is UV-finite: " <>
  "\[CapitalSigma]_iso = \[Zeta] \[Rho]_ss/(2 \[Pi] B \[Rho]0 lNoise^2). " <>
  "lNoise = 0 reproduces the grid-regulated sealed-wall result.  The kernel is " <>
  "built on the NEUMANN Laplacian, so it REFLECTS off the wall rather than " <>
  "vanishing there.\n\n" <>
  "  lyapunovParams : rules {Nint->..., fbox->...}; both positive integers, " <>
  "h = 2 fbox/Nint is kept exact so odd Nint puts the defect core on a node\n" <>
  "  modelParams    : same rules as SolveActiveNematicSteady, PLUS lNoise >= 0 " <>
  "(the noise correlation length, in the same length units as fbox; " <>
  "lNoise >= 5 h is required for the continuum value to within 1%)\n" <>
  "  steadyState    : <|\"rho\"->..., \"Q1\"->..., \"Q2\"->..., \"mesh\"->...|>, " <>
  "normally from SolveActiveNematicSteady[..., \"Neumann\"]\n\n" <>
  "Returns <|\"C\", \"sigmaRho\", \"sigmaRhoFn\", \"sigmaQ1\", \"sigmaQ2\", " <>
  "\"gridInt\", \"h\", \"Nint\", \"fbox\", \"bc\", \"lNoise\", \"lNoiseOverH\", " <>
  "\"Wsymbol\", \"eigenvalues\", \"maxReLambda\", " <>
  "\"nullIndex\", \"nullLambdaRel\", \"nullOverlap\", \"consA\", \"consAT\", \"consQ\", " <>
  "\"residual\", \"residualRefined\", \"modelParams\", \"lyapunovParams\", " <>
  "\"timings\"|>.";

SolveFluctuationsLyapunovNeumannCutoff::badgrid =
  "Nint and fbox must be positive integers (h = 2 fbox/Nint is kept exact); got `1`.";
SolveFluctuationsLyapunovNeumannCutoff::badss =
  "steadyState must be an association with keys \"rho\", \"Q1\", \"Q2\"; got `1`.";
SolveFluctuationsLyapunovNeumannCutoff::badlnoise =
  "modelParams must carry a numeric lNoise >= 0 (the noise correlation length); got `1`.";
SolveFluctuationsLyapunovNeumannCutoff::subgrid =
  "lNoise = `1` is only `2` grid spacings: the Gaussian is barely resolved and the \
result is grid-limited, not the continuum value.  Use lNoise >= 5 h.";
SolveFluctuationsLyapunovNeumannCutoff::notcons =
  "The sealed-wall conservation identity `1` failed: `2` against a matrix scale of `3`. \
The discrete divergence is leaking mass, so the conserved mode is driven and the \
Lyapunov equation has no solution.";
SolveFluctuationsLyapunovNeumannCutoff::nullmode =
  "Expected exactly one null eigenvalue (the conserved uniform mode) but found \
|\[Lambda]_null|/max|\[Lambda]| = `1`, second smallest = `2`, overlap with the uniform \
mode = `3`.";
SolveFluctuationsLyapunovNeumannCutoff::unstable =
  "Linearization is not strictly stable away from the conserved mode: \
max Re(\[Lambda]) = `1`; results may be unphysical.";
SolveFluctuationsLyapunovNeumannCutoff::notpsd =
  "Condition (ii) Q >= 0 is violated even with the exact moments at `1` grid points; \
check the closure assembly.";
SolveFluctuationsLyapunovNeumannCutoff::norefine =
  "Refinement did not improve the residual (`1` -> `2`); consider LyapunovSolve.";
SolveFluctuationsLyapunovNeumannCutoff::kappa =
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
   writes many runs into one directory and nothing may be overwritten.  Four
   decimals, never scientific, so 1.0000 / 1.7321 / 5.0000 are all distinct and
   sort sensibly.  The "." is harmless in a filename and to the ".bin" ->
   "_meta.m" replacement readers use. *)
lTag[l_] := "_lN" <> ToString@NumberForm[N[l], {8, 4},
   ExponentFunction -> (Null &), NumberPadding -> {"", "0"}];

(* "_Neu" separates sealed from reservoir, lTag separates the cutoff lengths. *)
runTag[corr_] := "box" <> ToString[corr["fbox"]] <> "_N" <> ToString[corr["Nint"]] <>
                 "_Neu" <> lTag[corr["lNoise"]];

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

(* One conservation gate.  These are the sealed-wall analogue of R1's
   lapIdentity: if any of them fails the scheme is void, so abort rather than
   produce a plausible-looking covariance. *)
consCheck[lbl_, val_, scale_] := (
  Print["    ", StringPadRight[lbl, 40], " = ", val,
        If[Abs[val] <= 1.*^-10 scale, "   (exact to round-off)",
           "   *** CONSERVATION VIOLATED ***"]];
  If[!(Abs[val] <= 1.*^-10 scale),
     Message[SolveFluctuationsLyapunovNeumannCutoff::notcons, lbl, val, scale];
     Print["    ABORT: the sealed-wall scheme is void."];
     Abort[]];
  val);


(* ::Subsection:: *)
(*Main solver*)

SolveFluctuationsLyapunovNeumannCutoff[
    lyapunovParams : {___Rule},
    modelParams : {___Rule},
    steadyState_Association] :=
Module[
  {
    nInt, fboxV, timings = {}, rec,
    Bv, xi0, xir, zet, rho0v, Lam, av, bv, Lv, cL,
    rhoSs, Q1Ss, Q2Ss,
    h, gridInt, iZero, xFace,
    D1, L1, L1q, Idn, Dx, Dy, Lx, Ly, Dxy, Lap, Lapq,
    G1, Pav, Gx, Gy, lapIdentity, anchorGap,
    xflat, yflat, sample, diag,
    xeX, xeY, yeX, yeY, sampleXE, sampleYE,
    rhoV, Q1V, Q2V, rhoLap, Q1x, Q1y, Q2x, Q2y,
    rhoXE, Q1XE, Q2XE, rhoYE, Q1YE, Q2YE,
    SXE, SYE, g2XE, g3XE, g2YE, g3YE, e2XE, e3XE, e2YE, e3YE,
    h2cXE, h2sXE, h3cXE, h3sXE, h2cYE, h2sYE, h3cYE, h3sYE,
    Svals, g2V, g3V, kapResid, tSamp, tClo, tKchk, tCond,
    minEigLoc, scaleLoc, minEigQ, tQpsd,
    A11, A12, A13, A21, A22, A23, A31, A32, A33, Atilde, tAsmA,
    kap, pf, locPf, divBlockNeu, h2c, h2s, h3c, h3s, e2, e3,
    Qrr, QQ1Q1, QQ2Q2, QQ1Q2, QQ2Q1, QQ1r, QQ2r, QrQ1, QrQ2, Qtilde, tAsmQ,
    vCons, scaleA, scaleQ, consA, consAT, consQ,
    Adense, Qdense, tDense, eigvals, eigvecs, tEs, maxReLam, Pev,
    lamMax, iNull, nullRel, lam2Rel, nullOverlap, vUnit,
    eigLyapSolve, lamSum, Cmat, tSol, resid, tRes1, tRef, residRef, tRes3,
    sigRhoVec, sigRhoMat, sigRhoFn, sigQ1, sigQ2, tExtr,
    lNoiseV, Wsym, Emat, Esym, Eone, smoothBlock, rng, tSmooth, Qasym
  },

  rec[lbl_, t_] := (AppendTo[timings, {lbl, t}]; t);

  (* --- unpack + validate ---------------------------------------------
     The keys on the right MUST be Global`: the caller builds the rule lists
     in Global`, and a bare name here would resolve into this private
     context instead and match nothing. *)
  {nInt, fboxV} = {Global`Nint, Global`fbox} /. lyapunovParams;
  If[!IntegerQ[nInt] || !IntegerQ[fboxV] || nInt <= 1 || fboxV <= 0,
    Message[SolveFluctuationsLyapunovNeumannCutoff::badgrid, {fboxV, nInt}]; Return[$Failed]];
  If[!AllTrue[{"rho", "Q1", "Q2"}, KeyExistsQ[steadyState, #] &],
    Message[SolveFluctuationsLyapunovNeumannCutoff::badss, Keys[steadyState]]; Return[$Failed]];

  {Bv, xi0, xir, zet, rho0v, Lam, av, bv, Lv} =
    {Global`B, Global`\[Xi]0, Global`\[Xi]r, Global`\[Zeta], Global`\[Rho]0,
     Global`\[CapitalLambda], Global`a, Global`b, Global`L} /. modelParams;
  cL = Lv + zet xir/(4 xi0);   (* elastic constant in \[Delta]H *)

  (* the cutoff parameter.  A missing key would silently stay the symbol lNoise
     and poison the mode weights, so check it here rather than failing deep in [5b]. *)
  lNoiseV = Global`lNoise /. modelParams;
  If[!NumericQ[lNoiseV] || N[lNoiseV] < 0,
    Message[SolveFluctuationsLyapunovNeumannCutoff::badlnoise, lNoiseV]; Return[$Failed]];
  lNoiseV = N[lNoiseV];

  {rhoSs, Q1Ss, Q2Ss} = steadyState /@ {"rho", "Q1", "Q2"};

  (* --- grid: CELL-CENTRED, so the wall is a face --------------------- *)
  (* x_i = -fbox + (i - 1/2) h needs i = (Nint+1)/2 to hit 0: ODD Nint puts the
     origin (defect core, S = 0) on the grid, EVEN Nint does not.  The Nint-1
     INTERIOR faces sit at -fbox + m h; the two wall faces carry no flux and so
     no coefficient. *)
  h = (2 fboxV)/nInt;
  gridInt = Table[-fboxV + (i - 1/2) h, {i, nInt}];
  xFace   = Table[-fboxV + m h, {m, 1, nInt - 1}];
  Print["[2] Grid: fbox=", fboxV, ", Nint=", nInt,
        ", h=", ToString @ NumberForm[N[h], 5], ", state-dim=", 3 nInt^2,
        "  (cell-centred, sealed walls)"];
  Print["    cells span [", N @ First @ gridInt, ", ", N @ Last @ gridInt,
        "],  ", nInt - 1, " interior faces in [", N @ First @ xFace, ", ",
        N @ Last @ xFace, "],  walls at +-", N @ fboxV];
  Print["    h = 2 fbox/Nint here, so this run matches Dirichlet Nint = ", nInt - 1,
        " (h = ", ToString @ NumberForm[N[2 fboxV/nInt], 5], ") exactly"];
  If[OddQ[nInt],
     iZero = (nInt + 1)/2;
     Print["    Nint odd  -> origin ON the grid at gridInt[[", iZero,
           "]] = ", gridInt[[iZero]], ", flat index k = ", iZero + (iZero - 1) nInt,
           " of ", nInt^2],
     iZero = None;
     Print["    Nint even -> origin NOT on the grid; nearest nodes at x = ",
           ToString @ NumberForm[N @ gridInt[[nInt/2]], 4], " and ",
           ToString @ NumberForm[N @ gridInt[[nInt/2 + 1]], 4]]];

  (* --- 1D FD matrices, built from G1 so the adjoint pair is exact -----
     G1: forward difference, the nInt cell centres -> the nInt-1 INTERIOR faces,
         (G1 u)_m = (u_{m+1} - u_m)/h.  Blocked flux = the two wall faces simply
         do not appear, which is exactly dean.tex eq:G-neumann.
     L1: DEFINED as -G1^T G1, hence tridiag(1,-2,1)/h^2 with -1/h^2 in the
         corners -- the Neumann Laplacian, adjoint pair by construction.
     D1: the node<-face average of the same gradient, P^T G1, i.e. the
         reflective central difference.  G1.1 = 0 and D1.1 = 0. *)
  G1 = SparseArray[Join[Table[{m, m + 1} ->  1/h, {m, nInt - 1}],
                        Table[{m, m}     -> -1/h, {m, nInt - 1}]], {nInt - 1, nInt}];
  L1 = SparseArray[-Transpose[G1] . G1];
  Pav = SparseArray[Join[Table[{m, m}     -> 1/2, {m, nInt - 1}],
                         Table[{m, m + 1} -> 1/2, {m, nInt - 1}]], {nInt - 1, nInt}];
  D1 = SparseArray[Transpose[Pav] . G1];

  (* State ordering k = i + (j-1) N, x fastest.  Hence Dx = I (x) D1, Dy = D1 (x) I. *)
  (* The ANCHORED Laplacian for the dQ sector: the same finite-volume gradient
     but with the two wall faces present, flux 2 u/h over half a cell.  Equal to
     L1 with its corners shifted by -2/h^2, i.e. tridiag(1,-2,1)/h^2 with -3/h^2
     in the corners.  Negative definite -- no null mode, which is the whole point:
     it is what confines the defect.  Used in A22/A33 and nowhere else. *)
  L1q = SparseArray[L1 + SparseArray[{{1, 1} -> -2/h^2, {nInt, nInt} -> -2/h^2},
                                     {nInt, nInt}]];

  Idn = IdentityMatrix[nInt, SparseArray];
  Gx = KroneckerProduct[Idn, G1];
  Gy = KroneckerProduct[G1, Idn];
  Dx = KroneckerProduct[Idn, D1];
  Dy = KroneckerProduct[D1, Idn];
  Lx = KroneckerProduct[Idn, L1];
  Ly = KroneckerProduct[L1, Idn];
  Lap = Lx + Ly;
  Lapq = KroneckerProduct[Idn, L1q] + KroneckerProduct[L1q, Idn];
  (* d_x d_y as a symmetrised DOUBLE DIVERGENCE, so its column sums vanish.
     Equals Kron[D1,D1] identically whenever D1^T = -D1, i.e. in the Dirichlet
     case -- see the header, point (3). *)
  Dxy = -(1/2) (KroneckerProduct[D1, Transpose[D1]] +
                KroneckerProduct[Transpose[D1], D1]);

  (* Trivially true given L1 := -G1^T G1, but kept as a tripwire against edits.
     Exact in rationals, so test == 0. *)
  lapIdentity = Max @ Abs @ Normal[
     Lap + (Transpose[Gx] . Gx + Transpose[Gy] . Gy)];
  Print["    |Lap + (Gx^T Gx + Gy^T Gy)| = ", lapIdentity,
        If[lapIdentity == 0, "   (exact -- adjoint pair)",
           "   *** NOT AN ADJOINT PAIR, R != 1 ***"]];
  If[lapIdentity != 0,
     Print["    ABORT: the staggered/compact identity failed; the scheme is void."];
     Abort[]];
  anchorGap = -Max @ Eigenvalues[N @ Normal @ L1q];
  Print["    dQ anchored: L1q = L1 with corners -3/h^2; its least-damped mode is ",
        -anchorGap, " (strictly negative -- this is what confines the defect)"];
  Print["    Gx.1 = ", Max @ Abs @ Normal[Gx . ConstantArray[1, nInt^2]],
        ",  Dx.1 = ", Max @ Abs @ Normal[Dx . ConstantArray[1, nInt^2]],
        ",  1^T.Dxy = ", Max @ Abs @ Normal[ConstantArray[1, nInt^2] . Dxy],
        "   (all exact 0: no-flux walls)"];

  (* --- the noise-correlation operator W = Kron[Emat, Emat] -------------
     W = MatrixExp[(lNoise^2/4) Lap] and Lap = Lx + Ly with Lx.Ly = Ly.Lx, so the
     exponential splits exactly into one nInt x nInt factor.  The cell-centred
     NEUMANN L1 is diagonalised ANALYTICALLY by the DCT-II basis
         v_p(i) = c_p Cos[Pi p (i - 1/2)/nInt],  p = 0..nInt-1,
         c_0 = Sqrt[1/nInt],  c_p = Sqrt[2/nInt],
         mu_p = -4 Sin[Pi p/(2 nInt)]^2/h^2,
     so Emat = V . diag(Wsym) . V^T with Wsym_p = Exp[-(lNoise/h)^2 Sin[...]^2].
     Closed form beats MatrixExp here: it is one symmetric matmul, it is exact,
     and when lNoise >> h the weights underflow to 0 -- which is the correct
     behaviour (those modes carry no noise) rather than an overflow.

     p = 0 is the CONSTANT mode with mu_0 = 0, so Wsym_0 = 1 exactly and E.1 = 1:
     the congruence Q -> W.Q.W cannot inject noise into the conserved mode.  The
     kernel therefore REFLECTS off the sealed wall instead of vanishing on it. *)
  If[lNoiseV == 0.,
     Wsym = ConstantArray[1., nInt];
     Emat = None;   (* short-circuit: W = I, so this run reproduces _neumann.m *)
     Print["    lNoise = 0 -> W = I; this run is the grid-regulated sealed scheme"],
     (* else *)
     Wsym = Exp[-(lNoiseV/N[h])^2 Sin[N[Pi] Range[0, nInt - 1]/(2 nInt)]^2];
     Emat = With[{V = Table[Sqrt[If[p == 0, 1., 2.]/nInt] *
                            Cos[N[Pi] p (i - 1/2)/nInt],
                            {i, nInt}, {p, 0, nInt - 1}]},
                 V . (Wsym * Transpose[V])];
     Esym = Max @ Abs[Emat - Transpose[Emat]];
     Eone = Max @ Abs[Emat . ConstantArray[1., nInt] - 1.];
     Print["    lNoise = ", lNoiseV, " = ", ToString@NumberForm[lNoiseV/N[h], 4],
           " h;  W = Kron[E,E] on the NEUMANN (reflecting) Laplacian"];
     Print["    E symmetry residual = ", Esym, ",  |E.1 - 1| = ", Eone,
           ",  Wsym[[1]] - 1 = ", Wsym[[1]] - 1., ",  max |E| = ", Max@Abs@Emat];
     Print["    W mode weights \[Element] ", MinMax@Wsym,
           "  (Exp[-(lNoise/h)^2 Sin[Pi p/(2 Nint)]^2], p = 0..Nint-1)"];
     If[Esym > 1.*^-12, Print["    WARNING: E is not symmetric to 1e-12."]];
     If[Eone > 1.*^-12,
        Print["    ABORT: E.1 != 1, so W would drive the conserved mode."]; Abort[]];
     If[lNoiseV < 5 N[h],
        Message[SolveFluctuationsLyapunovNeumannCutoff::subgrid, lNoiseV,
                ToString@NumberForm[lNoiseV/N[h], 3]]]];

  (* --- sample the steady state -------------------------------------- *)
  (* Batched sampling, ~200x faster than scalar calls.  N@ is load-bearing:
     exact rationals bypass the fast path and are 19x slower. *)
  xflat = Flatten @ ConstantArray[N @ gridInt, nInt];
  yflat = Flatten @ Transpose @ ConstantArray[N @ gridInt, nInt];
  sample[f_] := f[xflat, yflat];

  (* The diagonal flux coefficients live on the interior cell faces, so sample the
     FEM steady state THERE rather than averaging node values -- the interpolant is
     defined everywhere, so this costs nothing in accuracy.
       x-faces: (face_x, node_y), (nInt-1) x nInt, face index fastest
       y-faces: (node_x, face_y), nInt x (nInt-1), node index fastest
     The layout matches Gx = Idn (x) G1 and Gy = G1 (x) Idn respectively. *)
  xeX = Flatten @ ConstantArray[N @ xFace, nInt];
  xeY = Flatten @ Transpose @ ConstantArray[N @ gridInt, nInt - 1];
  yeX = Flatten @ ConstantArray[N @ gridInt, nInt - 1];
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
     (* the same three fields on the two face sets *)
     rhoXE = sampleXE[rhoSs];  Q1XE = sampleXE[Q1Ss];  Q2XE = sampleXE[Q2Ss];
     rhoYE = sampleYE[rhoSs];  Q1YE = sampleYE[Q1Ss];  Q2YE = sampleYE[Q2Ss];]];
  Print["    sampling (14 batched FEM calls: 8 nodal + 6 on the faces): ", tSamp, " s"];
  Print["    \[Rho]_ss \[Element] ", MinMax@rhoV,
        "   Q1_ss \[Element] ", MinMax@Q1V,
        "   Q2_ss \[Element] ", MinMax@Q2V];

  Svals = Sqrt[Q1V^2 + Q2V^2]/rhoV;
  SXE = Sqrt[Q1XE^2 + Q2XE^2]/rhoXE;
  SYE = Sqrt[Q1YE^2 + Q2YE^2]/rhoYE;
  (* The closure is needed wherever a coefficient is needed, so three times over.
     ~3 Nint^2 extra root-finds, against a multi-hour Eigensystem. *)
  {tClo, {g2V, g3V, g2XE, g3XE, g2YE, g3YE}} = AbsoluteTiming[
     {g2Exact /@ Svals, g3Exact /@ Svals,
      g2Exact /@ SXE,   g3Exact /@ SXE,
      g2Exact /@ SYE,   g3Exact /@ SYE}];
  rec["[3] kappa(S) inversion + g2,g3 (nodes + both face sets)", tClo];
  (* kappaOfS is not memoized, so this repeats all Nint^2 root-finds. *)
  {tKchk, kapResid} = AbsoluteTiming[
     Max @ Abs[(BesselI[1, #]/BesselI[0, #] & /@ (kappaOfS /@ Svals)) - Svals]];
  rec["[3] kappa(S) verification (diagnostic)", tKchk];
  Print["    S = |Q_ss|/\[Rho]_ss \[Element] ", MinMax@Svals, ",  mean ", Mean@Svals];
  Print["    S on x-faces \[Element] ", MinMax@SXE,
        ",  on y-faces \[Element] ", MinMax@SYE];
  Print["    g2 \[Element] ", MinMax@g2V, "   g3 \[Element] ", MinMax@g3V,
        "   (small-S limits 1/2, 1/6)"];
  Print["    \[Kappa](S) inversion: ", tClo, " s, max |I1/I0 - S| = ", kapResid];
  If[kapResid > 1.*^-8, Message[SolveFluctuationsLyapunovNeumannCutoff::kappa, kapResid]];

  (* --- drift matrix A ------------------------------------------------ *)
  (* Q = {{Q1,Q2},{Q2,-Q1}}.  The comoving term cancels the one from
     Div(dQ Grad rho_ss).  Sign a -> a + Lambda xir: rotational diffusion
     DESTROYS order.

     A11, A12 are exactly the staggered operators:  Lx = -Gx^T Gx and
     Ly = -Gy^T Gy separately, so
        A11 = (B/xi0) Lap        = -(B/xi0)(Gx^T Gx + Gy^T Gy)
        A12 = (zeta/xi0)(Lx-Ly)  = (zeta/xi0)(-Gx^T Gx + Gy^T Gy)
     both symmetric with the constant in their kernel.  A13 uses the
     double-divergence Dxy of the header, which is what makes the rho row of the
     drift conserve mass exactly.  All three stay NEUMANN even though dQ is
     anchored: A12 and A13 carry the flux of rho driven by grad Q, and that flux
     is blocked at the wall whatever dQ does.

     A22 and A33 are the only blocks that see the anchoring, through Lapq.

     A21, A31 keep the divergence form -(B/xi0) sum_i Gi^T diag(Q|face) Gi, which
     is exactly self-adjoint as Div(c Grad .) is in the continuum, and reuses the
     face samples the noise already needs. *)
  Print["[4] Assembling A ..."];
  diag[c_] := DiagonalMatrix[SparseArray @ c];
  tAsmA = rec["[4] assemble A", First @ AbsoluteTiming[
     A11 = (Bv/xi0) Lap;
     A12 = (zet/xi0) (Lx - Ly);
     A13 = (2 zet/xi0) Dxy;

     A21 = -(Bv/xi0) (Transpose[Gx] . (Q1XE * Gx) + Transpose[Gy] . (Q1YE * Gy));
     (* Lapq, not Lap: the wall anchors the director.  Without it the defect is
        not confined and Atilde has three unstable modes -- see header point (5). *)
     A22 = (Bv/xi0) diag[rhoLap] +
           (4/xir) (diag[-(av + Lam xir) - bv (6 Q1V^2 + 2 Q2V^2)] + cL Lapq);
     A23 = -(16 bv/xir) diag[Q1V Q2V];

     A31 = -(Bv/xi0) (Transpose[Gx] . (Q2XE * Gx) + Transpose[Gy] . (Q2YE * Gy));
     A32 = A23;
     A33 = (Bv/xi0) diag[rhoLap] +
           (4/xir) (diag[-(av + Lam xir) - bv (2 Q1V^2 + 6 Q2V^2)] + cL Lapq);

     Atilde = ArrayFlatten[{{A11, A12, A13},
                            {A21, A22, A23},
                            {A31, A32, A33}}];]];
  Print["    A dims = ", Dimensions[Atilde], ",  assembled in ", tAsmA, " s"];

  (* --- [4a] the sealed-wall conservation gate on A -------------------- *)
  vCons = Join[ConstantArray[1., nInt^2], ConstantArray[0., 2 nInt^2]];
  scaleA = Max @ Abs @ Atilde["NonzeroValues"];
  Print["[4a] Conservation of the uniform \[Delta]\[Rho] mode (scale ", scaleA, ") ..."];
  consA  = consCheck["A . v   (right null vector)",
                     Max @ Abs @ Normal[Atilde . vCons], scaleA];
  consAT = consCheck["v^T . A (left null vector)",
                     Max @ Abs @ Normal[vCons . Atilde], scaleA];

  (* --- noise matrix Q ------------------------------------------------ *)
  (* f = d_i F_i, <F_i F_j> = C_ij delta(r-r') -> sum_ij B_i diag(C_ij) B_j^T with
     B_i = -Gi^T (staggered, face-sampled) or -Di^T (collocated, node-sampled).
     Writing EVERY divergence as minus a transposed gradient is what makes
     1^T B_i = -(Gi.1)^T = 0, i.e. no noise reaches the conserved mode.
     EXACT max-entropy moments:
       <QQ> = rho J + (g2/2rho) U,  <Q'Q'> = 4 rho J - (2 g2/rho) U,
       <QQQ> = (1/2) T + (g3/rho^2) W;   g2 = g3 = 0 gives the old closure.

     Setting every coefficient to 1 makes the first two terms exactly -Lap, which
     is what forces R == 1.  The cross terms stay collocated (no corner
     interpolation), so no phase cos((q-p)/2) and no spurious anisotropic
     pedestal. *)
  Print["[5] Assembling Q ..."];
  kap = 2 zet/(xi0 rho0v);
  pf = kap/h^2;
  locPf = 2 Lam/(rho0v h^2);
  divBlockNeu[cxxF_, cyyF_, cxy_, cyx_] :=
    Transpose[Gx] . (cxxF * Gx) + Transpose[Gy] . (cyyF * Gy) +
    Transpose[Dx] . (cxy * Dy) + Transpose[Dy] . (cyx * Dx);

  tAsmQ = rec["[5] assemble Q (24 sparse products)", First @ AbsoluteTiming[
     h2c = Q1V^2 - Q2V^2;        h2s = 2 Q1V Q2V;
     h3c = Q1V^3 - 3 Q1V Q2V^2;  h3s = 3 Q1V^2 Q2V - Q2V^3;
     e2 = g2V/(2 rhoV);          (* U terms *)
     e3 = g3V/(4 rhoV^2);        (* W terms *)
     (* the same combinations on the two face sets, for the ii flux terms *)
     h2cXE = Q1XE^2 - Q2XE^2;        h2sXE = 2 Q1XE Q2XE;
     h3cXE = Q1XE^3 - 3 Q1XE Q2XE^2; h3sXE = 3 Q1XE^2 Q2XE - Q2XE^3;
     e2XE = g2XE/(2 rhoXE);          e3XE = g3XE/(4 rhoXE^2);
     h2cYE = Q1YE^2 - Q2YE^2;        h2sYE = 2 Q1YE Q2YE;
     h3cYE = Q1YE^3 - 3 Q1YE Q2YE^2; h3sYE = 3 Q1YE^2 Q2YE - Q2YE^3;
     e2YE = g2YE/(2 rhoYE);          e3YE = g3YE/(4 rhoYE^2);

     (* rho-rho: no closure enters *)
     Qrr = pf * divBlockNeu[Q1XE + rhoXE, -Q1YE + rhoYE, Q2V, Q2V];

     QQ1Q1 = (pf * divBlockNeu[3 Q1XE/4 + rhoXE/2 + e2XE h2cXE + e3XE h3cXE,
                               -3 Q1YE/4 + rhoYE/2 + e2YE h2cYE - e3YE h3cYE,
                               Q2V/4 + e3 h3s,
                               Q2V/4 + e3 h3s] +
              locPf * diag[2 rhoV - 4 e2 h2c]);

     QQ2Q2 = (pf * divBlockNeu[Q1XE/4 + rhoXE/2 - e2XE h2cXE - e3XE h3cXE,
                               -Q1YE/4 + rhoYE/2 - e2YE h2cYE + e3YE h3cYE,
                               3 Q2V/4 - e3 h3s,
                               3 Q2V/4 - e3 h3s] +
              locPf * diag[2 rhoV + 4 e2 h2c]);

     (* NB Q1-Q2 carries a LOCAL part: <Q'Q'> is no longer proportional to J. *)
     QQ1Q2 = (pf * divBlockNeu[Q2XE/4 + e2XE h2sXE + e3XE h3sXE,
                               -Q2YE/4 + e2YE h2sYE - e3YE h3sYE,
                               Q1V/4 - e3 h3c,
                               Q1V/4 - e3 h3c] +
              locPf * diag[-4 e2 h2s]);
     QQ2Q1 = Transpose[QQ1Q2];

     QQ1r = pf * divBlockNeu[Q1XE + rhoXE/2 + e2XE h2cXE,
                             Q1YE - rhoYE/2 - e2YE h2cYE,
                             e2 h2s,
                             e2 h2s];
     QQ2r = pf * divBlockNeu[Q2XE + e2XE h2sXE,
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
     all -- a fraction of a percent of the 27 nInt^6 Eigensystem below.

     IN PLACE, ONE BLOCK PAIR AT A TIME.  Qtilde is densified here (stage [6]
     needs it dense anyway, so this moves that cost rather than adding it) and each
     block is then overwritten via Part assignment.  Building all nine blocks and
     ArrayFlatten-ing them would instead hold 9 + 9 blocks of temporaries, and a
     whole-matrix (Q + Q^T)/2 would add two more full copies -- hence the pairwise
     symmetrisation below, which is block-sized and enforces exact global symmetry.
     Folding W into Gx/Dx instead would make them dense and cost O(nInt^6).

     Note this runs BEFORE the [5a] conservation gate, so the gate tests the matrix
     the Lyapunov solve actually uses.  It must still pass: E.1 = 1 exactly. *)
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
        Qtilde = Developer`ToPackedArray @ N @ Normal @ Qtilde;
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

  (* --- [5a] the sealed-wall conservation gate on Q -------------------- *)
  scaleQ = If[Head[Qtilde] === SparseArray,
               Max @ Abs @ Qtilde["NonzeroValues"], Max @ Abs @ Qtilde];
  Print["[5a] Noise into the conserved mode (scale ", scaleQ, ") ..."];
  consQ = consCheck["v^T . Q (no noise on the total)",
                    Max @ Abs @ Normal[vCons . Qtilde], scaleQ];

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
     Message[SolveFluctuationsLyapunovNeumannCutoff::notpsd,
             Count[minEigLoc, x_ /; x < -1.*^-10 scaleLoc]]];

  (* The node-wise 6x6 above is still the RIGHT physical condition -- at symbol
     level the ii terms give 4(cxx sp^2 + cyy sq^2) and the cross term
     8 cxy sp cp sq cq, and 4(cxx sp^2 + cyy sq^2) >= 8 Sqrt[cxx cyy] sp sq by
     AM-GM while cp cq <= 1, so PSD follows from Sqrt[cxx cyy] >= |cxy|, which is
     exactly pointwise PSD of the local flux covariance.  But the mixed scheme is
     no longer a single Gram matrix D M D^T, so check the assembled Qtilde too.
     Q is only PSD, never PD: the conserved mode is an exact null direction, so
     ask for FIVE eigenvalues and expect one at round-off.
     Arnoldi on the sparse matrix; skipped above n = 20000 where it gets slow. *)
  If[3 nInt^2 <= 20000,
     {tQpsd, minEigQ} = AbsoluteTiming[
        Quiet @ Check[
           Min @ Re @ Eigenvalues[N @ Qtilde, -5, Method -> "Arnoldi"],
           Indeterminate]];
     rec["[5] condition (ii): min eig of assembled Q (Arnoldi)", tQpsd];
     Print["    assembled Q min eig = ", minEigQ, "   [", tQpsd, " s]",
           If[NumericQ[minEigQ] && minEigQ < -1.*^-8 Max@Abs@Diagonal@Qtilde,
              "   *** Q IS NOT PSD ***", "   (PSD; one exact zero is the sealed mode)"]],
     minEigQ = Missing["skipped, n > 20000"];
     Print["    assembled Q min eig: skipped (n = ", 3 nInt^2, " > 20000)"]];

  (* --- Lyapunov solve ------------------------------------------------- *)
  (* LyapunovSolve is single-threaded; diagonalize instead: A = P D P^-1 =>
     C = P Y P^T, Y_ij = -(P^-1 Q P^-T)_ij/(lam_i + lam_j).  A is NOT self-adjoint
     (dropped Onsager partners), but the error is set by kappa(P)^2 eps ~ 1e-10;
     one refinement step then reaches ~1e-14.  The printed residual verifies it.

     SEALED WALLS: lam = 0 is an exact eigenvalue (the conserved uniform mode), so
     Y_{i0 i0} would be 0/0.  Zero the whole i0 row and column of the transformed
     right-hand side instead -- that is the canonical projection of eq:canonical,
     and doing it INSIDE eigLyapSolve makes the refinement step safe too. *)
  Print["[6] Lyapunov solve via eigendecomposition + refinement ..."];
  (* ToPackedArray is load-bearing.  G1, L1 and D1 are built from EXACT
     RATIONALS here so that the adjoint-pair identity can be tested with "== 0";
     if any of that exactness survives into the dense arrays they come back
     unpacked, every BLAS path is lost, and the back-substitute below runs
     interpreted -- measured at 600x slower than R1's. *)
  tDense = rec["[6] densify A, Q (Normal[])", First @ AbsoluteTiming[
     Adense = Developer`ToPackedArray @ N @ Normal @ Atilde;
     Qdense = Developer`ToPackedArray @ N @ Normal @ Qtilde;]];
  Print["    [0] densify A,Q        ", tDense, " s   (",
        ToString @ NumberForm[N[2 * 8 * Length[Adense]^2/2^30], 3], " GiB)",
        "   packed: A ", Developer`PackedArrayQ[Adense],
        ", Q ", Developer`PackedArrayQ[Qdense]];

  {tEs, {eigvals, eigvecs}} = AbsoluteTiming @ Eigensystem[Adense];
  rec["[6] Eigensystem (dense, nonsymmetric)", tEs];
  Print["    [a] Eigensystem        ", tEs, " s"];

  (* identify the conserved mode and check it is the only one *)
  lamMax = Max @ Abs @ eigvals;
  iNull = First @ Ordering[Abs[eigvals], 1];
  nullRel = Abs[eigvals[[iNull]]]/lamMax;
  lam2Rel = Sort[Abs @ eigvals][[2]]/lamMax;
  vUnit = Normalize @ vCons;
  nullOverlap = Abs[Normalize[eigvecs[[iNull]]] . vUnit];
  Print["    conserved mode: index ", iNull, ",  |\[Lambda]|/max|\[Lambda]| = ", nullRel,
        ",  2nd smallest = ", lam2Rel];
  Print["    overlap of its eigenvector with uniform \[Delta]\[Rho] = ", nullOverlap,
        If[nullRel < 1.*^-8 && lam2Rel > 1.*^-6 && nullOverlap > 0.99,
           "   (a single, clean null mode)", "   *** UNEXPECTED NULL STRUCTURE ***"]];
  If[!(nullRel < 1.*^-8 && lam2Rel > 1.*^-6 && nullOverlap > 0.99),
     Message[SolveFluctuationsLyapunovNeumannCutoff::nullmode, nullRel, lam2Rel, nullOverlap]];

  (* stability, EXCLUDING the conserved mode -- it sits exactly at 0 by design *)
  maxReLam = Max @ Re @ Delete[eigvals, iNull];
  Print["    max Re(\[Lambda]) off the null mode = ", maxReLam,
        If[maxReLam < 0, "  (stable)", "  (UNSTABLE!)"]];
  If[maxReLam >= 0, Message[SolveFluctuationsLyapunovNeumannCutoff::unstable, maxReLam]];
  Print["    min |\[Lambda]_i+\[Lambda]_j| off the null mode = ",
        2 Min @ Abs @ Re @ Delete[eigvals, iNull],
        ",  complex fraction = ",
        N[100 Count[eigvals, z_ /; Abs[Im[z]] > 1.*^-8 Abs[z]]/Length[eigvals]], "%"];

  Pev = Developer`ToPackedArray @ Transpose @ eigvecs;
  lamSum = Developer`ToPackedArray @ Outer[Plus, eigvals, eigvals];
  (* lam_null is 0 BY DESIGN, so lamSum[[iNull,iNull]] can come back EXACTLY 0., and
     then the 0/0 the comment above warns about is only half-fixed by zeroing the
     numerator: 0./0. is Indeterminate, not 0, and a single such entry is smeared
     over all of Cmat by the two matmuls in eigLyapSolve.  (Measured: 2 of 38
     set3_2 runs; the rest escaped only because Eigensystem returned |lam| ~ 1e-16
     instead of 0.)  Shift that one denominator off zero.  The numerator there is
     exactly 0., so the quotient is 0 either way and eq:canonical is untouched.
     ADD rather than assign, and add lamMax (a machine Real): the sum keeps the
     element type of the array, so a packed Complex stays packed Complex and a
     packed Real stays packed Real.  Assigning a bare 1. into a packed COMPLEX
     array unpacks it and loses every BLAS path below. *)
  lamSum[[iNull, iNull]] = lamSum[[iNull, iNull]] + lamMax;
  Print["    packed: eigvals ", Developer`PackedArrayQ[eigvals],
        ", Pev ", Developer`PackedArrayQ[Pev],
        ", lam_i+lam_j ", Developer`PackedArrayQ[lamSum]];
  (* Solve A X + X A^T = -Rhs in the eigenbasis; reused by the refinement.
     The conserved mode is projected out by zeroing row and column iNull of the
     transformed right-hand side -- eq:canonical.  Multiply by 0. rather than
     ASSIGNING 0.: assigning a real into a packed COMPLEX array unpacks it, and
     everything after that runs interpreted. *)
  eigLyapSolve[Rhs_] := Module[{num},
     num = -LinearSolve[Pev, Transpose @ LinearSolve[Pev, Transpose[Rhs]]];
     num[[iNull]] = num[[iNull]] 0.;
     num[[All, iNull]] = num[[All, iNull]] 0.;
     Re[Pev . (num / lamSum) . Transpose[Pev]]];

  {tSol, Cmat} = AbsoluteTiming @ eigLyapSolve[Qdense];
  rec["[6] back-substitute", tSol];
  Print["    [b] back-substitute    ", tSol, " s"];

  (* Each residual is two dense n x n matmuls, same O(n^3) as the Eigensystem. *)
  {tRes1, resid} = AbsoluteTiming[
     Max @ Abs[Adense . Cmat + Cmat . Transpose[Adense] + Qdense]];
  rec["[6] residual eval #1 (unrefined)", tRes1];
  Print["    residual (unrefined)   = ", resid, "   [", tRes1, " s]"];
  (* Non-numeric here means the solve produced Indeterminate/NaN entries.  Stop:
     the refinement step and exportData would otherwise write the whole
     unevaluated expression to disk (measured: 3.8 GB of .m and 12.3 GB of log).
     Every other gate in this module tests the INPUTS; this is the only one that
     tests that a number came out.  Max@Abs[] propagates Indeterminate from a
     single entry and NumberQ[Indeterminate] is False, so one test suffices. *)
  If[!NumberQ[resid],
     Print["    ABORT: residual is ", resid,
           " -- the Lyapunov solve returned non-numeric entries."];
     Abort[]];

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
     Message[SolveFluctuationsLyapunovNeumannCutoff::norefine, resid, residRef]];
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
  (* the canonical rank-one subtraction, eq:canonical: the total is fixed, so the
     sum of the covariance over BOTH indices is exactly zero *)
  Print["    canonical check: sum of the rho-rho block = ",
        Total @ Flatten @ Cmat[[1 ;; nInt^2, 1 ;; nInt^2]],
        "   (exactly 0 for a sealed box)"];

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
     "QasymBeforeSym" -> Qasym,
     "bc" -> "sealed + anchoring wall: \[Delta]\[Rho] Neumann, noise flux blocked \
(all sectors), \[Delta]Q strongly anchored (Dirichlet), canonical (uniform \
\[Delta]\[Rho] mode projected out)",
     "eigenvalues" -> eigvals,
     "maxReLambda" -> maxReLam,
     "nullIndex" -> iNull,
     "nullLambdaRel" -> nullRel,
     "nullLambda2Rel" -> lam2Rel,
     "nullOverlap" -> nullOverlap,
     "consA" -> consA,
     "consAT" -> consAT,
     "consQ" -> consQ,
     "residual" -> resid,
     "residualRefined" -> residRef,
     "minEigLocal" -> MinMax@minEigLoc,
     "minEigQ" -> minEigQ,
     "lapIdentity" -> lapIdentity,
     "scheme" -> "sealed/anchored: cell-centred grid; L1 = -G1^T G1 (Neumann \
Laplacian) for the rho sector and all noise; L1q = -G1q^T diag(1/2,1..1,1/2) G1q \
(anchored Dirichlet Laplacian, corners -3/h^2) in A22/A33 only; every divergence a \
transposed gradient; symmetrised double-divergence d_x d_y; canonical projection of \
the conserved uniform mode; PLUS a Gaussian noise correlation length lNoise \
imposed as Q -> W.Q.W with W = MatrixExp[(lNoise^2/4) Lap_Neumann] = Kron[E,E] \
(reflecting kernel)",
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
  (* the interpolant only spans the cell centres, +-(fbox - h/2), so clamp to
     that; with sealed walls the boundary is the point of the run, so plotBox is
     normally left at fbox rather than cropping. *)
  {xLo, xHi} = {Max[-plotBox, First@gridInt], Min[plotBox, Last@gridInt]};
  plotPts = Min[correlations["Nint"], 200];
  lbl = ("Var(\[Delta]\[Rho])(r)  [sealed wall: \[Rho] Neumann, Q anchored, h=" <>
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

  savePlot[dir, "sigma_rho_neumann_" <> tag <> ".png", plot2D];
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
       "bc" -> correlations["bc"],
       (* lNoise is in modelParams too, but every reader of a sweep wants it
          without unpacking a rule list, and lNoiseOverH is the number that says
          whether the run is grid-limited or converged. *)
       "lNoise" -> correlations["lNoise"],
       "lNoiseOverH" -> correlations["lNoiseOverH"],
       "sigmaRho" -> correlations["sigmaRho"],
       "sigmaQ1" -> correlations["sigmaQ1"],
       "sigmaQ2" -> correlations["sigmaQ2"],
       "maxReLambda" -> correlations["maxReLambda"],
       "nullLambdaRel" -> correlations["nullLambdaRel"],
       "nullOverlap" -> correlations["nullOverlap"],
       "consA" -> correlations["consA"],
       "consAT" -> correlations["consAT"],
       "consQ" -> correlations["consQ"],
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
          "bc" -> correlations["bc"],
          "lNoise" -> correlations["lNoise"],
          "lNoiseOverH" -> correlations["lNoiseOverH"],
          "blocks" -> {"rho", "Q1", "Q2"},
          "ordering" -> "row-major Real64; grid index k = i + (j-1) Nint, x fastest",
          "gridNote" -> "CELL-CENTRED: x_i = -fbox + (i - 1/2) h, h = 2 fbox/Nint; \
walls are faces, not nodes",
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
