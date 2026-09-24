(* ::Package:: *)

(* ::Section:: *)
(*Gaussian fluctuations around the comoving \[Rho]\[Dash]Q steady state: discrete Lyapunov solve.*)

(* psi = (d\[Rho], dQ1, dQ2)^T;  psi_t = A psi + noise;  Sigma solves
   A Sigma + Sigma A^T + Q = 0.

   ====================================================================
   WHAT THIS MODULE IS
   ====================================================================
   One solver, two boundary conditions, replacing the five
   ../lyapunov_solver_module*.m files.  Three things changed relative to them:

   (1) THE DRIFT IS THE COMOVING ONE.  See "The advection term" below.  This is
       a physics correction, not a refactor: the old modules are wrong.
   (2) lNoise IS CALIBRATED against the continuum theory by default, so a run at
       a given lNoise is directly comparable with dean.tex eq:iso-gauss without
       needing a fine grid.  See "Cutoff calibration".
   (3) BOTH consistent boundary conditions of dean.tex eq:bc-triple live here,
       selected by the bcType argument, instead of in separate files.

   The steady state MUST come from comoving_defect/comoving_steady_solver.m --
   the module aborts on one that carries no frame velocity "u".

   ====================================================================
   THE ADVECTION TERM  (dean.tex eq. linearized-Q-equation)
   ====================================================================
   dean.tex eq. Q-comoving is

       d_t Q = Div( Q (u + (B/xi0) Grad rho) ) + (4/xir) H_active

   with u a CONSTANT vector, obtained as the Lagrange multiplier conjugate to
   the pinning constraints Q1(0,0) = Q2(0,0) = 0 (dean.tex sec-comoving-solver).
   Linearising it gives

       d_t dQ = (B/xi0) Div(Q_ss Grad drho)
                + Div(dQ w_ss)                        <---- THIS TERM
                + (4/xir) dH_active + noise,
       w_ss   = u + (B/xi0) Grad rho_ss.

   THE OLD MODULES DROP Div(dQ w_ss) ENTIRELY.  Their A22/A33 carry no
   first-derivative operator at all, only (B/xi0) diag[Lap rho_ss].  The reason
   is archaeological: earlier versions of dean.tex wrote the dQ equation with
   +(B/xi0)(Grad rho_ss . Grad) dQ on the LEFT-hand side, which is -(u.Grad) dQ
   evaluated with the SUPERSEDED pointwise estimate u = -(B/xi0) Grad rho
   (FRIDGE/wrong_comoving.tex:45).  With that substitution the advection cancels
   identically.  With u a constant it does not: the u of the frame cancels the u
   inside w_ss, and the (B/xi0) Grad rho_ss advection SURVIVES.

   It is not small.  At the parameters of dean.tex sec-parameter-choice, at the
   defect core,

       (B/xi0) |Grad rho_ss| = 1.07 um/min     |u| = 0.28 um/min
       ->  |w_ss| = 0.78 um/min,  i.e. |w_ss| xir elld/(4K') = 0.15
           relative to the elastic relaxation of dQ.

   IMPLEMENTATION, and the trap.  Written in conservative (finite-volume) form,

       Adv = -( Gx^T diag(w_x|x-face) Px + Gy^T diag(w_y|y-face) Py )

   with Px, Py the node->face AVERAGING operators (the 1/2-valued twins of Gx,
   Gy).  Because u is constant,

       Div(dQ w_ss) = (w_ss.Grad) dQ + (Div w_ss) dQ,   Div w_ss = (B/xi0) Lap rho_ss

   so Adv ALREADY CONTAINS the (B/xi0) diag[rhoLap] term the old modules have.
   That term is therefore DELETED from A22/A33 here.  Keeping both double-counts
   it.  The identity is asserted at run time in [4b] as Adv.1 vs (B/xi0) rhoLap,
   on interior nodes only (the constant field is not representable under either
   wall condition, so the outermost rows differ legitimately).

   Two things this does NOT disturb:
     - MASS CONSERVATION.  Adv enters A22/A33 only.  The rho row block and the
       drho column block are untouched, so the sealed-wall gates [4a]/[5a] pass
       exactly as before.  dQ is not a conserved quantity, so nothing is owed.
     - THE NEED FOR UPWINDING.  The cell Peclet number is |w| h/(2 kappa) with
       kappa = 4K'/xir = 4.15; at |w| ~ 0.8 and h ~ 0.3 that is ~0.03.  Centred
       differencing is correct here by three orders of magnitude.  It is printed
       every run and warned on above 1.

   At a sealed wall dQ is anchored (dQ = 0), so the advective flux through the
   wall vanishes and the interior-face Gx/Px pair is already the right operator:
   no G1q analogue is needed for Adv.

   ====================================================================
   CUTOFF CALIBRATION  (dean.tex eq:lattice-bz, eq:ln-eff)
   ====================================================================
   The noise carries a Gaussian spatial correlation length, imposed as a
   congruence Q -> W.Q.W with W = MatrixExp[(l^2/4) Lap], so that the continuum
   equal-point variance is UV-finite,

       Sigma_iso(r) = zeta rho_ss(r)/(2 Pi B rho0 lNoise^2).         (eq:iso-gauss)

   But the LATTICE heat kernel decays more slowly than the continuum Gaussian
   (sin^2(phi/2) < phi^2/4), so it keeps too much short-wavelength noise and the
   measured pedestal OVERSHOOTS.  dean.tex eq:lattice-bz gives the factor in
   closed form: with s = l/h and F(s) = Exp[-s^2] BesselI[0, s^2],

       Sigma_lattice(l) = (zeta rho_ss/(B rho0 h^2)) F(s)^2.

   So a run at face value measures a cutoff that is not the one asked for.  By
   default this module fixes that: given the PHYSICAL lNoise it solves

       F(l/h) = 1/(Sqrt[2 Pi] t),        t = lNoise/h                 <--- [1b]

   for the SOLVER length l (reported as lNoiseEff) and builds W from l, so that
   Sigma_lattice(l) equals eq:iso-gauss at the physical lNoise.  Three checks
   that this is the right rule rather than a fitted one:

     - F decreases monotonically from F(0) = 1 to F(s) ~ 1/(Sqrt[2Pi] s), so
       l -> lNoise as h -> 0.  The correction vanishes on fine grids.
     - It is solvable exactly when t >= 1/Sqrt[2 Pi] = 0.3989, i.e.
       lNoise >= h/Sqrt[2 Pi].  That is eq:ln-eff recovered as the solvability
       edge: at equality l = 0 and the BARE GRID already is the correct cutoff.
       Below it no smoothing can help -- smoothing only removes variance -- and
       the module aborts with ::uncalibratable rather than return a number.
     - Expanding, l ~ lNoise + h^2/(8 lNoise), which inverts the +(1/4)(h/l)^2
       bias the old modules documented.  At lNoise = 1.21 h that is +8.5% in l,
       undoing the 17% variance error tabulated in
       ../jobscripts/set2_2_comove/fluctuations.wls:58-69.

   SCOPE, honestly.  This is ONE SCALAR, calibrated on the isotropic rho-rho
   pedestal.  It cannot simultaneously be exact for the anisotropic term
   (eq:aniso-kernel) or for the dQ sector, whose lattice correction factors are
   different functions.  It buys direct comparability of the pedestal, which is
   the quantity dean.tex sec-cutoff-numerics tests, and nothing more.

   rescaleCutoff -> False restores the old behaviour (W built from lNoise
   itself).  Output files are ALWAYS tagged by the physical lNoise, so a sweep
   is keyed to the physics; lNoiseEff is reported alongside.

   ====================================================================
   BOUNDARY CONDITIONS  (dean.tex eq:bc-triple)
   ====================================================================
   Posing the problem on a bounded domain is THREE choices, not one: on the
   steady state, on the fluctuating field, and on the noise flux.  Only two of
   the eight combinations are coherent, and both are implemented:

     bcType = "Dirichlet"  RESERVOIR.  rho_ss pinned, psi|_bdry = 0, flux FREE.
                           A window cut from a larger colony.
     bcType = "Neumann"    SEALED.  n.Grad rho_ss = 0, n.Grad psi = 0, flux
                           BLOCKED, dQ strongly ANCHORED.  A container wall.

   The pairing is not a convention.  Whether psi is pinned at the wall and
   whether F.n vanishes there are the same question -- can a bacterium cross --
   asked of the mean and of the fluctuation; answering it differently for the
   two puts the wall out of detailed balance by fiat.

   WHY dQ IS ANCHORED AT A SEALED WALL and not left Neumann: sealing is a
   statement about MASS, and only rho is conserved, so the flux argument fixes
   the rho sector alone.  Leaving dQ Neumann too leaves A with exactly three
   unstable eigenvalues that CONVERGE rather than refine away (measured at
   Nint = 21..61 in ../claude_experiments/fluctuations_neumann/): the two defect
   TRANSLATION modes and the global director ROTATION.  A free-anchoring wall
   images a +1/2 defect attractively, so the centred defect is a saddle; a wall
   that blocks mass but exerts no torque cannot confine the defect, and no
   equal-time covariance exists to compute.  The anchored Laplacian L1q enters
   A22/A33 and nothing else.

   ====================================================================
   INHERITED, AND NOT TO BE "SIMPLIFIED"
   ====================================================================
   Every item here cost real debugging time in the modules this one replaces.

   CONSISTENT STENCILS (R == 1).  The DIAGONAL noise flux terms use Gx, Gy with
   coefficients sampled on the corresponding cell faces; the CROSS terms keep
   the COLLOCATED differences with node-sampled coefficients.  That split is
   forced, not a compromise: a corner-interpolated cross term picks up a phase
   cos((q-p)/2) from the half-cell offset between an x-face and a y-face, which
   destroys the parity in p behind the continuum angular cancellation for
   traceless Q_ss and injects a SPURIOUS anisotropic pedestal
   0.087 Q2 zeta/(B rho0 h^2) -- an ~8% h^-2 modulation at S ~ 0.86.  Measured
   Brillouin-zone averages:

       scheme                             c_iso    c_Q1   c_Q2
       central everywhere                 0.36338    0      0
       fully staggered, corner cross      1          0      0.0870
       staggered diagonal + central cross 1          0      0

   Setting every coefficient to 1 makes the first two noise terms exactly -Lap,
   which is what forces R == 1 and the flat pedestal.  Asserted in [2] as
   lapIdentity == 0 -- EXACT in rationals, hence the "== 0" test, hence G1/L1/D1
   are built from exact rationals and only converted with
   Developer`ToPackedArray @ N @ Normal at [6].  Unpacked dense arrays lose
   every BLAS path and the back-substitute runs ~600x slower.

   DIVERGENCE FORM for A21/A31, -(B/xi0) sum_i Gi^T diag(Q_ss|face) Gi, rather
   than the product-rule expansion (Q1V Lap + Q1x Dx + Q1y Dy).  The two differ
   by 2.8% at Nint = 41 and only the divergence form is exactly self-adjoint, as
   Div(c Grad .) is in the continuum.

   EVERY DISCRETE DIVERGENCE IS MINUS A TRANSPOSED GRADIENT, B_i in {-Gi^T,
   -Di^T}, so 1^T B_i = -(Gi.1)^T = 0 identically.  Under Dirichlet this is
   invisible because D1^T = -D1; under Neumann it is the difference between
   conserving mass and creating it in the four corner cells at O(h^-2).  Same
   rule forces the symmetrised double divergence Dxy for the sealed case.

   Q -> W.Q.W IS APPLIED IN PLACE, ONE BLOCK PAIR AT A TIME, with pairwise
   symmetrisation.  Building all nine blocks and ArrayFlatten-ing them holds
   ~42 GiB of temporaries at Nint = 131; a whole-matrix (Q + Q^T)/2 adds two
   more 21 GiB copies.  Folding W into Gx/Dx instead would make them dense and
   cost O(Nint^6).

   THE SEALED NULL MODE.  v = (1,0,0) is an exact right AND left null vector of
   A, and Q.v = 0, so lambda = 0 is exact and the Lyapunov equation is solvable
   but NOT unique: C is fixed only up to c v v^T.  A sealed box conserves the
   total exactly, so the uniform mode has ZERO variance and the wanted solution
   is the canonical one.  eigLyapSolve zeroes row and column iNull of the
   transformed right-hand side, which makes the refinement step safe too, and
   shifts that one denominator off zero -- 0./0. is Indeterminate, not 0, and a
   single such entry is smeared over all of C by the two matmuls.

   ====================================================================
   NAMES
   ====================================================================
   Entry point SolveFluctuationsLyapunov.  This module opens its OWN private
   context, LyapunovSolver`Private`, and its public export names other than
   MergeParams/ScriptParamOverrides (which are identical in every module) do not
   collide with the old ones -- so it can be loaded alongside
   ../lyapunov_solver_module_cutoff.m, which verify.wls needs for its regression
   test.  The one name to watch is SolveFluctuationsLyapunov itself, which the
   superseded ../lyapunov_solver_module.m also defines; do not load that file.

   Grid ordering is k = i + (j-1) Nint -- x fastest -- matching the rest of the
   repo. *)

SetSystemOptions["ParallelOptions" -> "ParallelThreadNumber" -> $ProcessorCount];

Clear[SolveFluctuationsLyapunov, visualizeFluctuationsComoving,
      exportDataComoving, MergeParams, ScriptParamOverrides,
      LyapunovReferenceKernel];


SolveFluctuationsLyapunov::usage =
  "SolveFluctuationsLyapunov[lyapunovParams, modelParams, steadyState, bcType] " <>
  "solves A \[CapitalSigma] + \[CapitalSigma] A^T + Q = 0 for the equal-time " <>
  "covariance of (\[Delta]\[Rho], \[Delta]Q1, \[Delta]Q2) linearized about the " <>
  "COMOVING steady state, on a uniform Nint x Nint grid over [-fbox, fbox]^2.\n\n" <>
  "  lyapunovParams : {Nint->..., fbox->...} (positive integers), optionally " <>
  "rescaleCutoff->True, advection->True\n" <>
  "  modelParams    : {B, a, b, L, \[Zeta], \[Xi]0, \[Xi]r, \[Rho]0, " <>
  "\[CapitalLambda]} plus lNoise >= 0, the PHYSICAL noise correlation length\n" <>
  "  steadyState    : from SolveActiveNematicComovingSteady; must carry \"u\"\n" <>
  "  bcType         : \"Dirichlet\" (reservoir wall) or \"Neumann\" (sealed, " <>
  "anchoring wall) -- the two consistent rows of dean.tex eq:bc-triple\n\n" <>
  "The drift includes the comoving advection Div(\[Delta]Q w_ss) with " <>
  "w_ss = u + (B/\[Xi]0)\[Del]\[Rho]_ss (dean.tex eq. linearized-Q-equation), " <>
  "which the superseded ../lyapunov_solver_module*.m files omit.\n\n" <>
  "With rescaleCutoff->True (default) the smoothing length is calibrated from " <>
  "lNoise so that the lattice pedestal matches the continuum " <>
  "\[CapitalSigma]_iso = \[Zeta]\[Rho]_ss/(2\[Pi] B \[Rho]0 lNoise^2) of " <>
  "dean.tex eq:iso-gauss; the solved-for length is returned as \"lNoiseEff\". " <>
  "Requires lNoise >= h/Sqrt[2\[Pi]] (dean.tex eq:ln-eff); below that no " <>
  "smoothing can reach the target and the call aborts.\n\n" <>
  "Returns <|\"C\", \"sigmaRho\", \"sigmaRhoFn\", \"sigmaQ1\", \"sigmaQ2\", " <>
  "\"gridInt\", \"h\", \"Nint\", \"fbox\", \"bcType\", \"u\", \"lNoise\", " <>
  "\"lNoiseEff\", \"lNoiseOverH\", \"lNoiseEffOverH\", \"rescaleCutoff\", " <>
  "\"advection\", \"peclet\", \"advIdentity\", \"Wsymbol\", \"eigenvalues\", " <>
  "\"maxReLambda\", \"residual\", \"residualRefined\", \"minEigLocal\", " <>
  "\"minEigQ\", \"lapIdentity\", \"scheme\", \"modelParams\", \"lyapunovParams\", " <>
  "\"timings\"|>, plus the sealed-wall diagnostics \"nullIndex\", " <>
  "\"nullLambdaRel\", \"nullOverlap\", \"consA\", \"consAT\", \"consQ\".";

SolveFluctuationsLyapunov::badgrid =
  "Nint and fbox must be positive integers (h is kept exact); got `1`.";
SolveFluctuationsLyapunov::badbc =
  "Unknown bcType `1`; expected \"Dirichlet\" (reservoir) or \"Neumann\" (sealed).";
SolveFluctuationsLyapunov::badss =
  "steadyState must be an association with keys \"rho\", \"Q1\", \"Q2\"; got `1`.";
SolveFluctuationsLyapunov::nou =
  "steadyState carries no \"u\": it is not a comoving steady state.  This module \
linearizes dean.tex eq. Q-comoving, whose drift contains the frame velocity, so a \
lab-frame profile would be a state that is not its own fixed point.  Use \
comoving_defect/comoving_steady_solver.m.";
SolveFluctuationsLyapunov::badlnoise =
  "modelParams must carry a numeric lNoise >= 0 (the PHYSICAL noise correlation \
length); got `1`.";
SolveFluctuationsLyapunov::uncalibratable =
  "lNoise = `1` is below h/Sqrt[2 Pi] = `2`, so the bare grid already produces LESS \
variance than the continuum theory at this cutoff and no amount of smoothing can \
reach it (dean.tex eq:ln-eff).  Raise Nint to at least `3` at this fbox, or raise \
lNoise, or set rescaleCutoff -> False to report the grid-limited value instead.";
SolveFluctuationsLyapunov::marginal =
  "lNoise = `1` is only `2` grid spacings.  The pedestal is still calibrated, but \
the run is grid-dominated: the anisotropic term and the \[Delta]Q sector carry \
lattice corrections that this single scalar does not remove.";
SolveFluctuationsLyapunov::smallbox =
  "fbox = `1` exceeds the steady state's own box `2`, so the steady interpolant \
would be extrapolated near the wall.";
SolveFluctuationsLyapunov::peclet =
  "Cell Peclet number `1` exceeds 1: the advection Div(\[Delta]Q w_ss) is no longer \
resolved by centred differencing on this grid and the \[Delta]Q sector may show \
oscillations.  Refine.";
SolveFluctuationsLyapunov::advid =
  "The advection identity Adv.1 = (B/\[Xi]0) \[Del]^2\[Rho]_ss is off by `1` \
relative in the interior.  The two sides are different discretisations of \
\[Del].w_ss, so they disagree at O((h/\[ScriptL]_d)^2) and a few percent is \
normal on a coarse grid -- but this is large enough to suspect the face \
sampling of w_ss rather than the resolution.  It must fall like h^2: check the \
trend in Nint before acting on it.";
SolveFluctuationsLyapunov::notcons =
  "The sealed-wall conservation identity `1` failed: `2` against a matrix scale of \
`3`.  The discrete divergence is leaking mass, so the conserved mode is driven and \
the Lyapunov equation has no solution.";
SolveFluctuationsLyapunov::nullmode =
  "Expected exactly one null eigenvalue (the conserved uniform mode) but found \
|\[Lambda]_null|/max|\[Lambda]| = `1`, second smallest = `2`, overlap with the \
uniform mode = `3`.";
SolveFluctuationsLyapunov::unstable =
  "Linearization is not strictly stable: max Re(\[Lambda]) = `1`; results may be \
unphysical.  On coarse grids this is usually the defect core being under-resolved \
and falls like h^3 -- read the trend in Nint, not a single grid.";
SolveFluctuationsLyapunov::notpsd =
  "Condition (ii) Q >= 0 is violated even with the exact moments at `1` grid \
points; check the closure assembly.";
SolveFluctuationsLyapunov::norefine =
  "Refinement did not improve the residual (`1` -> `2`); consider LyapunovSolve.";
SolveFluctuationsLyapunov::kappa =
  "\[Kappa](S) root-find did not converge at some grid points (max |I1/I0 - S| = `1`).";
SolveFluctuationsLyapunov::badkernel =
  "The Lyapunov solve kernel returned `1` instead of an association carrying the \
keys `2`.  $LyapunovSolveKernel is set to a backend that does not honour the \
kernel contract; see LyapunovReferenceKernel::usage.";


(* ---------------------------------------------------------------------------
   The solve kernel seam.

   Block [6] -- the O(n^3) part, and ~98% of the runtime -- is reached through
   the function held in $LyapunovSolveKernel rather than inlined, so an
   alternative backend can replace it WITHOUT a second copy of the physics in
   blocks [1]-[5].  Everything upstream (stencils, closure, cutoff congruence,
   boundary conditions) and downstream (extraction, tagging, export) is shared.

   Loading lyapunov_solver.m always resets the kernel to the reference
   implementation, so a file that installs another backend must be loaded
   AFTER this one.
   --------------------------------------------------------------------------- *)

$LyapunovKernelWantsDense::usage =
  "$LyapunovKernelWantsDense -> False tells block [6] NOT to build the dense " <>
  "Adense/Qdense copies, and to hand the kernel the assembled Atilde (a " <>
  "SparseArray) and Qtilde directly instead.  An out-of-core backend that " <>
  "streams them to disk sets this to save 2 x 8 n^2 bytes.  Default True.";

$LyapunovSolveKernel::usage =
  "$LyapunovSolveKernel holds the function used for the dense Lyapunov solve, " <>
  "block [6] of SolveFluctuationsLyapunov.  Defaults to LyapunovReferenceKernel; " <>
  "set it to install another backend (see lyapunov_solver_optimized.m).  Reset " <>
  "to the default every time lyapunov_solver.m is loaded.";

LyapunovReferenceKernel::usage =
  "LyapunovReferenceKernel[Adense, Qdense, meta] solves " <>
  "A \[CapitalSigma] + \[CapitalSigma] A^T + Q = 0 by dense eigendecomposition " <>
  "plus one refinement step, and is the default value of $LyapunovSolveKernel.\n\n" <>
  "This is the KERNEL CONTRACT that any replacement backend must honour.\n" <>
  "  Adense, Qdense : packed Real64 n x n matrices, n = 3 Nint^2; Q symmetric\n" <>
  "  meta           : <|\"neumann\", \"nInt\", \"n\", \"vCons\", \"scratchDir\", " <>
  "\"bcType\"|>\n" <>
  "Returns <|\"C\" (Real n x n), \"eigenvalues\", \"maxReLambda\", \"residual\", " <>
  "\"residualRefined\", \"nullIndex\", \"nullLambdaRel\", \"nullLambda2Rel\", " <>
  "\"nullOverlap\", \"timings\" ({{label, seconds}, ...}), \"backend\" (a " <>
  "string)|>.\n\n" <>
  "Under Neumann the drift is singular by design (the conserved uniform mode), " <>
  "so the kernel must project that mode out of the transformed right-hand side " <>
  "(dean.tex eq:canonical) and report the four null-* diagnostics; under " <>
  "Dirichlet those four come back Missing[\"Dirichlet\"].";

visualizeFluctuationsComoving::usage =
  "visualizeFluctuationsComoving[correlations, outputDir, plotBox] returns " <>
  "{plot2D, plot3D} of Var(\[Delta]\[Rho])(r) over [-plotBox, plotBox]^2 and writes " <>
  "them as PNGs into outputDir.  outputDir = None plots without saving.";

exportDataComoving::usage =
  "exportDataComoving[correlations, outputDir, saveCovariance] writes " <>
  "sigma_rho_<tag>.m into outputDir and, when saveCovariance is True, the full " <>
  "covariance as row-major Real64 covariance_<tag>.bin plus a meta file.  The tag " <>
  "carries fbox, Nint, the wall, and the PHYSICAL lNoise.  Returns the list of " <>
  "written paths; outputDir = None writes nothing.";

MergeParams::usage =
  "MergeParams[defaults, overrides] returns defaults with the keys given in " <>
  "overrides replaced, keeping every key not mentioned. " <>
  "MergeParams[defaults, overrides, label] aborts on a key absent from defaults.";

ScriptParamOverrides::usage =
  "ScriptParamOverrides[str] parses \"name1={k->v,...}; name2={...}\" into " <>
  "<|\"name1\"->{k->v,...}, ...|> without evaluating the left-hand sides. " <>
  "ScriptParamOverrides[] reads $ScriptCommandLine[[2]], or returns <||> when no " <>
  "argument was given.";


Begin["LyapunovSolver`Private`"];

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
(*Cutoff calibration (dean.tex eq:lattice-bz)*)

(* F(s) = Exp[-s^2] BesselI[0, s^2], written in x = s^2.
   BesselI[0, x] overflows as a machine number above x ~ 715 while Exp[-x]
   underflows below -745, so the machine product is only safe on a limited
   range -- and x = (l/h)^2 reaches a few thousand on fine grids.  Above the
   safe range evaluate at 50 digits: the two factors are ~1e+/-1085 there and
   their product is O(1e-2), which arbitrary precision tracks exactly.  It is
   called a few dozen times per run, so the cost is irrelevant. *)
scaledI0[x_?NumericQ] := Which[
  x <= 0., 1.,
  x <= 50., Exp[-N[x]] BesselI[0, N[x]],
  True, With[{xp = SetPrecision[x, 50]}, N[Exp[-xp] BesselI[0, xp]]]];

(* Solve scaledI0[x] == c for x by bisection.  scaledI0 is smooth and strictly
   decreasing from 1 to 0, so bisection is bulletproof and 200 evaluations cost
   nothing; FindRoot on a piecewise NumericQ function is not worth the risk. *)
invScaledI0[c_?NumericQ] := Module[{lo = 0., hi = 1., mid, k = 0},
  If[c >= 1., Return[0.]];
  While[scaledI0[hi] > c && hi < 1.*^8, hi *= 2.];
  Do[mid = (lo + hi)/2;
     If[scaledI0[mid] > c, lo = mid, hi = mid];
     If[hi - lo <= 1.*^-14 (1. + hi), Break[]],
     {k, 200}];
  (lo + hi)/2];

(* The physical lNoise -> the length to build W from.  Returns $Failed when the
   target is unreachable, i.e. t < 1/Sqrt[2 Pi] (dean.tex eq:ln-eff). *)
cutoffLength[lPhys_?NumericQ, hV_?NumericQ] := Module[{t, c},
  If[lPhys == 0., Return[0.]];
  t = lPhys/hV;
  c = 1./(Sqrt[2. Pi] t);
  If[c >= 1., Return[$Failed]];
  hV Sqrt[invScaledI0[c]]];


(* ::Subsection:: *)
(*Output helpers*)

prepareDir[outputDir_] := Which[
  outputDir === None, None,
  !StringQ[outputDir],
    Print["ERROR: outputDir must be a string or None; got ", outputDir]; Abort[],
  True, Quiet @ CreateDirectory[outputDir, CreateIntermediateDirectories -> True];
        If[!DirectoryQ[outputDir],
           Print["ERROR: cannot create outputDir ", outputDir]; Abort[]];
        ExpandFileName[outputDir]];

savePlot[dir_, name_, plot_] :=
  If[dir === None, None,
     With[{p = FileNameJoin[{dir, name}]},
       Export[p, plot, ImageResolution -> 150]; Print["Saved: ", p]; p]];

(* Compact decimal: up to six places, never scientific, trailing zeros trimmed.
   1.73 -> "1.73", 8. -> "8", 0.8 -> "0.8", 15 -> "15".  ExponentFunction is
   load-bearing -- without it NumberForm switches to a superscripted "m x 10^e"
   for small values, which is not a filename.

   ONE regex, not two rules.  StringReplace applies a rule LIST in a single
   left-to-right pass, so stripping "0+$" and then "\\.$" as two rules leaves
   "8." -- the second rule never sees the string the first produced.
   "\\.?0+$" takes the trailing zeros and the now-bare point together.
   Six decimals is the resolution: two values closer than 1e-6 tag alike. *)
numTag[x_] := Module[
  {s = ToString@NumberForm[N[x], {12, 6}, ExponentFunction -> (Null &)]},
  If[StringContainsQ[s, "."],
     s = StringReplace[s, RegularExpression["\\.?0+$"] -> ""]];
  If[s === "" || s === "-", s = "0"];
  s];

(* The tag has to separate every run in a sweep, so it carries the five things
   this project varies: grid, box, wall, and the three physical scales.  The
   PHYSICAL lNoise goes in, not the calibrated lNoiseEff -- a sweep is keyed to
   the physics, not to whatever length the lattice correction happened to pick.

   Returned to the caller as correlations["runTag"] so that drivers writing
   their own side files (steady_*.m) use THIS string rather than rebuilding it.
   The previous set2_2_comove driver did rebuild it, with a comment reading
   "KEEP THE TWO IN STEP"; when they drift the steady dump silently stops
   pairing with its covariance. *)
(* B and ssBox are appended only when supplied, so a tag built without them is
   byte-identical to the pre-2026-09 form and old output stays parseable.  They
   MUST be supplied whenever either is swept: neither enters elld (K' absorbs
   zeta, and B enters no derived length at all), so without them a run at
   B = 10^4 and one at 3*10^5 collide on one filename and the second silently
   overwrites the first.  Same for the steady state's own box, now that it is
   varied independently of fbox. *)
runTag[corr_] := "box" <> numTag[corr["fbox"]] <>
                 "_N"  <> ToString[corr["Nint"]] <>
                 If[corr["bcType"] === "Neumann", "_Neu", "_Dir"] <>
                 "_lN" <> numTag[corr["lNoise"]] <>
                 "_z"  <> numTag[corr["zeta"]] <>
                 "_ld" <> numTag[corr["elld"]] <>
                 If[NumericQ @ Lookup[corr, "B", Missing[]],
                    "_B" <> numTag[corr["B"]], ""] <>
                 If[NumericQ @ Lookup[corr, "ssBox", Missing[]],
                    "_ss" <> numTag[corr["ssBox"]], ""] <>
                 (* Last, so the numeric fields above still parse by prefix.
                    A run with max Re(lambda) >= 0 is NOT a covariance: the
                    Lyapunov equation has no stable solution there and the
                    variances can come out negative.  Marking it in the FILENAME
                    means such a run cannot be swept up by a glob and averaged in
                    by accident -- the diagnostic is in the meta file either way,
                    but nobody reads that before plotting. *)
                 If[TrueQ @ Lookup[corr, "unstable", False], "_unstable", ""];

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

(* One conservation gate.  These are the sealed-wall analogue of lapIdentity: if
   any fails the scheme is void, so abort rather than produce a plausible-looking
   covariance. *)
consCheck[lbl_, val_, scale_] := (
  Print["    ", StringPadRight[lbl, 40], " = ", val,
        If[Abs[val] <= 1.*^-10 scale, "   (exact to round-off)",
           "   *** CONSERVATION VIOLATED ***"]];
  If[!(Abs[val] <= 1.*^-10 scale),
     Message[SolveFluctuationsLyapunov::notcons, lbl, val, scale];
     Print["    ABORT: the sealed-wall scheme is void."];
     Abort[]];
  val);


(* ::Subsection:: *)
(*Solve kernel -- reference (dense eigendecomposition) backend*)

(* LyapunovSolve is single-threaded; diagonalize instead: A = P D P^-1 =>
   C = P Y P^T, Y_ij = -(P^-1 Q P^-T)_ij/(lam_i + lam_j).  A is NOT self-adjoint
   (dropped Onsager partners, and now the advection too), but the error is set
   by kappa(P)^2 eps ~ 1e-10; one refinement step then reaches ~1e-14.

   Lifted verbatim out of block [6] when the kernel seam was introduced -- the
   arithmetic is unchanged, so this remains the reference every other backend is
   validated against. *)
LyapunovReferenceKernel[Adense_, Qdense_, meta_Association] :=
Module[
  {neuQ, vCons, tim, rc, eigvals, eigvecs, tEs, lamMax, iNull, nullRel,
   lam2Rel, vUnit, nullOverlap, maxReLam, Pev, lamSum, eigLyapSolve,
   Cmat, tSol, resid, tRes1, tRef, residRef, tRes3},

  neuQ  = TrueQ @ meta["neumann"];
  vCons = meta["vCons"];
  tim   = {};
  rc[lbl_, t_] := (AppendTo[tim, {lbl, t}]; t);

  {tEs, {eigvals, eigvecs}} = AbsoluteTiming @ Eigensystem[Adense];
  rc["[6] Eigensystem (dense, nonsymmetric)", tEs];
  Print["    [a] Eigensystem        ", tEs, " s"];

  lamMax = Max @ Abs @ eigvals;
  If[neuQ,
     (* identify the conserved mode and check it is the only one *)
     iNull = First @ Ordering[Abs[eigvals], 1];
     nullRel = Abs[eigvals[[iNull]]]/lamMax;
     lam2Rel = Sort[Abs @ eigvals][[2]]/lamMax;
     vUnit = Normalize @ vCons;
     nullOverlap = Abs[Normalize[eigvecs[[iNull]]] . vUnit];
     Print["    conserved mode: index ", iNull, ",  |\[Lambda]|/max|\[Lambda]| = ",
           nullRel, ",  2nd smallest = ", lam2Rel];
     Print["    overlap of its eigenvector with uniform \[Delta]\[Rho] = ", nullOverlap,
           If[nullRel < 1.*^-8 && lam2Rel > 1.*^-6 && nullOverlap > 0.99,
              "   (a single, clean null mode)", "   *** UNEXPECTED NULL STRUCTURE ***"]];
     If[!(nullRel < 1.*^-8 && lam2Rel > 1.*^-6 && nullOverlap > 0.99),
        Message[SolveFluctuationsLyapunov::nullmode, nullRel, lam2Rel, nullOverlap]];
     maxReLam = Max @ Re @ Delete[eigvals, iNull],
  (* else *)
     iNull = None; nullRel = lam2Rel = nullOverlap = Missing["Dirichlet"];
     maxReLam = Max @ Re @ eigvals];

  Print["    max Re(\[Lambda])", If[neuQ, " off the null mode", ""], " = ", maxReLam,
        If[maxReLam < 0, "  (stable)", "  (UNSTABLE!)"]];
  If[maxReLam >= 0, Message[SolveFluctuationsLyapunov::unstable, maxReLam]];
  Print["    min |\[Lambda]_i+\[Lambda]_j|",
        If[neuQ, " off the null mode", ""], " = ",
        2 Min @ Abs @ Re @ If[neuQ, Delete[eigvals, iNull], eigvals],
        ",  complex fraction = ",
        N[100 Count[eigvals, z_ /; Abs[Im[z]] > 1.*^-8 Abs[z]]/Length[eigvals]], "%"];

  Pev = Developer`ToPackedArray @ Transpose @ eigvecs;
  lamSum = Developer`ToPackedArray @ Outer[Plus, eigvals, eigvals];
  If[neuQ,
     (* lam_null is 0 BY DESIGN, so lamSum[[iNull,iNull]] can come back EXACTLY 0.,
        and then 0./0. is Indeterminate -- not 0 -- and a single such entry is
        smeared over all of Cmat by the two matmuls.  Shift that one denominator
        off zero; the numerator there is exactly 0., so the quotient is 0 either
        way and the canonical projection is untouched.  ADD rather than assign,
        and add lamMax (a machine Real): the sum keeps the element type, so a
        packed Complex stays packed Complex. *)
     lamSum[[iNull, iNull]] = lamSum[[iNull, iNull]] + lamMax];

  (* Solve A X + X A^T = -Rhs in the eigenbasis; reused by the refinement.  Under
     Neumann the conserved mode is projected out by zeroing row and column iNull
     of the transformed right-hand side (dean.tex eq:canonical); doing it INSIDE
     eigLyapSolve makes the refinement step safe too.  Multiply by 0. rather than
     ASSIGNING 0.: assigning a real into a packed COMPLEX array unpacks it. *)
  eigLyapSolve[Rhs_] := Module[{num},
     num = -LinearSolve[Pev, Transpose @ LinearSolve[Pev, Transpose[Rhs]]];
     If[neuQ,
        num[[iNull]] = num[[iNull]] 0.;
        num[[All, iNull]] = num[[All, iNull]] 0.];
     Re[Pev . (num / lamSum) . Transpose[Pev]]];

  {tSol, Cmat} = AbsoluteTiming @ eigLyapSolve[Qdense];
  rc["[6] back-substitute", tSol];
  Print["    [b] back-substitute    ", tSol, " s"];

  (* Each residual is two dense n x n matmuls, same O(n^3) as the Eigensystem. *)
  {tRes1, resid} = AbsoluteTiming[
     Max @ Abs[Adense . Cmat + Cmat . Transpose[Adense] + Qdense]];
  rc["[6] residual eval #1 (unrefined)", tRes1];
  Print["    residual (unrefined)   = ", resid, "   [", tRes1, " s]"];
  (* Non-numeric here means the solve produced Indeterminate/NaN entries.  Stop:
     the refinement step and exportDataComoving would otherwise write the whole
     unevaluated expression to disk (measured on the old modules: 3.8 GB of .m
     and 12.3 GB of log).  Every other gate tests the INPUTS; this is the only
     one that tests that a number came out. *)
  If[!NumberQ[resid],
     Print["    ABORT: residual is ", resid,
           " -- the Lyapunov solve returned non-numeric entries."];
     Abort[]];

  {tRef, Cmat} = AbsoluteTiming[
     Cmat + eigLyapSolve[Adense . Cmat + Cmat . Transpose[Adense] + Qdense]];
  rc["[6] refinement step (incl. residual eval #2)", tRef];
  {tRes3, residRef} = AbsoluteTiming[
     Max @ Abs[Adense . Cmat + Cmat . Transpose[Adense] + Qdense]];
  rc["[6] residual eval #3 (verification)", tRes3];
  Print["    [c] refinement step    ", tRef, " s"];
  Print["    residual (refined)     = ", residRef, "   [", tRes3, " s]",
        "   (improved ", ToString@NumberForm[N[resid/residRef], 4], "x)"];
  If[residRef > resid,
     Message[SolveFluctuationsLyapunov::norefine, resid, residRef]];
  Print["    total solve time = ", tEs + tSol + tRef, " s"];

  <| "C" -> Cmat,
     "eigenvalues" -> eigvals,
     "maxReLambda" -> maxReLam,
     "residual" -> resid,
     "residualRefined" -> residRef,
     "nullIndex" -> iNull,
     "nullLambdaRel" -> nullRel,
     "nullLambda2Rel" -> lam2Rel,
     "nullOverlap" -> nullOverlap,
     "timings" -> tim,
     "backend" -> "Mathematica/eigendecomposition" |>
];

$LyapunovSolveKernel = LyapunovReferenceKernel;
$LyapunovKernelWantsDense = True;


(* ::Subsection:: *)
(*Main solver*)

SolveFluctuationsLyapunov[
    lyapunovParams : {___Rule},
    modelParams : {___Rule},
    steadyState_Association,
    bcType_String : "Dirichlet"] :=
Module[
  {
    nInt, fboxV, neuQ, rescaleQ, advectQ, scratchV, timings = {}, rec,
    kernelOut, kernelKeys,
    Bv, xi0, xir, zet, rho0v, Lam, av, bv, Lv, cL, kapQ,
    aPrime, elldV, S0V, alphaV, derived, tagStr,
    rhoSs, Q1Ss, Q2Ss, uVec, ux, uy, ssBox,
    h, gridFull, gridInt, iZero, xFace, nFace,
    d1Full, l1Full, D1, L1, L1q, Pav, Idn, Dx, Dy, Lx, Ly, Dxy, Lap, Lapq,
    G1, P1, Gx, Gy, Px, Py, lapIdentity, anchorGap,
    lNoiseV, lEff, tRatio, nMin, Wsym, Emat, Esym, Eone,
    xflat, yflat, sample, diag,
    xeX, xeY, yeX, yeY, sampleXE, sampleYE,
    rhoV, Q1V, Q2V, rhoLap, rhoXd, rhoYd, wxE, wyE, peclet,
    pedPred, pedCont, boxFactor,
    rhoXE, Q1XE, Q2XE, rhoYE, Q1YE, Q2YE,
    SXE, SYE, g2XE, g3XE, g2YE, g3YE, e2XE, e3XE, e2YE, e3YE,
    h2cXE, h2sXE, h3cXE, h3sXE, h2cYE, h2sYE, h3cYE, h3sYE,
    Svals, g2V, g3V, kapResid, tSamp, tClo, tKchk, tCond,
    minEigLoc, scaleLoc, minEigQ, tQpsd,
    Adv, advIdent, intMask,
    A11, A12, A13, A21, A22, A23, A31, A32, A33, Atilde, tAsmA,
    vCons, scaleA, consA, consAT, consQ, scaleQ,
    kap, pf, locPf, divBlk, h2c, h2s, h3c, h3s, e2, e3,
    Qrr, QQ1Q1, QQ2Q2, QQ1Q2, QQ2Q1, QQ1r, QQ2r, QrQ1, QrQ2, Qtilde, tAsmQ,
    smoothBlock, rng, tSmooth, Qasym,
    Adense, Qdense, tDense, eigvals, eigvecs, tEs, maxReLam, Pev, lamSum,
    lamMax, iNull, nullRel, lam2Rel, vUnit, nullOverlap,
    eigLyapSolve, Cmat, tSol, resid, tRes1, tRef, residRef, tRes3,
    sigRhoVec, sigRhoMat, sigRhoFn, sigQ1, sigQ2, tExtr
  },

  rec[lbl_, t_] := (AppendTo[timings, {lbl, t}]; t);

  (* --- unpack + validate ---------------------------------------------
     The keys on the right MUST be Global`: the caller builds the rule lists in
     Global`, and a bare name here would resolve into this private context
     instead and match nothing.  Lookup against an Association, not the
     "{k,...} = {k,...} /. params" idiom: the replacement form self-assigns
     (fbox = fbox) when a key is absent, and a self-referential OwnValue blows
     the recursion limit the next time the symbol is touched. *)
  If[!MemberQ[{"Dirichlet", "Neumann"}, bcType],
    Message[SolveFluctuationsLyapunov::badbc, bcType]; Return[$Failed]];
  neuQ = (bcType === "Neumann");

  With[{lp = Association[lyapunovParams]},
    nInt     = Lookup[lp, Global`Nint, Missing[]];
    fboxV    = Lookup[lp, Global`fbox, Missing[]];
    rescaleQ = TrueQ @ Lookup[lp, Global`rescaleCutoff, True];
    advectQ  = TrueQ @ Lookup[lp, Global`advection, True];
    (* Only the out-of-core backends read this; the reference kernel ignores it.
       Default to the node-local scratch SLURM hands out (TMPDIR), which on these
       nodes is NVMe under /scratch -- NOT the shared GPFS the output lands on. *)
    scratchV = Lookup[lp, Global`scratchDir, Automatic];
  ];
  If[scratchV === Automatic,
     scratchV = SelectFirst[
        {Environment["TMPDIR"], "/scratch/" <> ToString @ Environment["SLURM_JOB_ID"],
         "/scratch", $TemporaryDirectory},
        StringQ[#] && DirectoryQ[#] &, $TemporaryDirectory]];
  If[!IntegerQ[nInt] || !IntegerQ[fboxV] || nInt <= 0 || fboxV <= 0,
    Message[SolveFluctuationsLyapunov::badgrid, {fboxV, nInt}]; Return[$Failed]];

  If[!AllTrue[{"rho", "Q1", "Q2"}, KeyExistsQ[steadyState, #] &],
    Message[SolveFluctuationsLyapunov::badss, Keys[steadyState]]; Return[$Failed]];
  (* The frame velocity is part of the drift, so a lab-frame steady state is not
     merely less accurate here -- it is a profile that is not its own fixed
     point.  Refuse it rather than silently reproducing the old inconsistency. *)
  If[!KeyExistsQ[steadyState, "u"],
    Message[SolveFluctuationsLyapunov::nou]; Return[$Failed]];
  uVec = N @ steadyState["u"];
  If[!(VectorQ[uVec, NumericQ] && Length[uVec] >= 2),
    Message[SolveFluctuationsLyapunov::nou]; Return[$Failed]];
  {ux, uy} = uVec[[1 ;; 2]];

  With[{mp = Association[modelParams]},
    Bv    = Lookup[mp, Global`B, Missing[]];
    av    = Lookup[mp, Global`a, Missing[]];
    bv    = Lookup[mp, Global`b, Missing[]];
    Lv    = Lookup[mp, Global`L, Missing[]];
    zet   = Lookup[mp, Global`\[Zeta], Missing[]];
    xi0   = Lookup[mp, Global`\[Xi]0, Missing[]];
    xir   = Lookup[mp, Global`\[Xi]r, Missing[]];
    rho0v = Lookup[mp, Global`\[Rho]0, Missing[]];
    Lam   = Lookup[mp, Global`\[CapitalLambda], Missing[]];
    lNoiseV = Lookup[mp, Global`lNoise, Missing[]];
  ];
  cL   = Lv + zet xir/(4 xi0);    (* K', the elastic constant in dH *)
  kapQ = 4 cL/xir;                (* the dQ diffusivity, for the Peclet number *)

  (* --- the activity-renormalised block (dean.tex primed-values, :665-672) ----
     a' = a + Lambda xir and K' = K + zeta xir/(4 xi0), so the defect core size
     and the far-field order are elld = Sqrt[K'/(-a')] and S0 = Sqrt[-a'/(2b)].
     Derived here rather than taken on trust: they are what the tag is keyed to,
     and a caller who edits zeta without re-deriving K would otherwise produce
     two different runs under one filename.  NB varying zeta at fixed elld moves
     K, not a' or K' -- K' absorbs the change exactly -- so what actually shifts
     is the activity alpha and the noise amplitude kapNoise. *)
  aPrime = av + Lam xir;
  elldV  = Sqrt[cL/(-aPrime)];
  S0V    = Sqrt[-aPrime/(2 bv)];
  alphaV = zet xir S0V/(4 xi0 cL);
  derived = <|"aPrime" -> aPrime, "KPrime" -> cL, "elld" -> elldV, "S0" -> S0V,
              "alpha" -> alphaV, "gamma0" -> Bv/xi0, "kappaQ" -> kapQ,
              "kappaNoise" -> 2 zet/(xi0 rho0v)|>;

  If[!NumericQ[lNoiseV] || N[lNoiseV] < 0,
    Message[SolveFluctuationsLyapunov::badlnoise, lNoiseV]; Return[$Failed]];
  lNoiseV = N[lNoiseV];

  {rhoSs, Q1Ss, Q2Ss} = steadyState /@ {"rho", "Q1", "Q2"};
  ssBox = Quiet @ Lookup[Lookup[steadyState, "mesh", <||>], "box", Missing[]];
  If[NumericQ[ssBox] && fboxV > ssBox,
     Message[SolveFluctuationsLyapunov::smallbox, fboxV, ssBox]];

  Print["[1] ", bcType, " wall (", If[neuQ, "sealed + anchoring", "reservoir"],
        "),  u = ", uVec[[1 ;; 2]],
        ",  advection ", If[advectQ, "ON", "OFF"],
        ",  cutoff rescaling ", If[rescaleQ, "ON", "OFF"]];

  (* --- grid -----------------------------------------------------------
     Dirichlet: NODE-centred, the wall is a node that is not in the state.
     Neumann:   CELL-centred, the wall is a FACE, so no coefficient is ever
                sampled on the boundary itself.
     Both keep h exact and put the defect core on a node for odd Nint.  Note
     Neumann Nint = M+1 has exactly the same h as Dirichlet Nint = M. *)
  If[neuQ,
     h = (2 fboxV)/nInt;
     gridInt = Table[-fboxV + (i - 1/2) h, {i, nInt}];
     xFace   = Table[-fboxV + m h, {m, 1, nInt - 1}],
  (* else *)
     h = (2 fboxV)/(nInt + 1);
     gridFull = Table[-fboxV + k h, {k, 0, nInt + 1}];
     gridInt  = gridFull[[2 ;; -2]];
     xFace    = Table[-fboxV + (m - 1/2) h, {m, nInt + 1}]];
  nFace = Length[xFace];

  Print["[2] Grid: fbox=", fboxV, ", Nint=", nInt,
        ", h=", ToString @ NumberForm[N[h], 5], ", state-dim=", 3 nInt^2,
        If[neuQ, "  (cell-centred, sealed walls)", "  (node-centred)"]];
  If[OddQ[nInt],
     iZero = (nInt + 1)/2;
     Print["    Nint odd  -> origin ON the grid at gridInt[[", iZero,
           "]] = ", gridInt[[iZero]], ", flat index k = ", iZero + (iZero - 1) nInt,
           " of ", nInt^2],
     iZero = None;
     Print["    Nint even -> origin NOT on the grid; nearest nodes at x = ",
           ToString @ NumberForm[N @ gridInt[[nInt/2]], 4], " and ",
           ToString @ NumberForm[N @ gridInt[[nInt/2 + 1]], 4]]];

  (* --- 1D operators, in EXACT rationals so lapIdentity can test "== 0" --- *)
  If[neuQ,
     (* G1: forward difference, Nint cell centres -> Nint-1 INTERIOR faces.  The
        two wall faces carry no flux and so do not appear (dean.tex eq:G-neumann).
        L1 is DEFINED as its Gram matrix, so the adjoint pair holds by
        construction.  D1 is the node<-face average of the same gradient, i.e.
        the reflective central difference; G1.1 = 0 and D1.1 = 0. *)
     G1 = SparseArray[Join[Table[{m, m + 1} ->  1/h, {m, nInt - 1}],
                           Table[{m, m}     -> -1/h, {m, nInt - 1}]], {nInt - 1, nInt}];
     L1 = SparseArray[-Transpose[G1] . G1];
     Pav = SparseArray[Join[Table[{m, m}     -> 1/2, {m, nInt - 1}],
                            Table[{m, m + 1} -> 1/2, {m, nInt - 1}]], {nInt - 1, nInt}];
     D1 = SparseArray[Transpose[Pav] . G1];
     P1 = Pav;
     (* The ANCHORED Laplacian for the dQ sector: the same finite-volume gradient
        with the two wall faces present, flux 2u/h over half a cell.  Equal to L1
        with its corners shifted by -2/h^2.  Negative definite, no null mode --
        which is the whole point: it is what confines the defect. *)
     L1q = SparseArray[L1 + SparseArray[{{1, 1} -> -2/h^2, {nInt, nInt} -> -2/h^2},
                                        {nInt, nInt}]],
  (* else: Dirichlet *)
     d1Full = NDSolve`FiniteDifferenceDerivative[Derivative[1], gridFull,
        "DifferenceOrder" -> 2]["DifferentiationMatrix"];
     l1Full = NDSolve`FiniteDifferenceDerivative[Derivative[2], gridFull,
        "DifferenceOrder" -> 2]["DifferentiationMatrix"];
     D1 = SparseArray @ d1Full[[2 ;; -2, 2 ;; -2]];
     L1 = SparseArray @ l1Full[[2 ;; -2, 2 ;; -2]];
     (* Forward difference, interior NODES -> the Nint+1 cell EDGES.  Edge m lies
        between nodes m-1 and m, with u_0 = u_{Nint+1} = 0.  The (Nint+1)-th row
        is what makes -G^T G = L1 hold EXACTLY: it carries the stencil entry an
        Nint x Nint collocated backward difference drops. *)
     G1 = SparseArray[Join[Table[{m, m} -> 1/h, {m, nInt}],
                           Table[{m + 1, m} -> -1/h, {m, nInt}]], {nInt + 1, nInt}];
     P1 = SparseArray[Join[Table[{m, m} -> 1/2, {m, nInt}],
                           Table[{m + 1, m} -> 1/2, {m, nInt}]], {nInt + 1, nInt}];
     L1q = L1];

  (* State ordering k = i + (j-1) N, x fastest.  Hence Dx = I (x) D1, Dy = D1 (x) I. *)
  Idn = IdentityMatrix[nInt, SparseArray];
  Dx = KroneckerProduct[Idn, D1];
  Dy = KroneckerProduct[D1, Idn];
  Lx = KroneckerProduct[Idn, L1];
  Ly = KroneckerProduct[L1, Idn];
  Gx = KroneckerProduct[Idn, G1];
  Gy = KroneckerProduct[G1, Idn];
  Px = KroneckerProduct[Idn, P1];
  Py = KroneckerProduct[P1, Idn];
  Lap = Lx + Ly;
  Lapq = If[neuQ, KroneckerProduct[Idn, L1q] + KroneckerProduct[L1q, Idn], Lap];
  (* d_x d_y.  Under Neumann it must be a symmetrised DOUBLE DIVERGENCE so its
     column sums vanish; that form equals Kron[D1,D1] identically whenever
     D1^T = -D1, which is what the Dirichlet central difference is.  So this is
     a generalisation, not a competing scheme. *)
  Dxy = If[neuQ,
           -(1/2) (KroneckerProduct[D1, Transpose[D1]] +
                   KroneckerProduct[Transpose[D1], D1]),
           KroneckerProduct[D1, D1]];

  (* The identity the whole scheme rests on.  Exact in rationals, so test == 0. *)
  lapIdentity = Max @ Abs @ Normal[
     Lap + (Transpose[Gx] . Gx + Transpose[Gy] . Gy)];
  Print["    |Lap + (Gx^T Gx + Gy^T Gy)| = ", lapIdentity,
        If[lapIdentity == 0, "   (exact -- adjoint pair)",
           "   *** NOT AN ADJOINT PAIR, R != 1 ***"]];
  If[lapIdentity != 0,
     Print["    ABORT: the staggered/compact identity failed; the scheme is void."];
     Abort[]];
  If[neuQ,
     anchorGap = -Max @ Eigenvalues[N @ Normal @ L1q];
     Print["    dQ anchored: L1q = L1 with corners -3/h^2; least-damped mode ",
           -anchorGap, " (strictly negative -- this is what confines the defect)"];
     Print["    Gx.1 = ", Max @ Abs @ Normal[Gx . ConstantArray[1, nInt^2]],
           ",  Dx.1 = ", Max @ Abs @ Normal[Dx . ConstantArray[1, nInt^2]],
           ",  1^T.Dxy = ", Max @ Abs @ Normal[ConstantArray[1, nInt^2] . Dxy],
           "   (all exact 0: no-flux walls)"]];

  (* --- [1b] calibrate the cutoff against the continuum theory ----------
     See "Cutoff calibration" in the header.  W is built from lEff; every output
     is tagged by the physical lNoise. *)
  tRatio = If[lNoiseV == 0., 0., lNoiseV/N[h]];
  If[lNoiseV == 0. || !rescaleQ,
     lEff = lNoiseV;
     If[lNoiseV > 0. && !rescaleQ,
        Print["[1b] rescaleCutoff -> False: W built from lNoise itself, so the ",
              "pedestal overshoots eq:iso-gauss by ",
              ToString@NumberForm[N[2 Pi tRatio^2 scaledI0[tRatio^2]^2], 6],
              "x (eq:lattice-bz)"]],
  (* else: calibrate *)
     lEff = cutoffLength[lNoiseV, N[h]];
     If[lEff === $Failed,
        nMin = Ceiling[2 fboxV/(Sqrt[2. Pi] lNoiseV)] - If[neuQ, 0, 1];
        Message[SolveFluctuationsLyapunov::uncalibratable,
                lNoiseV, N[h]/Sqrt[2. Pi], nMin];
        Return[$Failed]];
     Print["[1b] Cutoff calibrated (eq:lattice-bz): lNoise = ", lNoiseV,
           " = ", ToString@NumberForm[tRatio, 4], " h",
           "  ->  lNoiseEff = ", ToString@NumberForm[lEff, 8],
           " = ", ToString@NumberForm[lEff/N[h], 4], " h"];
     Print["     the lattice pedestal now equals ",
           "zeta rho_ss/(2 Pi B rho0 lNoise^2) of eq:iso-gauss; without this it ",
           "would overshoot by ",
           ToString@NumberForm[N[2 Pi tRatio^2 scaledI0[tRatio^2]^2], 6], "x"];
     If[tRatio < 0.6,
        Message[SolveFluctuationsLyapunov::marginal, lNoiseV,
                ToString@NumberForm[tRatio, 3]]]];

  (* --- the noise-correlation operator W = Kron[Emat, Emat] -------------
     W = MatrixExp[(lEff^2/4) Lap] and Lap = Lx + Ly with Lx.Ly = Ly.Lx, so the
     exponential splits EXACTLY into one Nint x Nint factor.  L1 is diagonalised
     ANALYTICALLY -- by the Dirichlet sine basis, or by the DCT-II basis for the
     cell-centred Neumann L1 -- so no numerical matrix exponential is ever
     needed and the mode weights underflow gracefully to 0 when lEff >> h.

     Neumann: p = 0 is the CONSTANT mode with mu_0 = 0, so Wsym[[1]] = 1 exactly
     and E.1 = 1 -- the congruence cannot inject noise into the conserved mode,
     and the kernel REFLECTS off the sealed wall.  Dirichlet: W is the Dirichlet
     heat kernel and its diagonal VANISHES at the wall, so the rods shrink as
     they approach it.  That contrast, and its sign, is the sharpest single
     check that the boundary condition is real. *)
  If[lEff == 0.,
     Wsym = ConstantArray[1., nInt];
     Emat = None;
     Print["    lNoise = 0 -> W = I; this run is the grid-regulated scheme"],
  (* else *)
     If[neuQ,
        Wsym = Exp[-(lEff/N[h])^2 Sin[N[Pi] Range[0, nInt - 1]/(2 nInt)]^2];
        Emat = With[{V = Table[Sqrt[If[p == 0, 1., 2.]/nInt] *
                               Cos[N[Pi] p (i - 1/2)/nInt],
                               {i, nInt}, {p, 0, nInt - 1}]},
                    V . (Wsym * Transpose[V])],
        Wsym = Exp[-(lEff/N[h])^2 Sin[N[Pi] Range[nInt]/(2 (nInt + 1))]^2];
        Emat = With[{V = Sqrt[2./(nInt + 1)] Table[
                          Sin[N[Pi] p i/(nInt + 1)], {i, nInt}, {p, nInt}]},
                    V . (Wsym * Transpose[V])]];
     Esym = Max @ Abs[Emat - Transpose[Emat]];
     Print["    W = Kron[E,E] on the ", If[neuQ, "NEUMANN (reflecting)", "DIRICHLET"],
           " Laplacian;  E symmetry residual = ", Esym, ", max |E| = ", Max@Abs@Emat];
     Print["    W mode weights \[Element] ", MinMax@Wsym];
     If[Esym > 1.*^-12, Print["    WARNING: E is not symmetric to 1e-12."]];
     If[neuQ,
        Eone = Max @ Abs[Emat . ConstantArray[1., nInt] - 1.];
        Print["    |E.1 - 1| = ", Eone, ",  Wsym[[1]] - 1 = ", Wsym[[1]] - 1.];
        If[Eone > 1.*^-12,
           Print["    ABORT: E.1 != 1, so W would drive the conserved mode."];
           Abort[]]]];

  (* --- sample the steady state ----------------------------------------
     Batched sampling, ~200x faster than scalar calls.  N@ is load-bearing:
     exact rationals bypass the fast path and are 19x slower. *)
  xflat = Flatten @ ConstantArray[N @ gridInt, nInt];
  yflat = Flatten @ Transpose @ ConstantArray[N @ gridInt, nInt];
  sample[f_] := f[xflat, yflat];

  (* The diagonal flux coefficients live on the cell faces, so sample the steady
     interpolant THERE rather than averaging node values -- it is defined
     everywhere, so this costs nothing in accuracy.
       x-faces: (face_x, node_y), nFace x Nint, face index fastest
       y-faces: (node_x, face_y), Nint x nFace, node index fastest
     matching Gx = Idn (x) G1 and Gy = G1 (x) Idn respectively. *)
  xeX = Flatten @ ConstantArray[N @ xFace, nInt];
  xeY = Flatten @ Transpose @ ConstantArray[N @ gridInt, nFace];
  yeX = Flatten @ ConstantArray[N @ gridInt, nFace];
  yeY = Flatten @ Transpose @ ConstantArray[N @ xFace, nInt];
  sampleXE[f_] := f[xeX, xeY];
  sampleYE[f_] := f[yeX, yeY];

  Print["[3] Sampling steady state ..."];
  tSamp = rec["[3] sample steady state", First @ AbsoluteTiming[
     rhoV = sample[rhoSs];
     Q1V = sample[Q1Ss];
     Q2V = sample[Q2Ss];
     (* rhoLap is no longer a term in A -- Adv subsumes it -- but it is still
        needed as the reference for the [4b] tripwire and for advection -> False. *)
     rhoLap = sample[Derivative[2, 0][rhoSs][##] + Derivative[0, 2][rhoSs][##] &];
     rhoXE = sampleXE[rhoSs];  Q1XE = sampleXE[Q1Ss];  Q2XE = sampleXE[Q2Ss];
     rhoYE = sampleYE[rhoSs];  Q1YE = sampleYE[Q1Ss];  Q2YE = sampleYE[Q2Ss];
     (* the comoving advection velocity, on the faces the noise already uses.
        NB the nodal Q gradients the product-rule A21/A31 needed are NOT sampled:
        the divergence form uses the face samples instead.  The superseded
        modules still take those four samples and never use them. *)
     rhoXd = sampleXE[Derivative[1, 0][rhoSs][##] &];
     rhoYd = sampleYE[Derivative[0, 1][rhoSs][##] &];]];
  Print["    sampling (12 batched calls): ", tSamp, " s"];
  Print["    \[Rho]_ss \[Element] ", MinMax@rhoV,
        "   Q1_ss \[Element] ", MinMax@Q1V,
        "   Q2_ss \[Element] ", MinMax@Q2V];

  wxE = ux + (Bv/xi0) rhoXd;
  wyE = uy + (Bv/xi0) rhoYd;
  peclet = N[Max[Max@Abs@wxE, Max@Abs@wyE] h/(2 kapQ)];
  Print["    w_ss = u + (B/\[Xi]0)\[Del]\[Rho]_ss:  w_x \[Element] ", MinMax@wxE,
        ",  w_y \[Element] ", MinMax@wyE];
  Print["    cell Peclet |w| h/(2 \[Kappa]) = ", ToString@NumberForm[peclet, 4],
        "   (\[Kappa] = 4K'/\[Xi]r = ", ToString@NumberForm[N@kapQ, 5], ")",
        If[peclet < 1, "   -- centred differencing is fine", "   *** TOO COARSE ***"]];
  If[peclet >= 1, Message[SolveFluctuationsLyapunov::peclet,
                          ToString@NumberForm[peclet, 4]]];

  (* --- [3b] what the pedestal is expected to come out at ----------------
     In the rho-only sector W commutes with the drift, so C = W C0 W with C0 the
     flat pedestal, and the finite-box value at the core follows in closed form
     from the very basis W was built in.  Worth reporting because for a SEALED
     box it is NOT the continuum eq:iso-gauss: the canonical projection removes
     the uniform mode, so a correlation length comparable with the box leaves
     genuinely less room for a density fluctuation.  Measured, the factor
     depends on lNoise/fbox ALONE -- identical at fbox = 10, 20 and 40 -- and is
     flat in h, so it is geometry and not discretisation:

         lNoise/fbox   0.05    0.1     0.2     0.3     0.4
         boxFactor     0.997   0.985   0.938   0.860   0.750

     The cutoff calibration neither should nor does remove it: it corrects the
     LATTICE, and this is the finite sealed box.  Under Dirichlet the factor is
     1 to six digits at every grid and every lNoise, so a boxFactor that is not
     ~1 there would mean something is wrong. *)
  If[OddQ[nInt] && lNoiseV > 0.,
     Module[{vc, dC, rhoC, sig0},
       vc = If[neuQ,
          Table[(If[p == 0, 1., 2.]/nInt) Cos[N[Pi] p (iZero - 1/2)/nInt]^2,
                {p, 0, nInt - 1}],
          Table[(2./(nInt + 1)) Sin[N[Pi] p iZero/(nInt + 1)]^2, {p, nInt}]];
       dC = Total[Wsym^2 vc];
       rhoC = rhoV[[iZero + (iZero - 1) nInt]];
       sig0 = zet rhoC/(Bv rho0v h^2);
       pedPred = sig0 (dC^2 - If[neuQ, 1./nInt^2, 0.]);
       pedCont = zet rhoC/(2 Pi Bv rho0v lNoiseV^2);
       boxFactor = pedPred/pedCont];
     Print["[3b] predicted \[Rho]-only pedestal at the core = ", pedPred,
           ",  eq:iso-gauss = ", pedCont,
           ",  ratio = ", ToString@NumberForm[boxFactor, 6]];
     If[neuQ && Abs[boxFactor - 1.] > 0.01,
        Print["     the sealed box suppresses it by ",
              ToString@NumberForm[N[100 (1 - boxFactor)], 3],
              "% at lNoise/fbox = ", ToString@NumberForm[lNoiseV/fboxV, 3],
              " -- finite-box physics, not a numerical error"]],
     pedPred = pedCont = boxFactor = Missing["needs odd Nint and lNoise > 0"]];

  Svals = Sqrt[Q1V^2 + Q2V^2]/rhoV;
  SXE = Sqrt[Q1XE^2 + Q2XE^2]/rhoXE;
  SYE = Sqrt[Q1YE^2 + Q2YE^2]/rhoYE;
  (* The closure is needed wherever a coefficient is needed, so three times over. *)
  {tClo, {g2V, g3V, g2XE, g3XE, g2YE, g3YE}} = AbsoluteTiming[
     {g2Exact /@ Svals, g3Exact /@ Svals,
      g2Exact /@ SXE,   g3Exact /@ SXE,
      g2Exact /@ SYE,   g3Exact /@ SYE}];
  rec["[3] kappa(S) inversion + g2,g3 (nodes + both face sets)", tClo];
  {tKchk, kapResid} = AbsoluteTiming[
     Max @ Abs[(BesselI[1, #]/BesselI[0, #] & /@ (kappaOfS /@ Svals)) - Svals]];
  rec["[3] kappa(S) verification (diagnostic)", tKchk];
  Print["    S = |Q_ss|/\[Rho]_ss \[Element] ", MinMax@Svals, ",  mean ", Mean@Svals];
  Print["    g2 \[Element] ", MinMax@g2V, "   g3 \[Element] ", MinMax@g3V,
        "   (small-S limits 1/2, 1/6)"];
  If[kapResid > 1.*^-8, Message[SolveFluctuationsLyapunov::kappa, kapResid]];

  (* --- drift matrix A --------------------------------------------------
     Q = {{Q1,Q2},{Q2,-Q1}}.  Sign a -> a + Lambda xir: rotational diffusion
     DESTROYS order.

     A11, A12, A13 are the rho row and stay put under both walls.  A12 and A13
     carry the flux of rho driven by grad Q, and THAT flux is blocked at a sealed
     wall whatever dQ does, so the Neumann divergence there is the correct
     statement and not a compromise -- it is what keeps mass conservation exact.

     A21, A31 keep the divergence form, self-adjoint as Div(c Grad .) is in the
     continuum, reusing the face samples the noise already needs.

     A22, A33 carry the comoving advection.  See "The advection term" in the
     header: Adv SUBSUMES the (B/xi0) diag[rhoLap] term that the superseded
     modules have, so that term must NOT also appear here. *)
  Print["[4] Assembling A ..."];
  diag[c_] := DiagonalMatrix[SparseArray @ c];
  tAsmA = rec["[4] assemble A", First @ AbsoluteTiming[
     Adv = If[advectQ,
              -(Transpose[Gx] . (wxE * Px) + Transpose[Gy] . (wyE * Py)),
              (* advection -> False reproduces the superseded modules *)
              (Bv/xi0) diag[rhoLap]];

     A11 = (Bv/xi0) Lap;
     A12 = (zet/xi0) (Lx - Ly);
     A13 = (2 zet/xi0) Dxy;

     A21 = -(Bv/xi0) (Transpose[Gx] . (Q1XE * Gx) + Transpose[Gy] . (Q1YE * Gy));
     A22 = Adv +
           (4/xir) (diag[-(av + Lam xir) - bv (6 Q1V^2 + 2 Q2V^2)] + cL Lapq);
     A23 = -(16 bv/xir) diag[Q1V Q2V];

     A31 = -(Bv/xi0) (Transpose[Gx] . (Q2XE * Gx) + Transpose[Gy] . (Q2YE * Gy));
     A32 = A23;
     A33 = Adv +
           (4/xir) (diag[-(av + Lam xir) - bv (2 Q1V^2 + 6 Q2V^2)] + cL Lapq);

     Atilde = ArrayFlatten[{{A11, A12, A13},
                            {A21, A22, A23},
                            {A31, A32, A33}}];]];
  Print["    A dims = ", Dimensions[Atilde], ",  assembled in ", tAsmA, " s"];

  (* --- [4b] the advection tripwire -------------------------------------
     Div(w_ss) = (B/xi0) Lap rho_ss because u is constant, so Adv applied to a
     uniform dQ must reproduce the term it replaced.  Checked on INTERIOR nodes
     only: the constant field is not representable under either wall condition,
     so the outermost two layers differ legitimately. *)
  If[advectQ && nInt >= 7,
     intMask = Flatten @ Table[
        If[2 < i < nInt - 1 && 2 < j < nInt - 1, 1, 0], {j, nInt}, {i, nInt}];
     advIdent = With[{d = Normal[Adv . ConstantArray[1., nInt^2]] - (Bv/xi0) rhoLap},
        Max @ Abs[Pick[d, intMask, 1]]/Max[1.*^-30, Max @ Abs @ Pick[
           N[(Bv/xi0) rhoLap], intMask, 1]]];
     (* CForm, not NumberForm: NumberForm renders a small Real as a
        SUPERSCRIPTED "m x 10^e", which OutputForm spreads over two lines and
        shreds the log -- the trap comoving_steady_solver.m:245-248 records. *)
     Print["[4b] advection identity  |Adv.1 - (B/\[Xi]0)\[Del]^2\[Rho]_ss|/scale = ",
           ToString[CForm[advIdent]], "   (interior nodes; ",
           "discretisation error only)"];
     (* Threshold 0.5, not a few percent.  Adv.1 is a finite-volume divergence
        of face-sampled w_ss; (B/xi0) rhoLap is a nodal Laplacian of the steady
        interpolant.  They agree only to O((h/ell_d)^2), which is 3.9% at
        Nint = 50, fbox = 15 and several times that on the coarsest grids of a
        box sweep -- all legitimate.  A real assembly error (a sign, a factor 2,
        a transposed operator) lands at ~1 or above, so 0.5 still catches one. *)
     If[advIdent > 0.5,
        Message[SolveFluctuationsLyapunov::advid, ToString[CForm[advIdent]]]],
     advIdent = Missing["not checked"]];

  (* --- [4a] the sealed-wall conservation gate on A --------------------- *)
  vCons = Join[ConstantArray[1., nInt^2], ConstantArray[0., 2 nInt^2]];
  If[neuQ,
     scaleA = Max @ Abs @ Atilde["NonzeroValues"];
     Print["[4a] Conservation of the uniform \[Delta]\[Rho] mode (scale ", scaleA, ") ..."];
     (* Adv lives in A22/A33 only, so neither of these sees it: the rho row block
        and the drho column block are untouched by the advection. *)
     consA  = consCheck["A . v   (right null vector)",
                        Max @ Abs @ Normal[Atilde . vCons], scaleA];
     consAT = consCheck["v^T . A (left null vector)",
                        Max @ Abs @ Normal[vCons . Atilde], scaleA],
     consA = consAT = Missing["Dirichlet"]];

  (* --- noise matrix Q --------------------------------------------------
     f = d_i F_i, <F_i F_j> = C_ij delta(r-r') -> sum_ij B_i diag(C_ij) B_j^T with
     B_i = -Gi^T (staggered, face-sampled) or -Di^T (collocated, node-sampled).
     Writing EVERY divergence as minus a transposed gradient is what makes
     1^T B_i = -(Gi.1)^T = 0, i.e. no noise reaches the conserved mode.
     EXACT max-entropy moments:
       <QQ> = rho J + (g2/2rho) U,  <Q'Q'> = 4 rho J - (2 g2/rho) U,
       <QQQ> = (1/2) T + (g3/rho^2) W;   g2 = g3 = 0 gives the old closure.
     Setting every coefficient to 1 makes the first two terms exactly -Lap, which
     is what forces R == 1.  The cross terms stay COLLOCATED -- see the header. *)
  Print["[5] Assembling Q ..."];
  kap = 2 zet/(xi0 rho0v);
  pf = kap/h^2;
  locPf = 2 Lam/(rho0v h^2);
  (* Dirichlet's Dx.(c * Dy^T) and Neumann's Dx^T.(c * Dy) agree identically when
     D1^T = -D1; the transposed-gradient form is the one that generalises. *)
  divBlk[cxxF_, cyyF_, cxy_, cyx_] :=
    Transpose[Gx] . (cxxF * Gx) + Transpose[Gy] . (cyyF * Gy) +
    Transpose[Dx] . (cxy * Dy) + Transpose[Dy] . (cyx * Dx);

  tAsmQ = rec["[5] assemble Q (24 sparse products)", First @ AbsoluteTiming[
     h2c = Q1V^2 - Q2V^2;        h2s = 2 Q1V Q2V;
     h3c = Q1V^3 - 3 Q1V Q2V^2;  h3s = 3 Q1V^2 Q2V - Q2V^3;
     e2 = g2V/(2 rhoV);          (* U terms *)
     e3 = g3V/(4 rhoV^2);        (* W terms *)
     h2cXE = Q1XE^2 - Q2XE^2;        h2sXE = 2 Q1XE Q2XE;
     h3cXE = Q1XE^3 - 3 Q1XE Q2XE^2; h3sXE = 3 Q1XE^2 Q2XE - Q2XE^3;
     e2XE = g2XE/(2 rhoXE);          e3XE = g3XE/(4 rhoXE^2);
     h2cYE = Q1YE^2 - Q2YE^2;        h2sYE = 2 Q1YE Q2YE;
     h3cYE = Q1YE^3 - 3 Q1YE Q2YE^2; h3sYE = 3 Q1YE^2 Q2YE - Q2YE^3;
     e2YE = g2YE/(2 rhoYE);          e3YE = g3YE/(4 rhoYE^2);

     (* rho-rho: no closure enters *)
     Qrr = pf * divBlk[Q1XE + rhoXE, -Q1YE + rhoYE, Q2V, Q2V];

     QQ1Q1 = (pf * divBlk[3 Q1XE/4 + rhoXE/2 + e2XE h2cXE + e3XE h3cXE,
                          -3 Q1YE/4 + rhoYE/2 + e2YE h2cYE - e3YE h3cYE,
                          Q2V/4 + e3 h3s,
                          Q2V/4 + e3 h3s] +
              locPf * diag[2 rhoV - 4 e2 h2c]);

     QQ2Q2 = (pf * divBlk[Q1XE/4 + rhoXE/2 - e2XE h2cXE - e3XE h3cXE,
                          -Q1YE/4 + rhoYE/2 - e2YE h2cYE + e3YE h3cYE,
                          3 Q2V/4 - e3 h3s,
                          3 Q2V/4 - e3 h3s] +
              locPf * diag[2 rhoV + 4 e2 h2c]);

     (* NB Q1-Q2 carries a LOCAL part: <Q'Q'> is no longer proportional to J. *)
     QQ1Q2 = (pf * divBlk[Q2XE/4 + e2XE h2sXE + e3XE h3sXE,
                          -Q2YE/4 + e2YE h2sYE - e3YE h3sYE,
                          Q1V/4 - e3 h3c,
                          Q1V/4 - e3 h3c] +
              locPf * diag[-4 e2 h2s]);
     QQ2Q1 = Transpose[QQ1Q2];

     QQ1r = pf * divBlk[Q1XE + rhoXE/2 + e2XE h2cXE,
                        Q1YE - rhoYE/2 - e2YE h2cYE,
                        e2 h2s,
                        e2 h2s];
     QQ2r = pf * divBlk[Q2XE + e2XE h2sXE,
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

  (* --- [5b] impose the cutoff: Q -> W . Q . W ---------------------------
     W = Kron[Emat, Emat] acting on the grid index of each of the 9 blocks.  Never
     form W: with k = i + (j-1) Nint an Nint^2 x Nint^2 block reshapes to the
     4-tensor X[[j,i,j',i']] and W.X.W is Emat contracted onto each of the four
     slots -- 36 Nint^5 flops in all, 0.04% of the 27 Nint^6 Eigensystem below.
     IN PLACE, ONE BLOCK PAIR AT A TIME; see the header for why. *)
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

  (* --- [5a] the sealed-wall conservation gate on Q --------------------- *)
  If[neuQ,
     scaleQ = If[Head[Qtilde] === SparseArray,
                 Max @ Abs @ Qtilde["NonzeroValues"], Max @ Abs @ Qtilde];
     Print["[5a] Noise into the conserved mode (scale ", scaleQ, ") ..."];
     consQ = consCheck["v^T . Q (no noise on the total)",
                       Max @ Abs @ Normal[vCons . Qtilde], scaleQ],
     consQ = Missing["Dirichlet"]];

  (* --- condition (ii): local PSD of the flux covariance ---------------- *)
  tCond = rec["[5] condition (ii): N^2 local 6x6 eigensolves", First @ AbsoluteTiming[
     minEigLoc = Min /@ Eigenvalues /@
        MapThread[localFluxCovariance, {rhoV, Q1V, Q2V, g2V, g3V}];]];
  scaleLoc = Max @ Abs @ Flatten @ {rhoV, Q1V, Q2V};
  Print["    local 6x6 min eig (exact closure)  \[Element] ", MinMax@minEigLoc,
        ";  negative at ", Count[minEigLoc, x_ /; x < -1.*^-10 scaleLoc],
        "/", Length@minEigLoc, " grid points   [", tCond, " s]"];
  If[Min@minEigLoc < -1.*^-10 scaleLoc,
     Message[SolveFluctuationsLyapunov::notpsd,
             Count[minEigLoc, x_ /; x < -1.*^-10 scaleLoc]]];

  (* The node-wise 6x6 above is still the RIGHT physical condition, but the mixed
     scheme is no longer a single Gram matrix D M D^T, so check the assembled
     Qtilde too.  Under Neumann Q is only PSD, never PD -- the conserved mode is
     an exact null direction -- so ask for one eigenvalue more and expect one at
     round-off.  Also note this runs AFTER [5b], i.e. it tests the SMOOTHED Q;
     a congruence preserves PSD, and W only shrinks the spectrum. *)
  If[3 nInt^2 <= 20000,
     {tQpsd, minEigQ} = AbsoluteTiming[
        Quiet @ Check[
           Min @ Re @ Eigenvalues[N @ Qtilde, If[neuQ, -5, -4], Method -> "Arnoldi"],
           Indeterminate]];
     rec["[5] condition (ii): min eig of assembled Q (Arnoldi)", tQpsd];
     Print["    assembled Q min eig = ", minEigQ, "   [", tQpsd, " s]",
           If[NumericQ[minEigQ] && minEigQ < -1.*^-8 Max@Abs@Diagonal@Qtilde,
              "   *** Q IS NOT PSD ***",
              If[neuQ, "   (PSD; one exact zero is the sealed mode)", "   (PSD)"]]],
     minEigQ = Missing["skipped, n > 20000"];
     Print["    assembled Q min eig: skipped (n = ", 3 nInt^2, " > 20000)"]];

  (* --- Lyapunov solve --------------------------------------------------
     Dispatched through $LyapunovSolveKernel (see the kernel contract at
     LyapunovReferenceKernel::usage) so a faster backend can be swapped in
     without duplicating blocks [1]-[5]. *)
  Print["[6] Lyapunov solve ..."];
  (* ToPackedArray is load-bearing.  G1, L1, D1 are built from EXACT RATIONALS so
     the adjoint-pair identity can be tested with "== 0"; if any of that
     exactness survives into the dense arrays they come back unpacked, every BLAS
     path is lost, and the back-substitute runs ~600x slower. *)
  (* An out-of-core backend streams A and Q to disk itself, so materialising
     dense copies here would cost 2 x 8 n^2 bytes for nothing -- and Qtilde is
     ALREADY dense (densified in [5b] for the cutoff congruence), so Qdense was
     a second full copy of it.  Hand the assembled matrices straight through. *)
  If[TrueQ[$LyapunovKernelWantsDense],
     tDense = rec["[6] densify A, Q (Normal[])", First @ AbsoluteTiming[
        Adense = Developer`ToPackedArray @ N @ Normal @ Atilde;
        Qdense = Developer`ToPackedArray @ N @ Normal @ Qtilde;]];
     Print["    [0] densify A,Q        ", tDense, " s   (",
           ToString @ NumberForm[N[2 * 8 * Length[Adense]^2/2^30], 3], " GiB)",
           "   packed: A ", Developer`PackedArrayQ[Adense],
           ", Q ", Developer`PackedArrayQ[Qdense]],
  (* else: pass through, no copy *)
     Adense = Atilde; Qdense = Qtilde;
     Print["    [0] densify A,Q        skipped (kernel streams them itself)"]];

  (* "release" lets an out-of-core kernel drop the inputs the moment they are on
     disk, BEFORE it allocates C.  Without it A, Q and C are all live at the read
     and the peak is 3 x 8 n^2 instead of 1.  It closes over this Module's locals,
     which is the only way to clear a caller's binding from inside the callee;
     calling it is optional and the reference kernel never does. *)
  kernelOut = $LyapunovSolveKernel[Adense, Qdense,
     <| "neumann" -> neuQ, "bcType" -> bcType, "nInt" -> nInt, "n" -> 3 nInt^2,
        "vCons" -> vCons, "scratchDir" -> scratchV,
        "release" -> Function[Null,
           Adense =.; Qdense =.; Atilde =.; Qtilde =.;
           ClearSystemCache[];] |>];

  (* A backend that returns anything else would otherwise surface far downstream
     as a cascade of Missing[] inside the export.  Fail here, where the cause is
     still legible. *)
  kernelKeys = {"C", "eigenvalues", "maxReLambda", "residual", "residualRefined",
                "nullIndex", "nullLambdaRel", "nullLambda2Rel", "nullOverlap",
                "timings", "backend"};
  If[!AssociationQ[kernelOut] || !AllTrue[kernelKeys, KeyExistsQ[kernelOut, #] &],
     Message[SolveFluctuationsLyapunov::badkernel,
             If[AssociationQ[kernelOut], Keys[kernelOut], Head[kernelOut]],
             kernelKeys];
     Abort[]];

  Cmat        = kernelOut["C"];
  eigvals     = kernelOut["eigenvalues"];
  maxReLam    = kernelOut["maxReLambda"];
  resid       = kernelOut["residual"];
  residRef    = kernelOut["residualRefined"];
  iNull       = kernelOut["nullIndex"];
  nullRel     = kernelOut["nullLambdaRel"];
  lam2Rel     = kernelOut["nullLambda2Rel"];
  nullOverlap = kernelOut["nullOverlap"];
  timings     = Join[timings, kernelOut["timings"]];
  Print["    backend = ", kernelOut["backend"]];

  (* --- extract --------------------------------------------------------- *)
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
  If[neuQ,
     (* the canonical rank-one subtraction: the total is fixed, so the sum of the
        covariance over BOTH indices is exactly zero *)
     Print["    canonical check: sum of the rho-rho block = ",
           Total @ Flatten @ Cmat[[1 ;; nInt^2, 1 ;; nInt^2]],
           "   (exactly 0 for a sealed box)"]];

  (* runTag needs only these six fields, so it can be built before the full
     result exists.  Returned below so drivers never rebuild it. *)
  tagStr = runTag[<|"fbox" -> fboxV, "Nint" -> nInt, "bcType" -> bcType,
                   "lNoise" -> lNoiseV, "zeta" -> zet, "elld" -> elldV,
                   "B" -> Bv,
                   "ssBox" -> If[NumericQ[ssBox], ssBox, fboxV],
                   (* same test as the ::unstable message above *)
                   "unstable" -> TrueQ[NumericQ[maxReLam] && maxReLam >= 0]|>];
  Print["    runTag = ", tagStr];

  printTimings[timings, nInt];

  <| "C" -> Cmat,
     "runTag" -> tagStr,
     "zeta" -> zet,
     "elld" -> elldV,
     "derived" -> derived,
     "sigmaRho" -> sigRhoMat,
     "sigmaRhoFn" -> sigRhoFn,
     "sigmaQ1" -> sigQ1,
     "sigmaQ2" -> sigQ2,
     "gridInt" -> gridInt,
     "h" -> N[h],
     "Nint" -> nInt,
     "fbox" -> fboxV,
     "bcType" -> bcType,
     "bc" -> If[neuQ,
        "sealed + anchoring wall: \[Delta]\[Rho] Neumann, noise flux blocked (all \
sectors), \[Delta]Q strongly anchored (Dirichlet), canonical (uniform \
\[Delta]\[Rho] mode projected out)",
        "reservoir wall: \[Rho]_ss pinned, \[Psi]|_bdry = 0, noise flux free"],
     "u" -> uVec[[1 ;; 2]],
     "lNoise" -> lNoiseV,
     "lNoiseEff" -> lEff,
     "lNoiseOverH" -> If[lNoiseV == 0., 0., lNoiseV/N[h]],
     "lNoiseEffOverH" -> If[lEff == 0., 0., lEff/N[h]],
     "rescaleCutoff" -> rescaleQ,
     "advection" -> advectQ,
     "peclet" -> peclet,
     "advIdentity" -> advIdent,
     (* the rho-only pedestal this grid/box/cutoff should produce at the core,
        the continuum eq:iso-gauss value, and their ratio.  Under Neumann the
        ratio is the finite sealed-box suppression and depends on lNoise/fbox
        alone; under Dirichlet it should be 1 to six digits. *)
     "pedestalPredicted" -> pedPred,
     "pedestalContinuum" -> pedCont,
     "boxFactor" -> boxFactor,
     "Wsymbol" -> Wsym,
     "QasymBeforeSym" -> Qasym,
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
     "scheme" -> ("comoving drift Div(dQ w_ss), w_ss = u + (B/xi0) Grad rho_ss, in \
conservative form -Gi^T diag(w|face) Pi (dean.tex eq. linearized-Q-equation); R1 \
consistent stencils (adjoint-pair diagonal noise flux, collocated cross terms, \
divergence-form A21/A31); every divergence a transposed gradient; " <>
        If[neuQ, "cell-centred sealed wall with anchored L1q in A22/A33, \
symmetrised double-divergence d_x d_y, canonical projection of the conserved \
uniform mode; reflecting (DCT-II) ", "node-centred reservoir wall; Dirichlet \
(sine) "] <> "noise kernel Q -> W.Q.W with W = MatrixExp[(lEff^2/4) Lap]" <>
        If[rescaleQ, "; lEff CALIBRATED from the physical lNoise against \
eq:lattice-bz so the isotropic pedestal matches eq:iso-gauss (pedestal channel \
only)", "; lEff = lNoise, uncalibrated"]),
     "modelParams" -> modelParams,
     "lyapunovParams" -> lyapunovParams,
     "timings" -> timings |>
];


(* ::Subsection:: *)
(*Visualization and export*)

visualizeFluctuationsComoving[correlations_Association, outputDir_,
                              plotBox_?NumericQ] :=
Module[{dir, gridInt, sigRhoFn, xLo, xHi, plotPts, lbl, plot2D, plot3D, tag, x, y},
  dir = prepareDir[outputDir];
  gridInt = correlations["gridInt"];
  sigRhoFn = correlations["sigmaRhoFn"];
  tag = Lookup[correlations, "runTag", runTag[correlations]];
  (* plotBox crops the wall boundary layer out of the figures; the interpolant
     only spans the interior nodes, so clamp to that. *)
  {xLo, xHi} = {Max[-plotBox, First@gridInt], Min[plotBox, Last@gridInt]};
  plotPts = Min[correlations["Nint"], 200];
  lbl = ("Var(\[Delta]\[Rho])(r)  [" <> correlations["bcType"] <> ", h=" <>
         ToString @ NumberForm[correlations["h"], 3] <>
         ", l_noise=" <> ToString @ NumberForm[correlations["lNoise"], 4] <>
         " (eff " <> ToString @ NumberForm[correlations["lNoiseEff"], 4] <> ")]");

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

  savePlot[dir, "sigma_rho_" <> tag <> ".png", plot2D];
  savePlot[dir, "sigma_rho_surface_" <> tag <> ".png", plot3D];
  {plot2D, plot3D}];

exportDataComoving[correlations_Association, outputDir_,
                   saveCovariance_ : False] :=
Module[{dir, tag, out, dumpPath, binPath, metaPath, Cmat, n, lyapTime, ssTime,
        prov, paths = {}},
  dir = prepareDir[outputDir];
  If[dir === None, Return[{}]];
  tag = Lookup[correlations, "runTag", runTag[correlations]];
  out[stem_, ext_] := FileNameJoin[{dir, stem <> tag <> ext}];

  (* EVERY parameter the run used, written into BOTH files so either one alone
     reconstructs the run.  solverParams (the COMOVING steady solver's) and
     otherParams are not arguments of this module -- the driver attaches them to
     `correlations` the same way it attaches steadyStateTime -- so they are read
     with Lookup and come back None if it did not.  The scalars duplicated out
     of modelParams are there so a reader does not have to unpack a rule list to
     filter a sweep. *)
  prov = <|
     "runTag" -> tag,
     "bcType" -> correlations["bcType"], "bc" -> correlations["bc"],
     "u" -> correlations["u"],
     "zeta" -> correlations["zeta"], "elld" -> correlations["elld"],
     (* BOTH cutoff lengths: lNoise is the physical one the run claims,
        lNoiseEff the one W was built from, and their ratio is the lattice
        correction that was removed. *)
     "lNoise" -> correlations["lNoise"],
     "lNoiseEff" -> correlations["lNoiseEff"],
     "lNoiseOverH" -> correlations["lNoiseOverH"],
     "lNoiseEffOverH" -> correlations["lNoiseEffOverH"],
     "rescaleCutoff" -> correlations["rescaleCutoff"],
     "advection" -> correlations["advection"],
     "peclet" -> correlations["peclet"],
     "advIdentity" -> correlations["advIdentity"],
     "boxFactor" -> correlations["boxFactor"],
     "pedestalPredicted" -> correlations["pedestalPredicted"],
     "pedestalContinuum" -> correlations["pedestalContinuum"],
     "scheme" -> correlations["scheme"],
     "modelParams" -> correlations["modelParams"],
     "lyapunovParams" -> correlations["lyapunovParams"],
     "derived" -> Lookup[correlations, "derived", None],
     "solverParams" -> Lookup[correlations, "solverParams", None],
     "otherParams" -> Lookup[correlations, "otherParams", None]|>;

  dumpPath = out["sigma_rho_", ".m"];
  Export[dumpPath,
     Join[prov,
       <|"fbox" -> correlations["fbox"], "Nint" -> correlations["Nint"],
         "h" -> correlations["h"], "gridInt" -> correlations["gridInt"],
         "sigmaRho" -> correlations["sigmaRho"],
         "sigmaQ1" -> correlations["sigmaQ1"],
         "sigmaQ2" -> correlations["sigmaQ2"],
         "maxReLambda" -> correlations["maxReLambda"],
         "residualRefined" -> correlations["residualRefined"],
         "minEigQ" -> correlations["minEigQ"],
         "lapIdentity" -> correlations["lapIdentity"]|>]];
  Print["Saved: ", dumpPath];
  AppendTo[paths, dumpPath];

  (* Full covariance, row-major Real64, n = 3 Nint^2 with blocks (drho, dQ1, dQ2).
     Read back with ArrayReshape[BinaryReadList[f, "Real64"], {n, n}]. *)
  If[TrueQ[saveCovariance],
     Cmat = correlations["C"];
     n = Length[Cmat];
     binPath = out["covariance_", ".bin"];
     metaPath = out["covariance_", "_meta.m"];
     Export[binPath, Cmat, {"Binary", "Real64"}];
     lyapTime = Total[correlations["timings"][[All, 2]]];
     ssTime = Lookup[correlations, "steadyStateTime", None];
     (* n, Nint, fbox, h and gridInt keep their names and meaning, so
        binaryReadTools.m's readCovRho/readCovQ/readCovRhoQ/readCorrQ read these
        files unchanged.  Everything else is additive. *)
     Export[metaPath,
        Join[prov,
          <|"n" -> n, "Nint" -> correlations["Nint"], "fbox" -> correlations["fbox"],
            "h" -> correlations["h"], "gridInt" -> correlations["gridInt"],
            "blocks" -> {"rho", "Q1", "Q2"},
            "ordering" -> "row-major Real64; grid index k = i + (j-1) Nint, x fastest",
            "lyapunovTime" -> lyapTime,
            "steadyStateTime" -> ssTime,
            "totalTime" -> If[NumericQ[ssTime], ssTime + lyapTime, lyapTime],
            "timings" -> correlations["timings"]|>]];
     Print["Saved: ", binPath, "  (",
           ToString@NumberForm[N[FileByteCount[binPath]/2^20], {6, 1}], " MiB)"];
     Print["Saved: ", metaPath];
     paths = Join[paths, {binPath, metaPath}]];
  paths];

End[];
