#!/usr/bin/env wolframscript
(* ::Package:: *)

(* set2_2_comove -- Gaussian fluctuations about the COMOVING +1/2 defect steady
   state, reservoir (Dirichlet) wall on both solvers.

   ----------------------------------------------------------------------------
   WHAT CHANGED, AND WHY THE OLD DATA IN THIS DIRECTORY IS SUPERSEDED
   ----------------------------------------------------------------------------
   This driver used to call lyapunov_solver_module_cutoff.m, whose drift
   generator DROPS the comoving advection entirely -- its own header recorded
   that as a known defect ("the correct generator carries an extra -u.Del").  It
   now calls lyapunov_solver/lyapunov_solver.m, where the term is present:

       d_t dQ = (B/xi0) Div(Q_ss Grad drho) + Div(dQ w_ss)
                + (4/xir) dH_active + noise,   w_ss = u + (B/xi0) Grad rho_ss

   (dean.tex eq. linearized-Q-equation).  The old cancellation assumed
   u = -(B/xi0) Grad rho POINTWISE, which is the superseded derivation in
   FRIDGE/wrong_comoving.tex:45; with u the constant Lagrange multiplier that
   comoving_defect/ actually solves for, what survives is the (B/xi0) Grad rho_ss
   advection, and at the core that is the LARGER piece -- 1.07 um/min against
   |u| = 0.28.  Numbers produced before this change are wrong in the dQ sector.

   Three other changes:

     zeta        10 -> 8.  This moves ONLY K, 0.05368 -> 0.06368: a = -0.172,
                 xi_r = 0.1, a' = -0.162, K' = 0.10368, ell_d = 0.8 and S0 = 0.9
                 are all unchanged, because K' = K + zeta xi_r/(4 xi_0) absorbs
                 it exactly.  What does move is the activity
                 alpha = zeta xi_r S0/(4 xi_0 K'), 0.434 -> 0.347, and the noise
                 amplitude 2 zeta/(xi_0 rho_0), so the pedestal drops 20%.
     cutoff      lNoise is now CALIBRATED against dean.tex eq:lattice-bz by
                 default: the module solves F(l/h) = 1/(Sqrt[2 Pi] lNoise/h) for
                 the length W is built from, so the lattice pedestal matches the
                 continuum eq:iso-gauss at the lNoise asked for.  The old
                 "lNoise >= 5 h" rule this file used to quote is superseded --
                 that bias is what the calibration removes.  Both lengths are
                 reported and exported; files are tagged by the PHYSICAL one.
     tag         now carries Nint, box, lNoise, zeta and ell_d, so runs that
                 differ in activity or core size can share a directory.  It
                 comes back as correlations["runTag"] and is NOT rebuilt here --
                 the previous version duplicated the module's private tag with a
                 comment reading "KEEP THE TWO IN STEP", which is exactly the
                 kind of thing that drifts and silently unpairs steady_* from
                 covariance_*.

   Output files therefore look like

       sigma_rho_box15_N101_Dir_lN1.73_z8_ld0.8.m
       steady_box15_N101_Dir_lN1.73_z8_ld0.8.m
       covariance_box15_N101_Dir_lN1.73_z8_ld0.8.bin   (+ _meta.m)

   which cannot collide with the old _lN1.7300 tags.  The same directory will
   hold both families; filter on the tag, not on a glob.

   ----------------------------------------------------------------------------
   TWO THINGS A READER OF THIS DATA STILL HAS TO KNOW
   ----------------------------------------------------------------------------
   (1) EXPECT ::unstable ON THE COARSE GRIDS.  maxReLambda comes out > 0 and the
       module warns and keeps going.  It is the defect core being under-resolved
       and it refines away fast -- measured at fbox = 10, ell_d = 0.8:

           Nint             15      21      27      33
           ell_d/h        0.64    0.88    1.12    1.36
           maxReLambda   0.498  0.0964  0.0342  0.0117

       falling faster than h^3 while Sigma_rho converges.  It is NOT the
       advection: the same sweep with advection -> False gives 0.423, 0.0789,
       0.0275, 0.0056 -- same order, same trend.  Read the trend in Nint, never
       a single grid.  In THIS sweep the worst corner is box = 40, Nint = 50,
       where ell_d/h = 0.51; it is a convergence-curve point, not a
       quantitative one.

       THE CONSEQUENCE, which is easy to mistake for a bug: where
       maxReLambda > 0 the Lyapunov solution is not a covariance, and the dQ
       variances come out NEGATIVE.  Measured at box = 15, Nint = 50
       (maxReLambda = 0.028): Sigma_rho at the core is a healthy 1.36e-4 while
       Sigma_Q1 is -3.26.  Sigma_rho is the robust channel and stays sensible
       well below the grid where the dQ sector does; do not read Sigma_Q1/Q2
       from a task whose maxReLambda is positive.

   (2) THE PEDESTAL CALIBRATION IS THE rho-rho CHANNEL ONLY.  It is one scalar,
       fitted to the isotropic pedestal, so it cannot also be exact for the
       anisotropic term or the dQ sector.  Under the Dirichlet wall used here
       correlations["boxFactor"] comes out 1 to six digits, which is the check
       that it did what it claims.

   See lyapunov_solver/README.md for the full account of both.               *)

thisDir = DirectoryName[$InputFileName];
repoRoot = ParentDirectory[ParentDirectory[thisDir]];

(*
visualizeSteadyState[steadyState, outputDir, plotBox, tag]
   The ONLY reason this file is loaded.  lyapunov_solver.m is self-contained --
   it carries its own prepareDir/savePlot and opens its own private context --
   so, unlike the superseded modules, it does not borrow anything from here.
   Loaded FIRST so that lyapunov_solver.m's definitions of the two shared
   helpers (MergeParams, ScriptParamOverrides) are the ones left standing; they
   are identical either way.
*)
Import[FileNameJoin[{repoRoot, "iterative_solver_module.m"}]]

(*
SolveFluctuationsLyapunov[lyapunovParams, modelParams, steadyState, bcType]
visualizeFluctuationsComoving[correlations, outputDir, plotBox]
exportDataComoving[correlations, outputDir, saveCovariance]
MergeParams, ScriptParamOverrides
   NB do NOT also load any lyapunov_solver_module*.m: the old family shares one
   private context and SolveFluctuationsLyapunov is a name the base module also
   defines, so the last one loaded would win.
*)
Import[FileNameJoin[{repoRoot, "lyapunov_solver", "lyapunov_solver.m"}]]

(* The Python solve backend.  MUST come after the base module, which resets
   $LyapunovSolveKernel to the reference implementation on every load.  Dirichlet
   is accelerated; Neumann would delegate back to the reference kernel. *)
Import[FileNameJoin[{repoRoot, "lyapunov_solver", "lyapunov_solver_optimized.m"}]]

(*
SolveActiveNematicComovingSteady[modelParams, solverParams, bcType]
   Self-Gets comoving_defect/comoving_core.m.  Loaded LAST: it opens no private
   context, so it defines many bare names in Global`, but Block localises all of
   them and it leaves nothing carrying a value, so it cannot disturb the two
   modules above.
*)
Import[FileNameJoin[{repoRoot, "comoving_defect", "comoving_steady_solver.m"}]]

(* --- Command-line overrides ------------------------------------------
   wolframscript -file fluctuations.wls \
     "solverParams={box->25.}; lyapunovSolverParams={Nint->101, fbox->25}"
   Each named list is UPDATED key-by-key: keys not mentioned keep their
   default.  An unknown set name or key aborts. *)

cliParams = ScriptParamOverrides[];
knownParamSets = {"modelParams", "solverParams",
                  "lyapunovSolverParams", "otherParams"};
If[Complement[Keys[cliParams], knownParamSets] =!= {},
   Print["ERROR: unknown parameter set(s) ",
         Complement[Keys[cliParams], knownParamSets],
         "; known: ", knownParamSets];
   Abort[]];
upd[name_String, defaults_] :=
  MergeParams[defaults, Lookup[cliParams, name, {}], name];

(* --- Model & solver parameters (defaults) ---------------------------
   The dean.tex parameter table (dean.tex:605-673) at zeta = 8.  a, b and L
   follow from ell_d, S0 and gamma through paramRelations, so changing zeta
   here re-derives L automatically -- which is the point: zeta and L are not
   independent once ell_d is fixed. *)
paramRelations={ell==Sqrt[(L+(\[Zeta] \[Xi]r)/(4 \[Xi]0))/(-a-\[CapitalLambda] \[Xi]r)],S0==Sqrt[(-a-\[CapitalLambda] \[Xi]r)/(2b)],\[Gamma]==8 b / \[Xi]r};
knownParams={ell->1.0,S0->0.9,\[Xi]0->5.,\[Zeta]->8.,\[Rho]0->0.3};

guessedParams={B->10000.,b->0.1,\[Gamma]->8.,\[CapitalLambda]->0.1};

inferedParams=First@Solve[paramRelations/.Join[knownParams,guessedParams],{\[Xi]r,a,L}];

(* At zeta = 8 and ell_d = 1.0 this gives
     xi_r = 0.1,  a = -0.172,  L (= K) = 0.122
     a' = a + Lam xi_r = -0.162,  K' = L + zeta xi_r/(4 xi_0) = 0.162
     -> ell_d = Sqrt[K'/(-a')] = 1.0,  S0 = Sqrt[-a'/(2b)] = 0.9,
        alpha = zeta xi_r S0/(4 xi_0 K') = 0.2222,
        kappa_Q = 4 K'/xi_r = 6.48.
   Raising ell_d from 0.8 moves ONLY K (0.06368 -> 0.122) and everything that
   follows from K': xi_r, a, a' and S0 are untouched, because they are fixed by
   S0 and gamma alone.  Two consequences worth knowing: the defect core is
   better resolved at the same h (ell_d/h up 25%), and the cell Peclet number
   FALLS, since kappa_Q rises from 4.1472 to 6.48.  The solver recomputes and
   exports all of it under "derived", so nothing here is taken on trust. *)

(* lNoise is NOT derived from the others -- it is the physical UV cutoff, and
   with rescaleCutoff the solver turns it into the lNoiseEff that W is built
   from.  One value only in this set. *)
derivedNames = {B,a,b,L,\[Zeta],\[Xi]0,\[Xi]r,\[Rho]0,\[CapitalLambda]};

modelParams = upd["modelParams",
  Append[
    Thread[derivedNames->(derivedNames/.Join[knownParams,guessedParams,inferedParams])],
    lNoise -> 0.5]];

(* The five required comoving keys plus its optional ones, all listed so the
   shell can override any of them: MergeParams aborts on a key absent from the
   defaults.  box is swept by the jobscripts and MUST equal fbox below. *)
solverParams = upd["solverParams",
  {\[Delta]mesh -> 0.04, hMax -> 2., box -> 15.,
   maxSteps -> 25, relativeDiff -> 1.*^-8,
   refine -> 1, rhoGauge -> "Point", zetaSteps -> 5, order -> 4,
   stepTolerance -> 1.*^-10, verbose -> True}];

(* fbox MUST be an Integer -- the module tests IntegerQ -- whereas solverParams
   box is a Real that only ever reaches N[box].  They must be the SAME LENGTH:
   the Lyapunov module samples the steady fields at its own nodes AND cell
   edges, out to +-(fbox - h/2), so box < fbox would extrapolate the steady
   interpolant. *)
(* scratchDir must be listed here even though only the Python backend reads it:
   MergeParams aborts on a key absent from the defaults, so a jobscript passing
   it would otherwise kill the run before any work happened.  Automatic ->
   $TMPDIR, then /scratch/$SLURM_JOB_ID, then /scratch. *)
lyapunovSolverParams = upd["lyapunovSolverParams",
  {Nint -> 21, fbox -> 15, scratchDir -> Automatic}];

otherParams = upd["otherParams",
  {outputDir -> "./output", saveData -> True, plotBox -> 12}];

Print["modelParams          = ", modelParams];
Print["solverParams         = ", solverParams];
Print["lyapunovSolverParams = ", lyapunovSolverParams];
Print["otherParams          = ", otherParams];

outDir  = outputDir /. otherParams;   (* None -> nothing is written *)
plotBoxV = plotBox  /. otherParams;
saveDataV = saveData /. otherParams;

boxV  = box  /. solverParams;
fboxV = fbox /. lyapunovSolverParams;
(* fbox <= box, NOT equality.  This set deliberately also runs the steady solver
   on a box 5 units larger than the Lyapunov box, to separate the effect of the
   steady state's own wall from the fluctuation wall.  A LARGER steady box is
   always safe -- the Lyapunov module samples the steady interpolant out to
   +-(fbox - h/2), so it stays inside.  A SMALLER one extrapolates and is still
   fatal. *)
If[N[fboxV] > N[boxV],
   Print["FATAL: lyapunovSolverParams fbox = ", fboxV,
         " exceeds solverParams box = ", boxV,
         ".  The steady interpolant would be extrapolated at the wall."];
   Abort[]];
If[N[boxV] != N[fboxV],
   Print["NOTE: steady box ", boxV, " > Lyapunov fbox ", fboxV,
         " -- the steady state's wall sits ", N[boxV - fboxV],
         " units outside the fluctuation domain (recorded in the tag as _ss",
         boxV, ")."]];

(* --- Steady state ---------------------------------------------------- *)

Print["[1] Comoving steady state, Dirichlet (reservoir) wall ..."];
{tSS, steadyState} = AbsoluteTiming[
  SolveActiveNematicComovingSteady[modelParams, solverParams, "Dirichlet"]];
If[steadyState === $Failed,
   Print["FATAL: SolveActiveNematicComovingSteady returned $Failed"]; Abort[]];
Print["    done in ", tSS, " s"];
Print["    u = ", steadyState["u"],
      "   core = ", steadyState["solver"]["core"],
      "   n = ", steadyState["mesh"]["n"],
      " (", steadyState["mesh"]["nx"], "^2)"];
If[!TrueQ[steadyState["solver"]["converged"]],
   Print["*** WARNING: the steady Newton did NOT converge.  residual = ",
         steadyState["solver"]["residual"],
         " -- the fluctuations below are linearised about an unconverged "  <>
         "profile. ***"]];

(* --- Fluctuations ---------------------------------------------------- *)

correlations = SolveFluctuationsLyapunov[lyapunovSolverParams, modelParams,
                                         steadyState, "Dirichlet"];
If[correlations === $Failed,
   Print["FATAL: SolveFluctuationsLyapunov returned $Failed"]; Abort[]];

(* Everything the module cannot see, attached here so exportDataComoving writes
   it into BOTH sigma_rho_*.m and covariance_*_meta.m.  It reads all four with
   Lookup, so this is additive and order-independent. *)
correlations = Join[correlations,
  <|"solverParams" -> solverParams,
    "otherParams" -> otherParams,
    "steadyStateTime" -> tSS|>];

runTagV = correlations["runTag"];
Print["[*] runTag = ", runTagV];

(* --- Steady-state dump ------------------------------------------------
   The module's exportDataComoving writes the fluctuation results; the steady
   state they were linearised about is reconstructed here, in the driver, so
   that no .m module has to know about it.  It lands under the SAME tag, taken
   from correlations["runTag"] rather than rebuilt, so one run's four outputs
   sort together and cannot drift apart. *)

(* Batched, not scalar: an InterpolatingFunction called on two equal-length
   packed vectors returns a packed vector, which is the fast path the Lyapunov
   module itself uses.  Result[[i,j]] = f(gx[[i]], gy[[j]]). *)
sampleOn[f_, gx_, gy_] := Transpose @
  Table[f[gx, ConstantArray[N[y], Length[gx]]], {y, gy}];

(* The module's own prepareDir is PRIVATE to LyapunovSolver`Private`, so calling
   it from this file would silently resolve to an undefined Global`prepareDir --
   FileNameJoin then stays unevaluated, Export writes nothing, and the run still
   exits 0.  (Measured: exactly that, losing the steady dump.)  Do the two lines
   here instead. *)
ensureDir[outputDir_] := Which[
  outputDir === None, None,
  True, Quiet @ CreateDirectory[outputDir, CreateIntermediateDirectories -> True];
        If[!DirectoryQ[outputDir],
           Print["FATAL: cannot create outputDir ", outputDir]; Abort[]];
        ExpandFileName[outputDir]];

exportSteadyState[ss_Association, corr_Association, outputDir_, tag_String] :=
Module[{dir, path, m, gx, gy, gi, fields},
  dir = ensureDir[outputDir];
  If[dir === None, Return[{}]];
  path = FileNameJoin[{dir, "steady_" <> tag <> ".m"}];
  m = ss["mesh"];
  {gx, gy} = N /@ {m["gx"], m["gy"]};
  gi = N @ corr["gridInt"];
  fields = ss /@ {"rho", "Q1", "Q2"};
  If[!StringQ[path],
     Print["FATAL: steady-state path did not resolve to a string: ", path];
     Abort[]];
  Export[path,
    <|"runTag" -> tag, "frame" -> "comoving", "bc" -> "Dirichlet",
      (* the point of this data set *)
      "u" -> ss["u"],
      (* every parameter the run used, so this file alone reconstructs it *)
      "modelParams" -> modelParams, "solverParams" -> solverParams,
      "lyapunovParams" -> lyapunovSolverParams, "otherParams" -> otherParams,
      "derived" -> corr["derived"],
      "zeta" -> corr["zeta"], "elld" -> corr["elld"],
      "lNoise" -> corr["lNoise"], "lNoiseEff" -> corr["lNoiseEff"],
      (* the full steady state on ITS OWN Pade-graded tensor grid: everything
         needed to re-interpolate it offline, independently of Nint *)
      "gx" -> gx, "gy" -> gy, "nx" -> m["nx"], "ny" -> m["ny"],
      "hMin" -> m["hMin"], "hMax" -> m["hMax"], "box" -> m["box"],
      "rhoGrid" -> sampleOn[fields[[1]], gx, gy],
      "Q1Grid"  -> sampleOn[fields[[2]], gx, gy],
      "Q2Grid"  -> sampleOn[fields[[3]], gx, gy],
      "gridNote" -> "rhoGrid[[i,j]] = rho(gx[[i]], gy[[j]]); non-uniform, graded " <>
                    "as 1/|S'(r)| by cdGridPade",
      (* and exactly what the Lyapunov module saw, on ITS uniform grid *)
      "gridInt" -> corr["gridInt"],
      "rhoInt" -> sampleOn[fields[[1]], gi, gi],
      "Q1Int"  -> sampleOn[fields[[2]], gi, gi],
      "Q2Int"  -> sampleOn[fields[[3]], gi, gi],
      (* "state" is dropped from ss["solver"]: the same numbers as *Grid, in a
         packing the reader would have to decode *)
      "solver" -> KeyDrop[ss["solver"], "state"],
      "steadyStateTime" -> Lookup[corr, "steadyStateTime", None]|>];
  Print["Saved: ", path];
  {path}];

(* --- Write ------------------------------------------------------------
   DATA FIRST, PLOTS LAST, and not as a matter of taste.  Rasterising a plot
   needs the front end, which dies with "Illegal instruction (core dumped)" on a
   headless node under some Mathematica builds -- observed here under 13.3.1
   even with DISPLAY unset, killing the kernel outright so that nothing after it
   runs.  With the old ordering that crash destroyed the entire run: at
   Nint = 151 that is ~19 h of eigensolve thrown away for a PNG.  The jobscripts
   load mathematica/14, where it does not happen, but the ordering is the thing
   that makes it survivable rather than the module version.

   Trailing ";" is load-bearing: a bare plot would make wolframscript print the
   whole Graphics expression -- every mesh vertex -- to stdout. *)
written = Join[
   exportDataComoving[correlations, outDir, saveDataV],
   exportSteadyState[steadyState, correlations, outDir, runTagV]];
If[outDir =!= None,
   If[!AllTrue[written, StringQ[#] && FileExistsQ[#] &],
      Print["FATAL: expected output files are missing: ",
            Select[written, !(StringQ[#] && FileExistsQ[#]) &]];
      Abort[]];
   Print["[*] wrote ", Length[written], " files:"];
   Do[Print["      ", FileNameTake[f], "  (",
            ToString@NumberForm[N[FileByteCount[f]/2^20], {8, 2}], " MiB)"],
      {f, written}]];

(* Plotting is best-effort and deliberately LAST, because the data above is the
   expensive product.  Quiet, because a kernel without a usable front end
   answers Rasterize with Rasterize::type, and that message ECHOES THE WHOLE
   Graphics EXPRESSION: measured 2.9 MB of log for one Nint = 50 surface plot,
   which would be gigabytes at Nint = 151.  Under other builds the front end
   simply dies ("Illegal instruction"), taking the kernel with it -- which is
   the other reason nothing that matters may run after this point. *)
Quiet @ Check[
   visualizeFluctuationsComoving[correlations, outDir, plotBoxV];
   visualizeSteadyState[steadyState, outDir, plotBoxV, "_" <> runTagV];,
   Print["WARNING: plotting failed (no usable front end?).  The data files ",
         "listed above are unaffected."]];
