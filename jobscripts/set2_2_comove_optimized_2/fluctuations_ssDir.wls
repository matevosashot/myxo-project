#!/usr/bin/env wolframscript
(* ::Package:: *)

(* set2_2_comove_optimized_2, ssDir variant -- Gaussian fluctuations about the
   COMOVING +1/2 defect steady state, with the parameter set of params.md.

   Identical to fluctuations.wls in this directory except for ONE line: the
   steady state uses the Dirichlet density wall (rho = 1) instead of RhoNeumann.
   The runTag cannot tell the two apart, so this variant MUST write to its own
   outputDir.

   Same pipeline as set2_2_comove_optimized/fluctuations.wls.  What differs:

     parameters  params.md: ell_d = 1, S0 = 0.9, xi0 = 5, rho0 = 0.25,
                 Lambda = 0.13, gamma = 10, xi_r = 3 xi0, zeta = 2.5 xi0.  This
                 gives xi_r = 15, zeta = 12.5, a' = -18.75, a = -20.7,
                 b = 11.574, K' = 18.75, K = 9.375, alpha = 0.45.  B is swept by
                 the jobscripts; the default is params.md's guess, 1.4e5.
     steady      comoving_steady_Qfree_solver.m, free (Neumann) Q wall with the
                 mirror line Q2(x,0) = 0, and the Dirichlet density wall
                 (rho = 1).  There is no mu for this wall; the solver reports
                 0, which is printed and saved.  "ReuseFactorization" -> False:
                 at xi_r = 3 xi0 the default modified Newton diverges at the first
                 zeta step.
     Lyapunov    Dirichlet (reservoir) wall, as before.  The steady state and the
                 fluctuations therefore use the same density wall; the Q walls
                 still differ (steady n.grad Q = 0, fluctuation dQ = 0).
     plots       only the fluctuation plots.  The steady-state plots needed
                 iterative_solver_module.m (now in FRIDGE/); the steady state is
                 saved in full to steady_<tag>.m anyway.

   The runTag carries box, Nint, wall, lNoise, zeta, ell_d, B and the steady box,
   but NOT rho0, Lambda or the steady walls.  Keep this set in its own outputDir.
   The steady walls are recorded in the meta files (otherParams "steadyBC",
   solverParams qWall).

   For the two caveats on reading this data (::unstable on coarse grids, and the
   pedestal calibration being the rho-rho channel only) see the header of
   set2_2_comove_optimized/fluctuations.wls and lyapunov_solver/README.md. *)

thisDir = DirectoryName[$InputFileName];
repoRoot = ParentDirectory[ParentDirectory[thisDir]];

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
SolveActiveNematicComovingSteadyQfree[modelParams, solverParams, bcType]
   Self-Gets comoving_defect/comoving_core.m.  Loaded LAST: it opens no private
   context, so it defines many bare names in Global`, but Block localises all of
   them and it leaves nothing carrying a value, so it cannot disturb the two
   modules above.
*)
Import[FileNameJoin[{repoRoot, "comoving_defect", "comoving_steady_Qfree_solver.m"}]]

steadyBC    = "Dirichlet";    (* density wall of the steady state *)
lyapunovBC  = "Dirichlet";    (* density wall of the fluctuations *)

(* --- Command-line overrides ------------------------------------------
   wolframscript -file fluctuations.wls \
     "modelParams={B->110000.}; solverParams={box->20.}; lyapunovSolverParams={Nint->101, fbox->20}"
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
   params.md.  a1 = a' and L1 = K' are the activity-renormalised values;
   gamma0 = B/xi0. *)
paramRelations = {ell == Sqrt[L1/(-a1)],
   S0 == Sqrt[(-a1)/(2 b)], \[Gamma] == -8 a1/\[Xi]r, \[Xi]r == 3 \[Xi]0,
   a1 == a + \[CapitalLambda] \[Xi]r,
   L1 == L + (\[Zeta] \[Xi]r)/(4 \[Xi]0), \[Gamma]0 == B/\[Xi]0,
   \[Zeta] == 2.5 \[Xi]0};
knownParams = {ell -> 1, S0 -> 0.9, \[Xi]0 -> 5, \[Rho]0 -> 0.25,
   \[CapitalLambda] -> 0.13};

(* A B given on the command line enters HERE, before Solve, so gamma0 = B/xi0
   follows it.  Merging it into modelParams afterwards would leave gamma0 at the
   default B. *)
guessedParams = {B -> Lookup[Association @ Lookup[cliParams, "modelParams", {}],
                             B, 140000.],
                 \[Gamma] -> 10};

inferedParams = First@Solve[paramRelations /. Join[knownParams, guessedParams],
    {b, a, L, a1, L1, \[Gamma]0, \[Xi]r, \[Zeta]}];

(* a1, L1 and gamma0 are carried for the record only: both solvers Lookup the
   keys they need and ignore the rest.  lNoise is NOT derived from the others --
   it is the physical UV cutoff, and with rescaleCutoff the solver turns it into
   the lNoiseEff that W is built from. *)
derivedNames = {B, a, b, L, \[Zeta], \[Xi]0, \[Xi]r, \[Rho]0, \[CapitalLambda],
   a1, L1, \[Gamma]0};

modelParams = upd["modelParams",
  Append[
    Thread[derivedNames -> (derivedNames /.
       Join[knownParams, guessedParams, inferedParams])],
    lNoise -> 0.5]];

(* The required comoving keys plus the optional ones, all listed so the shell
   can override any of them: MergeParams aborts on a key absent from the
   defaults.  qWall, symLine, symOps and "ReuseFactorization" are the Q-free
   solver's; box is swept by the jobscripts and must be >= fbox below. *)
solverParams = upd["solverParams",
  {\[Delta]mesh -> 0.04, hMax -> 2., box -> 20.,
   maxSteps -> 25, relativeDiff -> 1.*^-8,
   refine -> 1, rhoGauge -> "Point", zetaSteps -> 5, order -> 4,
   stepTolerance -> 1.*^-10, verbose -> True,
   qWall -> "Neumann", symLine -> True, symOps -> True,
   "ReuseFactorization" -> False}];

(* fbox MUST be an Integer -- the module tests IntegerQ -- whereas solverParams
   box is a Real that only ever reaches N[box]. *)
(* scratchDir must be listed here even though only the Python backend reads it:
   MergeParams aborts on a key absent from the defaults, so a jobscript passing
   it would otherwise kill the run before any work happened.  Automatic ->
   $TMPDIR, then /scratch/$SLURM_JOB_ID, then /scratch. *)
lyapunovSolverParams = upd["lyapunovSolverParams",
  {Nint -> 21, fbox -> 20, scratchDir -> Automatic}];

otherParams = upd["otherParams",
  {outputDir -> "./output", saveData -> True, plotBox -> 17}];

Print["modelParams          = ", modelParams];
Print["solverParams         = ", solverParams];
Print["lyapunovSolverParams = ", lyapunovSolverParams];
Print["otherParams          = ", otherParams];
Print["walls                = steady ", steadyBC, " (Q ", qWall /. solverParams,
      "),  Lyapunov ", lyapunovBC];

outDir  = outputDir /. otherParams;   (* None -> nothing is written *)
plotBoxV = plotBox  /. otherParams;
saveDataV = saveData /. otherParams;

boxV  = box  /. solverParams;
fboxV = fbox /. lyapunovSolverParams;
(* fbox <= box.  A LARGER steady box is safe -- the Lyapunov module samples the
   steady interpolant out to +-(fbox - h/2), so it stays inside.  A SMALLER one
   extrapolates and is fatal. *)
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

Print["[1] Comoving steady state, free Q wall, ", steadyBC, " density wall ..."];
{tSS, steadyState} = AbsoluteTiming[
  SolveActiveNematicComovingSteadyQfree[modelParams, solverParams, steadyBC]];
If[steadyState === $Failed,
   Print["FATAL: SolveActiveNematicComovingSteadyQfree returned $Failed"]; Abort[]];
Print["    done in ", tSS, " s"];
Print["    u = ", steadyState["u"],
      "   mu = ", steadyState["solver"]["mu"],
      "   lineForce = ", steadyState["solver"]["lineForce"],
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
                                         steadyState, lyapunovBC];
If[correlations === $Failed,
   Print["FATAL: SolveFluctuationsLyapunov returned $Failed"]; Abort[]];

(* Everything the module cannot see, attached here so exportDataComoving writes
   it into BOTH sigma_rho_*.m and covariance_*_meta.m.  It saves only a fixed
   set of keys, so the steady walls and mu ride inside "otherParams"; qWall is
   already in solverParams. *)
correlations = Join[correlations,
  <|"solverParams" -> solverParams,
    "otherParams" -> Join[otherParams,
       {"steadyBC" -> steadyBC, "lyapunovBC" -> lyapunovBC,
        "steadyMu" -> steadyState["solver"]["mu"]}],
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
   exits 0.  Do the two lines here instead. *)
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
    <|"runTag" -> tag, "frame" -> "comoving",
      "bc" -> steadyBC, "qWall" -> (qWall /. solverParams),
      "lyapunovBC" -> lyapunovBC,
      (* the point of this data set *)
      "u" -> ss["u"], "mu" -> ss["solver"]["mu"],
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
   DATA FIRST, PLOTS LAST.  Rasterising a plot needs the front end, which can
   kill the kernel outright on a headless node; the data above is the expensive
   product.  Trailing ";" is load-bearing: a bare plot would make wolframscript
   print the whole Graphics expression to stdout. *)
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

(* Plotting is best-effort and deliberately LAST.  Quiet, because a kernel
   without a usable front end answers Rasterize with Rasterize::type, and that
   message echoes the whole Graphics expression into the log. *)
Quiet @ Check[
   visualizeFluctuationsComoving[correlations, outDir, plotBoxV];,
   Print["WARNING: plotting failed (no usable front end?).  The data files ",
         "listed above are unaffected."]];
