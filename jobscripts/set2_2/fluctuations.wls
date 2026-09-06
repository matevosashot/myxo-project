#!/usr/bin/env wolframscript
(* ::Package:: *)

(* set2_2 -- IDENTICAL physics to set-2/set2_1, consistent stencils (R == 1), and a
   PHYSICAL ultraviolet cutoff: the noise carries a Gaussian spatial correlation of
   length lNoise instead of being delta-correlated in space.

       Sigma_iso(r) = zeta rho_ss(r)/(2 Pi B rho0 lNoise^2)

   which is dean.tex eq:iso-final with l_UV read as the noise correlation length.
   Unlike set-2 and set2_1 this has a finite h -> 0 limit, so a sweep over Nint at
   fixed lNoise measures convergence to a physical number rather than to h^-2.
   See dean.tex appendix sec-noise-cutoff.

   Every model and solver parameter apart from lNoise is byte-for-byte the set-2
   value, so the three data sets are directly comparable.

   NB lNoise must exceed ~5 h for the continuum value to within 1%; below lNoise ~ h
   the Gaussian is invisible to the grid and the run reverts to the grid-limited
   set2_1 answer.  The module warns when lNoise < 5 h.

   Output files carry BOTH the grid and an "_lN<value>" tag, so no two (Nint, lNoise)
   runs can collide. *)

thisDir = DirectoryName[$InputFileName];
repoRoot = ParentDirectory[ParentDirectory[thisDir]];

(*
SolveActiveNematicSteady[modelParams, solverParams, bcType]
visualizeSteadyState[steadyState, outputDir, plotBox]
*)
Import[FileNameJoin[{repoRoot, "iterative_solver_module.m"}]]

(*
SolveFluctuationsLyapunovCutoff[lyapunovParams, modelParams, steadyState]
visualizeFluctuations[correlations, outputDir, plotBox]
exportData[correlations, outputDir, saveData]
MergeParams, ScriptParamOverrides
   NB do NOT also load lyapunov_solver_module.m or _R1.m: all three share the
   private helpers (runTag, exportData, ...) and the last one loaded wins.
*)
Import[FileNameJoin[{repoRoot, "lyapunov_solver_module_cutoff.m"}]]

(* --- Command-line overrides ------------------------------------------
   wolframscript -file fluctuations.wls \
     "modelParams={lNoise->1.}; lyapunovSolverParams={Nint->51}"
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
   Identical to set-2 apart from lNoise.  a, b, L follow from ell, S0, gamma. *)
paramRelations={ell==Sqrt[(L+(\[Zeta] \[Xi]r)/(4 \[Xi]0))/(-a-\[CapitalLambda] \[Xi]r)],S0==Sqrt[(-a-\[CapitalLambda] \[Xi]r)/(2b)],\[Gamma]==8 b / \[Xi]r};
knownParams={ell->2.5,S0->0.9,\[Xi]0->5.,\[Zeta]->10,\[Rho]0->0.3};

guessedParams={B->1000.,b->0.1,\[Gamma]->8.,\[CapitalLambda]->0.1};

inferedParams=First@Solve[paramRelations/.Join[knownParams,guessedParams],{\[Xi]r,a,L}];

(* lNoise is NOT derived from the others -- it is the new physical input.  Default
   Sqrt[3.]: matching the second moment of a uniform rod of length 6 um gives
   lNoise^2 = L^2/12 = 3 um^2. *)
derivedNames = {B,a,b,L,\[Zeta],\[Xi]0,\[Xi]r,\[Rho]0,\[CapitalLambda]};

modelParams = upd["modelParams",
  Append[
    Thread[derivedNames->(derivedNames/.Join[knownParams,guessedParams,inferedParams])],
    lNoise -> Sqrt[3.]]];

solverParams = upd["solverParams",
  {\[Delta]mesh -> 0.061, hMax -> 3, box -> 20,
   maxSteps -> 8, relativeDiff -> 1.*^-5}];

lyapunovSolverParams = upd["lyapunovSolverParams",
  {Nint -> 21, fbox -> 10}];

otherParams = upd["otherParams",
  {outputDir -> "./output", saveData -> True, plotBox -> 9}];

Print["modelParams          = ", modelParams];
Print["solverParams         = ", solverParams];
Print["lyapunovSolverParams = ", lyapunovSolverParams];
Print["otherParams          = ", otherParams];

outDir  = outputDir /. otherParams;   (* None -> nothing is written *)
plotBoxV = plotBox  /. otherParams;
saveDataV = saveData /. otherParams;

(* --- Steady state ---------------------------------------------------- *)

Print["[1] Steady state with Dirichlet BC ..."];
{tSS, steadyState} = AbsoluteTiming[
  SolveActiveNematicSteady[modelParams, solverParams, "Dirichlet"]];
Print["    done in ", tSS, " s"];

(* --- Fluctuations ---------------------------------------------------- *)

correlations = SolveFluctuationsLyapunovCutoff[lyapunovSolverParams, modelParams,
                                              steadyState];

(* so the meta file can report a wall time that includes the FEM steady state *)
correlations = Append[correlations, "steadyStateTime" -> tSS];

(* visualizeSteadyState is deliberately NOT called: it writes UNTAGGED
   steady_rho.png / steady_Q.png, and this sweep sends 90 array tasks into one
   directory, so they would be clobbered 90 times.  The steady state does not
   depend on Nint or lNoise and is already saved under set2_1.

   Trailing ";" is load-bearing: a bare plot would make wolframscript print the
   whole Graphics expression -- every mesh vertex -- to stdout. *)
visualizeFluctuations[correlations, outDir, plotBoxV];
exportData[correlations, outDir, saveDataV];
