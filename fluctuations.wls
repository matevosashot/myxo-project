#!/usr/bin/env wolframscript
(* ::Package:: *)

repoRoot = DirectoryName[$InputFileName];

(*
SolveActiveNematicSteady[modelParams, solverParams, bcType]
visualizeSteadyState[steadyState, outputDir, plotBox]
*)
Import[FileNameJoin[{repoRoot, "iterative_solver_module.m"}]]

(*
SolveFluctuationsLyapunov[lyapunovParams, modelParams, steadyState]
visualizeFluctuations[correlations, outputDir, plotBox]
exportData[correlations, outputDir, saveData]
MergeParams, ScriptParamOverrides
*)
Import[FileNameJoin[{repoRoot, "lyapunov_solver_module.m"}]]

(* --- Command-line overrides ------------------------------------------
   wolframscript -file fluctuations.wls \
     "modelParams={B->2000., \[Zeta]->5.}; lyapunovSolverParams={Nint->51}"
   Each named list is UPDATED key-by-key: keys not mentioned keep their
   default.  preParams is applied before modelParams is derived from it,
   so an explicit modelParams key wins over the derived one.  An unknown
   set name or key aborts. *)

cliParams = ScriptParamOverrides[];
knownParamSets = {"preParams", "modelParams", "solverParams",
                  "lyapunovSolverParams", "otherParams"};
If[Complement[Keys[cliParams], knownParamSets] =!= {},
   Print["ERROR: unknown parameter set(s) ",
         Complement[Keys[cliParams], knownParamSets],
         "; known: ", knownParamSets];
   Abort[]];
upd[name_String, defaults_] :=
  MergeParams[defaults, Lookup[cliParams, name, {}], name];

(* --- Model & solver parameters (defaults) --------------------------- *)

preParams = upd["preParams",
  {ell -> 2.5, S0 -> 0.9, A -> -1}];

modelParams = upd["modelParams",
  {B -> 1000., a -> A, b -> -A/(2 S0^2), \[Xi]0 -> 5.,
   \[Xi]r -> 0.1, \[Zeta] -> 10., \[Rho]0 -> 1.,
   L -> ell^2 (-A), \[CapitalLambda] -> 0.1} /. preParams];

solverParams = upd["solverParams",
  {\[Delta]mesh -> 0.061, hMax -> 3, box -> 20,
   maxSteps -> 8, relativeDiff -> 1.*^-4}];

lyapunovSolverParams = upd["lyapunovSolverParams",
  {Nint -> 21, fbox -> 10}];

otherParams = upd["otherParams",
  {outputDir -> "./output", saveData -> True, plotBox -> 9}];

Print["preParams            = ", preParams];
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

correlations = SolveFluctuationsLyapunov[lyapunovSolverParams, modelParams,
                                         steadyState];

(* Trailing ";" is load-bearing: a bare plot would make wolframscript print
   the whole Graphics expression -- every mesh vertex -- to stdout. *)
visualizeSteadyState[steadyState, outDir, plotBoxV];
visualizeFluctuations[correlations, outDir, plotBoxV];
exportData[correlations, outDir, saveDataV];
