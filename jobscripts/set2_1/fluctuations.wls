#!/usr/bin/env wolframscript
(* ::Package:: *)

(* set2_1 -- IDENTICAL physics to set-2, CONSISTENT stencils.

   The only difference from jobscripts/set-2/fluctuations.wls is that this loads
   lyapunov_solver_module_R1.m instead of lyapunov_solver_module.m, so the
   divergence-form noise is discretized with the adjoint partner of the drift
   Laplacian and the UV pedestal comes out

       Sigma_iso(r) = zeta rho_ss(r)/(B rho0 h^2)      (R == 1)

   instead of <R> = 0.363380 times that.  Every model and solver parameter below
   is byte-for-byte the set-2 value, so the two runs are directly comparable.
   See dean.tex appendix "Linear OU field theory" and the R1 module header.

   Output files carry an "_R1" tag, so they cannot collide with set-2 even if the
   two ever share a directory. *)

thisDir = DirectoryName[$InputFileName];
repoRoot = ParentDirectory[ParentDirectory[thisDir]];

(*
SolveActiveNematicSteady[modelParams, solverParams, bcType]
visualizeSteadyState[steadyState, outputDir, plotBox]
*)
Import[FileNameJoin[{repoRoot, "iterative_solver_module.m"}]]

(*
SolveFluctuationsLyapunovR1[lyapunovParams, modelParams, steadyState]
visualizeFluctuations[correlations, outputDir, plotBox]
exportData[correlations, outputDir, saveData]
MergeParams, ScriptParamOverrides
   NB do NOT also load lyapunov_solver_module.m: the two share the private
   helpers (runTag, exportData, ...) and the last one loaded wins.
*)
Import[FileNameJoin[{repoRoot, "lyapunov_solver_module_R1.m"}]]

(* --- Command-line overrides ------------------------------------------
   wolframscript -file fluctuations.wls \
     "modelParams={B->2000., \[Zeta]->5.}; lyapunovSolverParams={Nint->51}"
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
   Identical to set-2.  a, b, L follow from ell, S0, gamma. *)
paramRelations={ell==Sqrt[(L+(\[Zeta] \[Xi]r)/(4 \[Xi]0))/(-a-\[CapitalLambda] \[Xi]r)],S0==Sqrt[(-a-\[CapitalLambda] \[Xi]r)/(2b)],\[Gamma]==8 b / \[Xi]r};
knownParams={ell->2.5,S0->0.9,\[Xi]0->5.,\[Zeta]->10,\[Rho]0->0.3};

guessedParams={B->1000.,b->0.1,\[Gamma]->8.,\[CapitalLambda]->0.1};

inferedParams=First@Solve[paramRelations/.Join[knownParams,guessedParams],{\[Xi]r,a,L}];

modelParamsNames={B,a,b,L,\[Zeta],\[Xi]0,\[Xi]r,\[Rho]0,\[CapitalLambda]};

(* set-2 assigned this directly, so a CLI "modelParams={...}" passed the whitelist
   and was then silently ignored.  Routed through upd here. *)
modelParams = upd["modelParams",
  Thread[modelParamsNames->(modelParamsNames/.Join[knownParams,guessedParams,inferedParams])]];

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

correlations = SolveFluctuationsLyapunovR1[lyapunovSolverParams, modelParams,
                                           steadyState];

(* so the meta file can report a wall time that includes the FEM steady state *)
correlations = Append[correlations, "steadyStateTime" -> tSS];

(* Trailing ";" is load-bearing: a bare plot would make wolframscript print
   the whole Graphics expression -- every mesh vertex -- to stdout. *)
(* TAGGED, because the steady state depends on solverParams box: an untagged
   steady_rho.png would be overwritten whenever a box-20 and a box-40 sweep share
   one outputDir.  Nint is in the tag too so no two concurrent array tasks Export
   to the same path. *)
visualizeSteadyState[steadyState, outDir, plotBoxV,
  "_box" <> ToString[box /. solverParams] <>
  "_N"   <> ToString[Nint /. lyapunovSolverParams]];
visualizeFluctuations[correlations, outDir, plotBoxV];
exportData[correlations, outDir, saveDataV];
