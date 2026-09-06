#!/usr/bin/env wolframscript
(* ::Package:: *)

(* set3_1 -- the SEALED-WALL counterpart of set2_1.

   set2_1 realises the RESERVOIR boundary condition: rho_ss pinned to 1 on
   \[PartialD]\[CapitalOmega], psi|_bdry = 0, noise flux free.  This realises the
   other consistent combination of dean.tex sec-bc, the SEALED one:

       n.Grad rho_ss = 0,    n.Grad d(rho) = 0,    noise flux BLOCKED,

   i.e. a container wall that nothing crosses, rather than a window cut out of a
   larger colony.

   The wall also ANCHORS THE DIRECTOR, dQ|_bdry = 0, exactly as the steady state
   does.  Sealing is a statement about mass and only rho is conserved, so the flux
   argument of sec-bc fixes the rho sector alone; the director condition is the
   separate physical question of whether the surface exerts a torque, and it must.
   Measured (claude_experiments/fluctuations_neumann/instability_probe.wls): with
   dQ Neumann as well, Atilde carries three unstable eigenvalues that CONVERGE
   under refinement (+0.104 -> +0.064 at Nint = 21 -> 61) -- two defect
   translations and the global director rotation.  A free-anchoring wall images a
   +1/2 defect attractively, so the centred defect is a saddle, there is no stable
   steady state, and no equal-time covariance exists.  Anchoring removes all three
   and leaves exactly one zero eigenvalue: the conserved total mass of
   eq:canonical.

   Two changes from set2_1 and nothing else:

     - SolveActiveNematicSteady[..., "Neumann"] instead of "Dirichlet";
     - lyapunov_solver_module_neumann.m instead of lyapunov_solver_module_R1.m.

   Every model and solver parameter is byte-for-byte the set-2 / set2_1 value, so
   the runs are directly comparable.  The Neumann grid is cell-centred with
   h = 2 fbox/Nint, so Nint here has EXACTLY the same h as Nint-1 in set2_1 --
   the shared NintValues list therefore supplies an h-matched Dirichlet partner
   for every odd run.

   What to expect, from dean.tex sec-bc: the pedestal is FLAT and identical to
   set2_1's in the bulk, because A = -M K and Q = 2 T M make the mobility M --
   which is where the boundary condition lives -- cancel out of C = T K^-1.  The
   two boundary conditions differ in exactly two places: the rank-one 1/|Omega|
   canonical subtraction of eq:canonical (the total is now conserved), and the
   sign of the anisotropic wall residue of eq:wall-kernel, since the Neumann mode
   weight at the wall favours SMALL normal wavenumber where the Dirichlet weight
   favoured large.

   Output files carry a "_Neu" tag, so they cannot collide with either earlier
   set even if they share a directory. *)

thisDir = DirectoryName[$InputFileName];
repoRoot = ParentDirectory[ParentDirectory[thisDir]];

(*
SolveActiveNematicSteady[modelParams, solverParams, bcType]
visualizeSteadyState[steadyState, outputDir, plotBox]
*)
Import[FileNameJoin[{repoRoot, "iterative_solver_module.m"}]]

(*
SolveFluctuationsLyapunovNeumann[lyapunovParams, modelParams, steadyState]
visualizeFluctuations[correlations, outputDir, plotBox]
exportData[correlations, outputDir, saveData]
MergeParams, ScriptParamOverrides
   NB do NOT also load lyapunov_solver_module.m or _R1.m: they share the private
   helpers (runTag, exportData, ...) and the last one loaded wins.
*)
Import[FileNameJoin[{repoRoot, "lyapunov_solver_module_neumann.m"}]]

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
   Identical to set-2 / set2_1.  a, b, L follow from ell, S0, gamma. *)
paramRelations={ell==Sqrt[(L+(\[Zeta] \[Xi]r)/(4 \[Xi]0))/(-a-\[CapitalLambda] \[Xi]r)],S0==Sqrt[(-a-\[CapitalLambda] \[Xi]r)/(2b)],\[Gamma]==8 b / \[Xi]r};
knownParams={ell->2.5,S0->0.9,\[Xi]0->5.,\[Zeta]->10,\[Rho]0->0.3};

guessedParams={B->1000.,b->0.1,\[Gamma]->8.,\[CapitalLambda]->0.1};

inferedParams=First@Solve[paramRelations/.Join[knownParams,guessedParams],{\[Xi]r,a,L}];

modelParamsNames={B,a,b,L,\[Zeta],\[Xi]0,\[Xi]r,\[Rho]0,\[CapitalLambda]};

modelParams = upd["modelParams",
  Thread[modelParamsNames->(modelParamsNames/.Join[knownParams,guessedParams,inferedParams])]];

solverParams = upd["solverParams",
  {\[Delta]mesh -> 0.061, hMax -> 3, box -> 20,
   maxSteps -> 8, relativeDiff -> 1.*^-5}];

lyapunovSolverParams = upd["lyapunovSolverParams",
  {Nint -> 21, fbox -> 10}];

(* plotBox defaults to fbox: with sealed walls the boundary IS the point of the
   run, so nothing is cropped out of the figures. *)
otherParams = upd["otherParams",
  {outputDir -> "./output", saveData -> True, plotBox -> 20}];

Print["modelParams          = ", modelParams];
Print["solverParams         = ", solverParams];
Print["lyapunovSolverParams = ", lyapunovSolverParams];
Print["otherParams          = ", otherParams];

outDir  = outputDir /. otherParams;   (* None -> nothing is written *)
plotBoxV = plotBox  /. otherParams;
saveDataV = saveData /. otherParams;

(* --- Steady state ---------------------------------------------------- *)

(* "Neumann": rho is free on the outer boundary and pinned to 1 only at the
   single boundary node x == 0, y = -box, which fixes the additive constant of an
   otherwise singular pure-Neumann problem.  Q keeps its Pade Dirichlet
   anchoring -- so this is "rho sealed, director anchored"; see the note at the
   top of lyapunov_solver_module_neumann.m. *)
Print["[1] Steady state with Neumann BC ..."];
{tSS, steadyState} = AbsoluteTiming[
  SolveActiveNematicSteady[modelParams, solverParams, "Neumann"]];
Print["    done in ", tSS, " s"];

(* --- Fluctuations ---------------------------------------------------- *)

correlations = SolveFluctuationsLyapunovNeumann[lyapunovSolverParams, modelParams,
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
