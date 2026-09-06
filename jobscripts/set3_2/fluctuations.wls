#!/usr/bin/env wolframscript
(* ::Package:: *)

(* set3_2 -- set3_1 plus a PHYSICAL ultraviolet cutoff.

   Same sealed, anchoring wall as set3_1, but the noise carries a Gaussian
   spatial correlation of length lNoise instead of being delta-correlated in
   space:

       Sigma_iso(r) = zeta rho_ss(r)/(2 Pi B rho0 lNoise^2)

   which is dean.tex eq:iso-final with l_UV read as the noise correlation
   length.  Unlike set3_1 this has a finite h -> 0 limit, so a sweep over Nint
   at fixed lNoise measures convergence to a physical number rather than to
   h^-2.  See dean.tex sec-noise-cutoff.

   The kernel is built on the NEUMANN Laplacian, so it REFLECTS off the sealed
   wall: a rod cannot straddle a wall it cannot cross.  Its diagonal is
   therefore ENHANCED near the wall, toward a factor 2 for lNoise >> h, where
   set2_2's Dirichlet kernel is SUPPRESSED to zero there.  That contrast is the
   sharpest check that the boundary condition is doing what it claims.

   NB lNoise must exceed ~5 h for the continuum value to within 1%; below
   lNoise ~ h the Gaussian is invisible to the grid and the run reverts to the
   grid-limited set3_1 answer.  The module warns when lNoise < 5 h.

   ---------------------------------------------------------------------
   set3_1 -- the SEALED-WALL counterpart of set2_1.

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
     - lyapunov_solver_module_neumann_cutoff.m instead of the R1 module.

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

   Output files carry BOTH the grid and a "_Neu_lN<value>" tag, so no two
   (Nint, lNoise) runs can collide. *)

thisDir = DirectoryName[$InputFileName];
repoRoot = ParentDirectory[ParentDirectory[thisDir]];

(*
SolveActiveNematicSteady[modelParams, solverParams, bcType]
visualizeSteadyState[steadyState, outputDir, plotBox]
*)
Import[FileNameJoin[{repoRoot, "iterative_solver_module.m"}]]

(*
SolveFluctuationsLyapunovNeumannCutoff[lyapunovParams, modelParams, steadyState]
visualizeFluctuations[correlations, outputDir, plotBox]
exportData[correlations, outputDir, saveData]
MergeParams, ScriptParamOverrides
   NB do NOT also load lyapunov_solver_module.m, _R1.m, _cutoff.m or _neumann.m:
   they share the private helpers (runTag, exportData, ...) and the last one
   loaded wins.
*)
Import[FileNameJoin[{repoRoot, "lyapunov_solver_module_neumann_cutoff.m"}]]

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
   Identical to set-2 / set2_1 apart from lNoise.  a, b, L follow from
   ell, S0, gamma. *)
paramRelations={ell==Sqrt[(L+(\[Zeta] \[Xi]r)/(4 \[Xi]0))/(-a-\[CapitalLambda] \[Xi]r)],S0==Sqrt[(-a-\[CapitalLambda] \[Xi]r)/(2b)],\[Gamma]==8 b / \[Xi]r};
knownParams={ell->2.5,S0->0.9,\[Xi]0->5.,\[Zeta]->10,\[Rho]0->0.3};

guessedParams={B->1000.,b->0.1,\[Gamma]->8.,\[CapitalLambda]->0.1};

inferedParams=First@Solve[paramRelations/.Join[knownParams,guessedParams],{\[Xi]r,a,L}];

(* lNoise is NOT derived from the others -- it is the new physical input.  Default
   Sqrt[3.]: matching the second moment of a uniform rod of length 6 um gives
   lNoise^2 = L^2/12 = 3 um^2. *)
modelParamsNames={B,a,b,L,\[Zeta],\[Xi]0,\[Xi]r,\[Rho]0,\[CapitalLambda]};

modelParams = upd["modelParams",
  Append[
    Thread[modelParamsNames->(modelParamsNames/.Join[knownParams,guessedParams,inferedParams])],
    lNoise -> Sqrt[3.]]];

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

correlations = SolveFluctuationsLyapunovNeumannCutoff[lyapunovSolverParams,
                                                      modelParams, steadyState];

(* so the meta file can report a wall time that includes the FEM steady state *)
correlations = Append[correlations, "steadyStateTime" -> tSS];

(* visualizeSteadyState is deliberately NOT called: it writes UNTAGGED
   steady_rho.png / steady_Q.png, and this sweep sends 54 array tasks into one
   directory, so they would be clobbered 54 times.  The steady state depends on
   neither Nint nor lNoise and is already saved under set3_1.

   Trailing ";" is load-bearing: a bare plot would make wolframscript print the
   whole Graphics expression -- every mesh vertex -- to stdout. *)
visualizeFluctuations[correlations, outDir, plotBoxV];
exportData[correlations, outDir, saveDataV];
