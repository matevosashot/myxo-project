(* ::Package:: *)

(* ::Section:: *)
(*Comoving-frame steady-state solver for the coupled \[Rho]\[Dash]Q active-nematic system.*)

(* The numerics live in comoving_core.m in this same directory; this file is a
   thin wrapper that puts them behind the repo's standard module interface
   (modelParams / solverParams / bcType -> <|"rho","Q1","Q2","mesh"|>).
   Guarded self-load, the idiom lyapunov_solver_module.m:21-22 uses. *)
If[Names["Global`cdSteady"] === {},
   Get[FileNameJoin[{DirectoryName[$InputFileName], "comoving_core.m"}]]];

Clear[SolveActiveNematicComovingSteady];

SolveActiveNematicComovingSteady::usage =
  "SolveActiveNematicComovingSteady[modelParams, solverParams, bcType] " <>
  "solves the coupled steady-state \[Rho]\[Dash]Q active-nematic system IN THE FRAME " <>
  "THAT MOVES WITH A +1/2 DEFECT, with the defect core held at the origin, and " <>
  "returns the fields together with the frame velocity u.\n\n" <>
  "  modelParams  : rules {B->..., a->..., b->..., \[Xi]0->..., \[Xi]r->..., " <>
  "\[Zeta]->..., \[Rho]0->..., L->..., \[CapitalLambda]->...}  (L is the Frank " <>
  "constant, as in iterative_solver_module.m)\n" <>
  "  solverParams : rules {\[Delta]mesh->..., hMax->..., box->..., " <>
  "maxSteps->..., relativeDiff->...}\n" <>
  "                 optional: refine->1 (scales every cell together), " <>
  "rhoGauge->\"Point\", zetaSteps->5, order->4, stepTolerance->1.*^-10, " <>
  "verbose->True\n" <>
  "  bcType       : \"Neumann\"   \[Rule] sealed wall, " <>
  "n.(B \[Del]\[Rho] + \[Zeta] \[Del].Q) = 0, plus a single-node density gauge;\n" <>
  "                 \"Dirichlet\" \[Rule] \[Rho] = 1 on the entire outer boundary.\n" <>
  "                 Q is held at the Pade +1/2 profile on the outer boundary in " <>
  "both cases.\n\n" <>
  "Returns <|\"rho\"->..., \"Q1\"->..., \"Q2\"->..., \"mesh\"->..., \"u\"->{ux,uy}, " <>
  "\"solver\"->...|>.  The three fields are 2-argument interpolations of (x,y) -- " <>
  "drop-in for visualizeSteadyState and the Lyapunov modules, exactly like " <>
  "SolveActiveNematicSteady.  \"mesh\" is the finite-difference tensor grid as an " <>
  "association (there is no ElementMesh here); it carries \"Coordinates\" in the " <>
  "same n x 2 shape an ElementMesh would, plus \"gx\", \"gy\", \"hMin\".  " <>
  "\"solver\" holds <|\"converged\",\"residual\",\"iterations\",\"core\",\"wall\", " <>
  "\"zetaPath\",\"gauge\"|>.\n\n" <>
  "u is NOT measured from the solution.  It is the Lagrange multiplier conjugate " <>
  "to the pinning constraints Q1(0,0) = Q2(0,0) = 0, so the constraint holds " <>
  "exactly and u falls out as a component of the solution vector.  See " <>
  "dean.tex, appendix \"Comoving-frame defect solver\".\n\n" <>
  "NB 1: this is NOT the lab-frame steady state.  The defect is pinned and the " <>
  "whole problem is posed in a frame translating at u, so the result is not " <>
  "comparable with SolveActiveNematicSteady except in the passive limit " <>
  "\[Zeta] -> 0, where u -> 0 and the two agree.\n\n" <>
  "NB 2: \"Neumann\" here is NOT the wall iterative_solver_module.m imposes.  " <>
  "That solver leaves \[Rho] free and pins it on the whole ray x==0, y<0 -- many " <>
  "nodes, over-constrained.  Here the wall is the full sealed condition and the " <>
  "residual constant null mode is removed by ONE node (rhoGauge->\"Point\", the " <>
  "default) or by the mean (\"Mean\") or by overwriting one near-wall PDE row " <>
  "(\"Replace\").  The three gauges differ only by an additive constant in " <>
  "\[Rho], to which u is blind; they agree on u to ten significant figures.  The " <>
  "density is shifted back to unit mean before it is returned.";

SolveActiveNematicComovingSteady::badbc =
  "Unknown boundary-condition type `1`; expected \"Neumann\" or \"Dirichlet\".";

SolveActiveNematicComovingSteady::badgauge =
  "Unknown rhoGauge `1`; expected \"Point\", \"Mean\" or \"Replace\".";

SolveActiveNematicComovingSteady::badmodel =
  "modelParams is missing or non-numeric for: `1`.  Expected rules for " <>
  "{B, a, b, L, \[Zeta], \[Xi]0, \[Xi]r, \[Rho]0, \[CapitalLambda]}.";

SolveActiveNematicComovingSteady::badsolver =
  "solverParams is missing or non-numeric for: `1`.  Expected rules for " <>
  "{\[Delta]mesh, hMax, box, maxSteps, relativeDiff}.";

SolveActiveNematicComovingSteady::noconv =
  "Bordered Newton did not converge at the target \[Zeta] within maxSteps = `1` " <>
  "iterations (final |step| tolerance `2`, residual = `3`).  u = `4` is reported " <>
  "anyway and should not be trusted.";

SolveActiveNematicComovingSteady::stalled =
  "\[Zeta]-continuation stalled at \[Zeta] = `1` (residual = `2`); continuing.";


SolveActiveNematicComovingSteady[
    modelParams : {___Rule},
    solverParams : {___Rule},
    bcType_String : "Neumann"] :=
Block[
  {
    (* solver knobs (unpacked from solverParams) *)
    \[Delta]mesh, hMax, box, maxSteps, relativeDiff,
    refine, rhoGauge, zetaSteps, order, stepTolerance, verbose,
    (* model parameters (unpacked from modelParams) *)
    B, a, b, L, \[Zeta], \[Xi]0, \[Xi]r, \[Rho]0, \[CapitalLambda],
    (* derived / working *)
    missing, mkModel, dpOf, dpTarget, dp, pf, gg, op, bcz, ne, n, k0,
    X, ex, st, zs, kk, wall0, rhoGaugeV
  },

  (* --- validate bcType up front ------------------------------------- *)
  If[!MemberQ[{"Neumann", "Dirichlet"}, bcType],
    Message[SolveActiveNematicComovingSteady::badbc, bcType];
    Return[$Failed]
  ];

  (* --- unpack solver parameters ------------------------------------- *)
  (* Lookup against an Association rather than the house "{k,...} = {k,...} /.
     solverParams" idiom.  The two agree whenever the key is present, but the
     replacement form self-assigns (box = box) when it is absent, and a
     self-referential OwnValue blows the recursion limit the next time the
     symbol is touched.  Lookup returns Missing instead, which the NumericQ
     screens below turn into a clean ::badsolver.  The keys are the caller's
     Global` symbols and Block localises their VALUES, not their identity, so
     they still match as keys. *)
  With[{sp = Association[solverParams]},
    \[Delta]mesh  = Lookup[sp, \[Delta]mesh,  Missing[]];
    hMax          = Lookup[sp, hMax,          Missing[]];
    box           = Lookup[sp, box,           Missing[]];
    maxSteps      = Lookup[sp, maxSteps,      Missing[]];
    relativeDiff  = Lookup[sp, relativeDiff,  Missing[]];
    (* optional, with their defaults *)
    refine        = Lookup[sp, refine,        1];
    rhoGauge      = Lookup[sp, rhoGauge,      "Point"];
    zetaSteps     = Lookup[sp, zetaSteps,     5];
    order         = Lookup[sp, order,         4];
    stepTolerance = Lookup[sp, stepTolerance, 1.*^-10];
    verbose       = TrueQ @ Lookup[sp, verbose, True];
  ];

  missing = Pick[{"\[Delta]mesh", "hMax", "box", "maxSteps", "relativeDiff"},
                 Not /@ NumericQ /@ {\[Delta]mesh, hMax, box, maxSteps,
                                     relativeDiff}];
  If[missing =!= {},
    Message[SolveActiveNematicComovingSteady::badsolver, missing];
    Return[$Failed]];

  If[!MemberQ[{"Point", "Mean", "Replace"}, rhoGauge],
    Message[SolveActiveNematicComovingSteady::badgauge, rhoGauge];
    Return[$Failed]];
  rhoGaugeV = rhoGauge;

  (* --- unpack model parameters -------------------------------------- *)
  With[{mpa = Association[modelParams]},
    B               = Lookup[mpa, B,               Missing[]];
    a               = Lookup[mpa, a,               Missing[]];
    b               = Lookup[mpa, b,               Missing[]];
    L               = Lookup[mpa, L,               Missing[]];
    \[Zeta]         = Lookup[mpa, \[Zeta],         Missing[]];
    \[Xi]0          = Lookup[mpa, \[Xi]0,          Missing[]];
    \[Xi]r          = Lookup[mpa, \[Xi]r,          Missing[]];
    \[Rho]0         = Lookup[mpa, \[Rho]0,         Missing[]];
    \[CapitalLambda]= Lookup[mpa, \[CapitalLambda], Missing[]];
  ];

  missing = Pick[{"B", "a", "b", "L", "\[Zeta]", "\[Xi]0", "\[Xi]r", "\[Rho]0",
                  "\[CapitalLambda]"},
                 Not /@ NumericQ /@ {B, a, b, L, \[Zeta], \[Xi]0, \[Xi]r,
                                     \[Rho]0, \[CapitalLambda]}];
  If[missing =!= {},
    Message[SolveActiveNematicComovingSteady::badmodel, missing];
    Return[$Failed]];

  (* comoving_core.m speaks string keys and calls the Frank constant "K" *)
  mkModel[z_] := <|"B" -> N[B], "a" -> N[a], "b" -> N[b], "K" -> N[L],
                   "zeta" -> N[z], "xi0" -> N[\[Xi]0], "xir" -> N[\[Xi]r],
                   "rho0" -> N[\[Rho]0], "Lambda" -> N[\[CapitalLambda]]|>;
  (* cdDerived returns a fresh association and does NOT carry rhoGauge across,
     so it has to be appended every time or the default silently applies. *)
  dpOf[z_] := Append[cdDerived[mkModel[z]], "rhoGauge" -> rhoGaugeV];

  dpTarget = dpOf[\[Zeta]];
  wall0 = AbsoluteTime[];

  (* --- grid, operators, boundary rows -------------------------------- *)
  (* The grid is built ONCE, at the target \[Zeta], and held fixed through the
     continuation: only the equations and the Pade far field move with \[Zeta].
     cdGridPade derives its point count from \[Delta]mesh and hMax by exactly
     the criterion iterative_solver_module.m:128-136 uses for targetCellSize,
     so those two knobs mean the same thing in both solvers. *)
  gg = cdGridPade[N[box], dpTarget, N[\[Delta]mesh], N[hMax], refine];
  op = cdOps[gg, gg, order];
  n  = op["n"]; k0 = op["k0"];
  ne = If[bcType === "Neumann" && rhoGaugeV =!= "Replace", 3, 2];

  If[verbose,
    Print["[comoving steady] bc = ", bcType, "   box = ", N[box],
          "   n = ", n, " (", op["nx"], "^2)",
          "   hMin = ", cdFmt@gg["hMin"]];
    Print["    ld = ", cdFmt@dpTarget["ld"], "   S0 = ", cdFmt@dpTarget["S0"],
          "   ld/hMin = ", cdFmt[dpTarget["ld"]/gg["hMin"]],
          "   alpha = ", cdFmt[\[Zeta] \[Xi]r dpTarget["S0"]/
                               (4 \[Xi]0 dpTarget["Kp"])],
          If[bcType === "Neumann", "   gauge = " <> rhoGaugeV, ""]]];

  (* --- \[Zeta]-continuation ------------------------------------------ *)
  (* A cold Newton from the Pade seed does NOT converge at the activities of
     interest (alpha ~ 0.43): measured, 40 wasted iterations and ~200 linear
     solves.  Ramping \[Zeta] up from zero keeps every Newton a small
     perturbation of the previous solution, ~3 iterations per step. *)
  zs = If[TrueQ[\[Zeta] == 0.], {0.},
          Prepend[Rest@Subdivide[0., N[\[Zeta]], zetaSteps], 0.]];

  dp  = dpOf[0.];
  pf  = cdPadeFuns[dp];
  bcz = cdBC[op, dp, pf, bcType];
  (* Uniform \[Rho] = 1 satisfies the reservoir wall but NOT the sealed one, and
     starting from it there makes the first Newton diverge; cdInitialRho solves
     the \[Rho] equation alone against the frozen Pade Q, one linear solve. *)
  X  = cdInitialRho[cdInitial[op, pf], op, dp, bcz, bcType];
  ex = ConstantArray[0., ne];

  Do[
    dp  = dpOf[zs[[kk]]];
    pf  = cdPadeFuns[dp];
    bcz = cdBC[op, dp, pf, bcType];
    st  = cdSteady[X, ex, op, dp, bcz, bcType,
            "MaxIterations" -> maxSteps,
            "StepTolerance" -> If[kk == Length[zs],
                                  N[stepTolerance], N[relativeDiff]],
            "Verbose" -> False];
    X = st["X"]; ex = st["ex"];
    If[verbose,
      Print["    zeta = ", StringPadRight[cdFmt[zs[[kk]]], 12],
            "  newton = ", StringPadRight[ToString@st["iterations"], 4],
            "  |R| = ", StringPadRight[cdFmt@st["residual"], 12],
            "  ux = ", cdFmt@ex[[1]]]];
    (* cdFmt the numerics: a Real reaching a Message template is formatted as a
       superscripted "m x 10^e", which the console renders across two lines and
       makes the message unreadable *)
    If[!TrueQ@st["converged"] && kk < Length[zs],
      Message[SolveActiveNematicComovingSteady::stalled,
              cdFmt@zs[[kk]], cdFmt@st["residual"]]],
    {kk, Length[zs]}];

  If[!TrueQ@st["converged"],
    Message[SolveActiveNematicComovingSteady::noconv,
            maxSteps, cdFmt@N[stepTolerance], cdFmt@st["residual"],
            cdFmt@ex[[1]]]];

  (* --- unit mean density for the sealed wall ------------------------- *)
  (* The gauges differ only by an additive constant, and u is blind to it (the
     Q equation sees \[Rho] through \[Del]\[Rho], the \[Rho] equation through
     \[Del]^2\[Rho]).  Normalising here makes the returned density comparable
     with the Dirichlet case, where the wall fixes \[Rho] = 1. *)
  If[bcType === "Neumann", X = cdShiftRhoMean[X, op]];

  If[verbose,
    (* DecimalForm, not NumberForm: NumberForm switches to a superscripted
       "m x 10^e" for small values, which OutputForm renders across TWO lines
       and shreds the log.  DecimalForm never leaves fixed point.  ToString@ is
       also load-bearing -- a bare form reaches Print as an unevaluated box. *)
    Print["    u = ", ToString@DecimalForm[ex[[1]], 12],
          "   uy = ", cdFmt@ex[[2]],
          "   core = ", cdFmt[{X[[n + k0]], X[[2 n + k0]]}],
          "   wall = ", Round[N[AbsoluteTime[] - wall0], 0.01], " s"]];

  (* --- package ------------------------------------------------------- *)
  (* cdInterp is ListInterpolation at InterpolationOrder -> 3 over the full
     tensor grid: vectorised calls and Derivative[2,0]/[0,2] both work, which
     is what the Lyapunov modules need of these three fields. *)
  <| "rho"  -> cdInterp[X, "rho", op],
     "Q1"   -> cdInterp[X, "Q1",  op],
     "Q2"   -> cdInterp[X, "Q2",  op],
     "mesh" -> <|"Coordinates" -> Transpose[{op["xs"], op["ys"]}],
                 "gx" -> op["gx"], "gy" -> op["gy"],
                 "n"  -> n, "nx" -> op["nx"], "ny" -> op["ny"],
                 "hMin" -> gg["hMin"], "hMax" -> gg["hMax"],
                 "box" -> N[box], "weights" -> op["w"], "grid" -> gg|>,
     "u"    -> ex[[1 ;; 2]],
     "solver" -> <|"converged"  -> TrueQ@st["converged"],
                   "residual"   -> st["residual"],
                   "iterations" -> st["iterations"],
                   "core"       -> {X[[n + k0]], X[[2 n + k0]]},
                   "zetaPath"   -> zs,
                   "gauge"      -> If[bcType === "Neumann", rhoGaugeV, None],
                   "mu"         -> If[Length[ex] >= 3, ex[[3]], 0.],
                   "state"      -> X,
                   "wall"       -> N[AbsoluteTime[] - wall0]|> |>
];
