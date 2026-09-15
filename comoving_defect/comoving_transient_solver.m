(* ::Package:: *)

(* ::Section:: *)
(*Comoving-frame transient (time-stepping) solver for the coupled \[Rho]\[Dash]Q active-nematic system.*)

(* The numerics live in comoving_core.m in this same directory; this file is a
   thin wrapper that puts them behind the repo's standard module interface.
   Guarded self-load, the idiom lyapunov_solver_module.m:21-22 uses. *)
If[Names["Global`cdEvolve"] === {},
   Get[FileNameJoin[{DirectoryName[$InputFileName], "comoving_core.m"}]]];

Clear[SolveActiveNematicComovingTransient];

SolveActiveNematicComovingTransient::usage =
  "SolveActiveNematicComovingTransient[modelParams, solverParams, bcType] " <>
  "integrates the coupled time-dependent \[Rho]\[Dash]Q active-nematic system IN " <>
  "THE FRAME THAT MOVES WITH A +1/2 DEFECT, with the core pinned at the origin, " <>
  "forward in blocks of length timeInterval, until the relative change of all " <>
  "three fields across a block drops below relativeDiff, or maxSteps blocks are " <>
  "reached.  All three fields and the frame velocity u are advanced together in " <>
  "one implicit-Euler solve per substep, so the coupling is fully implicit.\n\n" <>
  "  modelParams  : rules {B->..., a->..., b->..., \[Xi]0->..., \[Xi]r->..., " <>
  "\[Zeta]->..., \[Rho]0->..., L->..., \[CapitalLambda]->...}  (L is the Frank " <>
  "constant, as in transient_solver_module.m)\n" <>
  "  solverParams : rules {\[Delta]mesh->..., hMax->..., box->..., " <>
  "maxSteps->..., relativeDiff->..., timeInterval->...}\n" <>
  "                 optional: returnTransient->True (default False) -- also keep " <>
  "the whole approach to the steady state, see below;\n" <>
  "                 refine->1, rhoGauge->\"Point\", order->4, dt0->1.*^-5, " <>
  "dtMax->timeInterval, maxSubsteps->2000, verbose->True\n" <>
  "  bcType       : \"Neumann\"   \[Rule] sealed wall, " <>
  "n.(B \[Del]\[Rho] + \[Zeta] \[Del].Q) = 0, plus a single-node density gauge;\n" <>
  "                 \"Dirichlet\" \[Rule] \[Rho] = 1 on the entire outer boundary.\n" <>
  "                 Q is held at the Pade +1/2 profile on the outer boundary in " <>
  "both cases.\n\n" <>
  "Initial state: the Pade +1/2 defect core in Q, with \[Rho] = 1 (reservoir " <>
  "wall) or the quasi-static \[Rho] for the frozen Pade Q (sealed wall -- uniform " <>
  "\[Rho] does not satisfy that wall and makes the first implicit step diverge).\n\n" <>
  "Returns <|\"rho\"->..., \"Q1\"->..., \"Q2\"->..., \"mesh\"->..., \"u\"->{ux,uy}, " <>
  "\"solver\"->...|>, the fields as 2-argument interpolations of (x,y) -- drop-in " <>
  "for visualizeSteadyState and the Lyapunov modules, exactly like " <>
  "SolveActiveNematicComovingSteady.  \"solver\" carries the per-block history, " <>
  "so u(t) is available even without returnTransient.\n\n" <>
  "With returnTransient->True the same keys are returned unchanged (so the result " <>
  "stays drop-in) plus one extra key:\n" <>
  "  \"transient\" -> <|\"rho\"->f, \"Q1\"->f, \"Q2\"->f, \"u\"->g, " <>
  "\"tFinal\"->..., \"blocks\"->...|>\n" <>
  "where each f[t,x,y] takes GLOBAL time t \[Element] [0, tFinal] (scalar; x and y " <>
  "may be lists) and g[t] gives u_x, both clamped at the ends of the range -- so " <>
  "the whole approach to the steady state is one function, ready for Plot, " <>
  "DensityPlot or Animate.  NB \"blocks\" here holds the raw per-block states " <>
  "<|\"t\",\"X\",\"u\"|> on the finite-difference grid, NOT the NDSolveValue " <>
  "triples transient_solver_module.m returns; the fields between block ends are " <>
  "interpolated in t, not re-integrated.  Nothing heavy is retained when " <>
  "returnTransient is False.\n\n" <>
  "u is NOT measured from the solution.  It is the Lagrange multiplier conjugate " <>
  "to the pinning constraints Q1(0,0) = Q2(0,0) = 0, so the core sits at the " <>
  "origin exactly at every accepted step.  See dean.tex, appendix " <>
  "\"Comoving-frame defect solver\".\n\n" <>
  "NB 1: this is NOT the lab-frame transient.  The frame translates at u(t), so " <>
  "the result is not comparable with SolveActiveNematicTransient except in the " <>
  "passive limit \[Zeta] -> 0, where u -> 0.\n\n" <>
  "NB 2: \"Neumann\" is the same sealed wall transient_solver_module.m imposes, " <>
  "NOT the ray-pinned one of iterative_solver_module.m.  The residual constant " <>
  "null mode is removed by one node (rhoGauge->\"Point\", the default), and the " <>
  "density is shifted back to unit mean before it is returned.";

SolveActiveNematicComovingTransient::badbc =
  "Unknown boundary-condition type `1`; expected \"Neumann\" or \"Dirichlet\".";

SolveActiveNematicComovingTransient::badgauge =
  "Unknown rhoGauge `1`; expected \"Point\", \"Mean\" or \"Replace\".";

SolveActiveNematicComovingTransient::badtime =
  "timeInterval must be a positive number; got `1`.";

SolveActiveNematicComovingTransient::badmodel =
  "modelParams is missing or non-numeric for: `1`.  Expected rules for " <>
  "{B, a, b, L, \[Zeta], \[Xi]0, \[Xi]r, \[Rho]0, \[CapitalLambda]}.";

SolveActiveNematicComovingTransient::badsolver =
  "solverParams is missing or non-numeric for: `1`.  Expected rules for " <>
  "{\[Delta]mesh, hMax, box, maxSteps, relativeDiff, timeInterval}.";

SolveActiveNematicComovingTransient::noconv =
  "Transient solver did not reach relativeDiff = `1` within maxSteps = `2` " <>
  "blocks (final residual = `3`).";

SolveActiveNematicComovingTransient::stuck =
  "Block `1` did not reach t = `2` within maxSubsteps = `3` substeps; the time " <>
  "step has collapsed.  Stopping early at t = `4`.";


SolveActiveNematicComovingTransient[
    modelParams : {___Rule},
    solverParams : {___Rule},
    bcType_String : "Neumann"] :=
Block[
  {
    (* solver knobs (unpacked from solverParams) *)
    \[Delta]mesh, hMax, box, maxSteps, relativeDiff, timeInterval,
    returnTransient, refine, rhoGauge, order, dt0, dtMax, maxSubsteps, verbose,
    (* model parameters (unpacked from modelParams) *)
    B, a, b, L, \[Zeta], \[Xi]0, \[Xi]r, \[Rho]0, \[CapitalLambda],
    (* derived / working *)
    missing, mkModel, dp, pf, gg, op, bc, ne, n, k0, X, X0, ex, ev,
    blocks, hist, relChange, resid, iter, tGlobal, dtCur, wall0, elapsed,
    rhoGaugeV, done, prevFields, newFields, fieldsOf, mkTraj, tGrid, stack,
    nt, tFinal, uTraj
  },

  (* --- validate bcType up front ------------------------------------- *)
  If[!MemberQ[{"Neumann", "Dirichlet"}, bcType],
    Message[SolveActiveNematicComovingTransient::badbc, bcType];
    Return[$Failed]
  ];

  (* --- unpack solver parameters ------------------------------------- *)
  (* Lookup rather than the house "{k,...} = {k,...} /. solverParams" idiom:
     that form self-assigns (box = box) for an absent key, and a
     self-referential OwnValue blows the recursion limit later.  See the same
     note in comoving_steady_solver.m. *)
  With[{sp = Association[solverParams]},
    \[Delta]mesh    = Lookup[sp, \[Delta]mesh,    Missing[]];
    hMax            = Lookup[sp, hMax,            Missing[]];
    box             = Lookup[sp, box,             Missing[]];
    maxSteps        = Lookup[sp, maxSteps,        Missing[]];
    relativeDiff    = Lookup[sp, relativeDiff,    Missing[]];
    timeInterval    = Lookup[sp, timeInterval,    Missing[]];
    (* optional, with their defaults *)
    returnTransient = TrueQ @ Lookup[sp, returnTransient, False];
    refine          = Lookup[sp, refine,      1];
    rhoGauge        = Lookup[sp, rhoGauge,    "Point"];
    order           = Lookup[sp, order,       4];
    dt0             = Lookup[sp, dt0,         1.*^-5];
    dtMax           = Lookup[sp, dtMax,       Automatic];
    maxSubsteps     = Lookup[sp, maxSubsteps, 2000];
    verbose         = TrueQ @ Lookup[sp, verbose, True];
  ];

  missing = Pick[{"\[Delta]mesh", "hMax", "box", "maxSteps", "relativeDiff",
                  "timeInterval"},
                 Not /@ NumericQ /@ {\[Delta]mesh, hMax, box, maxSteps,
                                     relativeDiff, timeInterval}];
  If[missing =!= {},
    Message[SolveActiveNematicComovingTransient::badsolver, missing];
    Return[$Failed]];

  If[!TrueQ[timeInterval > 0],
    Message[SolveActiveNematicComovingTransient::badtime, timeInterval];
    Return[$Failed]];

  If[dtMax === Automatic, dtMax = N[timeInterval]];

  If[!MemberQ[{"Point", "Mean", "Replace"}, rhoGauge],
    Message[SolveActiveNematicComovingTransient::badgauge, rhoGauge];
    Return[$Failed]];
  rhoGaugeV = rhoGauge;

  (* --- unpack model parameters -------------------------------------- *)
  With[{mpa = Association[modelParams]},
    B                = Lookup[mpa, B,                Missing[]];
    a                = Lookup[mpa, a,                Missing[]];
    b                = Lookup[mpa, b,                Missing[]];
    L                = Lookup[mpa, L,                Missing[]];
    \[Zeta]          = Lookup[mpa, \[Zeta],          Missing[]];
    \[Xi]0           = Lookup[mpa, \[Xi]0,           Missing[]];
    \[Xi]r           = Lookup[mpa, \[Xi]r,           Missing[]];
    \[Rho]0          = Lookup[mpa, \[Rho]0,          Missing[]];
    \[CapitalLambda] = Lookup[mpa, \[CapitalLambda], Missing[]];
  ];

  missing = Pick[{"B", "a", "b", "L", "\[Zeta]", "\[Xi]0", "\[Xi]r", "\[Rho]0",
                  "\[CapitalLambda]"},
                 Not /@ NumericQ /@ {B, a, b, L, \[Zeta], \[Xi]0, \[Xi]r,
                                     \[Rho]0, \[CapitalLambda]}];
  If[missing =!= {},
    Message[SolveActiveNematicComovingTransient::badmodel, missing];
    Return[$Failed]];

  (* comoving_core.m speaks string keys and calls the Frank constant "K" *)
  mkModel = <|"B" -> N[B], "a" -> N[a], "b" -> N[b], "K" -> N[L],
              "zeta" -> N[\[Zeta]], "xi0" -> N[\[Xi]0], "xir" -> N[\[Xi]r],
              "rho0" -> N[\[Rho]0], "Lambda" -> N[\[CapitalLambda]]|>;
  (* cdDerived does NOT carry rhoGauge across, so it must be appended *)
  dp = Append[cdDerived[mkModel], "rhoGauge" -> rhoGaugeV];
  pf = cdPadeFuns[dp];
  wall0 = AbsoluteTime[];

  (* --- grid, operators, boundary rows -------------------------------- *)
  gg = cdGridPade[N[box], dp, N[\[Delta]mesh], N[hMax], refine];
  op = cdOps[gg, gg, order];
  n  = op["n"]; k0 = op["k0"];
  bc = cdBC[op, dp, pf, bcType];
  ne = If[bcType === "Neumann" && rhoGaugeV =!= "Replace", 3, 2];

  (* --- initial state ------------------------------------------------- *)
  (* Uniform \[Rho] = 1 satisfies the reservoir wall but NOT the sealed one:
     there \[Zeta] n.(\[Del].Q) != 0, so the first implicit step has to absorb a
     violent boundary layer and Newton diverges.  cdInitialRho holds Q frozen
     and solves the \[Rho] equation alone -- one linear solve. *)
  X0 = cdInitial[op, pf];
  If[bcType === "Neumann", X0 = cdInitialRho[X0, op, dp, bc, bcType]];
  X  = X0;
  ex = ConstantArray[0., ne];

  If[verbose,
    Print["[comoving transient] bc = ", bcType, "   box = ", N[box],
          "   n = ", n, " (", op["nx"], "^2)",
          "   hMin = ", cdFmt@gg["hMin"]];
    Print["    timeInterval = ", N[timeInterval], "   maxSteps = ", maxSteps,
          "   relativeDiff = ", N[relativeDiff],
          "   alpha = ", cdFmt[\[Zeta] \[Xi]r dp["S0"]/(4 \[Xi]0 dp["Kp"])],
          If[bcType === "Neumann", "   gauge = " <> rhoGaugeV, ""]]];

  (* relative max-norm change of one field across a block -- the same measure
     transient_solver_module.m:235-236 uses *)
  relChange[new_, old_] := Block[{scale = Max @ Abs @ new},
    If[scale > 0, Max @ Abs[new - old]/scale, Max @ Abs[new - old]]];
  fieldsOf[XX_] := {XX[[1 ;; n]], XX[[n + 1 ;; 2 n]], XX[[2 n + 1 ;; 3 n]]};

  (* --- march in blocks of timeInterval ------------------------------- *)
  blocks = If[returnTransient, {<|"t" -> 0., "X" -> X0, "u" -> {0., 0.}|>}, {}];
  hist    = {};
  iter    = 0;
  resid   = Infinity;
  tGlobal = 0.;
  dtCur   = N[dt0];
  elapsed = 0.;
  done    = False;
  prevFields = fieldsOf[X0];

  While[iter < maxSteps && resid > relativeDiff && !done,
    (* Tolerance -> 0 disables cdEvolve's own steady-state exit: the block must
       end at its boundary, and convergence is judged here, across the block,
       the way the FEM transient module judges it. *)
    ev = cdEvolve[X, ex, op, dp, bc, bcType,
           "dt0" -> dtCur, "dtMax" -> dtMax, "MaxSteps" -> maxSubsteps,
           "Tolerance" -> 0., "StopTime" -> N[timeInterval],
           "Verbose" -> False, "SnapshotEvery" -> 0];

    X = ev["X"]; ex = ev["ex"];
    (* hand the controller's next step back, or every block would restart the
       ramp from dt0 and crawl *)
    dtCur = ev["dtNext"];
    elapsed += ev["wall"];
    iter++;

    (* cdFmt the numerics: a Real reaching a Message template is formatted as a
       superscripted "m x 10^e", which the console renders across two lines and
       makes the message unreadable *)
    If[!TrueQ@ev["hitStopTime"],
      Message[SolveActiveNematicComovingTransient::stuck,
              iter, cdFmt@N[timeInterval], maxSubsteps, cdFmt@ev["t"]];
      done = True];

    tGlobal += ev["t"];
    newFields = fieldsOf[X];
    resid = Max @ MapThread[relChange, {newFields, prevFields}];
    prevFields = newFields;

    If[returnTransient,
      AppendTo[blocks, <|"t" -> tGlobal, "X" -> X, "u" -> ex[[1 ;; 2]]|>]];

    AppendTo[hist, <|"block" -> iter, "t" -> tGlobal,
                     "u" -> ex[[1 ;; 2]], "residual" -> resid,
                     "substeps" -> ev["steps"], "dtNext" -> dtCur,
                     "core" -> {X[[n + k0]], X[[2 n + k0]]},
                     "wall" -> ev["wall"]|>];

    If[verbose,
      Print["    step ", StringPadLeft[ToString[iter], 4],
            (* cdFmt / DecimalForm, never NumberForm: NumberForm switches to a
               superscripted "m x 10^e" for small values, which OutputForm
               renders across TWO lines and shreds the progress table *)
            "   t = ", StringPadRight[ToString@DecimalForm[tGlobal, 6], 12],
            "   ux = ", StringPadRight[cdFmt@ex[[1]], 12],
            "   uy = ", StringPadRight[cdFmt@ex[[2]], 12],
            "   residual = ", StringPadRight[cdFmt[resid], 12],
            "   substeps = ", StringPadRight[ToString@ev["steps"], 5],
            "   wall ", ToString@DecimalForm[N@ev["wall"], {12, 2}], " s",
            "   (total ", ToString@DecimalForm[N@elapsed, {12, 2}], " s)"]];
  ];

  If[resid > relativeDiff,
    Message[SolveActiveNematicComovingTransient::noconv,
            cdFmt@N[relativeDiff], maxSteps, cdFmt@resid]];

  (* --- unit mean density for the sealed wall ------------------------- *)
  (* The gauges differ only by an additive constant in \[Rho], to which u is
     blind.  Normalise the trajectory the same way as the final state, or the
     returned f[t,x,y] would sit at a different density level than "rho". *)
  If[bcType === "Neumann",
    X = cdShiftRhoMean[X, op];
    If[returnTransient,
      blocks = (Append[#, "X" -> cdShiftRhoMean[#["X"], op]] & /@ blocks)]];

  If[verbose,
    Print["    u = ", ToString@DecimalForm[ex[[1]], 12],
          "   uy = ", cdFmt@ex[[2]],
          "   core = ", cdFmt[{X[[n + k0]], X[[2 n + k0]]}],
          "   t = ", N[tGlobal],
          "   wall = ", Round[N[AbsoluteTime[] - wall0], 0.01], " s"]];

  (* --- global-time view of the trajectory ---------------------------- *)
  (* Like transient_solver_module.m:299-303, the returned object must be a
     SELF-CONTAINED closure -- mkTraj is a Block local and its definition dies
     with the Block -- so With inlines the interpolation and the clamp range
     into the Function body.  The Clip reproduces the clamping behaviour at
     t < 0 and t > tFinal that the FEM version has. *)
  If[returnTransient && Length[blocks] >= 2,
    tGrid  = #["t"] & /@ blocks;
    nt     = Length[tGrid];
    tFinal = Last[tGrid];
    stack[which_] := ListInterpolation[
       Transpose[cdFieldAt[#["X"], which, op]] & /@ blocks,
       {tGrid, op["gx"], op["gy"]},
       InterpolationOrder -> {Min[3, nt - 1], 3, 3}];
    mkTraj[which_] := With[{f = stack[which], t1 = tFinal},
       Function[{tt, xx, yy}, f[Clip[N[tt], {0., t1}], xx, yy]]];
    uTraj = With[{g = Interpolation[
         Transpose[{tGrid, #["u"][[1]] & /@ blocks}],
         InterpolationOrder -> Min[3, nt - 1]], t1 = tFinal},
       Function[tt, g[Clip[N[tt], {0., t1}]]]];
  ];

  (* --- package ------------------------------------------------------- *)
  Join[
    <| "rho"  -> cdInterp[X, "rho", op],
       "Q1"   -> cdInterp[X, "Q1",  op],
       "Q2"   -> cdInterp[X, "Q2",  op],
       "mesh" -> <|"Coordinates" -> Transpose[{op["xs"], op["ys"]}],
                   "gx" -> op["gx"], "gy" -> op["gy"],
                   "n"  -> n, "nx" -> op["nx"], "ny" -> op["ny"],
                   "hMin" -> gg["hMin"], "hMax" -> gg["hMax"],
                   "box" -> N[box], "weights" -> op["w"], "grid" -> gg|>,
       "u"    -> ex[[1 ;; 2]],
       "solver" -> <|"converged" -> TrueQ[resid <= relativeDiff],
                     "residual"  -> resid,
                     "blocks"    -> iter,
                     "tFinal"    -> N[tGlobal],
                     "core"      -> {X[[n + k0]], X[[2 n + k0]]},
                     "gauge"     -> If[bcType === "Neumann", rhoGaugeV, None],
                     "mu"        -> If[Length[ex] >= 3, ex[[3]], 0.],
                     "history"   -> hist,
                     "state"     -> X,
                     "wall"      -> N[AbsoluteTime[] - wall0]|> |>,
    If[returnTransient && Length[blocks] >= 2,
      <| "transient" -> <|
           "rho"    -> mkTraj["rho"],
           "Q1"     -> mkTraj["Q1"],
           "Q2"     -> mkTraj["Q2"],
           "u"      -> uTraj,
           "tFinal" -> N[tFinal],
           "blocks" -> blocks |> |>,
      <||>]]
];
