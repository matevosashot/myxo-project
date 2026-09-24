(* ::Package:: *)

(* ::Section:: *)
(*Transient (time-stepping) solver for the coupled \[Rho]\[Dash]Q active-nematic system.*)

Needs["NDSolve`FEM`"];

Clear[SolveActiveNematicTransient];

SolveActiveNematicTransient::usage =
  "SolveActiveNematicTransient[modelParams, solverParams, bcType] " <>
  "integrates the coupled time-dependent \[Rho]\[Dash]Q active-nematic system forward " <>
  "in blocks of length timeInterval, restarting each block from the previous " <>
  "block's final fields, until the relative change of all three fields across a " <>
  "block drops below relativeDiff, or maxSteps blocks are reached.  All three " <>
  "fields are advanced together in ONE NDSolveValue per block, so the coupling " <>
  "is fully implicit.\n\n" <>
  "  modelParams  : rules {B->..., a->..., b->..., \[Xi]0->..., \[Xi]r->..., " <>
  "\[Zeta]->..., \[Rho]0->..., L->..., \[CapitalLambda]->...}\n" <>
  "  solverParams : rules {\[Delta]mesh->..., hMax->..., box->..., " <>
  "maxSteps->..., relativeDiff->..., timeInterval->...}\n" <>
  "                 optional: returnTransient->True (default False) -- also " <>
  "keep the whole approach to the steady state, see below\n" <>
  "  bcType       : \"Neumann\"   \[Rule] sealed wall, " <>
  "n.(B \[Del]\[Rho] + \[Zeta] \[Del].Q) = 0, \[Rho] free (no pin -- the initial " <>
  "condition sets the level);\n" <>
  "                 \"Dirichlet\" \[Rule] \[Rho] = 1 on the entire outer boundary.\n" <>
  "                 Q is held at the Pade profile on the outer boundary in both " <>
  "cases.\n\n" <>
  "Initial state: uniform \[Rho] = 1 with the Pade defect core in Q.\n\n" <>
  "Returns <|\"rho\"->..., \"Q1\"->..., \"Q2\"->..., \"mesh\"->...|>, the fields as " <>
  "2-argument interpolations of (x,y) -- drop-in for visualizeSteadyState and the " <>
  "Lyapunov modules, exactly like SolveActiveNematicSteady.\n\n" <>
  "With returnTransient->True the same four keys are returned unchanged (so the " <>
  "result stays drop-in) plus one extra key:\n" <>
  "  \"transient\" -> <|\"rho\"->f, \"Q1\"->f, \"Q2\"->f, \"tFinal\"->..., " <>
  "\"blocks\"->...|>\n" <>
  "where each f[t,x,y] takes GLOBAL time t \[Element] [0, tFinal] (scalar; x and y " <>
  "may be lists) and dispatches internally to the block that covers it, so the " <>
  "whole approach to the steady state is one function -- ready for Plot, " <>
  "DensityPlot or Animate.  \"blocks\" holds the raw per-block NDSolveValue " <>
  "triples {\[Rho],Q1,Q2} in LOCAL time, block k covering global t in " <>
  "[(k-1) timeInterval, k timeInterval].  Nothing is retained when " <>
  "returnTransient is False, so the default memory cost is unchanged.\n\n" <>
  "NB the \"Neumann\" wall is NOT the one SolveActiveNematicSteady imposes.  That " <>
  "solver Activates the \[Zeta] terms, demoting them to a source, so its natural " <>
  "condition is n.\[Del]\[Rho] = 0.  Here Q1 and Q2 are unknowns, so the \[Zeta] " <>
  "terms must stay inside Inactive[Div] and they therefore carry their own " <>
  "boundary flux.  The two Neumann steady states differ; the Dirichlet ones agree.";

SolveActiveNematicTransient::noconv =
  "Transient solver did not reach relativeDiff = `1` within maxSteps = `2` " <>
  "blocks (final residual = `3`).";

SolveActiveNematicTransient::badbc =
  "Unknown boundary-condition type `1`; expected \"Neumann\" or \"Dirichlet\".";

SolveActiveNematicTransient::badtime =
  "timeInterval must be a positive number; got `1`.";

SolveActiveNematicTransient::nosol =
  "NDSolveValue failed to return three fields at block `1`.";


SolveActiveNematicTransient[
    modelParams : {___Rule},
    solverParams : {___Rule},
    bcType_String : "Neumann"] :=
Block[
  {
    (* solver knobs (unpacked from solverParams) *)
    \[Delta]mesh, hMax, box, maxSteps, relativeDiff, timeInterval,
    returnTransient, blocks, mkTraj,
    (* coordinate + field symbols *)
    t, x, y, \[Rho], Q1, Q2, Q, Hp1,
    (* tensor helpers *)
    tprod, contract, ddot, dm,
    (* Pade profile *)
    S0PadeOverR, Q0Pade, Q1Pade, Q2Pade, q1p, q2p, dS0, targetCellSize,
    (* PDE + BC bundle *)
    eq\[Rho], eqQ, qPenalty, system,
    (* mesh + nodal bookkeeping *)
    mesh, coords, xs, ys, tvec, nNode,
    \[Rho]Vals, Q1Vals, Q2Vals, \[Rho]Sol, Q1Sol, Q2Sol, newVals,
    relChange,
    (* iteration state *)
    iter, residual, elapsed, tm, sol,
    (* renormalized parameters *)
    a1, L1, b1
  },

  (* --- validate bcType up front ------------------------------------- *)
  If[!MemberQ[{"Neumann", "Dirichlet"}, bcType],
    Message[SolveActiveNematicTransient::badbc, bcType];
    Return[$Failed]
  ];

  (* --- unpack solver parameters ------------------------------------- *)
  {\[Delta]mesh, hMax, box, maxSteps, relativeDiff, timeInterval} =
    {\[Delta]mesh, hMax, box, maxSteps, relativeDiff, timeInterval} /. solverParams;

  If[!(NumericQ[timeInterval] && TrueQ[timeInterval > 0]),
    Message[SolveActiveNematicTransient::badtime, timeInterval];
    Return[$Failed]
  ];

  (* Optional, defaults to False: an absent key leaves the symbol unreplaced and
     TrueQ turns that into False. *)
  returnTransient = TrueQ[returnTransient /. solverParams];
  blocks = {};

  (* --- tensor-calculus helpers -------------------------------------- *)
  (* dm folds a scalar into an explicit 2x2 matrix.  Load-bearing: Dot binds
     TIGHTER than Times, so `-c IdentityMatrix[2] . Grad[u]` parses as a scalar
     times a vector and the FEM parser rejects it with "Inconsistent equation
     dimensions".  Every coefficient below is an explicit matrix with its sign
     already folded in. *)
  dm[s_] := {{s, 0}, {0, s}};
  tprod[A_, B_] := Outer[Times, A, B];
  contract[A_, B_, ind_] :=
    TensorContract[tprod[A, B],
      Map[{#[[1]], ArrayDepth[A] + #[[2]]} &, ind]];
  ddot[A_, B_] := contract[A, B, {{1, 1}, {2, 2}}];

  (* --- symmetric traceless order parameter Q(t,x,y) ----------------- *)
  Q[tt_, xx_, yy_] := {{Q1[tt, xx, yy],  Q2[tt, xx, yy]},
                       {Q2[tt, xx, yy], -Q1[tt, xx, yy]}};

  (* --- density equation --------------------------------------------- *)
  (* Canonical time-dependent FEM form: \[PartialD]_t u + Div[-c.Grad[u]] == f.
     Note the sign is OPPOSITE to iterative_solver_module.m, which writes
     Div[+c.Grad[\[Rho]]] -- harmless for a stationary solve, backward-heat once
     \[PartialD]_t is present.  The \[Zeta] terms stay inactive (Q1, Q2 are unknowns
     here, so a strong-form \[PartialD]^2 Q is not representable), which is why the
     Neumann wall below is the full sealed one. *)
  eq\[Rho] =
    D[\[Rho][t, x, y], t]
    + Inactive[Div][dm[-(B/\[Xi]0)] .
        Inactive[Grad][\[Rho][t, x, y], {x, y}], {x, y}]
    + Inactive[Div][{{-(\[Zeta]/\[Xi]0), 0}, {0, \[Zeta]/\[Xi]0}} .
        Inactive[Grad][Q1[t, x, y], {x, y}], {x, y}]
    + Inactive[Div][{{0, -(\[Zeta]/\[Xi]0)}, {-(\[Zeta]/\[Xi]0), 0}} .
        Inactive[Grad][Q2[t, x, y], {x, y}], {x, y}];

  (* --- local (relaxation) part of the Q potential ------------------- *)
  Hp1[tt_, xx_, yy_] := Block[{QM = Q[tt, xx, yy]},
    -a QM - b ddot[QM, QM] QM];

  (* --- Q equation (per-component tensor equation) ------------------- *)
  (* The first two terms are Div[(B/\[Xi]0) Q \[Del]\[Rho]] - (B/\[Xi]0) \[Del]\[Rho].\[Del]Q
     = (B/\[Xi]0) Q \[Del]^2 \[Rho], written so that no second derivative of an
     unknown ever appears outside a Div.  The cross term MUST be spelled with
     explicit D[...]; the Inactive[Grad].Inactive[Grad] spelling does not parse. *)
  eqQ =
    Table[
      D[Q[t, x, y][[i, j]], t]
      + Inactive[Div][dm[-(B/\[Xi]0) Q[t, x, y][[i, j]]] .
          Inactive[Grad][\[Rho][t, x, y], {x, y}], {x, y}]
      + (B/\[Xi]0) (D[\[Rho][t, x, y], x] D[Q[t, x, y][[i, j]], x]
                  + D[\[Rho][t, x, y], y] D[Q[t, x, y][[i, j]], y])
      + Inactive[Div][dm[-((4 L)/\[Xi]r + \[Zeta]/\[Xi]0)] .
          Inactive[Grad][Q[t, x, y][[i, j]], {x, y}], {x, y}],
      {i, 2}, {j, 2}
    ] - (4/\[Xi]r) Hp1[t, x, y] + 4 \[CapitalLambda] Q[t, x, y];

  (* --- Pade profile: sets the initial Q and its outer BC ------------ *)
  {a1, L1, b1} = {a + \[CapitalLambda] \[Xi]r,
                  L + (\[Zeta] \[Xi]r)/(4 \[Xi]0), b};

  S0PadeOverR[r2_] := Sqrt[-a1/(2 b1)] Sqrt[-a1/L1] Sqrt[
      (0.34 + 0.07 (-r2 a1/L1)) /
      (1 + 0.41 (-r2 a1/L1) + 0.07 (-r2 a1/L1)^2)
    ];
  Q0Pade[xv_, yv_] := S0PadeOverR[xv^2 + yv^2]  {{xv, yv}, {yv, -xv}};
  Q1Pade[xx_, yy_] = Q0Pade[xx, yy][[1, 1]];
  Q2Pade[xx_, yy_] = Q0Pade[xx, yy][[1, 2]];
  (* numeric, vectorizable versions for nodal work and for the BCs *)
  q1p[xx_, yy_] = Q1Pade[xx, yy] /. modelParams;
  q2p[xx_, yy_] = Q2Pade[xx, yy] /. modelParams;

  (* --- mesh: refine to resolve the Pade-derived core scale ---------- *)
  dS0[r_] = D[r S0PadeOverR[r^2], r];
  targetCellSize[xv_, yv_] := Block[{r = Sqrt[xv^2 + yv^2]/2},
    Min[\[Delta]mesh/dS0[r], hMax] /. modelParams];
  mesh = ToElementMesh[
    Rectangle[{-box, -box}, {box, box}],
    "MeshOrder"       -> 2,
    "IncludePoints"   -> {{0., 0.}},
    MeshRefinementFunction -> Function[{vertices, area},
      area > (targetCellSize @@ Mean[vertices])^2]
  ];
  Print["Mesh size: ", First@Dimensions[mesh["Coordinates"]]];

  (* --- boundary conditions ------------------------------------------ *)
  (* Mathematica 13.3.1 cannot time-step a FEM system of THREE OR MORE
     dependent variables in which only SOME carry a DirichletCondition: the
     assembly hands LinearSolve a singular matrix ("Zero pivot was detected").
     Verified in claude_experiments/probe_neumann{2,3,4,5}.wls -- it is neither
     the mesh, the coefficients, the coupling nor the constant null mode (a
     whole pinned edge does not help), and two variables mixed the same way are
     fine.  All-Dirichlet works and all-Neumann works.
       "Dirichlet" therefore uses native DirichletCondition on all three fields.
       "Neumann" leaves \[Rho] free, so Q's boundary data is imposed WEAKLY, as a
     Robin/penalty term: with n.(c \[Del]Q) = -qPenalty (Q - Q_Pade) the boundary
     error is n.(c \[Del]Q)/qPenalty.  At 10^6 c that is ~10^-7, well under the
     ~2*10^-5 floor set by the quadratic elements themselves. *)
  qPenalty = 1.*^6 ((4 L)/\[Xi]r + \[Zeta]/\[Xi]0) /. modelParams;

  system[ic_] := Join[
    Switch[bcType,
      "Dirichlet",
        {eq\[Rho]      == 0,
         eqQ[[1, 1]]   == 0,
         eqQ[[1, 2]]   == 0,
         DirichletCondition[\[Rho][t, x, y] == 1,          True],
         DirichletCondition[Q1[t, x, y]     == q1p[x, y],  True],
         DirichletCondition[Q2[t, x, y]     == q2p[x, y],  True]},
      "Neumann",
        {eq\[Rho]      == NeumannValue[0, True],
         eqQ[[1, 1]]   == NeumannValue[-qPenalty (Q1[t, x, y] - q1p[x, y]), True],
         eqQ[[1, 2]]   == NeumannValue[-qPenalty (Q2[t, x, y] - q2p[x, y]), True]}
    ] /. modelParams,
    ic];

  (* --- nodal bookkeeping -------------------------------------------- *)
  (* N@ is load-bearing for the batched InterpolatingFunction calls; exact
     rationals bypass the fast path. *)
  coords = N @ mesh["Coordinates"];
  xs     = coords[[All, 1]];
  ys     = coords[[All, 2]];
  nNode  = Length[xs];
  tvec   = ConstantArray[N @ timeInterval, nNode];

  (* relative max-norm change of one field across a block *)
  relChange[new_, old_] := Block[{scale = Max @ Abs @ new},
    If[scale > 0, Max @ Abs[new - old]/scale, Max @ Abs[new - old]]];

  (* --- initial state: uniform density, Pade defect core ------------- *)
  \[Rho]Vals = ConstantArray[1., nNode];
  Q1Vals     = q1p[xs, ys];
  Q2Vals     = q2p[xs, ys];
  {\[Rho]Sol, Q1Sol, Q2Sol} =
    ElementMeshInterpolation[{mesh}, #] & /@ {\[Rho]Vals, Q1Vals, Q2Vals};

  (* --- march in blocks of timeInterval ------------------------------ *)
  Print["[T] Transient solve: bcType = ", bcType,
        ", timeInterval = ", N[timeInterval],
        ", maxSteps = ", maxSteps,
        ", relativeDiff = ", N[relativeDiff]];

  iter     = 0;
  residual = Infinity;
  elapsed  = 0.;

  While[iter < maxSteps && residual > relativeDiff,
    {tm, sol} = AbsoluteTiming @ NDSolveValue[
      system[{\[Rho][0, x, y] == \[Rho]Sol[x, y],
              Q1[0, x, y]     == Q1Sol[x, y],
              Q2[0, x, y]     == Q2Sol[x, y]}],
      {\[Rho], Q1, Q2},
      {x, y} \[Element] mesh,
      {t, 0, timeInterval},
      Method -> {"MethodOfLines",
                 "SpatialDiscretization" -> {"FiniteElement"}}
    ];

    If[!(ListQ[sol] && Length[sol] == 3),
      Message[SolveActiveNematicTransient::nosol, iter + 1];
      Return[$Failed]
    ];

    If[returnTransient, AppendTo[blocks, sol]];

    newVals  = Through[sol[tvec, xs, ys]];
    residual = Max @ MapThread[relChange,
      {newVals, {\[Rho]Vals, Q1Vals, Q2Vals}}];

    {\[Rho]Vals, Q1Vals, Q2Vals} = newVals;
    {\[Rho]Sol, Q1Sol, Q2Sol} =
      ElementMeshInterpolation[{mesh}, #] & /@ newVals;

    iter++;
    elapsed += tm;
    Print["    step ", iter,
          "   t = ", N[iter timeInterval],
          "   wall ", Round[tm, 0.01], " s",
          "   (total ", Round[elapsed, 0.01], " s)",
          "   residual = ", residual];
  ];

  If[residual > relativeDiff,
    Message[SolveActiveNematicTransient::noconv,
            relativeDiff, maxSteps, residual]];

  (* --- global-time view of the trajectory --------------------------- *)
  (* The returned object must be a SELF-CONTAINED closure: mkTraj itself is a
     Block local and its definition dies with the Block, so With inlines the
     block list, the block length and the count into the Function body. *)
  mkTraj[idx_] := With[
    {bl = blocks[[All, idx]], dt = N[timeInterval], nb = Length[blocks]},
    Function[{tt, xx, yy},
      With[{k = Clip[Ceiling[N[tt]/dt], {1, nb}]},
        bl[[k]][Clip[N[tt] - (k - 1) dt, {0., dt}], xx, yy]]]];

  Join[
    <| "rho" -> \[Rho]Sol,
       "Q1"  -> Q1Sol,
       "Q2"  -> Q2Sol,
       "mesh"-> mesh |>,
    If[returnTransient && blocks =!= {},
      <| "transient" -> <|
           "rho"    -> mkTraj[1],
           "Q1"     -> mkTraj[2],
           "Q2"     -> mkTraj[3],
           "tFinal" -> N[iter timeInterval],
           "blocks" -> blocks |> |>,
      <||>]]
];
