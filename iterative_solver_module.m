(* ::Package:: *)

(* ::Section:: *)
(*Iterative steady-state solver for the coupled \[Rho]\[Dash]Q active-nematic system.*)

Needs["NDSolve`FEM`"];

Clear[SolveActiveNematicSteady];

SolveActiveNematicSteady::usage =
  "SolveActiveNematicSteady[modelParams, solverParams, bcType] " <>
  "iteratively solves the coupled steady-state \[Rho]\[Dash]Q active-nematic system, " <>
  "alternating between the density and order-parameter sub-problems until the " <>
  "max-norm change of all three fields drops below relativeDiff, or maxSteps " <>
  "iterations are reached.\n\n" <>
  "  modelParams  : rules {B->..., a->..., b->..., \[Xi]0->..., \[Xi]r->..., " <>
  "\[Zeta]->..., \[Rho]0->..., L->..., \[CapitalLambda]->...}\n" <>
  "  solverParams : rules {\[Delta]mesh->..., hMax->..., box->..., " <>
  "maxSteps->..., relativeDiff->...}\n" <>
  "  bcType       : \"Neumann\"   \[Rule] \[Rho] free on the outer boundary, " <>
  "pinned to 1 on the ray x==0, y<0;\n" <>
  "                 \"Dirichlet\" \[Rule] \[Rho] = 1 on the entire outer boundary.\n\n" <>
  "Returns <|\"rho\"->..., \"Q1\"->..., \"Q2\"->..., \"mesh\"->...|>.";

SolveActiveNematicSteady::noconv =
  "Iterative solver did not reach relativeDiff = `1` within maxSteps = `2` " <>
  "iterations (final residual = `3`).";

SolveActiveNematicSteady::badbc =
  "Unknown boundary-condition type `1`; expected \"Neumann\" or \"Dirichlet\".";


SolveActiveNematicSteady[
    modelParams : {___Rule},
    solverParams : {___Rule},
    bcType_String : "Neumann"] :=
Block[
  {
    (* solver knobs (unpacked from solverParams) *)
    \[Delta]mesh, hMax, box, maxSteps, relativeDiff,
    (* coordinate + field symbols *)
    x, y, \[Rho], Q1, Q2, Q, Hp1,
    (* tensor helpers *)
    tprod, contract, ddot,
    (* Pade profile *)
    S0PadeOverR, Q0Pade, Q1Pade, Q2Pade, dS0, targetCellSize,
    (* PDE + BC bundle *)
    eq\[Rho], eqQ, \[Rho]System,
    (* iteration state *)
    mesh, \[Rho]Sol, Q1Sol, Q2Sol,
    \[Rho]Step, QStep, iter, residual,
    (* renormalized parameters *)
    a1, L1
  },

  (* --- validate bcType up front ------------------------------------- *)
  If[!MemberQ[{"Neumann", "Dirichlet"}, bcType],
    Message[SolveActiveNematicSteady::badbc, bcType];
    Return[$Failed]
  ];

  (* --- unpack solver parameters ------------------------------------- *)
  {\[Delta]mesh, hMax, box, maxSteps, relativeDiff} =
    {\[Delta]mesh, hMax, box, maxSteps, relativeDiff} /. solverParams;

  (* --- tensor-calculus helpers -------------------------------------- *)
  tprod[A_, B_] := Outer[Times, A, B];
  contract[A_, B_, ind_] :=
    TensorContract[tprod[A, B],
      Map[{#[[1]], ArrayDepth[A] + #[[2]]} &, ind]];
  ddot[A_, B_] := contract[A, B, {{1, 1}, {2, 2}}];

  (* --- symmetric traceless order parameter Q(x,y) ------------------- *)
  Q[xx_, yy_] := {{Q1[xx, yy],  Q2[xx, yy]},
                  {Q2[xx, yy], -Q1[xx, yy]}};

  (* --- density equation --------------------------------------------- *)
  eq\[Rho] =
    Inactive[Div][
      {{B/\[Xi]0, 0}, {0, B/\[Xi]0}} .
        Inactive[Grad][\[Rho][x, y], {x, y}],
      {x, y}] +
    Activate@Inactive[Div][
      {{\[Zeta]/\[Xi]0, 0}, {0, -(\[Zeta]/\[Xi]0)}} .
        Inactive[Grad][Q1[x, y], {x, y}],
      {x, y}] +
    Activate@Inactive[Div][
      {{0, \[Zeta]/\[Xi]0}, {\[Zeta]/\[Xi]0, 0}} .
        Inactive[Grad][Q2[x, y], {x, y}],
      {x, y}];

  (* --- local (relaxation) part of the Q potential ------------------- *)
  Hp1[xx_, yy_] := Block[{QM = Q[xx, yy]},
    -a QM - b ddot[QM, QM] QM];

  (* --- Q equation (per-component tensor equation) ------------------- *)
  eqQ =
    Table[
      Inactive[Div][
        {{B/\[Xi]0 Q[x, y][[i, j]], 0},
         {0, B/\[Xi]0 Q[x, y][[i, j]]}} .
          Activate@Inactive[Grad][\[Rho][x, y], {x, y}],
        {x, y}] +
      Inactive[Div][
        {{(4 L)/\[Xi]r + \[Zeta]/\[Xi]0, 0},
         {0, (4 L)/\[Xi]r + \[Zeta]/\[Xi]0}} .
          Inactive[Grad][Q[x, y][[i, j]], {x, y}],
        {x, y}] -
      Activate@Inactive[Grad][B/\[Xi]0 \[Rho][x, y], {x, y}] .
        Inactive[Grad][Q[x, y][[i, j]], {x, y}],
      {i, 2}, {j, 2}
    ] + (4/\[Xi]r) Hp1[x, y] - 4 \[CapitalLambda] Q[x, y];

  (* --- Pade profile: sets initial Q and its Dirichlet outer BC ------ *)
  {a1, L1} = {a+\[CapitalLambda] \[Xi]r, L+(\[Zeta] \[Xi]r)/(4 \[Xi]0) };
  Print["a1 = ", a1/.modelParams, ", L1 = ", L1/.modelParams];
  S0PadeOverR[r2_] := Sqrt[-a1/(2 b)] Sqrt[-a1/L1] Sqrt[
      (0.34 + 0.07 (-r2 a1/L1)) /
      (1 + 0.41 (-r2 a1/L1) + 0.07 (-r2 a1/L1)^2)
    ];
  Q0Pade[xv_, yv_] := S0PadeOverR[xv^2 + yv^2]  {{xv, yv}, {yv, -xv}};
  Q1Pade[xx_, yy_] = Q0Pade[xx, yy][[1, 1]];
  Q2Pade[xx_, yy_] = Q0Pade[xx, yy][[1, 2]];

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

  (* --- boundary-condition bundle for the \[Rho] sub-problem --------- *)
  \[Rho]System = Switch[bcType,
    "Neumann",
      {eq\[Rho] == NeumannValue[0, True],
       DirichletCondition[\[Rho][x, y] == 1, x == 0 && y < 0]},
    "Dirichlet",
      {eq\[Rho] == 0,
       DirichletCondition[\[Rho][x, y] == 1, True]}
  ];

  (* --- density sub-step: freezes Q, solves for \[Rho] --------------- *)
  \[Rho]Step := Module[{\[Rho]SolNew, \[Delta]},
    \[Rho]SolNew = NDSolveValue[
      \[Rho]System /. {Q1 -> Q1Sol, Q2 -> Q2Sol} /. modelParams,
      \[Rho],
      {x, y} \[Element] mesh,
      InitialSeeding -> {\[Rho][x, y] == \[Rho]Sol[x, y]},
      Method -> {"FiniteElement"}
    ];
    \[Delta] = If[iter > 0,
      Max@Abs[\[Rho]Sol["ValuesOnGrid"] - \[Rho]SolNew["ValuesOnGrid"]],
      Infinity];
    \[Rho]Sol = \[Rho]SolNew;
    \[Delta]
  ];

  (* --- order-parameter sub-step: freezes \[Rho], solves for Q1, Q2 -- *)
  QStep := Module[{Q1SolNew, Q2SolNew, \[Delta]1, \[Delta]2},
    {Q1SolNew, Q2SolNew} = NDSolveValue[
      {eqQ[[1, 1]] == 0, eqQ[[1, 2]] == 0,
       DirichletCondition[Q1[x, y] == Q1Pade[x, y], True],
       DirichletCondition[Q2[x, y] == Q2Pade[x, y], True]} /.
        {\[Rho] -> \[Rho]Sol} /. modelParams,
      {Q1, Q2},
      {x, y} \[Element] mesh,
      InitialSeeding -> {Q1[x, y] == Q1Sol[x, y],
                         Q2[x, y] == Q2Sol[x, y]},
      Method -> {"FiniteElement"}
    ];
    \[Delta]1 = If[iter > 0,
      Max@Abs[Q1Sol["ValuesOnGrid"] - Q1SolNew["ValuesOnGrid"]],
      Infinity];
    \[Delta]2 = If[iter > 0,
      Max@Abs[Q2Sol["ValuesOnGrid"] - Q2SolNew["ValuesOnGrid"]],
      Infinity];
    Q1Sol = Q1SolNew;
    Q2Sol = Q2SolNew;
    Max[\[Delta]1, \[Delta]2]
  ];

  (* --- initial guess: uniform \[Rho], Pade Q ------------------------ *)
  \[Rho]Sol[xx_, yy_] = 1;
  Q1Sol[xx_, yy_]     = Q1Pade[xx, yy] /. modelParams;
  Q2Sol[xx_, yy_]     = Q2Pade[xx, yy] /. modelParams;

  (* --- iterate until convergence or maxSteps ------------------------ *)
  iter = 0;
  residual = Infinity;
  While[iter < maxSteps && residual > relativeDiff,
    residual = Max[\[Rho]Step, QStep];
    iter++
  ];
  If[residual > relativeDiff,
    Message[SolveActiveNematicSteady::noconv,
            relativeDiff, maxSteps, residual]];

  <| "rho" -> \[Rho]Sol,
     "Q1"  -> Q1Sol,
     "Q2"  -> Q2Sol,
     "mesh"-> mesh |>
];
