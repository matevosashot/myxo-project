(* ::Package:: *)

(* ::Section:: *)
(*Comoving-frame steady-state solver with a free Q wall.*)

(* Same problem and interface as comoving_steady_solver.m, but Q is not held
   at the Pade profile on the outer wall.  The numerics are comoving_core.m,
   used unchanged: only the boundary association handed to cdRJ / cdRes and
   one symmetrised operator are new.  Guarded self-load, as in
   comoving_steady_solver.m. *)
If[Names["Global`cdSteady"] === {},
   Get[FileNameJoin[{DirectoryName[$InputFileName], "comoving_core.m"}]]];

ClearAll[cdBCQfree, cdLineMask, cdMirrorOps, cdF2, cdLineForce,
         SolveActiveNematicComovingSteadyQfree];

(* -------------------------------------------------------------------- *)
(* What the Pade wall did, and what replaces it.                         *)
(*                                                                      *)
(* The Pade wall fixes three things: the x and y position of the defect *)
(* and its orientation.  With a free wall the position is still held by *)
(* the core pins Q1(0,0) = Q2(0,0) = 0, with u = (ux, uy) as their      *)
(* multipliers.  The orientation is held by the mirror line             *)
(*                                                                      *)
(*    Q2(x, 0) = 0   on every node of the row y = 0 except the core.     *)
(*                                                                      *)
(* The line is required.  At zeta = 0 the equations and the Neumann     *)
(* rows are invariant under a global rotation of (Q1, Q2), so without   *)
(* it the Jacobian is exactly singular and the first Newton step fails. *)
(* At zeta > 0 the rotation mode stays marginal (measured: slowest      *)
(* eigenvalue -2.5e-6 at box 10, zeta = 2.5 xi0), and Newton does not   *)
(* converge.                                                            *)
(*                                                                      *)
(* The mirror y -> -y (rho, Q1 even, Q2 odd, uy = 0) makes the Q2 PDE    *)
(* row vanish on y = 0, so the line replaces rows that a symmetric      *)
(* state satisfies anyway.  The core node is left out of the line: its  *)
(* Q2 PDE row then fixes uy, and uy = 0 checks the symmetry.  The       *)
(* returned "lineForce" is the Q2 PDE residual on the line nodes; it is *)
(* round-off when the line exerts no force.                             *)
(* -------------------------------------------------------------------- *)

(* 1 on the nodes of the row y = 0, except the core node k0 *)
cdLineMask[op_Association] := Module[{j0 = (op["ny"] + 1)/2},
  ReplacePart[Boole[# == j0] & /@ op["jj"], op["k0"] -> 0]];

(* NDSolve`FiniteDifferenceDerivative's 4th-order second derivative on a
   non-uniform grid uses a 6-point stencil, which cannot be centred, so Ly
   is NOT mirror-symmetric under y -> -y (measured: 2.4% of its largest
   entry at box 10).  Dx, Dy and Dxy are.  This is why comoving_steady_
   solver.m returns uy ~ 1e-5 rather than 0.  Averaging Ly with its mirror
   image keeps 4th order and makes the discrete problem exactly symmetric,
   so the line is then inactive to round-off.  Effect on u: 4e-7. *)
cdMirrorOps[op_Association] := Module[{n = op["n"], nx = op["nx"], perm, P, Ly},
  perm = Flatten@Reverse@Partition[Range[n], nx];
  P  = SparseArray[Thread[Transpose[{Range[n], perm}] -> 1.], {n, n}];
  Ly = (op["Ly"] + P . op["Ly"] . P)/2;
  Join[op, <|"Ly" -> Ly, "Lap" -> op["Lx"] + Ly|>]
];

(* -------------------------------------------------------------------- *)
(* Boundary rows.  The rho rows, the rho gauge and "imD" for rho come    *)
(* from cdBC unchanged.  The Q rows are rebuilt:                         *)
(*   qWall "Neumann" : n.grad Q1 = n.grad Q2 = 0, signed outward normal, *)
(*                     summed at corners as cdBC does for rho            *)
(*   qWall "Pade"    : Q = Pade, as in comoving_steady_solver.m          *)
(* With symLine the line rows REPLACE whatever row the node had, PDE or  *)
(* wall (the two wall nodes at y = 0).                                   *)
(* -------------------------------------------------------------------- *)
cdBCQfree[op_Association, dp_Association, pf_Association, bcType_String,
          qWall_String, symLine_] :=
Module[{n = op["n"], bm = op["bmask"], im = op["imask"], base, Z, Nq, Iq,
        line, Jq1, Jq2, rq1, rq2, Jbc, rhs, imQ2},
  base = cdBC[op, dp, pf, bcType];
  If[base === $Failed, Return[$Failed]];
  Z = SparseArray[{}, {n, n}];

  Switch[qWall,
    "Pade",
      Iq  = DiagonalMatrix[SparseArray[bm]];
      Jq1 = Iq; Jq2 = Iq;
      rq1 = bm pf["q1"][op["xs"], op["ys"]];
      rq2 = bm pf["q2"][op["xs"], op["ys"]];,
    "Neumann",
      Nq  = op["px"] op["Dx"] + op["py"] op["Dy"];
      Jq1 = Nq; Jq2 = Nq;
      rq1 = ConstantArray[0., n]; rq2 = ConstantArray[0., n];,
    _,
      Print["cdBCQfree: unknown qWall ", qWall]; Return[$Failed]];

  imQ2 = im;
  If[TrueQ[symLine],
    line = cdLineMask[op];
    Jq2  = (1 - line) Jq2 + DiagonalMatrix[SparseArray[line]];
    rq2  = (1 - line) rq2;
    imQ2 = (1 - line) im];

  Jbc = ArrayFlatten[{{base["Jbc"][[1 ;; n, 1 ;; n]],
                       base["Jbc"][[1 ;; n, n + 1 ;; 2 n]],
                       base["Jbc"][[1 ;; n, 2 n + 1 ;; 3 n]]},
                      {Z, Jq1, Z}, {Z, Z, Jq2}}];
  rhs = Join[base["rhs"][[1 ;; n]], rq1, rq2];

  <|"Jbc" -> SparseArray[Jbc], "rhs" -> rhs,
    "gaugeNode" -> base["gaugeNode"],
    "imD" -> Join[base["imD"][[1 ;; n]], im, imQ2],
    "qWall" -> qWall, "symLine" -> TrueQ[symLine]|>
];

(* the Q2 right-hand side F2 of comoving_core.m, on every node *)
cdF2[X_, ex_, op_Association, dp_Association] :=
Module[{n = op["n"], r, q1, q2, wx, wy, Bx, kap, Aq, bet},
  r = X[[1 ;; n]]; q1 = X[[n + 1 ;; 2 n]]; q2 = X[[2 n + 1 ;; 3 n]];
  {Bx, kap, Aq, bet} = Lookup[dp, {"Bx", "kap", "Aq", "bet"}];
  wx = ex[[1]] + Bx (op["Dx"] . r); wy = ex[[2]] + Bx (op["Dy"] . r);
  op["Dx"] . (q2 wx) + op["Dy"] . (q2 wy) + kap (op["Lap"] . q2)
    - Aq q2 - bet (q1^2 + q2^2) q2
];

(* Q2 PDE residual on the interior line nodes, relative to the largest
   diffusion term kap Lap Q2 anywhere.  Zero for a mirror-symmetric state. *)
cdLineForce[X_, ex_, op_Association, dp_Association] :=
Module[{n = op["n"], sel, scale},
  sel = Pick[Range[n], cdLineMask[op] op["imask"], 1];
  scale = Max@Abs[dp["kap"] (op["Lap"] . X[[2 n + 1 ;; 3 n]])];
  Max@Abs[cdF2[X, ex, op, dp][[sel]]]/scale
];


SolveActiveNematicComovingSteadyQfree::usage =
  "SolveActiveNematicComovingSteadyQfree[modelParams, solverParams, bcType] " <>
  "solves the same comoving steady problem as SolveActiveNematicComovingSteady, " <>
  "with the same arguments and return value, but with a free Q wall.\n\n" <>
  "  extra solverParams (all optional):\n" <>
  "    qWall   -> \"Neumann\" (default): n.\[Del]Q1 = n.\[Del]Q2 = 0 on the outer " <>
  "wall;  \"Pade\": Q held at the Pade profile, as in the original solver\n" <>
  "    symLine -> True (default): Q2(x,0) = 0 on the row y = 0, which fixes the " <>
  "orientation of the defect.  Required for qWall -> \"Neumann\".\n" <>
  "    symOps  -> True (default): mirror-symmetrise the y second derivative so " <>
  "the discrete problem is exactly symmetric under y -> -y.\n" <>
  "  bcType : the density wall, \"Neumann\" (sealed), \"Dirichlet\" or " <>
  "\"RhoNeumann\", exactly as in SolveActiveNematicComovingSteady.\n\n" <>
  "\"solver\" additionally carries \"qWall\", \"symLine\", \"symOps\", " <>
  "\"lineForce\" (Q2 PDE residual on the line, relative; round-off when the " <>
  "line exerts no force) and \"symmetry\" (cdSymmetry of the state).\n\n" <>
  "qWall -> \"Pade\", symLine -> True, symOps -> False reproduces " <>
  "SolveActiveNematicComovingSteady to 1e-9 in u.";

SolveActiveNematicComovingSteadyQfree::badbc =
  "Unknown boundary-condition type `1`; expected \"Neumann\", \"Dirichlet\" " <>
  "or \"RhoNeumann\".";
SolveActiveNematicComovingSteadyQfree::badqwall =
  "Unknown qWall `1`; expected \"Neumann\" or \"Pade\".";
SolveActiveNematicComovingSteadyQfree::noline =
  "qWall -> \"Neumann\" with symLine -> False: the orientation of the defect is " <>
  "not fixed.  The Jacobian is singular at \[Zeta] = 0 and Newton is not " <>
  "expected to converge.";
SolveActiveNematicComovingSteadyQfree::gaugeignored =
  "rhoGauge -> `1` is ignored for \"RhoNeumann\", which always pins " <>
  "\[Rho](0,0) = 1.";
SolveActiveNematicComovingSteadyQfree::badreuse =
  "\"ReuseFactorization\" -> `1`; expected True, False or Automatic.";
SolveActiveNematicComovingSteadyQfree::badgauge =
  "Unknown rhoGauge `1`; expected \"Point\", \"Mean\" or \"Replace\".";
SolveActiveNematicComovingSteadyQfree::badmodel =
  "modelParams is missing or non-numeric for: `1`.  Expected rules for " <>
  "{B, a, b, L, \[Zeta], \[Xi]0, \[Xi]r, \[Rho]0, \[CapitalLambda]}.";
SolveActiveNematicComovingSteadyQfree::badsolver =
  "solverParams is missing or non-numeric for: `1`.  Expected rules for " <>
  "{\[Delta]mesh, hMax, box, maxSteps, relativeDiff}.";
SolveActiveNematicComovingSteadyQfree::noconv =
  "Bordered Newton did not converge at the target \[Zeta] within maxSteps = `1` " <>
  "iterations (final |step| tolerance `2`, residual = `3`).  u = `4` is reported " <>
  "anyway and should not be trusted.";
SolveActiveNematicComovingSteadyQfree::stalled =
  "\[Zeta]-continuation stalled at \[Zeta] = `1` (residual = `2`); continuing.";


SolveActiveNematicComovingSteadyQfree[
    modelParams : {___Rule},
    solverParams : {___Rule},
    bcType_String : "Neumann"] :=
Block[
  {
    \[Delta]mesh, hMax, box, maxSteps, relativeDiff,
    refine, rhoGauge, zetaSteps, order, stepTolerance, verbose, reuseFac,
    qWallV, symLineV, symOpsV,
    B, a, b, L, \[Zeta], \[Xi]0, \[Xi]r, \[Rho]0, \[CapitalLambda],
    missing, mkModel, dpOf, dpTarget, dp, pf, gg, op, bcz, ne, n, k0,
    X, ex, st, zs, kk, wall0, rhoGaugeV, sym, lf
  },

  If[!MemberQ[{"Neumann", "Dirichlet", "RhoNeumann"}, bcType],
    Message[SolveActiveNematicComovingSteadyQfree::badbc, bcType];
    Return[$Failed]];

  (* Lookup against an Association, as in comoving_steady_solver.m: the
     "{k,...} = {k,...} /. solverParams" idiom self-assigns absent keys. *)
  With[{sp = Association[solverParams]},
    \[Delta]mesh  = Lookup[sp, \[Delta]mesh,  Missing[]];
    hMax          = Lookup[sp, hMax,          Missing[]];
    box           = Lookup[sp, box,           Missing[]];
    maxSteps      = Lookup[sp, maxSteps,      Missing[]];
    relativeDiff  = Lookup[sp, relativeDiff,  Missing[]];
    refine        = Lookup[sp, refine,        1];
    rhoGauge      = Lookup[sp, rhoGauge,      "Point"];
    zetaSteps     = Lookup[sp, zetaSteps,     5];
    order         = Lookup[sp, order,         4];
    stepTolerance = Lookup[sp, stepTolerance, 1.*^-10];
    verbose       = TrueQ @ Lookup[sp, verbose, True];
    reuseFac      = Lookup[sp, "ReuseFactorization", Automatic];
    qWallV        = Lookup[sp, qWall,         "Neumann"];
    symLineV      = TrueQ @ Lookup[sp, symLine, True];
    symOpsV       = TrueQ @ Lookup[sp, symOps,  True];
  ];

  missing = Pick[{"\[Delta]mesh", "hMax", "box", "maxSteps", "relativeDiff"},
                 Not /@ NumericQ /@ {\[Delta]mesh, hMax, box, maxSteps,
                                     relativeDiff}];
  If[missing =!= {},
    Message[SolveActiveNematicComovingSteadyQfree::badsolver, missing];
    Return[$Failed]];
  If[!MemberQ[{True, False, Automatic}, reuseFac],
    Message[SolveActiveNematicComovingSteadyQfree::badreuse, reuseFac];
    Return[$Failed]];
  If[!MemberQ[{"Point", "Mean", "Replace"}, rhoGauge],
    Message[SolveActiveNematicComovingSteadyQfree::badgauge, rhoGauge];
    Return[$Failed]];
  If[!MemberQ[{"Neumann", "Pade"}, qWallV],
    Message[SolveActiveNematicComovingSteadyQfree::badqwall, qWallV];
    Return[$Failed]];
  If[qWallV === "Neumann" && !symLineV,
    Message[SolveActiveNematicComovingSteadyQfree::noline]];
  If[bcType === "RhoNeumann" && rhoGauge =!= "Point",
    Message[SolveActiveNematicComovingSteadyQfree::gaugeignored, rhoGauge]];
  rhoGaugeV = rhoGauge;

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
    Message[SolveActiveNematicComovingSteadyQfree::badmodel, missing];
    Return[$Failed]];

  mkModel[z_] := <|"B" -> N[B], "a" -> N[a], "b" -> N[b], "K" -> N[L],
                   "zeta" -> N[z], "xi0" -> N[\[Xi]0], "xir" -> N[\[Xi]r],
                   "rho0" -> N[\[Rho]0], "Lambda" -> N[\[CapitalLambda]]|>;
  dpOf[z_] := Append[cdDerived[mkModel[z]], "rhoGauge" -> rhoGaugeV];

  dpTarget = dpOf[\[Zeta]];
  wall0 = AbsoluteTime[];

  (* --- grid, operators ---------------------------------------------- *)
  gg = cdGridPade[N[box], dpTarget, N[\[Delta]mesh], N[hMax], refine];
  op = cdOps[gg, gg, order];
  If[symOpsV, op = cdMirrorOps[op]];
  n  = op["n"]; k0 = op["k0"];
  ne = cdNe[dpTarget, bcType];

  If[verbose,
    Print["[comoving steady, free Q] bc = ", bcType, "   qWall = ", qWallV,
          "   symLine = ", symLineV, "   symOps = ", symOpsV,
          "   box = ", N[box], "   n = ", n, " (", op["nx"], "^2)",
          "   hMin = ", cdFmt@gg["hMin"]];
    Print["    ld = ", cdFmt@dpTarget["ld"], "   S0 = ", cdFmt@dpTarget["S0"],
          "   ld/hMin = ", cdFmt[dpTarget["ld"]/gg["hMin"]],
          "   alpha = ", cdFmt[\[Zeta] \[Xi]r dpTarget["S0"]/
                               (4 \[Xi]0 dpTarget["Kp"])],
          If[cdGauge[dpTarget, bcType] =!= None,
             "   gauge = " <> cdGauge[dpTarget, bcType], ""]]];

  (* --- zeta-continuation, as in comoving_steady_solver.m -------------- *)
  (* The zeta = 0 step also carries the Pade seed onto the free wall: the
     passive defect in a Neumann box differs from Pade by O(1) near the
     wall, and Newton takes 4-5 iterations for it. *)
  zs = If[TrueQ[\[Zeta] == 0.], {0.},
          Prepend[Rest@Subdivide[0., N[\[Zeta]], zetaSteps], 0.]];

  dp  = dpOf[0.];
  pf  = cdPadeFuns[dp];
  bcz = cdBCQfree[op, dp, pf, bcType, qWallV, symLineV];
  X   = cdInitialRho[cdInitial[op, pf], op, dp, bcz, bcType];
  ex  = ConstantArray[0., ne];

  Do[
    dp  = dpOf[zs[[kk]]];
    pf  = cdPadeFuns[dp];
    bcz = cdBCQfree[op, dp, pf, bcType, qWallV, symLineV];
    st  = cdSteady[X, ex, op, dp, bcz, bcType,
            "MaxIterations" -> maxSteps,
            "StepTolerance" -> If[kk == Length[zs],
                                  N[stepTolerance], N[relativeDiff]],
            "Verbose" -> False,
            "ReuseFactorization" -> reuseFac];
    X = st["X"]; ex = st["ex"];
    If[verbose,
      Print["    zeta = ", StringPadRight[cdFmt[zs[[kk]]], 12],
            "  newton = ", StringPadRight[ToString@st["iterations"], 4],
            "  |R| = ", StringPadRight[cdFmt@st["residual"], 12],
            "  ux = ", cdFmt@ex[[1]]]];
    If[!TrueQ@st["converged"] && kk < Length[zs],
      Message[SolveActiveNematicComovingSteadyQfree::stalled,
              cdFmt@zs[[kk]], cdFmt@st["residual"]]],
    {kk, Length[zs]}];

  If[!TrueQ@st["converged"],
    Message[SolveActiveNematicComovingSteadyQfree::noconv,
            maxSteps, cdFmt@N[stepTolerance], cdFmt@st["residual"],
            cdFmt@ex[[1]]]];

  (* unit mean density for the sealed wall, as in comoving_steady_solver.m *)
  If[bcType === "Neumann", X = cdShiftRhoMean[X, op]];

  sym = cdSymmetry[X, op];
  lf  = cdLineForce[X, ex, op, dp];

  If[verbose,
    Print["    u = ", ToString@DecimalForm[ex[[1]], 12],
          "   uy = ", cdFmt@ex[[2]],
          "   core = ", cdFmt[{X[[n + k0]], X[[2 n + k0]]}],
          If[ne >= 3, "   mu = " <> cdFmt@ex[[3]], ""],
          "   lineForce = ", cdFmt@lf,
          "   wall = ", Round[N[AbsoluteTime[] - wall0], 0.01], " s"]];

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
                   "gauge"      -> cdGauge[dpTarget, bcType],
                   "reuseFactorization" -> reuseFac,
                   "mu"         -> If[Length[ex] >= 3, ex[[3]], 0.],
                   "qWall"      -> qWallV,
                   "symLine"    -> symLineV,
                   "symOps"     -> symOpsV,
                   "lineForce"  -> lf,
                   "symmetry"   -> sym,
                   "state"      -> X,
                   "wall"       -> N[AbsoluteTime[] - wall0]|> |>
];
