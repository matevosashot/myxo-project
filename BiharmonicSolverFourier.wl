(* ::Package:: *)
(* ================================================================ *)
(*  BiharmonicSolver.wl                                             *)
(*  Solves: Lap_r Lap_{r'} f(r,r') = d_{r_i} d_{r'_j} G_{ij}(r,r')*)
(*  r, r' in R^2, f regular, f -> 0 at infinity                    *)
(* ================================================================ *)

BeginPackage["BiharmonicSolver`"]

(* Public symbols *)
solveBiharmonicFourier::usage =
  "solveBiharmonicFourier[gFunc, n, L] solves the biharmonic equation \
on an n^4 grid over [-L,L]^2 x [-L,L]^2. Returns an Association.";
buildInterpolation::usage =
  "buildInterpolation[sol] builds a 4D InterpolatingFunction from the solution.";
extractSlice::usage =
  "extractSlice[sol, x2, y2] extracts the 2D slice f(x,y,x2,y2) at fixed r'.";

Begin["`Private`"]

(* --- Compiled kernel for Fourier-space division --- *)
(* Avoids building full 4D wavevector arrays; saves memory *)
(* Inputs: real/imag parts of G_hat components + 1D wavevector array *)
biharmonicKernel = Compile[{
    {g11Re, _Real, 4}, {g11Im, _Real, 4},
    {g12Re, _Real, 4}, {g12Im, _Real, 4},
    {g21Re, _Real, 4}, {g21Im, _Real, 4},
    {g22Re, _Real, 4}, {g22Im, _Real, 4},
    {qv, _Real, 1}, {reg, _Real}},
  Module[{nn, qx, qy, qxp, qyp, q1s, q2s, dd, nR, nI},
    nn = Length[qv];
    Table[
      qx = qv[[a]]; qy = qv[[b]];
      qxp = qv[[c]]; qyp = qv[[e]];
      q1s = qx * qx + qy * qy;
      q2s = qxp * qxp + qyp * qyp;
      dd = q1s * q2s + reg;
      nR = -(qx*qxp*g11Re[[a,b,c,e]] + qx*qyp*g12Re[[a,b,c,e]] +
             qy*qxp*g21Re[[a,b,c,e]] + qy*qyp*g22Re[[a,b,c,e]]);
      nI = -(qx*qxp*g11Im[[a,b,c,e]] + qx*qyp*g12Im[[a,b,c,e]] +
             qy*qxp*g21Im[[a,b,c,e]] + qy*qyp*g22Im[[a,b,c,e]]);
      {nR / dd, nI / dd},
      {a, nn}, {b, nn}, {c, nn}, {e, nn}
    ]
  ],
  CompilationTarget -> "C",
  RuntimeOptions -> "Speed",
  Parallelization -> True
];

(* --- Main solver --- *)
Options[solveBiharmonicFourier] = {
  "Regularization" -> Automatic,
  "Method" -> "Compiled"  (* "Compiled" or "Vectorized" *)
};

solveBiharmonicFourier[gFunc_, n_Integer, L_?NumericQ,
    OptionsPattern[]] :=
Module[{dx, pts, qVals, gData,
    g11, g12, g21, g22,
    g11H, g12H, g21H, g22H,
    fHat, fGrid, eps, method, t0},

  method = OptionValue["Method"];
  Print["═══════════════════════════════════════════════"];
  Print["  Biharmonic Fourier Solver"];
  Print["═══════════════════════════════════════════════"];
  Print["  Grid size: ", n, "^4 = ", n^4, " points"];
  Print["  Domain: [-", L, ", ", L, "]^2 x [-", L, ", ", L, "]^2"];
  Print["  Method: ", method];
  Print["───────────────────────────────────────────────"];

  dx = 2.0 L / n;
  pts = N @ Range[-L, L - dx, dx];

  (* Wavevectors: 0, dq, 2dq, ..., n/2 dq, -(n/2-1)dq, ..., -dq *)
  qVals = N @ Table[
    If[k <= Quotient[n, 2], k, k - n],
    {k, 0, n - 1}
  ] * (Pi / L);

  (* --- Step 1: Evaluate G_{ij} on 4D grid (parallelized) --- *)
  Print["[Step 1/4] Evaluating G_{ij} on 4D grid (parallelized)..."];
  t0 = AbsoluteTime[];
  DistributeDefinitions[gFunc, pts, n];
  gData = ParallelTable[
    gFunc[pts[[i1]], pts[[j1]], pts[[i2]], pts[[j2]]],
    {i1, n}, {j1, n}, {i2, n}, {j2, n},
    Method -> "CoarsestGrained"
  ];

  (* Extract components as packed real arrays *)
  Print["  Extracting G components..."];
  g11 = Developer`ToPackedArray @ N @ gData[[All, All, All, All, 1, 1]];
  g12 = Developer`ToPackedArray @ N @ gData[[All, All, All, All, 1, 2]];
  g21 = Developer`ToPackedArray @ N @ gData[[All, All, All, All, 2, 1]];
  g22 = Developer`ToPackedArray @ N @ gData[[All, All, All, All, 2, 2]];
  gData = Null;
  Print["  Done. (", Round[AbsoluteTime[] - t0, 0.01], " s)"];
  Print["  Memory for G components: ~",
    Round[4 * n^4 * 8 / 1024.^2, 0.1], " MB"];

  (* --- Step 2: 4D FFT of each component (parallel) --- *)
  Print["[Step 2/4] Computing 4D FFT of G components (4 parallel FFTs)..."];
  t0 = AbsoluteTime[];
  {g11H, g12H, g21H, g22H} = WaitAll[{
    ParallelSubmit[{g11}, Fourier[g11, FourierParameters -> {1, -1}]],
    ParallelSubmit[{g12}, Fourier[g12, FourierParameters -> {1, -1}]],
    ParallelSubmit[{g21}, Fourier[g21, FourierParameters -> {1, -1}]],
    ParallelSubmit[{g22}, Fourier[g22, FourierParameters -> {1, -1}]]
  }];
  g11 = Null; g12 = Null; g21 = Null; g22 = Null;
  Print["  Done. (", Round[AbsoluteTime[] - t0, 0.01], " s)"];

  (* Regularization parameter *)
  eps = If[OptionValue["Regularization"] === Automatic,
    (Pi / (L * n))^4,
    OptionValue["Regularization"]
  ];

  (* --- Step 3: Pointwise division in Fourier space --- *)
  Print["[Step 3/4] Pointwise division in Fourier space..."];
  Print["  Regularization eps = ", ScientificForm[eps, 3]];
  t0 = AbsoluteTime[];
  fHat = Switch[method,

    (* --- Method A: Compiled C kernel (memory-efficient) --- *)
    "Compiled",
    Module[{raw},
      raw = biharmonicKernel[
        Re[g11H], Im[g11H], Re[g12H], Im[g12H],
        Re[g21H], Im[g21H], Re[g22H], Im[g22H],
        qVals, eps
      ];
      g11H = Null; g12H = Null;
      g21H = Null; g22H = Null;
      (* raw is {n,n,n,n,2} -> convert to complex {n,n,n,n} *)
      raw[[All, All, All, All, 1]] +
        I * raw[[All, All, All, All, 2]]
    ],

    (* --- Method B: Vectorized array ops (faster, more memory) --- *)
    "Vectorized",
    Module[{one, qx, qy, qxp, qyp, q1Sq, q2Sq, num, den},
      one = ConstantArray[1., n];
      qx  = Outer[Times, qVals, one, one, one];
      qy  = Outer[Times, one, qVals, one, one];
      qxp = Outer[Times, one, one, qVals, one];
      qyp = Outer[Times, one, one, one, qVals];
      one = Null;

      num = -(qx*qxp*g11H + qx*qyp*g12H +
              qy*qxp*g21H + qy*qyp*g22H);
      g11H = Null; g12H = Null;
      g21H = Null; g22H = Null;

      q1Sq = qx^2 + qy^2;  qx = Null; qy = Null;
      q2Sq = qxp^2 + qyp^2; qxp = Null; qyp = Null;
      den = q1Sq * q2Sq + eps;
      q1Sq = Null; q2Sq = Null;

      num / den
    ]
  ];
  Print["  Done. (", Round[AbsoluteTime[] - t0, 0.01], " s)"];

  (* Zero mode: enforce f -> 0 at infinity *)
  fHat[[1, 1, 1, 1]] = 0. + 0. I;

  (* --- Step 4: Inverse 4D FFT --- *)
  Print["[Step 4/4] Computing inverse 4D FFT..."];
  t0 = AbsoluteTime[];
  fGrid = Re @ InverseFourier[fHat, FourierParameters -> {1, -1}];
  fHat = Null;
  Print["  Done. (", Round[AbsoluteTime[] - t0, 0.01], " s)"];

  Print["───────────────────────────────────────────────"];
  Print["  Max |f| = ", ScientificForm[Max[Abs[fGrid]], 4]];
  Print["  Solution complete."];
  Print["═══════════════════════════════════════════════"];

  (* --- Return solution --- *)
  <|
    "f" -> fGrid,
    "grid" -> pts,
    "dx" -> dx,
    "L" -> N[L],
    "n" -> n
  |>
];

(* --- Utility: extract f at specific (r, r') via interpolation --- *)
buildInterpolation[sol_Association] :=
Module[{pts, fData, coords, interpData},
  pts = sol["grid"];
  fData = sol["f"];
  (* Build 4D InterpolatingFunction *)
  ListInterpolation[fData, {pts, pts, pts, pts}]
];

(* --- Utility: extract 2D slice f(r, r0) for fixed r' = r0 --- *)
extractSlice[sol_Association, x2_?NumericQ, y2_?NumericQ] :=
Module[{pts, n, i2, j2},
  pts = sol["grid"];
  n = sol["n"];
  (* Find nearest grid indices *)
  i2 = Nearest[pts -> "Index", x2][[1]];
  j2 = Nearest[pts -> "Index", y2][[1]];
  (* Return 2D array f(x, y, x2, y2) and grid *)
  <|"slice" -> sol["f"][[All, All, i2, j2]], "grid" -> pts|>
];


(* ================================================================ *)
End[]       (* `Private` *)
EndPackage[]

(* ================================================================ *)
(*  TEST EXAMPLE                                                     *)
(*  Isotropic Gaussian: G_{ij} = delta_{ij} * exp(-|r-r'|^2/2s^2)  *)
(*  Analytic check: f should be ~log-like at large |r-r'|           *)
(* ================================================================ *)

(* Compiled test G function *)
gTest = Compile[{{x1, _Real}, {y1, _Real}, {x2, _Real}, {y2, _Real}},
  Module[{s2 = 1.0, amp = 1.0, r2},
    r2 = (x1 - x2)^2 + (y1 - y2)^2;
    amp * Exp[-r2 / (2.0 s2)] * {{1., 0.}, {0., 1.}}
  ],
  RuntimeOptions -> "Speed"
];

(* Run solver *)
(* Uncomment to execute:

LaunchKernels[];  (* start parallel kernels if not already running *)

sol = solveBiharmonicFourier[gTest, 32, 8.0,
  "Method" -> "Compiled"];

(* Visualize diagonal slice f(x, y, x, y) *)
fDiag = Table[sol["f"][[i, j, i, j]], {i, sol["n"]}, {j, sol["n"]}];
ListPlot3D[fDiag, PlotRange -> All, PlotLabel -> "f(r, r)"]

(* Visualize f(x, y, 0, 0) *)
slice = extractSlice[sol, 0., 0.];
ListPlot3D[slice["slice"], DataRange -> {{-8, 8}, {-8, 8}},
  PlotRange -> All, PlotLabel -> "f(r, 0)"]

(* Build interpolating function for arbitrary evaluation *)
fInterp = buildInterpolation[sol];
fInterp[1.0, 0.5, -0.3, 0.2]  (* evaluate at any point *)

*)