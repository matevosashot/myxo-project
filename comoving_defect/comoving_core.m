(* ::Package:: *)

(* ==================================================================== *)
(*  comoving_core.m                                                     *)
(*                                                                      *)
(*  Coupled rho-Q active-nematic system in the comoving frame of a      *)
(*  +1/2 defect held at the origin.  Equations (dean.tex 720-731):      *)
(*                                                                      *)
(*    d_t rho = (B/xi0) Lap rho + (zeta/xi0) d_a d_b Q_ab               *)
(*    d_t Q   = Div( Q (u + (B/xi0) Grad rho) ) + (4/xir) H_active      *)
(*    H_active = -a' Q - b (Q:Q) Q + K' Lap Q                           *)
(*                                                                      *)
(*  with a' = a + Lambda xir,  K' = K + zeta xir/(4 xi0).               *)
(*                                                                      *)
(*  u is NOT measured at the core.  The semi-discrete system is treated *)
(*  as a DAE:  u is the multiplier conjugate to the pinning constraints *)
(*  Q1(0,0) = Q2(0,0) = 0, so the constraint holds exactly at every     *)
(*  step and no second derivative at the core is ever needed.           *)
(*                                                                      *)
(*  Q = {{Q1, Q2}, {Q2, -Q1}}.                                          *)
(*  Grid ordering k = i + (j-1) nx  (x fastest), matching the repo.     *)
(* ==================================================================== *)

ClearAll[cdGrid, cdOps, cdDerived, cdPadeFuns, cdInitial, cdBC, cdRJ,
         cdNewton, cdEvolve, cdSteady, cdFieldAt, cdCoreG, cdInterp,
         cdRes, cdRates, cdBorder, cdShiftPos, cdFmt, cdRadialBVP,
         cdGridPade, cdCoreEstimates, cdSymmetry, cdInitialRho, cdBScale];

(* Compact SINGLE-LINE numeric form.  ScientificForm renders as a
   two-line superscript under wolframscript's OutputForm, which breaks up
   every progress line into unreadable fragments, so the mantissa and
   exponent are assembled by hand as "1.056e-7". *)
cdFmt[x_?NumericQ] := Module[{v = N[x], e, m},
  Which[
    !NumericQ[v] || v === Indeterminate, ToString[v],
    v == 0., "0",
    !(Abs[v] < Infinity), ToString[v],
    True,
      e = Floor[Log10[Abs[v]]];
      m = v/10.^e;
      If[Abs[m] >= 9.9995, m = m/10.; e = e + 1];
      ToString[NumberForm[m, {5, 3}]] <> "e" <> ToString[e]]];
cdFmt[l_List] := "{" <> StringRiffle[cdFmt /@ l, ", "] <> "}";
cdFmt[x_] := ToString[x];

(* -------------------------------------------------------------------- *)
(* 1D stretched grid.  Odd point count so the origin is a node; sinh     *)
(* stretch is fine at the core and coarse at the wall.  s -> 0 is        *)
(* uniform.  w = trapezoid quadrature weights on the same grid.          *)
(* -------------------------------------------------------------------- *)
cdGrid[L_?NumericQ, n_Integer, s_?NumericQ] := Module[{xi, g, w},
  If[EvenQ[n] || n < 5,
    Print["cdGrid: n must be odd and >= 5, got ", n]; Return[$Failed]];
  xi = N @ Table[-1 + 2 (k - 1)/(n - 1), {k, n}];
  g  = N @ If[Abs[s] < 1.*^-8, L xi, L Sinh[s xi]/Sinh[s]];
  g[[(n + 1)/2]] = 0.;                       (* exact zero at the centre *)
  w = Table[
        Which[
          k == 1, (g[[2]] - g[[1]])/2,
          k == n, (g[[n]] - g[[n - 1]])/2,
          True,   (g[[k + 1]] - g[[k - 1]])/2],
        {k, n}];
  <|"g" -> g, "w" -> w, "n" -> n, "L" -> L, "s" -> s,
    "hMin" -> Min@Differences[g], "hMax" -> Max@Differences[g]|>
];

(* -------------------------------------------------------------------- *)
(* Pade-graded 1D grid -- the tensor-product analogue of the FEM         *)
(* solvers' refinement criterion (iterative_solver_module.m:126-138,     *)
(* transient_solver_module.m:182-191):                                   *)
(*                                                                      *)
(*    dS0[r]           = d/dr ( r S0PadeOverR[r^2] )   = S'(r)          *)
(*    targetCellSize   = Min[dmesh/dS0[r], hMax],  r = |x|/2            *)
(*                                                                      *)
(* i.e. cell size ~ 1/|S'(r)|: fine where the defect amplitude turns     *)
(* over, coarse in the far field where S -> S0 and S' -> 0.  dmesh is    *)
(* the change in S allowed across one cell.                             *)
(*                                                                      *)
(* targetCellSize is a 2D area criterion, so for a tensor grid we build  *)
(* the 1D version by counting cells: xi(x) = Int_0^x dx'/h(x') is the    *)
(* cell index, so sampling x at uniform xi gives spacing h(x) by         *)
(* construction.  n = 2m+1 is DERIVED from dmesh/hMax, as in the FEM     *)
(* solver, rather than dialled in.  refine (1,2,4,...) scales every cell *)
(* together, so grading is preserved under refinement and convergence    *)
(* orders stay meaningful.                                              *)
(* -------------------------------------------------------------------- *)
cdGridPade[LL_?NumericQ, dp_Association, dmesh_?NumericQ, hMaxv_?NumericQ,
           refine_: 1] :=
Module[{S0, ld, rr, sExpr, dsE, dS, hOf, nsub, xs, dx, invh, cum, tot, m,
        xOf, gpos, g, w, n},
  {S0, ld} = Lookup[dp, {"S0", "ld"}];

  sExpr = rr (S0/ld) Sqrt[(0.34 + 0.07 (rr^2/ld^2))/
                          (1 + 0.41 (rr^2/ld^2) + 0.07 (rr^2/ld^2)^2)];
  dsE   = D[sExpr, rr];
  dS[rv_] = dsE /. rr -> rv;                     (* Set: evaluates now *)

  hOf[x_] := Min[dmesh/Max[dS[Abs[x]/2], 1.*^-12], hMaxv];

  (* xi(x) = Int_0^x dx'/h  on a fine subgrid, by trapezoid *)
  nsub = 20000;
  xs   = N@Subdivide[0., LL, nsub];
  dx   = LL/nsub;
  invh = 1./(hOf /@ xs);
  cum  = Prepend[Accumulate[(Most[invh] + Rest[invh]) dx/2], 0.];
  tot  = Last[cum];
  m    = Max[3, Ceiling[refine tot]];

  (* invert xi -> x and sample at uniform xi *)
  xOf  = Interpolation[Transpose@{cum, xs}, InterpolationOrder -> 3];
  gpos = xOf /@ (Range[0, m] tot/m);
  gpos[[1]] = 0.; gpos[[-1]] = N[LL];

  g = Join[-Reverse[Rest[gpos]], gpos];          (* symmetric, 0 is a node *)
  n = Length[g];
  w = Table[
        Which[k == 1, (g[[2]] - g[[1]])/2,
              k == n, (g[[n]] - g[[n - 1]])/2,
              True,   (g[[k + 1]] - g[[k - 1]])/2],
        {k, n}];
  <|"g" -> g, "w" -> w, "n" -> n, "L" -> LL, "s" -> Missing["Pade"],
    "hMin" -> Min@Differences[g], "hMax" -> Max@Differences[g],
    "dmesh" -> dmesh, "hMaxTarget" -> hMaxv, "refine" -> refine,
    "cells" -> tot, "hOf" -> hOf|>
];

(* -------------------------------------------------------------------- *)
(* 2D tensor-product operators by Kronecker lift.  Same construction as  *)
(* lyapunov_solver_module.m:229-262, but on a non-uniform grid --        *)
(* FiniteDifferenceDerivative handles that with no change.               *)
(* -------------------------------------------------------------------- *)
cdOps[gx_Association, gy_Association, order_Integer: 4] :=
Module[{nx, ny, n, gxv, gyv, d1x, l1x, d1y, l1y, Ix, Iy,
        Dx, Dy, Lx, Ly, Dxy, xs, ys, i0, j0, k0, ii, jj,
        bmask, imask, px, py, wgt},
  gxv = gx["g"]; gyv = gy["g"];
  nx = Length[gxv]; ny = Length[gyv]; n = nx ny;

  d1x = SparseArray @ NDSolve`FiniteDifferenceDerivative[
          Derivative[1], gxv, "DifferenceOrder" -> order]["DifferentiationMatrix"];
  l1x = SparseArray @ NDSolve`FiniteDifferenceDerivative[
          Derivative[2], gxv, "DifferenceOrder" -> order]["DifferentiationMatrix"];
  d1y = SparseArray @ NDSolve`FiniteDifferenceDerivative[
          Derivative[1], gyv, "DifferenceOrder" -> order]["DifferentiationMatrix"];
  l1y = SparseArray @ NDSolve`FiniteDifferenceDerivative[
          Derivative[2], gyv, "DifferenceOrder" -> order]["DifferentiationMatrix"];

  Ix = IdentityMatrix[nx, SparseArray];
  Iy = IdentityMatrix[ny, SparseArray];

  (* k = i + (j-1) nx  =>  x index is the FAST one  =>  KroneckerProduct[y, x] *)
  Dx  = KroneckerProduct[Iy,  d1x];
  Dy  = KroneckerProduct[d1y, Ix ];
  Lx  = KroneckerProduct[Iy,  l1x];
  Ly  = KroneckerProduct[l1y, Ix ];
  Dxy = KroneckerProduct[d1y, d1x];

  xs = Flatten @ ConstantArray[gxv, ny];
  ys = Flatten @ Table[ConstantArray[gyv[[j]], nx], {j, ny}];
  ii = Flatten @ ConstantArray[Range[nx], ny];
  jj = Flatten @ Table[ConstantArray[j, nx], {j, ny}];

  i0 = (nx + 1)/2; j0 = (ny + 1)/2;
  k0 = i0 + (j0 - 1) nx;

  bmask = MapThread[If[#1 == 1 || #1 == nx || #2 == 1 || #2 == ny, 1, 0] &, {ii, jj}];
  imask = 1 - bmask;                                 (* 1 on interior nodes *)
  (* SIGNED outward normal components.  Using unsigned 0/1 masks here is
     wrong at the two corners where the normal is (+1,-1) or (-1,+1):
     the condition there is (d_x - d_y) rho = 0, not (d_x + d_y) rho = 0.
     With the unsigned version those two corner rows become exact linear
     combinations of the neighbouring edge rows, and the whole bordered
     matrix drops rank by exactly 2 at every grid size and stencil order. *)
  px    = Map[Which[# == 1, -1, # == nx, 1, True, 0] &, ii];
  py    = Map[Which[# == 1, -1, # == ny, 1, True, 0] &, jj];
  wgt   = Flatten @ Outer[Times, gy["w"], gx["w"]];  (* quadrature, k-ordered *)

  <|"nx" -> nx, "ny" -> ny, "n" -> n, "gx" -> gxv, "gy" -> gyv,
    "Dx" -> Dx, "Dy" -> Dy, "Lx" -> Lx, "Ly" -> Ly, "Lap" -> Lx + Ly,
    "Dxy" -> Dxy, "xs" -> xs, "ys" -> ys, "ii" -> ii, "jj" -> jj,
    "k0" -> k0, "bmask" -> bmask, "imask" -> imask,
    "px" -> px, "py" -> py, "w" -> wgt, "area" -> Total[wgt],
    "order" -> order|>
];

(* -------------------------------------------------------------------- *)
(* Model parameters -> the coefficients that actually appear.            *)
(* modelParams uses the repo convention {B,a,b,L,zeta,xi0,xir,rho0,Lam}, *)
(* where L is the Frank constant K.                                      *)
(* -------------------------------------------------------------------- *)
cdDerived[mp_Association] := Module[{B, a, b, K, ze, x0, xr, Lam, ap, Kp},
  {B, a, b, K, ze, x0, xr, Lam} =
    Lookup[mp, {"B", "a", "b", "K", "zeta", "xi0", "xir", "Lambda"}];
  ap = a + Lam xr;                     (* a' *)
  Kp = K + ze xr/(4 x0);               (* K' *)
  <|"Bx"  -> B/x0,                     (* B/xi0            *)
    "lam" -> ze/x0,                    (* zeta/xi0         *)
    "kap" -> 4 Kp/xr,                  (* 4K'/xir          *)
    "Aq"  -> 4 ap/xr,                  (* 4a'/xir          *)
    "bet" -> 8 b/xr,                   (* 8b/xir           *)
    "ap" -> ap, "Kp" -> Kp,
    "ld"  -> Sqrt[Kp/(-ap)],           (* defect core size *)
    "S0"  -> Sqrt[-ap/(2 b)]|>
];

(* -------------------------------------------------------------------- *)
(* Pade +1/2 defect: Q1 = S(r) cos(phi), Q2 = S(r) sin(phi), with        *)
(* S = S0 f(r/ld).  Identical to transient_solver_module.m:170-176.      *)
(* -------------------------------------------------------------------- *)
cdPadeFuns[dp_Association] := Module[{S0, ld, sOverR},
  {S0, ld} = Lookup[dp, {"S0", "ld"}];
  sOverR[r2_] := (S0/ld) Sqrt[
     (0.34 + 0.07 (r2/ld^2)) /
     (1 + 0.41 (r2/ld^2) + 0.07 (r2/ld^2)^2)];
  <|"q1" -> Function[{x, y}, sOverR[x^2 + y^2] x],
    "q2" -> Function[{x, y}, sOverR[x^2 + y^2] y],
    "sOverR" -> sOverR|>
];

(* -------------------------------------------------------------------- *)
(* Initial state: uniform rho = 1 with the Pade core in Q.               *)
(* -------------------------------------------------------------------- *)
cdInitial[op_Association, pf_Association] := Module[{xs, ys},
  xs = op["xs"]; ys = op["ys"];
  Join[ConstantArray[1., op["n"]], pf["q1"][xs, ys], pf["q2"][xs, ys]]
];

(* -------------------------------------------------------------------- *)
(* Quasi-static rho for a FROZEN Q: solve                                *)
(*    (B/xi0) Lap rho + (zeta/xi0) d_a d_b Q_ab + mu = 0                 *)
(* with the chosen wall condition (+ the mean row for "Neumann").        *)
(* The Q equation is untouched, so this is a single LINEAR solve.        *)
(*                                                                      *)
(* This matters: rho == 1 satisfies the Dirichlet wall but NOT the       *)
(* sealed one (zeta n.(div Q) != 0 there), so starting from it forces a  *)
(* violent boundary layer in the first implicit step.  With B/xi0 = 200  *)
(* the resulting grad rho feeds straight into the advection velocity     *)
(* w = u + (B/xi0) grad rho and Newton diverges.                         *)
(* -------------------------------------------------------------------- *)
cdInitialRho[X_?VectorQ, op_Association, dp_Association, bc_Association,
             bcType_String] :=
Module[{n = op["n"], q1, q2, Bx, lam, im, Arho, Aq, src, bq, rhs, sol, M, b},
  q1 = X[[n + 1 ;; 2 n]]; q2 = X[[2 n + 1 ;; 3 n]];
  {Bx, lam} = Lookup[dp, {"Bx", "lam"}];
  im = op["imask"];

  Arho = im (Bx op["Lap"]) + bc["Jbc"][[1 ;; n, 1 ;; n]];
  src  = im (lam ((op["Lx"] . q1) - (op["Ly"] . q1) + 2 (op["Dxy"] . q2)));
  bq   = bc["Jbc"][[1 ;; n, n + 1 ;; 3 n]] . Join[q1, q2];
  rhs  = -src - bq + bc["rhs"][[1 ;; n]];

  If[bcType === "Neumann",
    (* border with the mu column and the gauge row (see cdRhoGauge) *)
    M = cdBorder[SparseArray[Arho], SparseArray[Transpose[{im}]],
                 SparseArray[If[cdRhoGauge[dp] === "Mean",
                    {op["w"]},
                    {SparseArray[{op["k0"] -> 1.}, {n}]}]], 1];
    b = Append[rhs, If[cdRhoGauge[dp] === "Mean", op["area"], 1.]];
    sol = LinearSolve[M, b];
    Join[sol[[1 ;; n]], q1, q2],
  (* else Dirichlet: no multiplier *)
    sol = LinearSolve[Arho, rhs];
    Join[sol, q1, q2]]
];

(* -------------------------------------------------------------------- *)
(* Boundary rows, as a LINEAR pair (Jbc, rhsbc): the residual on a       *)
(* boundary row is  Jbc.X - rhsbc.  Jbc is zero on every interior row.   *)
(*                                                                      *)
(*   Q      : Dirichlet on the Pade profile, both bcType values.         *)
(*   rho    : "Dirichlet" -> rho = 1                                     *)
(*            "Neumann"   -> sealed wall n.(B grad rho + zeta div Q) = 0,*)
(*                           scaled by 1/xi0 for conditioning.  Corners  *)
(*                           get the sum of the x and y conditions.      *)
(* -------------------------------------------------------------------- *)
cdBC[op_Association, dp_Association, pf_Association, bcType_String] :=
Module[{n, bm, px, py, Dx, Dy, Bx, lam, Z, Ibc, Jrr, Jrq1, Jrq2, Jbc, rhs,
        qb1, qb2, kg, imRho},
  n = op["n"]; bm = op["bmask"]; px = op["px"]; py = op["py"];
  Dx = op["Dx"]; Dy = op["Dy"];
  {Bx, lam} = Lookup[dp, {"Bx", "lam"}];
  Z   = SparseArray[{}, {n, n}];
  Ibc = DiagonalMatrix[SparseArray[bm]];         (* identity on boundary rows *)

  Switch[bcType,
    "Dirichlet",
      Jrr = Ibc; Jrq1 = Z; Jrq2 = Z;,
    "Neumann",
      (* n.(Bx grad rho + lam div Q) = 0 with n = (px, py) the SIGNED
         outward normal (unnormalised -- the condition is homogeneous):
           n.div Q = px (d_x Q1 + d_y Q2) + py (d_x Q2 - d_y Q1)        *)
      Jrr  = Bx  (px Dx + py Dy);
      Jrq1 = lam (px Dx - py Dy);
      Jrq2 = lam (px Dy + py Dx);,
    _,
      Print["cdBC: unknown bcType ", bcType]; Return[$Failed]
  ];

  (* "Replace" gauge: overwrite the PDE row at one near-wall node with
     rho = 1, so no extra row and no multiplier are needed. *)
  kg = None; imRho = op["imask"];
  If[bcType === "Neumann" && cdRhoGauge[dp] === "Replace",
    kg = 2 + (2 - 1) op["nx"];            (* (i,j) = (2,2): just inside a wall *)
    Jrr = Jrr + SparseArray[{{kg, kg} -> 1.}, {n, n}];
    imRho = ReplacePart[imRho, kg -> 0]];

  Jbc = ArrayFlatten[{{Jrr, Jrq1, Jrq2}, {Z, Ibc, Z}, {Z, Z, Ibc}}];

  qb1 = bm pf["q1"][op["xs"], op["ys"]];
  qb2 = bm pf["q2"][op["xs"], op["ys"]];
  rhs = Join[If[bcType === "Dirichlet", bm 1., ConstantArray[0., n]], qb1, qb2];
  If[kg =!= None, rhs[[kg]] = 1.];

  <|"Jbc" -> SparseArray[Jbc], "rhs" -> rhs, "gaugeNode" -> kg,
    "imD" -> Join[imRho, op["imask"], op["imask"]]|>
];

(* -------------------------------------------------------------------- *)
(* Residual and bordered Jacobian.                                       *)
(*                                                                      *)
(*   X  = Join[rho, q1, q2]            (3n)                             *)
(*   ex = {ux, uy}                     ("Dirichlet")                    *)
(*      = {ux, uy, mu}                 ("Neumann"; mu is the multiplier *)
(*                                      of the mean-density row)        *)
(*                                                                      *)
(*   theta = 1 -> implicit Euler step of size dt from Xold              *)
(*   theta = 0 -> steady state (dt, Xold ignored)                       *)
(*                                                                      *)
(*  Interior rows:  theta (X - Xold)/dt - F(X, ex)                      *)
(*  Boundary rows:  Jbc.X - rhs                                          *)
(*  Border rows  :  q1[k0], q2[k0], and (w.rho - area) when "Neumann"    *)
(* -------------------------------------------------------------------- *)
(* Assemble [[Mtop, Cmat], [Esp, 0]] by index offset.
   ArrayFlatten CANNOT be used here: with ragged blocks (3n x 3n beside
   3n x ne) it leaves its sparse path and materialises the whole thing
   densely -- 9.3 GB and 22 s at n = 81^2.  Working from the packed
   NonzeroPositions arrays keeps it O(nnz) and vectorised. *)
(* NB: pos + {di,dj} would thread over ROWS (Plus is Listable over the
   outer dimension), silently producing garbage.  Shift each column. *)
cdShiftPos[pos_, di_Integer, dj_Integer] :=
  If[Length[pos] == 0, {},
     Transpose[{pos[[All, 1]] + di, pos[[All, 2]] + dj}]];

cdBorder[Mtop_SparseArray, Cmat_SparseArray, Esp_SparseArray, ne_Integer] :=
Module[{n3 = Length[Mtop], nt, pos, val},
  nt  = n3 + ne;
  pos = Join[cdShiftPos[Mtop["NonzeroPositions"], 0, 0],
             cdShiftPos[Cmat["NonzeroPositions"], 0, n3],
             cdShiftPos[Esp["NonzeroPositions"],  n3, 0]];
  val = Join[Mtop["NonzeroValues"], Cmat["NonzeroValues"], Esp["NonzeroValues"]];
  SparseArray[pos -> val, {nt, nt}]
];

(* Row scaling for the pinning/mean constraints.
   The constraint rows carry no 1/dt while the differential rows do, so
   the bordered Schur complement behaves like dt*G and the matrix becomes
   singular as dt -> 0 -- the usual index-2 DAE conditioning blow-up.
   (Measured: smallest singular value 2.2e-5 at dt=1e-4, 1.2e-1 at dt=1,
   i.e. exactly proportional to dt.)  Scaling the constraint rows by 1/dt
   is an exact row scaling -- it leaves the Newton step unchanged in exact
   arithmetic -- and makes the conditioning dt-independent. *)
cdBScale[theta_, dt_] := If[TrueQ[theta == 0], 1., 1./dt];

(* --------------------------------------------------------------------
   Gauge row for the sealed ("Neumann") density problem, which is singular
   through the CONSTANT null vector only.  Two rank-1 conditions kill it:

     "Mean"    : Sum_k w_k rho_k = |Omega|   (physical: mass fixed)
     "Point"   : rho(0,0) = 1                (extra row + multiplier mu)
     "Replace" : rho(kg)  = 1  REPLACING the PDE row at one node, no mu

   They differ only by an additive constant in rho, and u is blind to that:
   the Q equation sees rho through grad rho and the rho equation through
   lap rho.  Integrating the rho equation over the sealed box gives
   mu |Omega| = 0 either way.  So both give the same u, Q and grad rho.

   "Mean" puts a DENSE ROW across the rho block (factorisation 1.90 s vs
   0.99 s).  "Point" fixes that but still carries a dense mu COLUMN --
   dF_rho/dmu = 1 on every interior rho row, 5329 entries -- which is what
   actually makes Neumann ~7x slower end to end than Dirichlet at equal
   nnz and equal Newton count (59 iterations either way).

   "Replace" removes both: the gauge condition overwrites the PDE row at a
   single node, so there is no extra row and no multiplier, and the matrix
   is structurally the same as the Dirichlet one.  The sealed system's
   discrete equations are rank-deficient by exactly one (the compatibility
   condition mirroring the constant null vector), so dropping one of them
   loses nothing to leading order.  The node is taken NEAR A WALL, never at
   the core, so the row that is sacrificed is far from where accuracy
   matters.  cdShiftRhoMean restores the mean-1 normalisation afterwards.

   NB this is not the ray pin of iterative_solver_module.m, which fixes
   rho on a whole ray -- many nodes, over-constrained, hence its local
   distortion.  One node with one multiplier is exactly rank-1.
   -------------------------------------------------------------------- *)
cdRhoGauge[dp_Association] := Lookup[dp, "rhoGauge", "Point"];

(* shift rho by a constant so that the weighted mean is 1 *)
cdShiftRhoMean[X_, op_Association] := Module[{n = op["n"], sh},
  sh = 1. - op["w"] . X[[1 ;; n]]/op["area"];
  Join[X[[1 ;; n]] + sh, X[[n + 1 ;; 3 n]]]
];

(* Residual only -- no Jacobian assembly.  Used by cdRates and by the
   convergence check inside cdNewton, where building J is pure waste. *)
cdRes[X_?VectorQ, ex_?VectorQ, Xold_, dt_, theta_,
      op_Association, dp_Association, bc_Association, bcType_String] :=
Module[{n = op["n"], r, q1, q2, ux, uy, mu, wx, wy, q2sum, Fr, F1, F2,
        Bx, lam, kap, Aq, bet, imD, Rb, k0 = op["k0"]},
  r = X[[1 ;; n]]; q1 = X[[n + 1 ;; 2 n]]; q2 = X[[2 n + 1 ;; 3 n]];
  {ux, uy} = ex[[1 ;; 2]];
  mu = If[Length[ex] >= 3, ex[[3]], 0.];
  {Bx, lam, kap, Aq, bet} = Lookup[dp, {"Bx", "lam", "kap", "Aq", "bet"}];
  wx = ux + Bx (op["Dx"] . r); wy = uy + Bx (op["Dy"] . r);
  q2sum = q1^2 + q2^2;
  Fr = Bx (op["Lap"] . r)
       + lam ((op["Lx"] . q1) - (op["Ly"] . q1) + 2 (op["Dxy"] . q2)) + mu
       + If[TrueQ@Lookup[dp, "ugr", False],
            ux (op["Dx"] . r) + uy (op["Dy"] . r), 0.];
  F1 = op["Dx"] . (q1 wx) + op["Dy"] . (q1 wy) + kap (op["Lap"] . q1)
       - Aq q1 - bet q2sum q1;
  F2 = op["Dx"] . (q2 wx) + op["Dy"] . (q2 wy) + kap (op["Lap"] . q2)
       - Aq q2 - bet q2sum q2;
  imD = Lookup[bc, "imD", Join[op["imask"], op["imask"], op["imask"]]];
  Rb = {q1[[k0]], q2[[k0]]};
  If[bcType === "Neumann" && cdRhoGauge[dp] =!= "Replace",
    AppendTo[Rb, If[cdRhoGauge[dp] === "Mean",
                    op["w"] . r - op["area"], r[[k0]] - 1.]]];
  Join[imD (theta (X - Xold)/dt - Join[Fr, F1, F2])
       + bc["Jbc"] . X - bc["rhs"], cdBScale[theta, dt] Rb]
];

cdRJ[X_?VectorQ, ex_?VectorQ, Xold_, dt_, theta_,
     op_Association, dp_Association, bc_Association, bcType_String] :=
Module[{n, r, q1, q2, ux, uy, mu, Dx, Dy, Lx, Ly, Lap, Dxy, im, k0,
        Bx, lam, kap, Aq, bet, wx, wy, q2sum, F1, F2, Fr, Fvec,
        Jrr, Jrq1, Jrq2, J1r, J11, J12, J2r, J21, J22, Jf,
        Mtop, Cmat, Erows, Rtop, Rb, ne, imD, cols, R, M, bsc, ugr, rs},

  n = op["n"]; k0 = op["k0"];
  r  = X[[1 ;; n]]; q1 = X[[n + 1 ;; 2 n]]; q2 = X[[2 n + 1 ;; 3 n]];
  {ux, uy} = ex[[1 ;; 2]];
  ne = Length[ex];
  mu = If[ne >= 3, ex[[3]], 0.];

  Dx = op["Dx"]; Dy = op["Dy"]; Lx = op["Lx"]; Ly = op["Ly"];
  Lap = op["Lap"]; Dxy = op["Dxy"]; im = op["imask"];
  {Bx, lam, kap, Aq, bet} = Lookup[dp, {"Bx", "lam", "kap", "Aq", "bet"}];

  (* ---- right-hand side ------------------------------------------- *)
  wx = ux + Bx (Dx . r);
  wy = uy + Bx (Dy . r);
  q2sum = q1^2 + q2^2;

  (* Optional exact-frame term +u.grad rho in the rho equation.  dean.tex
     drops it as O((grad rho)^2); keeping it is a switch so its size can
     be measured instead of assumed. *)
  ugr = TrueQ@Lookup[dp, "ugr", False];
  Fr = Bx (Lap . r) + lam ((Lx . q1) - (Ly . q1) + 2 (Dxy . q2)) + mu
       + If[ugr, ux (Dx . r) + uy (Dy . r), 0.];
  F1 = Dx . (q1 wx) + Dy . (q1 wy) + kap (Lap . q1) - Aq q1 - bet q2sum q1;
  F2 = Dx . (q2 wx) + Dy . (q2 wy) + kap (Lap . q2) - Aq q2 - bet q2sum q2;
  Fvec = Join[Fr, F1, F2];

  (* ---- dF/dX ------------------------------------------------------ *)
  (* v M  scales the ROWS of M by v;  M.DiagonalMatrix scales COLUMNS.  *)
  Jrr  = Bx Lap + If[ugr, ux Dx + uy Dy, 0];
  Jrq1 = lam (Lx - Ly);
  Jrq2 = 2 lam Dxy;

  J1r  = Bx (Dx . (q1 Dx) + Dy . (q1 Dy));
  J2r  = Bx (Dx . (q2 Dx) + Dy . (q2 Dy));

  J11 = Dx . DiagonalMatrix[SparseArray[wx]] + Dy . DiagonalMatrix[SparseArray[wy]]
        + kap Lap - DiagonalMatrix[SparseArray[Aq + bet (3 q1^2 + q2^2)]];
  J22 = Dx . DiagonalMatrix[SparseArray[wx]] + Dy . DiagonalMatrix[SparseArray[wy]]
        + kap Lap - DiagonalMatrix[SparseArray[Aq + bet (q1^2 + 3 q2^2)]];
  J12 = DiagonalMatrix[SparseArray[-2 bet q1 q2]];
  J21 = J12;

  Jf = ArrayFlatten[{{Jrr, Jrq1, Jrq2}, {J1r, J11, J12}, {J2r, J21, J22}}];

  (* ---- assemble interior + boundary rows --------------------------- *)
  imD = Lookup[bc, "imD", Join[im, im, im]];    (* 1 on differential rows *)
  Rtop = imD (theta (X - Xold)/dt - Fvec) + bc["Jbc"] . X - bc["rhs"];
  Mtop = imD (theta DiagonalMatrix[SparseArray[ConstantArray[1./dt, 3 n]]] - Jf)
         + bc["Jbc"];

  (* ---- border columns  dR/dex ------------------------------------- *)
  cols = {Join[If[ugr, Dx . r, ConstantArray[0., n]], Dx . q1, Dx . q2],
          Join[If[ugr, Dy . r, ConstantArray[0., n]], Dy . q1, Dy . q2]};
  If[ne >= 3,   (* the mu column: dense over the rho block, hence "Replace" *)
    AppendTo[cols, Join[ConstantArray[1., n], ConstantArray[0., 2 n]]]];
  Cmat = Transpose[SparseArray[-imD # & /@ cols]];

  (* ---- border rows ------------------------------------------------ *)
  bsc   = cdBScale[theta, dt];
  Erows = {SparseArray[{n + k0 -> bsc}, {3 n}],
           SparseArray[{2 n + k0 -> bsc}, {3 n}]};
  Rb    = {q1[[k0]], q2[[k0]]};
  If[bcType === "Neumann" && ne >= 3,
    If[cdRhoGauge[dp] === "Mean",
      AppendTo[Erows, SparseArray@Join[bsc op["w"], ConstantArray[0., 2 n]]];
      AppendTo[Rb, op["w"] . r - op["area"]],
    (* else "Point": one entry, keeps the gauge ROW sparse *)
      AppendTo[Erows, SparseArray[{k0 -> bsc}, {3 n}]];
      AppendTo[Rb, r[[k0]] - 1.]]];
  Rb = bsc Rb;

  R = Join[Rtop, Rb];
  M = cdBorder[SparseArray[Mtop], SparseArray[Cmat], SparseArray[Erows], ne];

  (* ---- row equilibration ------------------------------------------- *)
  (* The rho rows carry B/xi0 = 2000 and the Q rows 4K'/xir = 4.15, so the
     blocks differ by ~500x before 1/h^2 multiplies both.  At h = 0.027
     that puts entries at 2.7e6 against 2.9e3 and the factorisation stops
     being useful: measured, the passive solve simply failed to move
     (|R| pinned at its initial 0.174).  Scaling every row by its own
     largest entry is an exact row scaling -- the Newton step is unchanged
     in exact arithmetic -- and it is what makes the fine grids solvable.
     The scale vector is returned so the caller can divide matching
     residuals from cdRes by it. *)
  {R, M, rs} = cdEquilibrate[R, M];
  {R, M, rs}
];

(* row scale factors of a sparse matrix, via the packed CSR arrays *)
cdEquilibrate[R_, M_SparseArray] :=
Module[{rp, vals, nr, rs},
  rp   = M["RowPointers"];
  vals = Abs@M["NonzeroValues"];
  nr   = Length[M];
  rs = Table[
     If[rp[[i + 1]] > rp[[i]], Max[vals[[rp[[i]] + 1 ;; rp[[i + 1]]]]], 1.],
     {i, nr}];
  rs = Map[If[# > 0., #, 1.] &, rs];
  {R/rs, (1./rs) M, rs}
];

(* -------------------------------------------------------------------- *)
(* Newton on the bordered system.  Returns {X, ex, converged, nIter,     *)
(* residualNorm}.                                                        *)
(* -------------------------------------------------------------------- *)
(* Convergence is judged on the Newton STEP, which is scale-free, rather
   than on |R|: the implicit-Euler rows carry a 1/dt factor, so an
   absolute residual tolerance would mean something different at every dt. *)
(* The Jacobian factorisation is by far the most expensive thing here
   (1.3 s vs 0.08 s to assemble at n = 81^2), so it is reused across
   Newton iterations: a modified Newton converges linearly instead of
   quadratically but each extra iteration costs only a back-substitution.
   The residual is always evaluated exactly, so the converged answer is
   unaffected -- only the iteration count is.  We refactor when the step
   stops shrinking by at least "RefactorRatio" per iteration. *)
(* StepTolerance must sit ABOVE the round-off floor of the residual, which
   scales with the matrix entries: with B/xi0 = 2000 and 1/h^2 ~ 180 those
   reach ~4e5, so |R| bottoms out near 1e-10 and the Newton step cannot go
   below ~1e-11 no matter how many iterations are spent.  A tolerance below
   that floor silently burns every iteration.  "StagnationTolerance" is the
   safety net: once the step is small AND has stopped halving, we are at the
   floor and further iterations are waste. *)
Options[cdNewton] = {"MaxIterations" -> 10, "StepTolerance" -> 1.*^-10,
                     "StagnationTolerance" -> 1.*^-8,
                     "Damping" -> 1., "Verbose" -> False,
                     "ReuseFactorization" -> True, "RefactorRatio" -> 0.25,
                     "MaxFactorizations" -> 4};

cdNewton[X0_, ex0_, Xold_, dt_, theta_,
         op_Association, dp_Association, bc_Association, bcType_String,
         opts : OptionsPattern[]] :=
Module[{X = X0, ex = ex0, n = op["n"], ne = Length[ex0], R, M, del, nrm,
        maxIt, stol, stag, damp, verb, reuse, ratio, maxFac, it = 0,
        ok = False, sc, step = Infinity, prev = Infinity, lsf = None,
        nFac = 0, needFac = True, rscale = 1.},
  {maxIt, stol, stag, damp, verb, reuse, ratio, maxFac} = OptionValue[
     {"MaxIterations", "StepTolerance", "StagnationTolerance", "Damping",
      "Verbose", "ReuseFactorization", "RefactorRatio", "MaxFactorizations"}];

  While[it < maxIt,
    If[needFac || !reuse || lsf === None,
      M = cdRJ[X, ex, Xold, dt, theta, op, dp, bc, bcType];
      rscale = M[[3]]; M = M[[2]];
      lsf = Quiet@Check[LinearSolve[M], $Failed];
      If[lsf === $Failed,
        If[verb, Print["      cdNewton: factorisation failed at it ", it]];
        Break[]];
      nFac++; needFac = False];

    (* cdRes is unscaled; apply the same row scaling the matrix carries *)
    R   = cdRes[X, ex, Xold, dt, theta, op, dp, bc, bcType]/rscale;
    del = Quiet@Check[lsf[-R], $Failed];
    If[del === $Failed,
      If[verb, Print["      cdNewton: solve failed at it ", it]]; Break[]];

    X  += damp del[[1 ;; 3 n]];
    ex += damp del[[3 n + 1 ;; 3 n + ne]];
    it++;
    sc   = Max[1., Max@Abs@X];
    step = Max@Abs@del/sc;
    If[verb, Print["      newton ", it, "  |R| = ", Max@Abs@R,
                   "  |step| = ", step, "  fac = ", nFac]];
    If[!(NumericQ[step] && step < 1.*^8), Break[]];         (* diverged *)
    If[step < stol, ok = True; Break[]];
    (* at the round-off floor: small and no longer halving -> done *)
    If[step < stag && step > 0.5 prev, ok = True; Break[]];
    (* stalling on a stale Jacobian -> refresh it *)
    If[reuse && step > ratio prev && nFac < maxFac, needFac = True];
    prev = step;
  ];
  nrm = Max@Abs@cdRes[X, ex, Xold, dt, theta, op, dp, bc, bcType];
  {X, ex, ok, it, nrm, step, nFac}
];

(* -------------------------------------------------------------------- *)
(* Steady-state residual norms of a state, split by field, so that the   *)
(* rho and Q sectors (which carry very different scales) can be judged   *)
(* separately.  These are the |d_t X| of the differential rows.          *)
(* -------------------------------------------------------------------- *)
cdRates[X_, ex_, op_Association, dp_Association, bc_Association,
        bcType_String] :=
Module[{n = op["n"], R},
  R = cdRes[X, ex, X, 1., 0, op, dp, bc, bcType];
  <|"rho" -> Max@Abs@R[[1 ;; n]],
    "Q1"  -> Max@Abs@R[[n + 1 ;; 2 n]],
    "Q2"  -> Max@Abs@R[[2 n + 1 ;; 3 n]],
    "all" -> Max@Abs@R[[1 ;; 3 n]]|>
];

(* -------------------------------------------------------------------- *)
(* PHASE 1: implicit Euler to the comoving steady state.                 *)
(*   - dt grows while Newton is comfortable, is cut on failure           *)
(*   - the pinning constraint is satisfied exactly at every accepted step *)
(* -------------------------------------------------------------------- *)
Options[cdEvolve] = {"dt0" -> 1.*^-4, "dtMax" -> 1., "dtGrow" -> 1.6,
                     "MaxSteps" -> 400, "Tolerance" -> 1.*^-6,
                     "Verbose" -> True, "PrintEvery" -> 10,
                     (* keep the fields every k-th accepted step (0 = off);
                        a snapshot is 3n doubles, ~135 kB at n = 75^2 *)
                     "SnapshotEvery" -> 0,
                     (* called as f[k, t, X, ex, hist] every k-th step, for
                        live checkpointing from the driver *)
                     "Checkpoint" -> None, "CheckpointEvery" -> 20,
                     (* stop at exactly this t, in addition to the usual
                        tolerance/MaxSteps exits.  Infinity (the default) is
                        inert and reproduces the old behaviour bit for bit.
                        comoving_transient_solver.m uses it to march in exact
                        blocks of timeInterval. *)
                     "StopTime" -> Infinity};

cdEvolve[X0_, ex0_, op_Association, dp_Association, bc_Association,
         bcType_String, opts : OptionsPattern[]] :=
Module[{X = X0, ex = ex0, n = op["n"], t = 0., dt, dtNat, dtMax, grow, maxSteps,
        tol, verb, every, hist = {}, Xn, exn, ok, it, nrm, step, nFac, rate,
        k = 0, nFail = 0, rates, G, tm, elapsed = 0., done = False,
        snapEvery, snaps = {}, ckpt, ckptEvery, tStop, hitStop = False},
  {dt, dtMax, grow, maxSteps, tol, verb, every, snapEvery, ckpt, ckptEvery,
   tStop} =
    OptionValue[{"dt0", "dtMax", "dtGrow", "MaxSteps", "Tolerance", "Verbose",
                 "PrintEvery", "SnapshotEvery", "Checkpoint", "CheckpointEvery",
                 "StopTime"}];
  (* dtNat is the step the adaptive controller wants; dt is what is actually
     applied.  They differ only on the step that lands on tStop.  Keeping them
     apart matters: clipping the controller state would make the next block
     restart from the (possibly tiny) clipped value and crawl. *)
  dtNat = dt;

  If[verb, Print["[evolve] bc = ", bcType, "  n = ", n,
                 "  dt0 = ", dt, "  dtMax = ", dtMax,
                 "  maxSteps = ", maxSteps]];

  While[k < maxSteps && !done && !hitStop,
    dt = If[tStop < Infinity, Min[dtNat, tStop - t], dtNat];
    {tm, {Xn, exn, ok, it, nrm, step, nFac}} = AbsoluteTiming@
      cdNewton[X, ex, X, dt, 1, op, dp, bc, bcType, "MaxIterations" -> 8];
    elapsed += tm;
    If[!ok,
      nFail++; dtNat /= 4.;
      If[verb, Print["   step ", k, " rejected (newton it = ", it,
                     "), dt -> ", dtNat]];
      If[dtNat < 1.*^-12, Print["[evolve] dt underflow, aborting"]; Break[]];
      Continue[]];

    rate = Max@Abs[Xn - X]/dt/Max[1., Max@Abs@Xn];
    X = Xn; ex = exn; t += dt; k++;
    (* snap t onto tStop so the caller's block boundaries are exact rather than
       one rounding away.  This is tracked SEPARATELY from done: "converged"
       must keep meaning "the rate fell below Tolerance", or a caller marching
       in blocks would read every block boundary as convergence. *)
    If[tStop < Infinity && t >= tStop - 1.*^-12 Max[1., Abs[tStop]],
      t = tStop; hitStop = True];

    rates = cdRates[X, ex, op, dp, bc, bcType];
    G     = cdCoreG[X, op];
    If[snapEvery > 0 && (Mod[k, snapEvery] == 0 || k == 1),
      AppendTo[snaps, <|"k" -> k, "t" -> t, "X" -> X, "u" -> ex[[1 ;; 2]]|>]];

    AppendTo[hist, <|"k" -> k, "t" -> t, "dt" -> dt, "u" -> ex[[1 ;; 2]],
                     (* ex has length 2 under the "Replace" gauge *)
                     "mu" -> If[Length[ex] >= 3, ex[[3]], 0.],
                     "rate" -> rate, "newton" -> it,
                     "rateQ" -> Max[rates["Q1"], rates["Q2"]],
                     "rateRho" -> rates["rho"],
                     "detG" -> Det[G],
                     "core" -> {X[[n + op["k0"]]], X[[2 n + op["k0"]]]}|>];

    If[verb && (Mod[k, every] == 0 || k == 1),
      Print["   k = ", StringPadLeft[ToString[k], 4],
            "  t = ", StringPadRight[cdFmt[t], 12],
            "  dt = ", StringPadRight[cdFmt[dt], 12],
            "  ux = ", StringPadRight[ToString@NumberForm[ex[[1]], 8], 12],
            "  uy = ", StringPadRight[cdFmt[ex[[2]]], 12],
            "  |dX/dt| = ", StringPadRight[cdFmt[rate], 12],
            "  nwt = ", it,
            "  ", Round[elapsed, 0.1], "s"]];

    If[ckpt =!= None && Mod[k, ckptEvery] == 0,
      ckpt[k, t, X, ex, hist, snaps]];

    If[rate < tol, done = True];
    If[it <= 3, dtNat = Min[dtNat grow, dtMax]];
  ];

  If[verb, Print["[evolve] ", If[done, "converged", "stopped"],
                 " at k = ", k, ", t = ", t,
                 ", u = ", ex[[1 ;; 2]],
                 ", wall = ", Round[elapsed, 0.01], " s",
                 ", rejected = ", nFail]];

  <|"X" -> X, "ex" -> ex, "t" -> t, "converged" -> done,
    "hitStopTime" -> hitStop,
    (* the step the controller wants NEXT, so a caller marching in blocks can
       hand it back as "dt0" instead of restarting the ramp from dt0 each time *)
    "dtNext" -> dtNat,
    "history" -> hist, "snapshots" -> snaps, "steps" -> k, "wall" -> elapsed|>
];

(* -------------------------------------------------------------------- *)
(* PHASE 2: steady state by bordered Newton (theta = 0), seeded from     *)
(* phase 1.  Same matrix, 1/dt -> 0.                                     *)
(* -------------------------------------------------------------------- *)
(* 1e-10, not 1e-12: see the note on cdNewton's StepTolerance -- with
   B/xi0 = 2000 the residual round-off floor is ~1e-10, so a tighter
   tolerance can never be met and every iteration past it is wasted. *)
Options[cdSteady] = {"MaxIterations" -> 25, "StepTolerance" -> 1.*^-10,
                     "Verbose" -> True};

cdSteady[X0_, ex0_, op_Association, dp_Association, bc_Association,
         bcType_String, opts : OptionsPattern[]] :=
Module[{X, ex, ok, it, nrm, step, nFac, tm, rates},
  {tm, {X, ex, ok, it, nrm, step, nFac}} = AbsoluteTiming@
    cdNewton[X0, ex0, X0, 1., 0, op, dp, bc, bcType,
             "MaxIterations" -> OptionValue["MaxIterations"],
             "StepTolerance" -> OptionValue["StepTolerance"],
             "Verbose" -> OptionValue["Verbose"]];
  rates = cdRates[X, ex, op, dp, bc, bcType];
  If[OptionValue["Verbose"],
    Print["[steady] ", If[ok, "converged", "FAILED"], " in ", it,
          " newton steps, wall = ", Round[tm, 0.01], " s"];
    Print["         u = ", ex[[1 ;; 2]],
          "   |R| = ", nrm, "   |step| = ", step];
    Print["         residual by field: ", rates]];
  <|"X" -> X, "ex" -> ex, "converged" -> ok, "iterations" -> it,
    "residual" -> nrm, "rates" -> rates, "wall" -> tm|>
];

(* -------------------------------------------------------------------- *)
(* Diagnostics                                                           *)
(* -------------------------------------------------------------------- *)
(* 2x2 core gradient matrix G = {{d_x Q1, d_y Q1}, {d_x Q2, d_y Q2}}|_0  *)
cdCoreG[X_, op_Association] := Module[{n = op["n"], q1, q2, k0 = op["k0"]},
  q1 = X[[n + 1 ;; 2 n]]; q2 = X[[2 n + 1 ;; 3 n]];
  {{(op["Dx"] . q1)[[k0]], (op["Dy"] . q1)[[k0]]},
   {(op["Dx"] . q2)[[k0]], (op["Dy"] . q2)[[k0]]}}
];

(* The two closed-form core estimates of u, for comparison with the
   converged multiplier.  Both come from d_t Q_i(0,0) = 0 with Q(0) = 0:
     naive  : u = -(B/xi0) grad rho(0)             [dean.tex:709]
     exact  : u = -(B/xi0) grad rho(0)
                 - (4K'/xir) G^-1 . (lap Q1(0), lap Q2(0))
   The correction term is amplified by 4K'/(xir |G|) ~ 190, which is why
   it is a diagnostic here and not the determination. *)
cdCoreEstimates[X_, op_Association, dp_Association] :=
Module[{n = op["n"], k0 = op["k0"], r, q1, q2, Bx, kap, G, gr, lapQ},
  r = X[[1 ;; n]]; q1 = X[[n + 1 ;; 2 n]]; q2 = X[[2 n + 1 ;; 3 n]];
  {Bx, kap} = Lookup[dp, {"Bx", "kap"}];
  G  = cdCoreG[X, op];
  gr = {(op["Dx"] . r)[[k0]], (op["Dy"] . r)[[k0]]};
  lapQ = {(op["Lap"] . q1)[[k0]], (op["Lap"] . q2)[[k0]]};
  <|"naive" -> -Bx gr,
    "exact" -> -Bx gr - kap Inverse[G] . lapQ,
    "gradRho0" -> gr, "lapQ0" -> lapQ, "G" -> G, "detG" -> Det[G]|>
];

(* Mirror symmetry about y = 0: rho, Q1 even and Q2 odd.  Nothing in the
   discretisation imposes this, so it is a free check on the advection
   signs. *)
cdSymmetry[X_, op_Association] :=
Module[{n = op["n"], nx = op["nx"], ny = op["ny"], m, flip},
  flip[v_] := Flatten@Reverse@Partition[v, nx];
  m = Max@Abs@X;
  <|"rhoEven" -> Max@Abs[X[[1 ;; n]] - flip[X[[1 ;; n]]]]/m,
    "Q1even"  -> Max@Abs[X[[n + 1 ;; 2 n]] - flip[X[[n + 1 ;; 2 n]]]]/m,
    "Q2odd"   -> Max@Abs[X[[2 n + 1 ;; 3 n]] + flip[X[[2 n + 1 ;; 3 n]]]]/m|>
];

(* field block as an nx*ny matrix, rows = y, cols = x *)
cdFieldAt[X_, which_String, op_Association] := Module[{n = op["n"], v},
  v = Switch[which,
        "rho", X[[1 ;; n]], "Q1", X[[n + 1 ;; 2 n]], "Q2", X[[2 n + 1 ;; 3 n]]];
  Partition[v, op["nx"]]
];

cdInterp[X_, which_String, op_Association] :=
  ListInterpolation[Transpose@cdFieldAt[X, which, op],
                    {op["gx"], op["gy"]}, InterpolationOrder -> 3];

(* -------------------------------------------------------------------- *)
(* Independent reference: the 1D radial defect profile                   *)
(*     S'' + S'/r - S/r^2 = (1/ld^2)(S^2/S0^2 - 1) S,                    *)
(*     S(0) = 0,  S(R) = S0.                                             *)
(* Plain second-order FD + Newton on a uniform radial grid.  Shares no   *)
(* code with the 2D solver, which is the point -- NDSolve's shooting     *)
(* method fails on this BVP (step size underflows near r ~ 2.6).         *)
(* -------------------------------------------------------------------- *)
cdRadialBVP[ld_, S0_, R_, m_Integer] :=
Module[{h, r, S, D2, D1, Id, res, jac, del, it = 0, w},
  h = R/m;
  r = Table[k h, {k, 0, m}];                      (* r[[1]] = 0 *)
  S = S0 Tanh[r/ld];                              (* seed *)
  S[[1]] = 0.; S[[-1]] = S0;
  w = Range[2, m];                                (* unknown interior indices *)
  D2 = SparseArray[{Band[{1, 1}] -> -2./h^2, Band[{1, 2}] -> 1./h^2,
                    Band[{2, 1}] -> 1./h^2}, {m + 1, m + 1}];
  D1 = SparseArray[{Band[{1, 2}] -> 0.5/h, Band[{2, 1}] -> -0.5/h},
                   {m + 1, m + 1}];
  (* evaluate only at interior nodes: r[[1]] = 0 would divide by zero *)
  res[s_] := (D2 . s)[[w]] + (D1 . s)[[w]]/r[[w]] - s[[w]]/r[[w]]^2
             - (1/ld^2) (s[[w]]^2/S0^2 - 1) s[[w]];
  jac[s_] := (D2[[w, w]]
              + SparseArray[Band[{1, 1}] -> 1./r[[w]]] . D1[[w, w]]
              - SparseArray[Band[{1, 1}] -> 1./r[[w]]^2]
              - SparseArray[Band[{1, 1}] ->
                  (1/ld^2) (3 s[[w]]^2/S0^2 - 1)]);
  While[it < 60,
    del = LinearSolve[jac[S], -res[S]];
    S[[w]] += del; it++;
    If[Max@Abs@del < 1.*^-13 Max[1., S0], Break[]]];
  <|"r" -> r, "S" -> S, "iterations" -> it,
    "residual" -> Max@Abs@res[S],
    "f" -> Interpolation[Transpose@{r, S}, InterpolationOrder -> 3]|>
];
