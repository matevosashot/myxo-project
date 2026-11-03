(* ::Package:: *)

(* ::Section:: *)
(*Gaussian density fluctuations of the FULL model about its comoving steady state*)

(* psi = (drho, dQ1, dQ2), psi_t = A psi + xi, <xi xi> = W N W,
   A C + C A^T + W N W = 0.

   This is ../lyapunov_solver/lyapunov_solver.m (Dirichlet reservoir wall only)
   with two extensions; that module is loaded read-only and its helpers are
   called by full name.

   (1) THE DRIFT IS THE LINEARISED FULL MODEL of
       ../comoving_defect_full/comoving_full_core.m, with the same weights
       {wT, wR, wH, wE}.  All weights 0 gives the paper's drift block for block
       (the blocks below marked "paper" are copied from lyapunov_solver.m:1165).
       In complex form, q = Q1 + i Q2, h = H1 + i H2, dbar = dx - i dy,
       C = g2(|q|^2) q^2 (rho = 1 in the moments, as in the steady state):
         rho row  -(2 wH/xi0) Re Div(dq Grad hbar + q Grad dhbar)
         Q rows   (zeta/xi0)(wT/2) dbar^2 dC
                  -(4 wR/xir)(dC hbar + C dhbar)
                  -(wH/xi0)[Lap dh + wE (Div(dC Grad hbar) + Div(C Grad dhbar))]
         dh = -(a + 4 b|q|^2) dq - 2 b q^2 dqbar + K Lap dq     (bare a, b, K)
       The fluctuation wall is the module's psi = 0; dh is formed with the same
       Dirichlet Lap, so the wall also has dH = 0.

   (2) THE NOISE IS SPLIT BY SOURCE.  N is the module's matrix (Dean noise with
       the exact max-entropy moments), written as the sum of five symmetric
       pieces, each assembled exactly as in lyapunov_solver.m:1265-1299:
         "rhoIso"   pf divBlk[rho, rho, 0, 0]         in rho-rho  (pedestal)
         "rhoAniso" pf divBlk[Q1, -Q1, Q2, Q2]        in rho-rho
         "cross"    the rho-Q blocks QQ1r, QQ2r and transposes
         "QQrev"    the pf divBlk parts of the Q-Q blocks (reversal flux noise)
         "QQrot"    the locPf diag parts of the Q-Q blocks (rotational noise)
       The first four come from the reversal noise k, the last from the
       rotational noise lambda.  The smoothing W = Kron[E, E] is NOT applied
       here: the backend lyap_multi.py applies it to the Schur vectors instead,
       so no dense n x n noise matrix is ever formed.

   Main entry points
     lfSteadyFromRun[cfRun, k]           the steady state of milestone k
     lfAssemble[lp, mp, ss]              grid, A, the five N_k, W, diagnostics
     lfWriteScratch[dir, asm, extra]     CSR files + meta.json for lyap_multi.py *)

lfRepoRoot = ParentDirectory[DirectoryName[$InputFileName]];
If[!ValueQ[$LyapunovSolveKernel],
   Get[FileNameJoin[{lfRepoRoot, "lyapunov_solver", "lyapunov_solver.m"}]]];
If[Length[DownValues[cfG2Funs]] == 0,
   Get[FileNameJoin[{lfRepoRoot, "comoving_defect_full", "comoving_full_core.m"}]]];

ClearAll[lfSteadyFromRun, lfGridOps, lfAssemble, lfWriteCSR, lfWriteScratch,
         lfSourceNames, lfReal2];

lfSourceNames = {"rhoIso", "rhoAniso", "cross", "QQrev", "QQrot"};

(* complex linear map dq -> P dq + R dqbar as the real 2x2 block on (dQ1, dQ2) *)
lfReal2[P_, R_] := {{Re[P + R], -Im[P - R]}, {Im[P + R], Re[P - R]}};

(* ------------------------------------------------------------------ *)
(* steady state from an out/<tag>.mx of run_full.wls                   *)
(* ------------------------------------------------------------------ *)
lfSteadyFromRun[run_Association, k_Integer] :=
Module[{ms = run["milestones"][[k]], gx, gy, nx, n, X, field},
  gx = run["mesh"]["gx"]; gy = run["mesh"]["gy"]; nx = run["mesh"]["nx"];
  n = run["mesh"]["n"]; X = ms["X"];
  If[run["stateLayout"] =!= "Join[rho - 1, Q1, Q2, H1, H2]",
     Print["lfSteadyFromRun: unexpected state layout ", run["stateLayout"]]; Abort[]];
  field[j_, off_ : 0.] := ListInterpolation[
     Transpose @ Partition[off + X[[(j - 1) n + 1 ;; j n]], nx], {gx, gy},
     InterpolationOrder -> 3];
  <|"rho" -> field[1, 1.], "Q1" -> field[2], "Q2" -> field[3],
    "H1" -> field[4], "H2" -> field[5],
    "u" -> ms["u"], "w" -> ms["w"], "mesh" -> <|"box" -> run["mesh"]["box"]|>,
    "source" -> <|"tag" -> run["tag"], "milestone" -> k, "g2Mode" -> run["g2Mode"]|>|>
];

(* ------------------------------------------------------------------ *)
(* the module's Dirichlet grid and operators (lyapunov_solver.m:874-950) *)
(* ------------------------------------------------------------------ *)
lfGridOps[nInt_Integer, fboxV_Integer] :=
Module[{h, gridFull, gridInt, xFace, d1Full, l1Full, D1, L1, G1, P1, Idn,
        Dx, Dy, Lx, Ly, Gx, Gy, Px, Py, Lap, Dxy, lapIdentity},
  h = (2 fboxV)/(nInt + 1);
  gridFull = Table[-fboxV + k h, {k, 0, nInt + 1}];
  gridInt  = gridFull[[2 ;; -2]];
  xFace    = Table[-fboxV + (m - 1/2) h, {m, nInt + 1}];
  d1Full = NDSolve`FiniteDifferenceDerivative[Derivative[1], gridFull,
     "DifferenceOrder" -> 2]["DifferentiationMatrix"];
  l1Full = NDSolve`FiniteDifferenceDerivative[Derivative[2], gridFull,
     "DifferenceOrder" -> 2]["DifferentiationMatrix"];
  D1 = SparseArray @ d1Full[[2 ;; -2, 2 ;; -2]];
  L1 = SparseArray @ l1Full[[2 ;; -2, 2 ;; -2]];
  G1 = SparseArray[Join[Table[{m, m} -> 1/h, {m, nInt}],
                        Table[{m + 1, m} -> -1/h, {m, nInt}]], {nInt + 1, nInt}];
  P1 = SparseArray[Join[Table[{m, m} -> 1/2, {m, nInt}],
                        Table[{m + 1, m} -> 1/2, {m, nInt}]], {nInt + 1, nInt}];
  Idn = IdentityMatrix[nInt, SparseArray];
  Dx = KroneckerProduct[Idn, D1];  Dy = KroneckerProduct[D1, Idn];
  Lx = KroneckerProduct[Idn, L1];  Ly = KroneckerProduct[L1, Idn];
  Gx = KroneckerProduct[Idn, G1];  Gy = KroneckerProduct[G1, Idn];
  Px = KroneckerProduct[Idn, P1];  Py = KroneckerProduct[P1, Idn];
  Lap = Lx + Ly;
  Dxy = KroneckerProduct[D1, D1];
  lapIdentity = Max @ Abs @ Normal[Lap + (Transpose[Gx] . Gx + Transpose[Gy] . Gy)];
  If[lapIdentity != 0, Print["lfGridOps: adjoint-pair identity failed"]; Abort[]];
  <|"h" -> h, "gridInt" -> gridInt, "xFace" -> xFace, "nInt" -> nInt, "fbox" -> fboxV,
    "Dx" -> N@Dx, "Dy" -> N@Dy, "Lx" -> N@Lx, "Ly" -> N@Ly, "Gx" -> N@Gx, "Gy" -> N@Gy,
    "Px" -> N@Px, "Py" -> N@Py, "Lap" -> N@Lap, "Dxy" -> N@Dxy,
    "lapIdentity" -> lapIdentity|>
];

(* ------------------------------------------------------------------ *)
(* assembly                                                            *)
(*   lp: <|"Nint", "fbox", "rescaleCutoff" (True)|>                     *)
(*   mp: the module's modelParams rules (Global` B, a, b, L, zeta, ...)  *)
(*   ss: lfSteadyFromRun[...]  (weights taken from ss["w"] unless       *)
(*       lp["w"] is given)                                              *)
(* ------------------------------------------------------------------ *)
lfAssemble[lp_Association, mpRules_List, ss_Association] :=
Module[{mp = Association[mpRules], Bv, av, bv, Lv, zet, xi0, xir, rho0v, Lam, lNoiseV,
        nInt, fboxV, ops, h, gridInt, xFace, nFace, Dx, Dy, Lx, Ly, Gx, Gy, Px, Py,
        Lap, Dxy, w, wT, wR, wH, wE, ux, uy, cL, lEff, Wsym, Emat,
        xflat, yflat, sample, xeX, xeY, yeX, yeY, sampleXE, sampleYE, diag,
        rhoSs, Q1Ss, Q2Ss, H1Ss, H2Ss,
        rhoV, Q1V, Q2V, H1V, H2V, rhoXE, Q1XE, Q2XE, rhoYE, Q1YE, Q2YE,
        rhoXd, rhoYd, wxE, wyE, H1xXE, H2xXE, H1yYE, H2yYE,
        Svals, SXE, SYE, g2V, g3V, g2XE, g3XE, g2YE, g3YE,
        Adv, A11, A12, A13, A21, A22, A23, A31, A32, A33,
        q, qXE, qYE, hb, hbxXE, hbyYE, g2f, g2pf, cfc, Cn, CXE, CYE, Cq, Cqb,
        Ah, Bh, KL, Ph, Rh, Phb, Rhb, divC, divV, dbar2, Pq, Rq, Pr, Rr,
        blkQ, rowR, Atilde,
        kap, pf, locPf, divBlk, h2c, h2s, h3c, h3s, e2, e3,
        h2cXE, h2sXE, h3cXE, h3sXE, e2XE, e3XE, h2cYE, h2sYE, h3cYE, h3sYE, e2YE, e3YE,
        Z, sa, flat3, QQ11, QQ22, QQ12, sym, NrrIso, NrrAni, Ncross, NQQrev, NQQrot, QQ1r, QQ2r,
        minEigLoc, timings = {}, rec, t0},
  rec[lbl_, t_] := (AppendTo[timings, {lbl, t}]; t);
  {Bv, av, bv, Lv, zet, xi0, xir, rho0v, Lam, lNoiseV} = N @ Lookup[mp,
     {Global`B, Global`a, Global`b, Global`L, Global`\[Zeta], Global`\[Xi]0,
      Global`\[Xi]r, Global`\[Rho]0, Global`\[CapitalLambda], Global`lNoise}];
  nInt = lp["Nint"]; fboxV = lp["fbox"];
  w = N @ Lookup[lp, "w", ss["w"]];
  {wT, wR, wH, wE} = w;
  {ux, uy} = N @ ss["u"][[1 ;; 2]];
  cL = Lv + zet xir/(4 xi0);

  rec["[2] grid + operators", First @ AbsoluteTiming[
     ops = lfGridOps[nInt, fboxV]]];
  {h, gridInt, xFace} = Lookup[ops, {"h", "gridInt", "xFace"}];
  nFace = Length[xFace];
  {Dx, Dy, Lx, Ly, Gx, Gy, Px, Py, Lap, Dxy} =
     Lookup[ops, {"Dx", "Dy", "Lx", "Ly", "Gx", "Gy", "Px", "Py", "Lap", "Dxy"}];

  (* calibrated cutoff and the Dirichlet heat kernel E (lyapunov_solver.m:973-1027) *)
  lEff = If[TrueQ @ Lookup[lp, "rescaleCutoff", True],
            LyapunovSolver`Private`cutoffLength[lNoiseV, N[h]], lNoiseV];
  If[lEff === $Failed, Print["lfAssemble: lNoise below the solvability edge"]; Abort[]];
  Wsym = Exp[-(lEff/N[h])^2 Sin[N[Pi] Range[nInt]/(2 (nInt + 1))]^2];
  Emat = With[{V = Sqrt[2./(nInt + 1)] Table[Sin[N[Pi] p i/(nInt + 1)],
                                             {i, nInt}, {p, nInt}]},
              V . (Wsym * Transpose[V])];
  Emat = (Emat + Transpose[Emat])/2;

  (* sampling, as lyapunov_solver.m:1043-1076 *)
  xflat = Flatten @ ConstantArray[N @ gridInt, nInt];
  yflat = Flatten @ Transpose @ ConstantArray[N @ gridInt, nInt];
  sample[f_] := f[xflat, yflat];
  xeX = Flatten @ ConstantArray[N @ xFace, nInt];
  xeY = Flatten @ Transpose @ ConstantArray[N @ gridInt, nFace];
  yeX = Flatten @ ConstantArray[N @ gridInt, nFace];
  yeY = Flatten @ Transpose @ ConstantArray[N @ xFace, nInt];
  sampleXE[f_] := f[xeX, xeY];
  sampleYE[f_] := f[yeX, yeY];
  diag[c_] := DiagonalMatrix[SparseArray @ c];
  {rhoSs, Q1Ss, Q2Ss, H1Ss, H2Ss} = ss /@ {"rho", "Q1", "Q2", "H1", "H2"};

  rec["[3] sample steady state", First @ AbsoluteTiming[
     rhoV = sample[rhoSs]; Q1V = sample[Q1Ss]; Q2V = sample[Q2Ss];
     H1V = sample[H1Ss]; H2V = sample[H2Ss];
     rhoXE = sampleXE[rhoSs]; Q1XE = sampleXE[Q1Ss]; Q2XE = sampleXE[Q2Ss];
     rhoYE = sampleYE[rhoSs]; Q1YE = sampleYE[Q1Ss]; Q2YE = sampleYE[Q2Ss];
     rhoXd = sampleXE[Derivative[1, 0][rhoSs][##] &];
     rhoYd = sampleYE[Derivative[0, 1][rhoSs][##] &];
     H1xXE = sampleXE[Derivative[1, 0][H1Ss][##] &];
     H2xXE = sampleXE[Derivative[1, 0][H2Ss][##] &];
     H1yYE = sampleYE[Derivative[0, 1][H1Ss][##] &];
     H2yYE = sampleYE[Derivative[0, 1][H2Ss][##] &];]];
  wxE = ux + (Bv/xi0) rhoXd;
  wyE = uy + (Bv/xi0) rhoYd;

  (* noise closure: the module's exact moments at S = |Q|/rho *)
  Svals = Sqrt[Q1V^2 + Q2V^2]/rhoV;
  SXE = Sqrt[Q1XE^2 + Q2XE^2]/rhoXE;
  SYE = Sqrt[Q1YE^2 + Q2YE^2]/rhoYE;
  rec["[3] kappa(S) inversion, g2 and g3", First @ AbsoluteTiming[
     {g2V, g3V, g2XE, g3XE, g2YE, g3YE} = {
        LyapunovSolver`Private`g2Exact /@ Svals, LyapunovSolver`Private`g3Exact /@ Svals,
        LyapunovSolver`Private`g2Exact /@ SXE,   LyapunovSolver`Private`g3Exact /@ SXE,
        LyapunovSolver`Private`g2Exact /@ SYE,   LyapunovSolver`Private`g3Exact /@ SYE}]];

  (* every block with background 0.: ArrayFlatten of blocks with different
     backgrounds returns a dense List (and a 0 / 0. mix crashed the kernel) *)
  sa[M_] := SparseArray[M, {nInt^2, nInt^2}, 0.];
  Z = sa[{}];
  flat3[b_] := SparseArray[ArrayFlatten[Map[sa, b, {2}]], {3 nInt^2, 3 nInt^2}, 0.];

  (* --- drift: paper blocks (lyapunov_solver.m:1165-1187, advection on) --- *)
  t0 = AbsoluteTime[];
  Adv = -(Transpose[Gx] . (wxE * Px) + Transpose[Gy] . (wyE * Py));
  A11 = (Bv/xi0) Lap;
  A12 = (zet/xi0) (Lx - Ly);
  A13 = (2 zet/xi0) Dxy;
  A21 = -(Bv/xi0) (Transpose[Gx] . (Q1XE * Gx) + Transpose[Gy] . (Q1YE * Gy));
  A22 = Adv + (4/xir) (diag[-(av + Lam xir) - bv (6 Q1V^2 + 2 Q2V^2)] + cL Lap);
  A23 = -(16 bv/xir) diag[Q1V Q2V];
  A31 = -(Bv/xi0) (Transpose[Gx] . (Q2XE * Gx) + Transpose[Gy] . (Q2YE * Gy));
  A32 = A23;
  A33 = Adv + (4/xir) (diag[-(av + Lam xir) - bv (2 Q1V^2 + 6 Q2V^2)] + cL Lap);

  (* --- drift: the full-model terms, in Wirtinger form ---------------- *)
  {g2f, g2pf} = cfG2Funs[Lookup[ss["source"], "g2Mode", "exact"]];
  cfc[qq_] := Module[{s2 = Re[qq Conjugate[qq]], sc, G, Gp},
     sc = cfClipS2[s2]; G = g2f[sc]; Gp = g2pf[sc] UnitStep[0.99 - s2];
     {G qq^2, qq (Gp s2 + 2 G), Gp qq^3}];
  q = Q1V + I Q2V; qXE = Q1XE + I Q2XE; qYE = Q1YE + I Q2YE;
  hb = H1V - I H2V; hbxXE = H1xXE - I H2xXE; hbyYE = H1yYE - I H2yYE;
  {Cn, Cq, Cqb} = cfc[q];
  CXE = First @ cfc[qXE]; CYE = First @ cfc[qYE];
  (* dh = Ph dq + Rh dqbar, dhbar = conj(Rh) dq + conj(Ph) dqbar *)
  Ah = -av - 4 bv Re[q Conjugate[q]];
  Bh = -2 bv q^2;
  KL = Lv Lap;
  Ph = diag[Ah] + KL;  Rh = diag[Bh];
  Phb = diag[Conjugate[Bh]];  Rhb = diag[Ah] + KL;
  divC[cX_, cY_] := -(Transpose[Gx] . (cX * Gx) + Transpose[Gy] . (cY * Gy));  (* Div(c Grad .) *)
  divV[gX_, gY_] := -(Transpose[Gx] . (gX * Px) + Transpose[Gy] . (gY * Py));  (* Div(. g)    *)
  dbar2 = Lx - Ly - 2 I Dxy;
  Pq = SparseArray[{}, {nInt^2, nInt^2}];  Rq = Pq;
  If[wT != 0.,
     Pq += (zet/xi0) (wT/2) dbar2 . diag[Cq];
     Rq += (zet/xi0) (wT/2) dbar2 . diag[Cqb]];
  If[wR != 0.,
     Pq += -(4 wR/xir) (diag[hb Cq] + diag[Cn] . Phb);
     Rq += -(4 wR/xir) (diag[hb Cqb] + diag[Cn] . Rhb)];
  If[wH != 0.,
     Pq += -(wH/xi0) Lap . Ph;
     Rq += -(wH/xi0) Lap . Rh;
     If[wE != 0.,
        Pq += -(wH wE/xi0) (divV[hbxXE, hbyYE] . diag[Cq] + divC[CXE, CYE] . Phb);
        Rq += -(wH wE/xi0) (divV[hbxXE, hbyYE] . diag[Cqb] + divC[CXE, CYE] . Rhb)]];
  blkQ = lfReal2[Pq, Rq];
  A22 += blkQ[[1, 1]]; A23 += blkQ[[1, 2]]; A32 += blkQ[[2, 1]]; A33 += blkQ[[2, 2]];
  If[wH != 0.,
     Pr = -(2 wH/xi0) (divV[hbxXE, hbyYE] + divC[qXE, qYE] . Phb);
     Rr = -(2 wH/xi0) (divC[qXE, qYE] . Rhb);
     rowR = lfReal2[Pr, Rr][[1]];
     A12 += rowR[[1]]; A13 += rowR[[2]]];
  Atilde = flat3[{{A11, A12, A13}, {A21, A22, A23}, {A31, A32, A33}}];
  If[Head[Atilde] =!= SparseArray, Print["lfAssemble: A is not sparse"]; Abort[]];
  rec["[4] assemble A", AbsoluteTime[] - t0];

  (* --- noise, by source (lyapunov_solver.m:1240-1299) ----------------- *)
  t0 = AbsoluteTime[];
  kap = 2 zet/(xi0 rho0v); pf = kap/h^2; locPf = 2 Lam/(rho0v h^2);
  divBlk[cxxF_, cyyF_, cxy_, cyx_] :=
    Transpose[Gx] . (cxxF * Gx) + Transpose[Gy] . (cyyF * Gy) +
    Transpose[Dx] . (cxy * Dy) + Transpose[Dy] . (cyx * Dx);
  h2c = Q1V^2 - Q2V^2;        h2s = 2 Q1V Q2V;
  h3c = Q1V^3 - 3 Q1V Q2V^2;  h3s = 3 Q1V^2 Q2V - Q2V^3;
  e2 = g2V/(2 rhoV);          e3 = g3V/(4 rhoV^2);
  h2cXE = Q1XE^2 - Q2XE^2;        h2sXE = 2 Q1XE Q2XE;
  h3cXE = Q1XE^3 - 3 Q1XE Q2XE^2; h3sXE = 3 Q1XE^2 Q2XE - Q2XE^3;
  e2XE = g2XE/(2 rhoXE);          e3XE = g3XE/(4 rhoXE^2);
  h2cYE = Q1YE^2 - Q2YE^2;        h2sYE = 2 Q1YE Q2YE;
  h3cYE = Q1YE^3 - 3 Q1YE Q2YE^2; h3sYE = 3 Q1YE^2 Q2YE - Q2YE^3;
  e2YE = g2YE/(2 rhoYE);          e3YE = g3YE/(4 rhoYE^2);
  (* every block with background 0.: ArrayFlatten of blocks with different
     backgrounds returns a dense List (and a 0 / 0. mix crashed the kernel) *)
  sym[M_] := SparseArray[(M + Transpose[M])/2, Dimensions[M], 0.];
  QQ1r = pf divBlk[Q1XE + rhoXE/2 + e2XE h2cXE, Q1YE - rhoYE/2 - e2YE h2cYE,
                   e2 h2s, e2 h2s];
  QQ2r = pf divBlk[Q2XE + e2XE h2sXE, Q2YE - e2YE h2sYE,
                   rhoV/2 - e2 h2c, rhoV/2 - e2 h2c];
  QQ11 = pf divBlk[3 Q1XE/4 + rhoXE/2 + e2XE h2cXE + e3XE h3cXE,
                   -3 Q1YE/4 + rhoYE/2 + e2YE h2cYE - e3YE h3cYE,
                   Q2V/4 + e3 h3s, Q2V/4 + e3 h3s];
  QQ22 = pf divBlk[Q1XE/4 + rhoXE/2 - e2XE h2cXE - e3XE h3cXE,
                   -Q1YE/4 + rhoYE/2 - e2YE h2cYE + e3YE h3cYE,
                   3 Q2V/4 - e3 h3s, 3 Q2V/4 - e3 h3s];
  QQ12 = pf divBlk[Q2XE/4 + e2XE h2sXE + e3XE h3sXE,
                   -Q2YE/4 + e2YE h2sYE - e3YE h3sYE,
                   Q1V/4 - e3 h3c, Q1V/4 - e3 h3c];
  NrrIso = sym @ flat3[{{pf divBlk[rhoXE, rhoYE, 0. Q2V, 0. Q2V], Z, Z},
                        {Z, Z, Z}, {Z, Z, Z}}];
  NrrAni = sym @ flat3[{{pf divBlk[Q1XE, -Q1YE, Q2V, Q2V], Z, Z},
                        {Z, Z, Z}, {Z, Z, Z}}];
  Ncross = sym @ flat3[{{Z, Transpose[QQ1r], Transpose[QQ2r]},
                        {QQ1r, Z, Z}, {QQ2r, Z, Z}}];
  NQQrev = sym @ flat3[{{Z, Z, Z}, {Z, QQ11, QQ12}, {Z, Transpose[QQ12], QQ22}}];
  NQQrot = sym @ flat3[{{Z, Z, Z},
     {Z, locPf diag[2 rhoV - 4 e2 h2c], locPf diag[-4 e2 h2s]},
     {Z, locPf diag[-4 e2 h2s], locPf diag[2 rhoV + 4 e2 h2c]}}];
  If[!AllTrue[{NrrIso, NrrAni, Ncross, NQQrev, NQQrot}, Head[#] === SparseArray &],
     Print["lfAssemble: a noise source is not sparse"]; Abort[]];
  rec["[5] assemble the five noise sources", AbsoluteTime[] - t0];

  (* condition (ii), local: the module's 6x6 flux covariance at every node *)
  rec["[5] condition (ii), local 6x6", First @ AbsoluteTiming[
     minEigLoc = Min /@ Eigenvalues /@ MapThread[
        LyapunovSolver`Private`localFluxCovariance, {rhoV, Q1V, Q2V, g2V, g3V}]]];

  <|"A" -> Atilde,
    "N" -> AssociationThread[lfSourceNames -> {NrrIso, NrrAni, Ncross, NQQrev, NQQrot}],
    "Emat" -> Emat, "Wsym" -> Wsym, "lEff" -> lEff, "lNoise" -> lNoiseV,
    "h" -> N[h], "gridInt" -> N[gridInt], "nInt" -> nInt, "fbox" -> fboxV,
    "w" -> w, "u" -> {ux, uy}, "gamma0" -> Bv/xi0,
    "fields" -> <|"rho" -> rhoV, "Q1" -> Q1V, "Q2" -> Q2V, "H1" -> H1V, "H2" -> H2V,
                  "g2" -> g2V, "g3" -> g3V|>,
    "pedestal" -> zet/(2 Pi Bv rho0v lNoiseV^2),
    "minEigLocal" -> MinMax[minEigLoc],
    "peclet" -> N[Max[Max@Abs@wxE, Max@Abs@wyE] h/(2 (4 cL/xir))],
    "lapIdentity" -> ops["lapIdentity"],
    "params" -> <|"B" -> Bv, "a" -> av, "b" -> bv, "K" -> Lv, "zeta" -> zet,
                  "xi0" -> xi0, "xir" -> xir, "rho0" -> rho0v, "Lambda" -> Lam,
                  "lNoise" -> lNoiseV|>,
    "ops" -> ops, "timings" -> timings|>
];

(* ------------------------------------------------------------------ *)
(* CSR export (0-based), as lyapunov_solver_optimized.m:writeA          *)
(* ------------------------------------------------------------------ *)
lfWriteCSR[dir_, name_, M_SparseArray] := Module[{n = Length[M], pos, vals, indptr},
  pos  = M["NonzeroPositions"];
  vals = M["NonzeroValues"];
  indptr = Prepend[Accumulate @ BinCounts[pos[[All, 1]], {1, n + 1, 1}], 0];
  BinaryWrite[FileNameJoin[{dir, name <> "_indptr.bin"}], indptr, "Integer64"];
  BinaryWrite[FileNameJoin[{dir, name <> "_indices.bin"}], pos[[All, 2]] - 1, "Integer64"];
  BinaryWrite[FileNameJoin[{dir, name <> "_data.bin"}], N[vals], "Real64"];
  Close /@ (FileNameJoin[{dir, name <> #}] & /@
            {"_indptr.bin", "_indices.bin", "_data.bin"});
];

(* extra: <|"threads", "points" -> {{x,y},..}, "saveC" -> bool, "refine" -> bool|> *)
lfWriteScratch[dir_String, asm_Association, extra_Association] := Module[{meta},
  Quiet @ CreateDirectory[dir, CreateIntermediateDirectories -> True];
  lfWriteCSR[dir, "A", asm["A"]];
  KeyValueMap[lfWriteCSR[dir, "N_" <> #1, #2] &, asm["N"]];
  BinaryWrite[FileNameJoin[{dir, "E.bin"}], Flatten @ asm["Emat"], "Real64"];
  Close[FileNameJoin[{dir, "E.bin"}]];
  BinaryWrite[FileNameJoin[{dir, "Wsym.bin"}], asm["Wsym"], "Real64"];
  Close[FileNameJoin[{dir, "Wsym.bin"}]];
  meta = Join[<|"n" -> 3 asm["nInt"]^2, "nInt" -> asm["nInt"], "h" -> asm["h"],
     "fbox" -> asm["fbox"], "gamma0" -> asm["gamma0"], "pedestal" -> asm["pedestal"],
     "sources" -> lfSourceNames, "gridInt" -> asm["gridInt"]|>, extra];
  Export[FileNameJoin[{dir, "meta.json"}], meta, "RawJSON"];
  dir
];

(* ------------------------------------------------------------------ *)
(* run lyap_multi.py on a scratch dir; the interpreter (python/3.14.7  *)
(* module build, OpenBLAS) is the one lyapunov_solver_optimized.m picks *)
(* ------------------------------------------------------------------ *)
If[Length[DownValues[LyapunovSolver`Private`resolvedSpec]] == 0,
   Module[{k = $LyapunovSolveKernel, w = $LyapunovKernelWantsDense},
     Get[FileNameJoin[{lfRepoRoot, "lyapunov_solver", "lyapunov_solver_optimized.m"}]];
     (* loading it switches the module's kernel; this file does not use the module's
        solve, so restore whatever was set *)
     $LyapunovSolveKernel = k; $LyapunovKernelWantsDense = w]];

lfScript = FileNameJoin[{DirectoryName[$InputFileName], "lyap_multi.py"}];

lfRunBackend[dir_String] := Module[{spec, res, t},
  spec = LyapunovSolver`Private`resolvedSpec[];
  If[spec === $Failed, Print["lfRunBackend: no usable python3"]; Abort[]];
  Print["    python: ", spec["path"]];
  {t, res} = AbsoluteTiming @ RunProcess[{spec["path"], lfScript, dir},
     ProcessEnvironment -> LyapunovSolver`Private`fullEnv[spec["env"]]];
  If[StringLength[res["StandardError"]] > 0, Print[res["StandardError"]]];
  Print["    backend exit ", res["ExitCode"], " after ", Round[t, 0.1], " s"];
  Import[FileNameJoin[{dir, "result.json"}], "RawJSON"]
];
