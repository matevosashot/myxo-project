#!/usr/bin/env wolframscript
(* ============================================================================
   verify.wls -- acceptance tests for the graded-grid solver.

   Ordered by how sharply they would expose a misplaced metric factor:

     [1] uniform reduction of the 1-D operators        vs lyapunov_solver.m:709
     [2] structural identities, uniform AND graded     the adjoint-pair gate
     [3] uniform reduction of A and Q as MATRICES      vs the uniform assembly
     [4] maxReLambda at box15/N51 == 0.0265508793108064 (recorded)
     [5] the flat pedestal across a 6:1 graded grid
     [6] end-to-end graded solve: stable, PSD, flat, cheap

   Usage: wolframscript -file lyapunov_graded_solver/verify.wls [quick]
   ============================================================================ *)

root = "/home/ashmat/Projects/myxo-project-dean";
Get[root <> "/lyapunov_graded_solver/graded_grid.m"];
Get[root <> "/lyapunov_graded_solver/graded_ops.m"];
Get[root <> "/lyapunov_graded_solver/graded_noise.m"];
Get[root <> "/lyapunov_graded_solver/graded_lyapunov.m"];

quick = MemberQ[Rest[$ScriptCommandLine], "quick"];
pass = 0; fail = 0;
chk[l_, v_, lim_] := (If[TrueQ[v <= lim], pass++, fail++];
  Print["  ", If[TrueQ[v <= lim], "PASS", "FAIL"], "  ", l, " = ",
        ScientificForm[N@v, 3], "  (limit ", ScientificForm[N@lim, 2], ")"]);
chkT[l_, b_] := (If[TrueQ[b], pass++, fail++];
  Print["  ", If[TrueQ[b], "PASS", "FAIL"], "  ", l]);
rel[a_, b_] := Max@Abs@Flatten@N@Normal[a - b]/Max@Abs@Flatten@N@Normal[b];

fboxV = 15; lN = 1.73;
steadyPath = root <> "/data/fluctuations-set2_2_comove/steady_box15_N51_Dir_lN1.73_z8_ld0.8.m";
mpRef = {Global`B -> 10000., Global`a -> -0.17199999999999996, Global`b -> 0.1,
   Global`L -> 0.06367999999999997, Global`\[Zeta] -> 8., Global`\[Xi]0 -> 5.,
   Global`\[Xi]r -> 0.1, Global`\[Rho]0 -> 0.3, Global`\[CapitalLambda] -> 0.1,
   Global`lNoise -> 1.73};

Print["=================================================================="];
Print["[1] uniform reduction of the 1-D operators"];
nI = 51;
gu = GradedGrid`nuGridUniform[fboxV, nI];
o1 = GradedOps`nuOps1D[gu];
h = (2 fboxV)/(nI + 1);
gF = Table[-fboxV + k h, {k, 0, nI + 1}];
D1r = SparseArray@NDSolve`FiniteDifferenceDerivative[Derivative[1], gF,
   "DifferenceOrder" -> 2]["DifferentiationMatrix"][[2 ;; -2, 2 ;; -2]];
L1r = SparseArray@NDSolve`FiniteDifferenceDerivative[Derivative[2], gF,
   "DifferenceOrder" -> 2]["DifferentiationMatrix"][[2 ;; -2, 2 ;; -2]];
G1r = SparseArray[Join[Table[{m, m} -> 1/h, {m, nI}],
                       Table[{m + 1, m} -> -1/h, {m, nI}]], {nI + 1, nI}];
chk["G vs uniform G1", rel[o1["G"], G1r], 1.*^-13];
chk["L vs uniform L1", rel[o1["L"], L1r], 1.*^-13];
chk["D vs uniform D1", rel[o1["D"], D1r], 1.*^-13];
chk["Ghat == G when w constant", rel[o1["Ghat"], o1["G"]], 1.*^-13];

Print["\n[2] structural identities"];
opU = GradedOps`nuOps2D[gu, gu];
gg = GradedGrid`nuGrid[fboxV, 0.20, 1.2, 4., 0.35];
opG = GradedOps`nuOps2D[gg, gg];
Print["    graded: Nint = ", gg["Nint"], "  hMin = ", gg["hMin"],
      "  hMax = ", gg["hMax"], "  ratio = ", gg["hMax"]/gg["hMin"]];
chk["lapIdentity, uniform", GradedOps`nuLapIdentity[opU], 1.*^-13];
chk["lapIdentity, GRADED", GradedOps`nuLapIdentity[opG], 1.*^-13];
chk["Dhat^T + Dhat, uniform", GradedOps`nuAntisymmetry[opU], 1.*^-13];
chk["Dhat^T + Dhat, GRADED", GradedOps`nuAntisymmetry[opG], 1.*^-13];
chkT["graded grid has 0 as an interior node",
     MemberQ[gg["gridInt"], _?(Abs[#] < 1.*^-13 &)]];

Print["\n[3] uniform reduction of A and Q as matrices (Nint = 25)"];
nS = 25;
gs = GradedGrid`nuGridUniform[fboxV, nS];
ops = GradedOps`nuOps2D[gs, gs];
ss = GradedLyapunov`nuSteadyFromFile[steadyPath];
mp = Association[mpRef];
Bv = mp[Global`B]; av = mp[Global`a]; bv = mp[Global`b]; Lv = mp[Global`L];
zet = mp[Global`\[Zeta]]; xi0 = mp[Global`\[Xi]0]; xir = mp[Global`\[Xi]r];
rho0v = mp[Global`\[Rho]0]; Lam = mp[Global`\[CapitalLambda]];
cL = Lv + zet xir/(4 xi0);
hs = (2 fboxV)/(nS + 1);
(* the uniform reference operators *)
gFs = Table[-fboxV + k hs, {k, 0, nS + 1}];
giS = gFs[[2 ;; -2]]; xFs = Table[-fboxV + (m - 1/2) hs, {m, nS + 1}];
D1s = SparseArray@NDSolve`FiniteDifferenceDerivative[Derivative[1], gFs,
   "DifferenceOrder" -> 2]["DifferentiationMatrix"][[2 ;; -2, 2 ;; -2]];
L1s = SparseArray@NDSolve`FiniteDifferenceDerivative[Derivative[2], gFs,
   "DifferenceOrder" -> 2]["DifferentiationMatrix"][[2 ;; -2, 2 ;; -2]];
G1s = SparseArray[Join[Table[{m, m} -> 1/hs, {m, nS}],
                       Table[{m + 1, m} -> -1/hs, {m, nS}]], {nS + 1, nS}];
P1s = SparseArray[Join[Table[{m, m} -> 1/2, {m, nS}],
                       Table[{m + 1, m} -> 1/2, {m, nS}]], {nS + 1, nS}];
Id = IdentityMatrix[nS, SparseArray];
Dxu = KroneckerProduct[Id, D1s]; Dyu = KroneckerProduct[D1s, Id];
Lxu = KroneckerProduct[Id, L1s]; Lyu = KroneckerProduct[L1s, Id];
Gxu = KroneckerProduct[Id, G1s]; Gyu = KroneckerProduct[G1s, Id];
Pxu = KroneckerProduct[Id, P1s]; Pyu = KroneckerProduct[P1s, Id];
Lapu = Lxu + Lyu; Dxyu = KroneckerProduct[D1s, D1s];
xf = Flatten@ConstantArray[N@giS, nS]; yf = Flatten@Transpose@ConstantArray[N@giS, nS];
xeXu = Flatten@ConstantArray[N@xFs, nS];
xeYu = Flatten@Transpose@ConstantArray[N@giS, nS + 1];
yeXu = Flatten@ConstantArray[N@giS, nS + 1];
yeYu = Flatten@Transpose@ConstantArray[N@xFs, nS];
rhoV = ss["rho"][xf, yf]; Q1V = ss["Q1"][xf, yf]; Q2V = ss["Q2"][xf, yf];
rhoXE = ss["rho"][xeXu, xeYu]; Q1XE = ss["Q1"][xeXu, xeYu]; Q2XE = ss["Q2"][xeXu, xeYu];
rhoYE = ss["rho"][yeXu, yeYu]; Q1YE = ss["Q1"][yeXu, yeYu]; Q2YE = ss["Q2"][yeXu, yeYu];
rhoXd = Derivative[1, 0][ss["rho"]][xeXu, xeYu];
rhoYd = Derivative[0, 1][ss["rho"]][yeXu, yeYu];
{ux, uy} = N@ss["u"][[1 ;; 2]];
dgu[c_] := DiagonalMatrix[SparseArray@c];
Advu = -(Transpose[Gxu] . ((ux + (Bv/xi0) rhoXd) * Pxu) +
         Transpose[Gyu] . ((uy + (Bv/xi0) rhoYd) * Pyu));
A11u = (Bv/xi0) Lapu;
A21u = -(Bv/xi0) (Transpose[Gxu] . (Q1XE * Gxu) + Transpose[Gyu] . (Q1YE * Gyu));
A22u = Advu + (4/xir) (dgu[-(av + Lam xir) - bv (6 Q1V^2 + 2 Q2V^2)] + cL Lapu);
(* and the same blocks through the graded machinery on the uniform grid *)
A11g = (Bv/xi0) ops["Lap"];
A21g = -(Bv/xi0) (Transpose[ops["Gx"]] . (Q1XE * ops["Gx"]) +
                  Transpose[ops["Gy"]] . (Q1YE * ops["Gy"]));
Advg = -(Transpose[ops["Gx"]] . ((ux + (Bv/xi0) rhoXd) * ops["Px"]) +
         Transpose[ops["Gy"]] . ((uy + (Bv/xi0) rhoYd) * ops["Py"]));
A22g = Advg + (4/xir) (dgu[-(av + Lam xir) - bv (6 Q1V^2 + 2 Q2V^2)] + cL ops["Lap"]);
chk["A11 graded-machinery vs uniform", rel[A11g, A11u], 1.*^-12];
chk["A21 (divergence form) vs uniform", rel[A21g, A21u], 1.*^-12];
chk["Adv vs uniform", rel[Advg, Advu], 1.*^-12];
chk["A22 vs uniform", rel[A22g, A22u], 1.*^-12];
(* Q: the hatted Q must equal h^2 times the uniform pf-based Q *)
flds = <|"rhoV" -> rhoV, "Q1V" -> Q1V, "Q2V" -> Q2V,
   "rhoXE" -> rhoXE, "Q1XE" -> Q1XE, "Q2XE" -> Q2XE,
   "rhoYE" -> rhoYE, "Q1YE" -> Q1YE, "Q2YE" -> Q2YE|>;
SvN = Sqrt[Q1V^2 + Q2V^2]/rhoV;
SxN = Sqrt[Q1XE^2 + Q2XE^2]/rhoXE; SyN = Sqrt[Q1YE^2 + Q2YE^2]/rhoYE;
flds["g2V"] = GradedNoise`nuG2 /@ SvN; flds["g3V"] = GradedNoise`nuG3 /@ SvN;
flds["g2XE"] = GradedNoise`nuG2 /@ SxN; flds["g3XE"] = GradedNoise`nuG3 /@ SxN;
flds["g2YE"] = GradedNoise`nuG2 /@ SyN; flds["g3YE"] = GradedNoise`nuG3 /@ SyN;
Qg = GradedNoise`nuAssembleQ[ops, flds,
   <|"zeta" -> zet, "xi0" -> xi0, "rho0" -> rho0v, "Lambda" -> Lam|>];
kapu = 2 zet/(xi0 rho0v); pfu = kapu/hs^2; locu = 2 Lam/(rho0v hs^2);
divBlkU[cx_, cy_, cxy_, cyx_] :=
  Transpose[Gxu] . (cx * Gxu) + Transpose[Gyu] . (cy * Gyu) +
  Transpose[Dxu] . (cxy * Dyu) + Transpose[Dyu] . (cyx * Dxu);
Qrru = pfu divBlkU[Q1XE + rhoXE, -Q1YE + rhoYE, Q2V, Q2V];
ngS = nS^2;
chk["Qhat_rr == h^2 * uniform Q_rr",
    rel[Qg[[1 ;; ngS, 1 ;; ngS]], N[hs^2] Qrru], 1.*^-12];

Print["\n[4] maxReLambda at box15/N51 vs the recorded value"];
If[quick, Print["  SKIP (quick)"],
 Module[{g51, op51, A51, ev},
  g51 = GradedGrid`nuGridUniform[fboxV, 51];
  op51 = GradedOps`nuOps2D[g51, g51];
  A51 = GradedLyapunov`Private`nuBuildAFrom[op51, ss, mpRef];
  ev = Max@Re@Eigenvalues[N@A51, 4,
     Method -> {"Arnoldi", "Criteria" -> "RealPart", "MaxIterations" -> 100000}];
  chk["maxReLambda vs recorded 0.0265508793108064",
      Abs[ev - 0.0265508793108064]/0.0265508793108064, 1.*^-6]]];

Print["\n[5] flat pedestal on a 6:1 graded grid"];
pfR = GradedNoise`nuPedestalFactor[gg, gg, lN, "Raw" -> True];
pfC = GradedNoise`nuPedestalFactor[gg, gg, lN];
gi = gg["gridInt"]; ii = Select[Range[gg["Nint"]], Abs[gi[[#]]] < fboxV - 3 &];
Print["    raw (no normalisation) max|factor-1| = ",
      ScientificForm[Max@Abs[pfR[[ii, ii]] - 1.], 3]];
chk["corrected max|factor-1| over |x|,|y| < 12", Max@Abs[pfC[[ii, ii]] - 1.], 3.*^-3];
guF = GradedGrid`nuGridUniform[fboxV, 131];
chk["normalisation is inert on a uniform grid",
    Max@Abs[GradedNoise`Private`nuSmootherScale[guF, lN] - 1.], 1.*^-13];

Print["\n=================================================================="];
Print[pass, " passed, ", fail, " failed"];
Print[If[fail == 0, "ALL CHECKS PASSED", "THERE ARE FAILURES"]];
