#!/usr/bin/env wolframscript
(* conditioning.wls -- smallest singular value of the bordered Jacobian on the
   production grid (box 20, dmesh 0.04, hMax 2, refine 1), at the Pade seed,
   by inverse power iteration on (M^T M)^-1.  Compares
     plain    : rows equilibrated only (rho as unknown column, no column scale)
     scaled   : rho columns scaled by xi0/B, then rows equilibrated (the solver)
     noH2line : scaled, but without H2 = 0 on the mirror line
   at weights 0 (paper), {1,1,0,0} and full.  Writes data/conditioning.json.
   Supports the statements of report.tex Sec. 4 (scaling, mirror line). *)
$HistoryLength = 0;
thisDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{ParentDirectory[thisDir], "comoving_defect_full", "comoving_full_core.m"}]];
mp = <|"B" -> 9.*^4, "a" -> -20.7, "b" -> 11.574074074074078, "K" -> 9.375,
       "zeta" -> 12.5, "xi0" -> 5., "xir" -> 15., "rho0" -> 0.25, "Lambda" -> 0.13|>;
dp0 = cfModel[mp, {0, 0, 0, 0}, "exact"];
gg = cdGridPade[20., dp0, 0.04, 2., 1];
op = cfOps[cdMirrorOps[cdOps[gg, gg, 4]]];
n = op["n"];
X = cfInitial[op, dp0]; ex = {0., 0.};
smallest[M_] := Module[{lsf, lsfT, v, w, s},
  lsf = Quiet@LinearSolve[M]; lsfT = Quiet@LinearSolve[Transpose[M]];
  SeedRandom[1]; v = RandomReal[{-1, 1}, Length[M]]; v /= Norm[v];
  Do[w = lsfT[lsf[v]]; s = Norm[w]; v = w/s, {40}];
  {1/Sqrt[s], v}];
share[v_] := Append[Table[Norm[v[[(k - 1) n + 1 ;; k n]]], {k, 5}], Norm[v[[5 n + 1 ;;]]]];
res = Table[
  Module[{dp = cfModel[mp, w, "exact"], out = <|"w" -> w|>, R, M, cs, s, v},
    Do[
      Module[{bc = cfBC[op, dp, "Dirichlet", h2], Mx},
        {R, M} = cfRJraw[X, ex, op, dp, bc];
        If[!h2 && !colsc, Continue[]];
        Mx = If[colsc, cdEquilibrate[R, M . DiagonalMatrix[SparseArray[cfColScale[op, dp]]]][[2]],
                       cdEquilibrate[R, M][[2]]];
        {s, v} = smallest[Mx];
        out[Which[!h2, "noH2line", colsc, "scaled", True, "plain"]] =
          <|"sigmaMin" -> s, "share(rho,Q1,Q2,H1,H2,u)" -> share[v]|>],
      {h2, {True, False}}, {colsc, {False, True}}];
    Print[w, "  ", out];
    out],
  {w, {{0, 0, 0, 0}, {1, 1, 0, 0}, {1, 1, 1, 1}}}];
Export[FileNameJoin[{thisDir, "data", "conditioning.json"}], res, "RawJSON"];
Print["wrote data/conditioning.json"];
