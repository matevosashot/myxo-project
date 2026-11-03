#!/usr/bin/env wolframscript
(* check_linearization.wls -- the full-model drift terms of lfAssemble against the
   continuum Gateaux derivative, computed symbolically.

   Base state: a smooth analytic texture q0 (a +1 vortex with |q0| < 0.9) and its
   molecular field h0 = -a q0 - 2b|q0|^2 q0 + K Lap q0; the base need not be a
   steady state.  Direction v: two shifted Gaussians, negligible at the wall.
   For each weight set, (A(w) - A(0)).v must converge to dF[q0; v] at O(h^2);
   printed for Nint = 41 and 81 on the interior |x|, |y| < 12.
   Usage: wolframscript -file check_linearization.wls *)

$HistoryLength = 0;
thisDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{thisDir, "lyapunov_full_solver.m"}]];

modelParams = {B -> 90000., a -> -20.7, b -> 11.574074074074078, L -> 9.375,
   \[Zeta] -> 12.5, \[Xi]0 -> 5., \[Xi]r -> 15., \[Rho]0 -> 0.25,
   \[CapitalLambda] -> 0.13, lNoise -> 0.5};
{av, bv, Kv, zet, xi0, xir} = {-20.7, 11.574074074074078, 9.375, 12.5, 5., 15.};
{g2f, g2pf} = cfG2Funs["exact"];

Clear[x, y, eps];
q10 = 0.9 x/Sqrt[1 + x^2 + y^2] + 0.1 Exp[-(x^2 + y^2)/30];
q20 = 0.9 y/Sqrt[1 + x^2 + y^2];
v1 = Exp[-((x - 1)^2 + (y + 0.5)^2)/4] (1 + 0.5 x);
v2 = Exp[-((x + 0.7)^2 + (y - 1)^2)/4] (0.3 - 0.5 y);
fn[e_] := Function @@ {{x, y}, e + 0. x};   (* + 0. x: a constant still returns a vector *)
lap[f_] := D[f, {x, 2}] + D[f, {y, 2}];
hOf[q_, s2_] := -av q - 2 bv s2 q + Kv lap[q];

F[w_] := Module[{Q1 = q10 + eps v1, Q2 = q20 + eps v2, q, qb, s2, h, hb, C, Fq, Fr},
  q = Q1 + I Q2; qb = Q1 - I Q2; s2 = Q1^2 + Q2^2;
  h = hOf[q, s2]; hb = hOf[qb, s2];
  C = g2f[s2] q^2;
  Fq = (zet/xi0) (w[[1]]/2) (D[C, {x, 2}] - D[C, {y, 2}] - 2 I D[C, x, y]) -
       (4 w[[2]]/xir) C hb -
       (w[[3]]/xi0) (lap[h] + w[[4]] (D[C D[hb, x], x] + D[C D[hb, y], y]));
  Fr = -(2 w[[3]]/xi0) (1/2) (D[q D[hb, x], x] + D[q D[hb, y], y] +
                              D[qb D[h, x], x] + D[qb D[h, y], y]);
  {Fr, Fq}];
dF[w_] := D[F[w], eps] /. eps -> 0;

h10 = ComplexExpand[Re[hOf[q10 + I q20, q10^2 + q20^2]]];
h20 = ComplexExpand[Im[hOf[q10 + I q20, q10^2 + q20^2]]];
ss = <|"rho" -> Function[{x, y}, 1. + 0. x], "Q1" -> Function @@ {{x, y}, q10},
       "Q2" -> Function @@ {{x, y}, q20}, "H1" -> Function @@ {{x, y}, h10},
       "H2" -> Function @@ {{x, y}, h20},
       "u" -> {0.1, 0.}, "w" -> {0., 0., 0., 0.}, "mesh" -> <|"box" -> 20.|>,
       "source" -> <|"g2Mode" -> "exact"|>|>;

wSets = If[Length[Rest[$ScriptCommandLine]] > 0, ToExpression[$ScriptCommandLine[[2]]],
   {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 1, 1}, {1, 1, 1, 1}}];
Do[
  Module[{asm0, ops, g, xs, ys, m, vvec, mask, ref, num, Fr, Fq, err},
    asm0 = lfAssemble[<|"Nint" -> nI, "fbox" -> 20, "w" -> {0., 0., 0., 0.}|>, modelParams, ss];
    g = asm0["gridInt"]; m = nI^2;
    xs = Flatten @ ConstantArray[g, nI]; ys = Flatten @ Transpose @ ConstantArray[g, nI];
    vvec = Join[ConstantArray[0., m], fn[v1][xs, ys],
                fn[v2][xs, ys]];
    mask = Flatten @ Table[If[Abs[g[[i]]] < 12 && Abs[g[[j]]] < 12, 1, 0], {j, nI}, {i, nI}];
    Do[
      Module[{asm = lfAssemble[<|"Nint" -> nI, "fbox" -> 20, "w" -> N[w]|>, modelParams, ss]},
        num = (asm["A"] - asm0["A"]) . vvec;
        {Fr, Fq} = dF[w];
        ref = With[{zq = fn[Fq][xs, ys], zr = fn[Fr][xs, ys]},
           Join[If[ListQ[zr], Re[zr], ConstantArray[0., m]], Re[zq], Im[zq]]];
        err = Table[With[{r = ref[[(c - 1) m + 1 ;; c m]], z = num[[(c - 1) m + 1 ;; c m]]},
                If[Max@Abs@Pick[r, mask, 1] < 1.*^-12, Max@Abs@Pick[z, mask, 1],
                   Max@Abs@Pick[z - r, mask, 1]/Max@Abs@Pick[r, mask, 1]]], {c, 3}];
        Print["Nint ", nI, "  w = ", w, "   rel.err (rho, Q1, Q2) = ", err]],
      {w, wSets}]],
  {nI, {41, 81, 161}}];
