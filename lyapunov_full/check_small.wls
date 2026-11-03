#!/usr/bin/env wolframscript
(* check_small.wls -- regression of lyapunov_full_solver.m against the paper's module.

   Usage: wolframscript -file check_small.wls [Nint] [milestone] [scratchDir]

   (1) at milestone 1 (paper model, weights 0) the assembled A and the smoothed
       sum of the five noise sources equal the A and W N W that
       SolveFluctuationsLyapunov hands to its kernel;
   (2) lyap_multi.py on these files reproduces the reference kernel's C;
   the scratch dir is left for check_multi.py (sources, pathway).
   With milestone 4 (full model) only (2) is run, against the reference kernel
   applied to the new A. *)

$HistoryLength = 0;
thisDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{thisDir, "lyapunov_full_solver.m"}]];
args = Rest[$ScriptCommandLine];
nInt = If[Length[args] >= 1, ToExpression[args[[1]]], 41];
kMs  = If[Length[args] >= 2, ToExpression[args[[2]]], 1];
scr  = If[Length[args] >= 3, args[[3]],
          FileNameJoin[{$TemporaryDirectory, "lf_check_N" <> ToString[nInt] <> "_m" <> ToString[kMs]}]];

modelParams = {B -> 90000., a -> -20.7, b -> 11.574074074074078, L -> 9.375,
   \[Zeta] -> 12.5, \[Xi]0 -> 5., \[Xi]r -> 15., \[Rho]0 -> 0.25,
   \[CapitalLambda] -> 0.13, lNoise -> 0.5};

Get[FileNameJoin[{lfRepoRoot, "comoving_defect_full", "out", "A1.mx"}]];
ss = lfSteadyFromRun[cfRun, kMs];
Print["steady state: milestone ", kMs, "  w = ", ss["w"], "  u = ", ss["u"]];

asm = lfAssemble[<|"Nint" -> nInt, "fbox" -> 20|>, modelParams, ss];
Print["assembled; timings ", asm["timings"]];
n = 3 nInt^2;
Wb = With[{W2 = KroneckerProduct[asm["Emat"], asm["Emat"]]},
          ArrayFlatten[{{W2, 0, 0}, {0, W2, 0}, {0, 0, W2}}]];
Nsm = Wb . Normal[Total[Values @ asm["N"]]] . Wb;
Ad = Normal @ asm["A"];
Print["Wb, Nsm, Ad ready"];

(* the reference: the module, with a kernel that records what it is handed *)
capture = <||>;
$LyapunovKernelWantsDense = True;
$LyapunovSolveKernel = Function[{A, Q, meta},
   capture = <|"A" -> A, "Q" -> Q|>;
   LyapunovReferenceKernel[A, Q, meta]];
If[kMs == 1,
  Print["calling the module ..."];
  ref = SolveFluctuationsLyapunov[{Nint -> nInt, fbox -> 20}, modelParams, ss, "Dirichlet"];
  Print["T1  |A - A_module|/|A|      = ", Max@Abs[Ad - capture["A"]]/Max@Abs[capture["A"]]];
  Print["T1  |WNW - Q_module|/|Q|    = ", Max@Abs[Nsm - capture["Q"]]/Max@Abs[capture["Q"]]];
  Cref = ref["C"],
  (* full model: the reference kernel on the new drift and the same noise *)
  Cref = LyapunovReferenceKernel[Ad, Nsm,
     <|"neumann" -> False, "bcType" -> "Dirichlet", "nInt" -> nInt, "n" -> n,
       "vCons" -> ConstantArray[0., n]|>]["C"]];

Print["reference done"];
lfWriteScratch[scr, asm, <|"threads" -> 8, "points" -> {{0., 0.}, {-10., 0.}},
                          "saveC" -> True, "refine" -> True, "nModes" -> 6|>];
res = lfRunBackend[scr];
Print["backend: ", KeyDrop[res, {"timings"}]];
Cpy = ArrayReshape[BinaryReadList[FileNameJoin[{scr, "C.bin"}], "Real64"], {n, n}];
Print["T5  |C_py - C_ref|/max|C_ref| = ", Max@Abs[Cpy - Cref]/Max@Abs[Cref]];
Print["T5  rho-rho diagonal rel.err  = ",
      Max@Abs[Diagonal[Cpy][[;; nInt^2]] - Diagonal[Cref][[;; nInt^2]]]/
      Max@Abs[Diagonal[Cref][[;; nInt^2]]]];
Print["max Re lambda = ", res["max_re_lambda"], "   pedestal = ", asm["pedestal"],
      "   Var drho at the core node / pedestal = ",
      Diagonal[Cpy][[(nInt + 1)/2 + ((nInt + 1)/2 - 1) nInt]]/asm["pedestal"]];
Print["scratch: ", scr];
