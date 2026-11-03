#!/usr/bin/env wolframscript
(* run_fluct.wls -- equal-time fluctuations about one steady state of the full-model
   ladder, with the noise split by source (lyapunov_full_solver.m, lyap_multi.py).

   Usage:
     wolframscript -file run_fluct.wls <Nint> <milestone> [runTag=A1] [saveC=0] [outRoot]
       milestone k of ../comoving_defect_full/out/<runTag>.mx; for A1:
         1 paper (0,0,0,0)   2 +wT (1,0,0,0)   3 +wR (1,1,0,0)   4 full (1,1,1,1)
   Parameters: those of the paper's fluctuation run
     data/fluctuations-set2_2_comove_optimized_2_ssDir/
       covariance_box20_N180_Dir_lN0.5_z12.5_ld1_B90000_ss20_meta.m
   (B = 9e4, lNoise = 0.5 calibrated, Dirichlet reservoir wall); fbox is the steady
   state's own box (20 for A1, 10/15/30/40 for L10..L40).
   Writes <outRoot>/<tag>/{out.npz, result.json, steady.json[, C.bin]},
   tag = box<fbox>_N<Nint>_lN0.5_B90000_<runTag>m<k>_w<wT><wR><wH><wE>. *)

$HistoryLength = 0;
thisDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{thisDir, "lyapunov_full_solver.m"}]];
args = Rest[$ScriptCommandLine];
If[Length[args] < 2, Print["usage: run_fluct.wls Nint milestone [runTag] [saveC] [outRoot]"]; Exit[1]];
nInt  = ToExpression[args[[1]]];
kMs   = ToExpression[args[[2]]];
ssTag = If[Length[args] >= 3, args[[3]], "A1"];
saveC = Length[args] >= 4 && args[[4]] === "1";
outRoot = If[Length[args] >= 5, args[[5]],
             FileNameJoin[{lfRepoRoot, "data", "fluctuations-full"}]];

modelParams = {B -> 90000., a -> -20.7, b -> 11.574074074074078, L -> 9.375,
   \[Zeta] -> 12.5, \[Xi]0 -> 5., \[Xi]r -> 15., \[Rho]0 -> 0.25,
   \[CapitalLambda] -> 0.13, lNoise -> 0.5};

Get[FileNameJoin[{lfRepoRoot, "comoving_defect_full", "out", ssTag <> ".mx"}]];
ss = lfSteadyFromRun[cfRun, kMs];
fboxV = Round[ss["mesh"]["box"]];
wTag = StringJoin[ToString /@ Round[ss["w"]]];
tag = "box" <> ToString[fboxV] <> "_N" <> ToString[nInt] <> "_lN0.5_B90000_" <> ssTag <> "m" <> ToString[kMs] <>
      "_w" <> wTag;
outDir = FileNameJoin[{outRoot, tag}];
Quiet @ CreateDirectory[outDir, CreateIntermediateDirectories -> True];
Print["tag = ", tag, "   host = ", $MachineName, "   w = ", ss["w"], "   u = ", ss["u"]];

{tAsm, asm} = AbsoluteTiming @ lfAssemble[<|"Nint" -> nInt, "fbox" -> fboxV|>, modelParams, ss];
Print["assembled in ", Round[tAsm, 0.1], " s: h = ", asm["h"], ", lEff = ", asm["lEff"],
      ", pedestal = ", asm["pedestal"], ", Peclet = ", asm["peclet"],
      ", local 6x6 min eig in ", asm["minEigLocal"]];
Print["timings ", asm["timings"]];

scr = SelectFirst[{Environment["TMPDIR"],
                   "/scratch/" <> ToString @ Environment["SLURM_JOB_ID"], "/scratch",
                   $TemporaryDirectory}, StringQ[#] && DirectoryQ[#] &];
scr = FileNameJoin[{scr, "lf_" <> tag <> "_" <> ToString[$ProcessID]}];
threads = With[{s = Environment["SLURM_CPUS_PER_TASK"]},
               If[StringQ[s], ToExpression[s], $ProcessorCount]];
hh = asm["h"]/2;
lfWriteScratch[scr, asm, <|"threads" -> threads, "saveC" -> saveC, "refine" -> True,
   "nModes" -> 8,
   "points" -> {{-hh, -hh}, {hh, -hh}, {-hh, hh}, {hh, hh},
                {-5., 0.}, {5., 0.}, {0., 5.}, {-10., 0.}, {-1.5, 1.5}}|>];
Print["scratch: ", scr];
res = lfRunBackend[scr];
Print["backend: ", KeyDrop[res, {"traceback"}]];
If[!TrueQ[res["ok"]], Print["FAILED; scratch kept"]; Exit[1]];

Do[CopyFile[FileNameJoin[{scr, f}], FileNameJoin[{outDir, f}], OverwriteTarget -> True],
   {f, Join[{"out.npz", "result.json", "meta.json"}, If[saveC, {"C.bin"}, {}]]}];
Export[FileNameJoin[{outDir, "steady.json"}],
  <|"tag" -> tag, "source" -> ss["source"], "w" -> ss["w"], "u" -> ss["u"],
    "gridInt" -> asm["gridInt"], "h" -> asm["h"], "nInt" -> nInt, "fbox" -> fboxV,
    "lNoise" -> asm["lNoise"], "lEff" -> asm["lEff"], "pedestal" -> asm["pedestal"],
    "peclet" -> asm["peclet"], "minEigLocal" -> asm["minEigLocal"],
    "params" -> asm["params"], "fields" -> asm["fields"],
    "assemblyTimings" -> asm["timings"], "host" -> $MachineName,
    "date" -> DateString["ISODateTime"]|>, "RawJSON"];
DeleteDirectory[scr, DeleteContents -> True];
Print["wrote ", outDir];
