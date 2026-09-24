(* ::Package:: *)

(* ::Section:: *)
(*Python-accelerated Lyapunov solve kernel*)

(* ---------------------------------------------------------------------------
   WHAT THIS IS

   A drop-in replacement for block [6] of lyapunov_solver.m -- the dense
   O(n^3) Lyapunov solve, which is ~98% of the runtime -- and nothing else.
   Every stencil, the closure, the cutoff congruence, the boundary conditions,
   the extraction and the export stay in lyapunov_solver.m and are shared, not
   copied.  Load order therefore MATTERS:

       Get["lyapunov_solver/lyapunov_solver.m"];            (* first  *)
       Get["lyapunov_solver/lyapunov_solver_optimized.m"];  (* second *)

   because loading the base module resets $LyapunovSolveKernel to the reference
   implementation.  Nothing else in a driver changes: the same
   SolveFluctuationsLyapunov call, the same returned association, the same
   exports, the same tags.

   WHY IT IS FASTER

   The reference kernel diagonalizes A (nonsymmetric, so Complex128 throughout)
   and back-substitutes in the eigenbasis.  This one takes the REAL Schur form
   and solves the triangular Lyapunov equation by a recursive Bartels-Stewart
   whose flops land almost entirely in GEMM.  Measured on an EPYC 9655
   (gaspra02, 192 cores), same node, both with a working BLAS:

       n      Mathematica eigen route     this backend
       2000          6.41 s                  1.71 s       3.7x
       3000         12.10 s                  3.46 s       3.5x
       4000         23.53 s                  5.99 s       3.9x

   Memory matters more than the speed here.  The reference route holds ~15.5
   dense n x n Real64-equivalents at peak, because eigvecs, Pev, lamSum, num,
   num/lamSum and Transpose[Pev] are all Complex128 and all live at once; this
   route never forms a complex matrix.  At Nint = 151 (n = 68403, one dense
   n x n = 37.4 GiB) that is the difference between ~580 GiB and the ~112 GiB
   of A + Q + C -- i.e. between the 2000G queue class, which reserves for days,
   and something that starts now.

   SCOPE: DIRICHLET ONLY

   Under Neumann the drift is singular BY DESIGN (the conserved uniform mode),
   so Bartels-Stewart would need an explicit null-mode deflation in the Schur
   basis that is not implemented.  Neumann therefore delegates to
   LyapunovReferenceKernel -- unchanged, still correct, just not accelerated.
   This is scope, not a failure path; it is announced but never silent.

   FAILURE POLICY: ABORT

   Any backend failure -- interpreter missing, numpy/scipy missing, scratch
   full, a non-finite residual -- aborts the run.  It does NOT silently fall
   back to the reference kernel, so every file in a sweep is known to have come
   from one numerical path.  The abort names the cause and leaves the scratch
   directory in place for inspection.
   --------------------------------------------------------------------------- *)

If[!ValueQ[$LyapunovSolveKernel],
   Print["ERROR: load lyapunov_solver.m BEFORE lyapunov_solver_optimized.m."];
   Abort[]];

Clear[LyapunovPythonKernel];

$LyapunovPython::usage =
  "$LyapunovPython is the python3 interpreter used by LyapunovPythonKernel.  It " <>
  "must have numpy and scipy with a working (threaded) BLAS.  Defaults to the " <>
  "python/3.14.7 module's interpreter, then to any python3 on PATH.  On these " <>
  "nodes do NOT use anaconda3: it resolves to /usr/bin/python3, whose numpy is " <>
  "linked against reference netlib BLAS and runs GEMM at 6.4 GFLOPS against " <>
  "4000 for the module build.";

$LyapunovPythonScript::usage =
  "$LyapunovPythonScript is the path to lyap_solve.py; defaults to the copy " <>
  "sitting next to this file.";

$LyapunovPythonThreads::usage =
  "$LyapunovPythonThreads is the BLAS thread count handed to the backend.  " <>
  "Defaults to SLURM_CPUS_PER_TASK when set, else $ProcessorCount -- the " <>
  "allocation, not the node, is what the job is entitled to.";

$LyapunovKeepScratch::usage =
  "$LyapunovKeepScratch -> True leaves the whole per-call scratch directory " <>
  "in place after a successful solve.  Default False, which removes it: at " <>
  "Nint = 150 the matrices alone are 102 GiB per run, and the directories " <>
  "accumulate one per solve.  Scratch is never deleted after a FAILED solve, " <>
  "so the inputs remain available for inspection.";

$LyapunovPythonRefine::usage =
  "$LyapunovPythonRefine -> False drops the iterative-refinement pass in the " <>
  "backend, saving two n x n buffers at the memory peak.  Default True, and " <>
  "you almost certainly want it: Bartels-Stewart is backward stable, but the " <>
  "FORWARD error scales with the conditioning of the Lyapunov operator, which " <>
  "is poor for this drift.  Measured on the real problem, one pass gives a " <>
  "residual of 5e-11 against 1e-14 refined -- 3-4 orders worse.  A random " <>
  "stable A does NOT show this and must not be used to justify turning it off.";

LyapunovPythonKernel::usage =
  "LyapunovPythonKernel[Adense, Qdense, meta] honours the kernel contract of " <>
  "LyapunovReferenceKernel, solving A C + C A^T + Q = 0 out of core via " <>
  "lyap_solve.py (real Schur + recursive Bartels-Stewart).  Neumann delegates " <>
  "to the reference kernel.";

LyapunovPythonKernel::nopython =
  "No usable python3 found (tried `1`).  Set $LyapunovPython to an interpreter \
with numpy and scipy.";
LyapunovPythonKernel::noscript =
  "The backend script is missing: `1`.  Set $LyapunovPythonScript.";
LyapunovPythonKernel::noscratch =
  "Scratch directory `1` is not writable.  Pass scratchDir -> \"...\" in \
lyapunovParams; on SLURM the node-local NVMe is $TMPDIR.";
LyapunovPythonKernel::pyfail =
  "The Python backend failed (exit code `1`): `2`";
LyapunovPythonKernel::badsize =
  "The backend returned `1` doubles for C, expected `2`.";
LyapunovPythonKernel::space =
  "Scratch `1` has `2` GiB free but this solve needs about `3` GiB (A + Q + C \
at n = `4`).";


Begin["LyapunovSolver`Private`"];

$LyapunovPythonScript = FileNameJoin[{DirectoryName[$InputFileName], "lyap_solve.py"}];

$LyapunovKeepScratch = False;

$LyapunovPythonThreads := With[
  {s = Environment["SLURM_CPUS_PER_TASK"]},
  If[StringQ[s] && StringMatchQ[s, DigitCharacter ..], ToExpression[s], $ProcessorCount]];

$LyapunovPython = Automatic;

(* Running one of the python/* module interpreters by absolute path is NOT enough:
   the interpreter is linked against libpython in its own lib64 and the modulefile
   is what supplies PYTHONHOME, LD_LIBRARY_PATH and PYTHONPATH.  Without them it
   dies with "error while loading shared libraries: libpython3.14.so.1.0".  So
   reconstruct that environment here rather than requiring every driver and
   jobscript to "module load" first -- though if one has, the bare "python3"
   candidate picks it up and inherits a correct environment anyway. *)
moduleSpec[root_String] := Module[{ver},
  ver = StringRiffle[Take[StringSplit[Last @ StringSplit[root, "Python-"], "."], 2], "."];
  <|"path" -> FileNameJoin[{root, "bin", "python3"}],
    "env" -> <|
       "PYTHONHOME" -> root,
       "LD_LIBRARY_PATH" -> StringRiffle[
          {FileNameJoin[{root, "lib"}], FileNameJoin[{root, "lib64"}],
           Replace[Environment["LD_LIBRARY_PATH"], Except[_String] -> Nothing]}, ":"],
       "PYTHONPATH" -> StringRiffle[
          {FileNameJoin[{root, "lib64", "python" <> ver, "lib-dynload"}],
           FileNameJoin[{root, "lib", "python" <> ver, "site-packages"}]}, ":"]|>|>];

fullEnv[over_Association] := Normal @ Join[Association @ GetEnvironment[], over];

(* The python/3.14.7 module first: measured numpy 2.5.2 / scipy 1.18.1 on
   OpenBLAS dispatching to its SkylakeX kernel, ~3900 GFLOPS on Zen5 -- identical
   to a hand-built venv.  Every candidate is probed for numpy AND scipy before
   use, which is also what rejects /usr/bin/python3: its numpy is linked against
   reference netlib BLAS and runs GEMM at 6.4 GFLOPS, ~600x slower. *)
pythonSpec[] := pythonSpec[] = Module[{cands, ok},
  cands = Flatten @ {
     If[StringQ @ Environment["LYAPUNOV_PYTHON"],
        <|"path" -> Environment["LYAPUNOV_PYTHON"], "env" -> <||>|>, Nothing],
     moduleSpec /@ {"/usr/local/python/Python-3.14.7",
                    "/usr/local/python/Python-3.13.9",
                    "/usr/local/python/Python-3.12.11"},
     <|"path" -> "python3", "env" -> <||>|>};
  ok[s_] := TrueQ @ Quiet @ Check[
     0 === RunProcess[{s["path"], "-c", "import numpy, scipy.linalg"}, "ExitCode",
                      ProcessEnvironment -> fullEnv[s["env"]]], False];
  SelectFirst[cands, ok, $Failed]];

(* A user-supplied $LyapunovPython (a bare path) runs with the ambient
   environment: whoever sets it is taking responsibility for it working. *)
resolvedSpec[] := If[StringQ[$LyapunovPython],
   <|"path" -> $LyapunovPython, "env" -> <||>|>, pythonSpec[]];


$LyapunovPythonRefine = True;

(* A is assembled sparse (~13 nonzeros per row).  Writing CSR instead of a dense
   dump costs ~10 MB instead of 34 GiB at Nint = 150 AND spares Mathematica the
   Normal[] copy; lyap_solve.py densifies on its side, where the dense A has to
   exist anyway for LAPACK.  A dense Qtilde is written as-is: the cutoff
   congruence already made it dense in [5b]. *)
writeA[dir_, A_, n_] := Module[{pos, vals, rows, indptr, cols},
  If[! MatrixQ[A, NumericQ] && Head[A] =!= SparseArray,
     Export[FileNameJoin[{dir, "A.bin"}], A, {"Binary", "Real64"}]; Return[Null]];
  If[Head[A] =!= SparseArray,
     Export[FileNameJoin[{dir, "A.bin"}], A, {"Binary", "Real64"}]; Return[Null]];
  pos  = A["NonzeroPositions"];
  vals = A["NonzeroValues"];
  rows = pos[[All, 1]];
  cols = pos[[All, 2]] - 1;                       (* python is 0-based *)
  (* BinCounts over rows -> CSR indptr.  NonzeroPositions is row-major sorted,
     which is exactly CSR order, so no re-sort is needed. *)
  indptr = Prepend[Accumulate @ BinCounts[rows, {1, n + 1, 1}], 0];
  Export[FileNameJoin[{dir, "A_indptr.bin"}], indptr, {"Binary", "Integer64"}];
  Export[FileNameJoin[{dir, "A_indices.bin"}], cols, {"Binary", "Integer64"}];
  Export[FileNameJoin[{dir, "A_data.bin"}], N[vals], {"Binary", "Real64"}];
];


(* --- the kernel ---------------------------------------------------------- *)

LyapunovPythonKernel[Adense_, Qdense_, meta_Association] :=
Module[
  {n, dir, spec, py, script, tim, rc, tWrite, tRun, tRead, proc, res, resJson,
   cFlat, Cmat, eigFlat, eigvals, maxReLam, needGiB, freeGiB, cPath},

  n = meta["n"];

  (* Neumann: scope, not failure.  Announce and hand back to the reference. *)
  If[TrueQ @ meta["neumann"],
     Print["    [python backend] Neumann -> delegating to the reference kernel ",
           "(A is singular by design; the null-mode deflation is not ",
           "implemented in the Schur backend)."];
     Return @ LyapunovReferenceKernel[Adense, Qdense, meta]];

  spec   = resolvedSpec[];
  script = $LyapunovPythonScript;
  If[spec === $Failed,
     Message[LyapunovPythonKernel::nopython,
             "LYAPUNOV_PYTHON, python/3.14.7, python/3.13.9, python/3.12.11, python3"];
     Abort[]];
  py = spec["path"];
  If[!FileExistsQ[script], Message[LyapunovPythonKernel::noscript, script]; Abort[]];

  (* One directory per call, so two solves in a session cannot collide and a
     failed run leaves its inputs behind under a name that identifies it. *)
  dir = FileNameJoin[{meta["scratchDir"],
                      "lyap_" <> ToString[$ProcessID] <> "_" <>
                      ToString @ Round[10^3 AbsoluteTime[]]}];
  Quiet @ CreateDirectory[dir, CreateIntermediateDirectories -> True];
  If[!DirectoryQ[dir], Message[LyapunovPythonKernel::noscratch, meta["scratchDir"]];
     Abort[]];

  (* 3 n^2 Real64 -- A and Q in, C out.  Checked up front: running out of scratch
     three hours in, after the solve, is the expensive way to discover this. *)
  needGiB = N[3 * 8 n^2/2^30];
  (* df, not FileSystemInformation: the latter is version-dependent and a silent
     non-evaluation here would skip the check rather than fail it. *)
  freeGiB = Quiet @ Check[
     N[ToExpression @ StringTrim @ Last @ StringSplit[
          RunProcess[{"df", "-B1", "--output=avail", dir}, "StandardOutput"],
          "\n"] / 2^30],
     Missing[]];
  If[NumericQ[freeGiB] && freeGiB < 1.05 needGiB,
     Message[LyapunovPythonKernel::space, dir,
             ToString@NumberForm[freeGiB, 4], ToString@NumberForm[needGiB, 4], n];
     Abort[]];

  tim = {};
  rc[lbl_, t_] := (AppendTo[tim, {lbl, t}]; t);
  Print["    [python backend] ", py];
  Print["    scratch = ", dir, "   (", ToString@NumberForm[needGiB, 4],
        " GiB of A + Q + C)"];

  (* --- write ------------------------------------------------------------
     Export[..., {"Binary","Real64"}] streams a packed array row-major without
     materialising a Flatten copy -- the same call exportDataComoving already
     uses for covariance_*.bin. *)
  tWrite = rc["[6] write A, Q to scratch", First @ AbsoluteTiming[
     writeA[dir, Adense, n];
     Export[FileNameJoin[{dir, "Q.bin"}], Qdense, {"Binary", "Real64"}];
     Export[FileNameJoin[{dir, "meta.json"}],
            <|"n" -> n, "neumann" -> False,
              (* Refinement stays ON.  It costs two more n x n buffers at the
                 peak, but without it the residual on the REAL drift is 5e-11
                 against 1e-14 with -- Bartels-Stewart bounds the backward error,
                 not the forward one, and this operator is ill conditioned. *)
              "refine" -> TrueQ[$LyapunovPythonRefine],
              "threads" -> $LyapunovPythonThreads|>, "JSON"];]];
  Print["    [w] write A,Q          ", tWrite, " s"];

  (* Everything the solve needs is now on disk.  Drop the caller's A and Q
     BEFORE python allocates, so the two peaks never overlap. *)
  If[MatchQ[meta["release"], _Function], meta["release"][]];

  (* --- solve ------------------------------------------------------------- *)
  tRun = rc["[6] python schur + recursive Bartels-Stewart", First @ AbsoluteTiming[
     proc = RunProcess[{py, script, dir}, All, "",
                       ProcessEnvironment -> fullEnv[spec["env"]]];]];

  If[proc["ExitCode"] =!= 0,
     Message[LyapunovPythonKernel::pyfail, proc["ExitCode"],
             StringTake[proc["StandardError"] <> proc["StandardOutput"],
                        UpTo[4000]]];
     Print["    scratch KEPT for inspection: ", dir];
     Abort[]];

  resJson = FileNameJoin[{dir, "result.json"}];
  res = Quiet @ Check[Association @ Import[resJson, "JSON"], $Failed];
  If[res === $Failed || !TrueQ[res["ok"]],
     Message[LyapunovPythonKernel::pyfail, proc["ExitCode"],
             If[AssociationQ[res], res["error"], "unreadable " <> resJson]];
     Print["    scratch KEPT for inspection: ", dir];
     Abort[]];

  (* --- read -------------------------------------------------------------- *)
  cPath = FileNameJoin[{dir, "C.bin"}];
  tRead = rc["[6] read C from scratch", First @ AbsoluteTiming[
     cFlat = BinaryReadList[cPath, "Real64"];
     If[Length[cFlat] =!= n^2,
        Message[LyapunovPythonKernel::badsize, Length[cFlat], n^2]; Abort[]];
     Cmat = Developer`ToPackedArray @ ArrayReshape[cFlat, {n, n}];
     cFlat =.;
     eigFlat = BinaryReadList[FileNameJoin[{dir, "eig.bin"}], "Real64"];
     eigvals = Complex @@@ Partition[eigFlat, 2];]];
  Print["    [r] read C             ", tRead, " s"];

  maxReLam = Max @ Re @ eigvals;
  Print["    max Re(\[Lambda]) = ", maxReLam,
        If[maxReLam < 0, "  (stable)", "  (UNSTABLE!)"]];
  If[maxReLam >= 0, Message[SolveFluctuationsLyapunov::unstable, maxReLam]];
  Print["    min |\[Lambda]_i+\[Lambda]_j| = ", 2 Min @ Abs @ Re @ eigvals,
        ",  complex fraction = ",
        N[100 Count[eigvals, z_ /; Abs[Im[z]] > 1.*^-8 Abs[z]]/Length[eigvals]], "%"];
  Print["    residual (unrefined)   = ", res["residual"]];
  Print["    residual (refined)     = ", res["residual_refined"],
        "   (improved ",
        ToString @ NumberForm[N[res["residual"]/res["residual_refined"]], 4], "x)"];
  If[res["residual_refined"] > res["residual"],
     Message[SolveFluctuationsLyapunov::norefine,
             res["residual"], res["residual_refined"]]];
  Print["    total solve time = ", tWrite + tRun + tRead, " s",
        "   (python internal: ",
        ToString @ NumberForm[N @ Total @ Values @ Association[res["timings"]], 4],
        " s)"];

  (* Remove the whole per-call directory, not just the big .bin files: eig.bin,
     meta.json and result.json are small but the DIRECTORIES accumulate, one per
     solve, and a sweep makes hundreds.  Guarded on the "lyap_" prefix this
     kernel itself generated, so a mistyped scratchDir can never turn this into
     a recursive delete of something the caller cares about.  On failure nothing
     is deleted -- the abort paths above return before reaching this. *)
  If[!TrueQ[$LyapunovKeepScratch],
     If[StringStartsQ[FileNameTake[dir], "lyap_"] && DirectoryQ[dir],
        Quiet @ DeleteDirectory[dir, DeleteContents -> True],
        Quiet @ DeleteFile @ Select[
           FileNameJoin[{dir, #}] & /@ {"A.bin", "Q.bin", "C.bin", "A_indptr.bin",
                                        "A_indices.bin", "A_data.bin"}, FileExistsQ]]];

  <| "C" -> Cmat,
     "eigenvalues" -> eigvals,
     "maxReLambda" -> maxReLam,
     "residual" -> res["residual"],
     "residualRefined" -> res["residual_refined"],
     "nullIndex" -> None,
     "nullLambdaRel" -> Missing["Dirichlet"],
     "nullLambda2Rel" -> Missing["Dirichlet"],
     "nullOverlap" -> Missing["Dirichlet"],
     "timings" -> tim,
     "backend" -> "python/schur+recursive-bartels-stewart (" <> py <> ")" |>
];

End[];

$LyapunovSolveKernel = LyapunovPythonKernel;
(* Take Atilde (SparseArray) and Qtilde (already dense from [5b]) as assembled:
   block [6] then builds no dense copies at all, saving 2 x 8 n^2. *)
$LyapunovKernelWantsDense = False;

Print["lyapunov_solver_optimized.m loaded: $LyapunovSolveKernel -> ",
      "LyapunovPythonKernel  (Dirichlet accelerated, Neumann delegates)"];
