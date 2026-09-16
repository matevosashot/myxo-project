(* ::Package:: *)

covStem[file_String] := If[StringEndsQ[file, ".bin"], StringDrop[file, -4], file];

readCovMeta[file_String] := Import[covStem[file] <> "_meta.m"];

(* Sigma[[k, k + shift]] for k = k1..k2.  shift = 0 is the main diagonal (the
   per-block variances); shift = ng picks the same-point cross-covariance
   between block b and block b+1, since blocks are ng = Nint^2 wide. *)
covDiag[stem_String, n_Integer, k1_Integer, k2_Integer, shift_Integer : 0] :=
Module[{str = OpenRead[stem <> ".bin", BinaryFormat -> True]},
  WithCleanup[
     Table[SetStreamPosition[str, 8 ((k - 1) (n + 1) + shift)];
           BinaryRead[str, "Real64"], {k, k1, k2}],
     Close[str]]];

(* Full block (a,b) of Sigma: rows (a-1)ng+1..a ng, cols (b-1)ng+1..b ng.
   Each block row is ng CONTIGUOUS Real64, so this is ng seeks + ng bulk reads,
   not the ng^2 seeks a naive extension of covDiag would do.
   ToPackedArray is not cosmetic: BinaryReadList returns a packed row, but Table
   wraps them in an unpacked outer list and ArrayReshape inherits that, costing
   3.3x the memory (7.6 vs 2.3 GB per block at Nint = 130) and making corrArr
   ~36x slower.  Packing costs 0.2 ms. *)
covBlock[stem_String, n_Integer, ng_Integer, a_Integer, b_Integer] :=
Module[{str = OpenRead[stem <> ".bin", BinaryFormat -> True]},
  WithCleanup[
     Developer`ToPackedArray @
        Table[SetStreamPosition[str, 8 (((a - 1) ng + k - 1) n + (b - 1) ng)];
              BinaryReadList[str, "Real64", ng], {k, ng}],
     Close[str]]];

covMat[vec_, Nint_] := Transpose[Partition[vec, Nint]];

covFn[mat_, grid_] :=
   ListInterpolation[mat, {grid, grid}, InterpolationOrder -> 2];

(* B[[k,k']] with k = i + (j-1) Nint  ->  A[[i,j,i',j']].  ArrayReshape splits
   each flat index slow-digit-first, giving T[[j,i,j',i']]; {2,1,4,3} swaps each
   pair back to x-first.  The 4-D analogue of covMat. *)
corrArr[b_, Nint_Integer] :=
   Transpose[ArrayReshape[b, {Nint, Nint, Nint, Nint}], {2, 1, 4, 3}];

corrFn[arr_, grid_] :=
   ListInterpolation[arr, {grid, grid, grid, grid}, InterpolationOrder -> 2];

covBase[meta_] :=
   <|"grid" -> N[meta["gridInt"]], "h" -> meta["h"],
     "fBox" -> meta["fbox"], "Nint" -> meta["Nint"]|>;

readCovRho[file_String] :=
Module[{stem = covStem[file], meta, Nint, grid, m},
  meta = readCovMeta[stem];
  Nint = meta["Nint"];
  grid = N[meta["gridInt"]];
  m = covMat[covDiag[stem, meta["n"], 1, Nint^2], Nint];
  Join[covBase[meta], <|"rhoValues" -> m, "rho" -> covFn[m, grid]|>]];

readCovQ[file_String] :=
Module[{stem = covStem[file], meta, Nint, ng, n, grid, v, c, m1, m2, m12},
  meta = readCovMeta[stem];
  Nint = meta["Nint"];
  ng = Nint^2;
  n = meta["n"];
  grid = N[meta["gridInt"]];
  v = covDiag[stem, n, ng + 1, 3 ng];
  c = covDiag[stem, n, ng + 1, 2 ng, ng];   (* Sigma[[ng+k, 2ng+k]] *)
  m1  = covMat[v[[1 ;; ng]], Nint];
  m2  = covMat[v[[ng + 1 ;; 2 ng]], Nint];
  m12 = covMat[c, Nint];
  Join[covBase[meta],
     <|"Q1Values" -> m1, "Q1" -> covFn[m1, grid],
       "Q2Values" -> m2, "Q2" -> covFn[m2, grid],
       "Q1Q2Values" -> m12, "Q1Q2" -> covFn[m12, grid]|>]];

readCovRhoQ[file_String] :=
Module[{stem = covStem[file], meta, Nint, ng, n, grid, v, c, m0, m1, m2, m12},
  meta = readCovMeta[stem];
  Nint = meta["Nint"];
  ng = Nint^2;
  n = meta["n"];
  grid = N[meta["gridInt"]];
  v = covDiag[stem, n, 1, 3 ng];
  c = covDiag[stem, n, ng + 1, 2 ng, ng];   (* Sigma[[ng+k, 2ng+k]] *)
  m0  = covMat[v[[1 ;; ng]], Nint];
  m1  = covMat[v[[ng + 1 ;; 2 ng]], Nint];
  m2  = covMat[v[[2 ng + 1 ;; 3 ng]], Nint];
  m12 = covMat[c, Nint];
  Join[covBase[meta],
     <|"rhoValues" -> m0, "rho" -> covFn[m0, grid],
       "Q1Values" -> m1, "Q1" -> covFn[m1, grid],
       "Q2Values" -> m2, "Q2" -> covFn[m2, grid],
       "Q1Q2Values" -> m12, "Q1Q2" -> covFn[m12, grid]|>]];

(* Two-point correlators Q1Q1(r,r') = <dQ1(r) dQ1(r')> etc., as 4-argument
   interpolants in (x, y, x', y').  Unlike the readers above these need FULL
   ng x ng blocks -- Nint^4 reals each, 2.9 GiB at Nint = 141 -- and an
   interpolant costs about 2x its array, so the returned association is roughly
   12 blocks' worth.  Q2Q1(r,r') = Q1Q2(r',r), so it is obtained by transposing
   Q1Q2 rather than re-reading the disk. *)
readCorrQ[file_String] :=
Module[{stem = covStem[file], meta, Nint, ng, n, grid, a11, a22, a12, a21},
  meta = readCovMeta[stem];
  Nint = meta["Nint"];
  ng = Nint^2;
  n = meta["n"];
  grid = N[meta["gridInt"]];
  a11 = corrArr[covBlock[stem, n, ng, 2, 2], Nint];
  a22 = corrArr[covBlock[stem, n, ng, 3, 3], Nint];
  a12 = corrArr[covBlock[stem, n, ng, 2, 3], Nint];
  a21 = Transpose[a12, {3, 4, 1, 2}];   (* Q2Q1(r,r') = Q1Q2(r',r) *)
  Join[covBase[meta],
     <|"Q1Q1Values" -> a11, "Q1Q1" -> corrFn[a11, grid],
       "Q2Q2Values" -> a22, "Q2Q2" -> corrFn[a22, grid],
       "Q1Q2Values" -> a12, "Q1Q2" -> corrFn[a12, grid],
       "Q2Q1Values" -> a21, "Q2Q1" -> corrFn[a21, grid]|>]];


(* ==================================================================== *)
(*  Tags, steady states, and whole runs                                 *)
(*                                                                      *)
(*  Everything above reads a covariance .bin given its path and is      *)
(*  UNCHANGED.  The set2_2_comove outputs written by                    *)
(*  lyapunov_solver/lyapunov_solver.m keep n, Nint, fbox, h and gridInt  *)
(*  under those names in covariance_*_meta.m, so readCovRho, readCovQ,  *)
(*  readCovRhoQ and readCorrQ read them with no change at all.  What    *)
(*  follows is additive: the run tag, the steady-state dump, and the    *)
(*  two together.                                                       *)
(* ==================================================================== *)

(* Two tag generations live side by side in fluctuations-set2_2_comove/:

     old   box15_N151_lN1.7300                 (zeta = 10, advection dropped)
     new   box15_N151_Dir_lN1.73_z8_ld0.8      (zeta = 8, advection present)

   parseTag reads both.  Fields the old form does not carry come back Missing,
   which is the point: a run whose "zeta" is Missing is one of the superseded
   ones, and that is worth being able to filter on.  Accepts a bare tag, a file
   name, or a full path with any of the four stems. *)
tagString[s_String] := Module[{t = FileNameTake[s]},
  (* FileNameTake + an EXPLICIT extension list, never FileBaseName: the tag ends
     in "_ld0.8", and FileBaseName reads that trailing ".8" as the extension and
     silently returns "..._ld0".  That turned elld into 0 and sent readSteady
     looking for a file that does not exist. *)
  t = StringReplace[t,
     ("." ~~ ("bin" | "m" | "png" | "pdf" | "csv" | "mx")) ~~ EndOfString -> ""];
  t = StringReplace[t,
     (* LONGEST first: StringReplace takes the first alternative that matches
        at a position, so "sigma_rho_" listed ahead of "sigma_rho_surface_"
        would leave a stray "surface_" behind. *)
     StartOfString ~~ ("sigma_rho_surface_" | "sigma_rho_dirichlet_" |
                       "sigma_rho_" | "covariance_" | "steady_") -> ""];
  StringReplace[t, "_meta" ~~ EndOfString -> ""]];

parseTag[s_String] := Module[{t = tagString[s], g},
  g[pat_] := Module[{m = StringCases[t, pat, 1]},
     If[m === {}, Missing["absent"], ToExpression @ First @ m]];
  <|"tag"    -> t,
    "box"    -> g["box" ~~ x : NumberString :> x],
    "Nint"   -> g["_N" ~~ x : DigitCharacter .. :> x],
    "bc"     -> Which[StringContainsQ[t, "_Neu"], "Neumann",
                      StringContainsQ[t, "_Dir"], "Dirichlet",
                      True, Missing["absent"]],
    "lNoise" -> g["_lN" ~~ x : NumberString :> x],
    "zeta"   -> g["_z" ~~ x : NumberString :> x],
    "elld"   -> g["_ld" ~~ x : NumberString :> x]|>];

(* The steady state the fluctuations were linearised about.  It is written by
   the DRIVER (jobscripts/set2_2_comove/fluctuations.wls), not by a module, and
   carries the fields twice over: on the Pade-graded tensor grid the comoving
   solver used, and resampled on the uniform Lyapunov grid.  Both are returned
   as interpolants; *Grid[[i,j]] = f(gx[[i]], gy[[j]]) in each case. *)
steadyPath[file_String] := Module[{dir = DirectoryName[ExpandFileName[file]]},
  FileNameJoin[{dir, "steady_" <> tagString[file] <> ".m"}]];

readSteady[file_String] := Module[{path, d, gx, gy, gi, mk},
  path = If[StringEndsQ[file, ".m"] && StringContainsQ[FileBaseName[file], "steady_"],
            ExpandFileName[file], steadyPath[file]];
  If[!FileExistsQ[path],
     Message[readSteady::nofile, path]; Return[$Failed]];
  d = Import[path];
  {gx, gy} = N /@ {d["gx"], d["gy"]};
  gi = N @ d["gridInt"];
  (* InterpolationOrder 3 so Derivative[2,0] survives, matching what the
     comoving solver itself hands out. *)
  mk[arr_, ax_, ay_] := ListInterpolation[arr, {ax, ay}, InterpolationOrder -> 3];
  Join[d,
     <|"rho" -> mk[d["rhoGrid"], gx, gy],
       "Q1"  -> mk[d["Q1Grid"],  gx, gy],
       "Q2"  -> mk[d["Q2Grid"],  gx, gy],
       "rhoOnInt" -> mk[d["rhoInt"], gi, gi],
       "Q1OnInt"  -> mk[d["Q1Int"],  gi, gi],
       "Q2OnInt"  -> mk[d["Q2Int"],  gi, gi]|>]];

readSteady::nofile = "No steady-state dump at `1`.";

(* Every run in a directory, newest tag scheme and old alike, as a list of
   associations ready for Dataset/Select.  Keyed off the covariance meta files,
   so a run whose .bin was never written does not appear. *)
listRuns[dir_String] := Module[{files},
  files = FileNames["covariance_*_meta.m", ExpandFileName[dir]];
  SortBy[
    Table[Join[parseTag[f],
       <|"meta" -> f,
         "bin" -> StringReplace[f, "_meta.m" ~~ EndOfString -> ".bin"],
         "sigmaRho" -> FileNameJoin[{DirectoryName[f],
            "sigma_rho_" <> tagString[f] <> ".m"}],
         "steady" -> steadyPath[f],
         "hasSteady" -> FileExistsQ[steadyPath[f]],
         "binMiB" -> If[FileExistsQ[StringReplace[f, "_meta.m" ~~ EndOfString -> ".bin"]],
            N[FileByteCount[StringReplace[f, "_meta.m" ~~ EndOfString -> ".bin"]]/2^20],
            Missing["absent"]]|>],
      {f, files}],
    {Lookup[#, "box", 0] &, Lookup[#, "Nint", 0] &}]];

(* Parameters + steady state + the covariance readers, in one object.  The
   covariance itself is NOT read here -- at Nint = 151 the .bin is 35 GiB.  The
   four covariance keys are stored as DELAYED rules, so the disk read happens on
   access and only for the ones you touch (and again on each access -- assign it
   if you need it twice).
     r = readRun["/data/.../covariance_box15_N50_Dir_lN1.73_z8_ld0.8.bin"];
     r["u"]      r["derived"]["alpha"]      r["steady"]["rho"][0.3, 0.1]
     sig = r["rhoQ"];    (* reads the .bin: the three variance fields *)     *)
readRun[file_String] := Module[{stem = covStem[file], meta, st},
  meta = readCovMeta[stem];
  st = Quiet @ readSteady[stem];
  Join[parseTag[stem], covBase[meta],
     <|"meta" -> meta,
       "n" -> meta["n"],
       "u" -> Lookup[meta, "u", Missing["absent"]],
       "derived" -> Lookup[meta, "derived", Missing["absent"]],
       "modelParams" -> Lookup[meta, "modelParams", Missing["absent"]],
       "solverParams" -> Lookup[meta, "solverParams", Missing["absent"]],
       "lyapunovParams" -> Lookup[meta, "lyapunovParams", Missing["absent"]],
       "otherParams" -> Lookup[meta, "otherParams", Missing["absent"]],
       "lNoiseEff" -> Lookup[meta, "lNoiseEff", Missing["absent"]],
       "boxFactor" -> Lookup[meta, "boxFactor", Missing["absent"]],
       "steady" -> st,
       "rho"   :> readCovRho[stem],
       "rhoQ"  :> readCovRhoQ[stem],
       "Q"     :> readCovQ[stem],
       "corrQ" :> readCorrQ[stem]|>]];
