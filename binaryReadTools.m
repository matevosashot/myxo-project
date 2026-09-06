(* ::Package:: *)

covStem[file_String] := If[StringEndsQ[file, ".bin"], StringDrop[file, -4], file];

readCovMeta[file_String] := Import[covStem[file] <> "_meta.m"];

covDiag[stem_String, n_Integer, k1_Integer, k2_Integer] :=
Module[{str = OpenRead[stem <> ".bin", BinaryFormat -> True]},
  WithCleanup[
     Table[SetStreamPosition[str, 8 (k - 1) (n + 1)];
           BinaryRead[str, "Real64"], {k, k1, k2}],
     Close[str]]];

covMat[vec_, Nint_] := Transpose[Partition[vec, Nint]];

covFn[mat_, grid_] :=
   ListInterpolation[mat, {grid, grid}, InterpolationOrder -> 2];

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
Module[{stem = covStem[file], meta, Nint, ng, grid, v, m1, m2},
  meta = readCovMeta[stem];
  Nint = meta["Nint"];
  ng = Nint^2;
  grid = N[meta["gridInt"]];
  v = covDiag[stem, meta["n"], ng + 1, 3 ng];
  m1 = covMat[v[[1 ;; ng]], Nint];
  m2 = covMat[v[[ng + 1 ;; 2 ng]], Nint];
  Join[covBase[meta],
     <|"Q1Values" -> m1, "Q1" -> covFn[m1, grid],
       "Q2Values" -> m2, "Q2" -> covFn[m2, grid]|>]];

readCovRhoQ[file_String] :=
Module[{stem = covStem[file], meta, Nint, ng, grid, v, m0, m1, m2},
  meta = readCovMeta[stem];
  Nint = meta["Nint"];
  ng = Nint^2;
  grid = N[meta["gridInt"]];
  v = covDiag[stem, meta["n"], 1, 3 ng];
  m0 = covMat[v[[1 ;; ng]], Nint];
  m1 = covMat[v[[ng + 1 ;; 2 ng]], Nint];
  m2 = covMat[v[[2 ng + 1 ;; 3 ng]], Nint];
  Join[covBase[meta],
     <|"rhoValues" -> m0, "rho" -> covFn[m0, grid],
       "Q1Values" -> m1, "Q1" -> covFn[m1, grid],
       "Q2Values" -> m2, "Q2" -> covFn[m2, grid]|>]];
