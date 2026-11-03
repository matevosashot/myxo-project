#!/usr/bin/env wolframscript
(* export_fluct_ref.wls -- the paper's fluctuation runs (ssDir, B = 9e4) as JSON for
   the comparison in Sec. "Density fluctuations": data/fluct_ref_N<Nint>.json with
   gridInt, sigmaRho, sigmaQ1, sigmaQ2 ([j][i] = (x_i, y_j)), maxReLambda, lNoiseEff, u. *)
dir = FileNameJoin[{DirectoryName[$InputFileName], "..", "data",
                    "fluctuations-set2_2_comove_optimized_2_ssDir"}];
Do[
  Module[{f = FileNameJoin[{dir, "sigma_rho_box20_N" <> ToString[nI] <>
                                 "_Dir_lN0.5_z12.5_ld1_B90000_ss20.m"}], a},
    a = Get[f];
    Export[FileNameJoin[{DirectoryName[$InputFileName], "data",
                         "fluct_ref_N" <> ToString[nI] <> ".json"}],
      <|"file" -> f, "gridInt" -> N@a["gridInt"], "sigmaRho" -> a["sigmaRho"],
        "sigmaQ1" -> a["sigmaQ1"], "sigmaQ2" -> a["sigmaQ2"],
        "maxReLambda" -> a["maxReLambda"], "lNoiseEff" -> a["lNoiseEff"],
        "residualRefined" -> a["residualRefined"], "u" -> N@a["u"]|>, "RawJSON"];
    Print["exported N", nI]],
  {nI, {150, 180}}];
