(* ::Package:: *)

(* ============================================================================
   graded_lyapunov.m -- the top-level driver of the fluctuation solver.

   Same physics and the same scheme as lyapunov_solver.m; the only difference is
   that every operator carries the metric of a graded grid.  Because the
   symmetrising transform of graded_ops.m puts the problem back into the
   plain-transpose form, EVERY assembly formula below is character-for-character
   the uniform one with G -> Ghat, D -> Dhat, and the 1/h^2 prefactors gone:

       A11 = (B/xi0) Lap                      lyapunov_solver.m:966
       A21 = -(B/xi0) Sum_i Gi^T diag(Q|face) Gi          :969
       Adv = -(Gx^T diag(w|face) Px + ...)               :961
       Q   = kap Sum_i Gi^T diag(M|face) Gi + ...        :1040

   Solve       Ahat Chat + Chat Ahat^T = -Qhat
   then        C = M^(-1/2) Chat M^(-1/2),  Sigma(r_k) = Chat_kk / w_k.

   The 1/w_k there is the whole story of the pedestal: delta(r-r') -> delta/w_k,
   so the BARE variance is mesh-dependent by construction and it is the smoother
   of graded_noise.m -- whose quadrature weight cancels it -- that makes the
   answer physical.

   SCOPE.  Dirichlet (reservoir wall) only; the steady state comes from
   comoving_steady_solver.m unchanged.  Neumann needs the sealed finite-volume
   gradient and the anchored L1q corner re-derived with local cell sizes, and is
   deliberately not attempted here.
   ============================================================================ *)

BeginPackage["GradedLyapunov`",
  {"GradedGrid`", "GradedOps`", "GradedNoise`"}];

nuSolve::usage =
  "nuSolve[gridSpec, modelParams, steadyState] solves the graded-grid Lyapunov \
problem.  gridSpec is either a grid association from nuGrid/nuGridUniform (used \
for both axes) or {gridX, gridY}.  modelParams is the usual rule list.  \
steadyState must supply callable \"rho\", \"Q1\", \"Q2\" and a comoving \"u\".";

nuSteadyFromFile::usage =
  "nuSteadyFromFile[path] rebuilds a callable steady state from a steady_*.m \
dump written by jobscripts/set2_2_comove/fluctuations.wls (keys gx, gy, \
rhoGrid, Q1Grid, Q2Grid, u).";

nuExport::usage =
  "nuExport[res, outputDir, tag] writes covariance_<tag>.bin and \
covariance_<tag>_meta.m in the layout binaryReadTools.m already reads, with the \
graded axis grids and the cell weights added.";

Begin["`Private`"];

nuSteadyFromFile[path_String] :=
Module[{s = Import[path], gx, gy, ip},
  gx = N@s["gx"]; gy = N@s["gy"];
  ip[gr_] := ListInterpolation[gr, {gx, gy}, InterpolationOrder -> 3];
  <|"rho" -> ip[s["rhoGrid"]], "Q1" -> ip[s["Q1Grid"]], "Q2" -> ip[s["Q2Grid"]],
    "u" -> N@s["u"], "modelParams" -> s["modelParams"], "source" -> path|>];

(* ---- reusable pieces, so verify.wls exercises the SAME code path ---- *)

nuPars[modelParams_] :=
Module[{mp = Association[modelParams], cL, aPrime},
  cL = mp[Global`L] + mp[Global`\[Zeta]] mp[Global`\[Xi]r]/(4 mp[Global`\[Xi]0]);
  aPrime = mp[Global`a] + mp[Global`\[CapitalLambda]] mp[Global`\[Xi]r];
  <|"B" -> mp[Global`B], "a" -> mp[Global`a], "b" -> mp[Global`b],
    "zeta" -> mp[Global`\[Zeta]], "xi0" -> mp[Global`\[Xi]0],
    "xir" -> mp[Global`\[Xi]r], "rho0" -> mp[Global`\[Rho]0],
    "Lambda" -> mp[Global`\[CapitalLambda]], "lNoise" -> N@mp[Global`lNoise],
    "KPrime" -> cL, "aPrime" -> aPrime,
    "elld" -> Sqrt[cL/(-aPrime)], "S0" -> Sqrt[-aPrime/(2 mp[Global`b])],
    "alpha" -> mp[Global`\[Zeta]] mp[Global`\[Xi]r] Sqrt[-aPrime/(2 mp[Global`b])]/
               (4 mp[Global`\[Xi]0] cL),
    "gamma0" -> mp[Global`B]/mp[Global`\[Xi]0],
    "kappaQ" -> 4 cL/mp[Global`\[Xi]r],
    "kappaNoise" -> 2 mp[Global`\[Zeta]]/(mp[Global`\[Xi]0] mp[Global`\[Rho]0])|>];

(* node and face samples plus the exact von Mises closure *)
nuSampleFields[op_Association, ss_Association] :=
Module[{a = <||>, rho = ss["rho"], q1 = ss["Q1"], q2 = ss["Q2"], Sv, Sx, Sy},
  a["rhoV"] = rho[op["xflat"], op["yflat"]];
  a["Q1V"] = q1[op["xflat"], op["yflat"]];
  a["Q2V"] = q2[op["xflat"], op["yflat"]];
  a["rhoXE"] = rho[op["xeX"], op["xeY"]];
  a["Q1XE"] = q1[op["xeX"], op["xeY"]];
  a["Q2XE"] = q2[op["xeX"], op["xeY"]];
  a["rhoYE"] = rho[op["yeX"], op["yeY"]];
  a["Q1YE"] = q1[op["yeX"], op["yeY"]];
  a["Q2YE"] = q2[op["yeX"], op["yeY"]];
  a["rhoXd"] = Derivative[1, 0][rho][op["xeX"], op["xeY"]];
  a["rhoYd"] = Derivative[0, 1][rho][op["yeX"], op["yeY"]];
  Sv = Sqrt[a["Q1V"]^2 + a["Q2V"]^2]/a["rhoV"];
  Sx = Sqrt[a["Q1XE"]^2 + a["Q2XE"]^2]/a["rhoXE"];
  Sy = Sqrt[a["Q1YE"]^2 + a["Q2YE"]^2]/a["rhoYE"];
  a["S"] = Sv;
  a["g2V"] = nuG2 /@ Sv;   a["g3V"] = nuG3 /@ Sv;
  a["g2XE"] = nuG2 /@ Sx;  a["g3XE"] = nuG3 /@ Sx;
  a["g2YE"] = nuG2 /@ Sy;  a["g3YE"] = nuG3 /@ Sy;
  a];

(* the drift, in hatted operators.  Formula-for-formula lyapunov_solver.m:960-982. *)
nuBuildA[op_Association, f_Association, p_Association, uVec_List] :=
Module[{Bv = p["B"], av = p["a"], bv = p["b"], zet = p["zeta"], xi0 = p["xi0"],
        xir = p["xir"], Lam = p["Lambda"], cL = p["KPrime"], dg, Adv,
        A11, A12, A13, A21, A22, A23, A31, A33, wxE, wyE},
  dg[c_] := DiagonalMatrix[SparseArray@c];
  wxE = uVec[[1]] + (Bv/xi0) f["rhoXd"];
  wyE = uVec[[2]] + (Bv/xi0) f["rhoYd"];
  Adv = -(Transpose[op["Gx"]] . (wxE * op["Px"]) +
          Transpose[op["Gy"]] . (wyE * op["Py"]));
  A11 = (Bv/xi0) op["Lap"];
  A12 = (zet/xi0) (op["Lx"] - op["Ly"]);
  A13 = (2 zet/xi0) op["Dxy"];
  A21 = -(Bv/xi0) (Transpose[op["Gx"]] . (f["Q1XE"] * op["Gx"]) +
                   Transpose[op["Gy"]] . (f["Q1YE"] * op["Gy"]));
  A22 = Adv + (4/xir) (dg[-(av + Lam xir) - bv (6 f["Q1V"]^2 + 2 f["Q2V"]^2)] +
                       cL op["Lap"]);
  A23 = -(16 bv/xir) dg[f["Q1V"] f["Q2V"]];
  A31 = -(Bv/xi0) (Transpose[op["Gx"]] . (f["Q2XE"] * op["Gx"]) +
                   Transpose[op["Gy"]] . (f["Q2YE"] * op["Gy"]));
  A33 = Adv + (4/xir) (dg[-(av + Lam xir) - bv (2 f["Q1V"]^2 + 6 f["Q2V"]^2)] +
                       cL op["Lap"]);
  SparseArray@ArrayFlatten[{{A11, A12, A13}, {A21, A22, A23}, {A31, A23, A33}}]];

(* convenience: sample and assemble in one call, for the tests *)
nuBuildAFrom[op_Association, ss_Association, modelParams_] :=
   nuBuildA[op, nuSampleFields[op, ss], nuPars[modelParams],
            N@ss["u"][[1 ;; 2]]];

nuSolve[gridSpec_, modelParams_, steadyState_Association] :=
Module[{grX, grY, op, mp, Bv, av, bv, Lv, zet, xi0, xir, rho0v, Lam, lNoiseV,
        cL, aPrime, elldV, S0V, derived, rhoSs, Q1Ss, Q2Ss, ux, uy,
        n, nx, ny, mass, misqrt, f, Svals, SXE, SYE, dg,
        Adv, A11, A12, A13, A21, A22, A23, A31, A32, A33, Ahat,
        Qhat, Ex, Ey, Adense, Qdense, eigvals, eigvecs, Pev, lamSum,
        maxReLam, Chat, Cmat, resid, minEigQ, lapId, antisym, pedFac,
        sigRho, sigQ1, sigQ2, kmat, t, dummy, timings = {}, rec},
  rec[lbl_, v_] := (AppendTo[timings, {lbl, v}]; v);

  (* ---------- [1] parameters ---------- *)
  {grX, grY} = If[ListQ[gridSpec], gridSpec, {gridSpec, gridSpec}];
  mp = Association[modelParams];
  Bv = mp[Global`B]; av = mp[Global`a]; bv = mp[Global`b]; Lv = mp[Global`L];
  zet = mp[Global`\[Zeta]]; xi0 = mp[Global`\[Xi]0]; xir = mp[Global`\[Xi]r];
  rho0v = mp[Global`\[Rho]0]; Lam = mp[Global`\[CapitalLambda]];
  lNoiseV = N@mp[Global`lNoise];
  cL = Lv + zet xir/(4 xi0);
  aPrime = av + Lam xir;
  elldV = Sqrt[cL/(-aPrime)];  S0V = Sqrt[-aPrime/(2 bv)];
  derived = <|"aPrime" -> aPrime, "KPrime" -> cL, "elld" -> elldV, "S0" -> S0V,
              "alpha" -> zet xir S0V/(4 xi0 cL), "gamma0" -> Bv/xi0,
              "kappaQ" -> 4 cL/xir, "kappaNoise" -> 2 zet/(xi0 rho0v)|>;
  {rhoSs, Q1Ss, Q2Ss} = steadyState /@ {"rho", "Q1", "Q2"};
  {ux, uy} = N@steadyState["u"][[1 ;; 2]];

  (* ---------- [2] grid, operators, the adjoint-pair gate ---------- *)
  Print["[2] Graded operators ..."];
  {t, op} = AbsoluteTiming@nuOps2D[grX, grY];
  rec["[2] operators", t];
  nx = op["nx"]; ny = op["ny"]; n = op["n"];
  mass = op["mass"]; misqrt = op["misqrt"];
  lapId = nuLapIdentity[op];
  antisym = nuAntisymmetry[op];
  Print["    Nint = ", nx, " x ", ny, "   n = ", 3 n,
        "   h in [", grX["hMin"], ", ", grX["hMax"], "]"];
  Print["    lapIdentity = ", ScientificForm[lapId, 3],
        If[lapId <= 1.*^-12, "   (adjoint pair)", "   *** NOT AN ADJOINT PAIR ***"]];
  Print["    |Dhat^T + Dhat| = ", ScientificForm[antisym, 3]];
  If[lapId > 1.*^-12,
     Print["    ABORT: the scheme is void without the adjoint pair."]; Abort[]];

  (* ---------- [3] sample the steady state, node and face ---------- *)
  Print["[3] Sampling steady state and closure ..."];
  {t, f} = AbsoluteTiming@nuSampleFields[op, steadyState];
  rec["[3] sample + exact closure", t];
  Print["    S in ", MinMax@f["S"], "   (exact von Mises closure)"];

  (* ---------- [4] drift ---------- *)
  Print["[4] Assembling A ..."];
  {t, Ahat} = AbsoluteTiming@nuBuildA[op, f,
     nuPars[modelParams], {ux, uy}];
  rec["[4] assemble A", t];

  (* ---------- [5] noise, then the physical cutoff ---------- *)
  Print["[5] Assembling Q ..."];
  {t, Qhat} = AbsoluteTiming@nuAssembleQ[op, f,
     <|"zeta" -> zet, "xi0" -> xi0, "rho0" -> rho0v, "Lambda" -> Lam|>];
  rec["[5] assemble Q", t];
  Print["[5b] Cutoff: Gaussian quadrature smoother, lNoise = ", lNoiseV, " ..."];
  {t, dummy} = AbsoluteTiming[
     Ex = nuSmoother1D[grX, lNoiseV];
     Ey = nuSmoother1D[grY, lNoiseV]];
  rec["[5b] build W", t];
  pedFac = nuPedestalFactor[grX, grY, lNoiseV];
  Print["    local pedestal factor over the interior: ",
        MinMax@Flatten@pedFac[[3 ;; -3, 3 ;; -3]], "   (1 = eq:iso-gauss)"];
  If[lNoiseV > 0.,
     {t, Qhat} = AbsoluteTiming@nuSmoothQ[Qhat, Ex, Ey, nx, ny];
     rec["[5b] apply W.Q.W", t]];

  Qdense = Developer`ToPackedArray@N@Normal@Qhat;
  minEigQ = If[3 n <= 20000,
     Quiet@Min@Re@Eigenvalues[Qdense, -4, Method -> "Arnoldi"],
     Missing["skipped, n > 20000"]];
  Print["    assembled Q min eig = ", minEigQ,
        If[NumericQ[minEigQ] && minEigQ < -1.*^-8 Max@Abs@Diagonal@Qdense,
           "   *** Q IS NOT PSD ***", "   (PSD)"]];

  (* ---------- [6] solve ---------- *)
  Print["[6] Lyapunov solve, ", 3 n, " x ", 3 n, " ..."];
  Adense = Developer`ToPackedArray@N@Normal@Ahat;
  {t, {eigvals, eigvecs}} = AbsoluteTiming@Eigensystem[Adense];
  rec["[6] Eigensystem", t];
  maxReLam = Max@Re@eigvals;
  Print["    max Re(lambda) = ", maxReLam,
        If[maxReLam < 0, "   (stable)", "   *** UNSTABLE: not a covariance ***"]];
  Pev = Developer`ToPackedArray@Transpose@eigvecs;
  lamSum = Developer`ToPackedArray@Outer[Plus, eigvals, eigvals];
  {t, Chat} = AbsoluteTiming@Re[
     Pev . ((-LinearSolve[Pev, Transpose@LinearSolve[Pev, Transpose[Qdense]]])/
            lamSum) . Transpose[Pev]];
  rec["[6] back-substitute", t];
  resid = Max@Abs[Adense . Chat + Chat . Transpose[Adense] + Qdense]/
          Max@Abs@Qdense;
  Print["    relative residual = ", ScientificForm[resid, 3]];
  Clear[Adense, Qdense, Pev, lamSum];

  (* ---------- [7] undo the metric ---------- *)
  (* C = M^(-1/2) Chat M^(-1/2); the variance at node k is Chat_kk / w_k. *)
  kmat[b_] := (b - 1) n + 1 ;; b n;
  sigRho = Transpose@Partition[Diagonal[Chat[[kmat[1], kmat[1]]]]/mass, nx];
  sigQ1 = Transpose@Partition[Diagonal[Chat[[kmat[2], kmat[2]]]]/mass, nx];
  sigQ2 = Transpose@Partition[Diagonal[Chat[[kmat[3], kmat[3]]]]/mass, nx];
  Cmat = Developer`ToPackedArray[
     Chat * Outer[Times, Join[misqrt, misqrt, misqrt], Join[misqrt, misqrt, misqrt]]];

  Print["    Sigma_rho in ", MinMax@sigRho,
        "   continuum pedestal = ", zet/(2 Pi Bv rho0v lNoiseV^2)];

  <|"C" -> Cmat, "Chat" -> Chat, "sigmaRho" -> sigRho,
    "sigmaQ1" -> sigQ1, "sigmaQ2" -> sigQ2,
    "gridX" -> grX["gridInt"], "gridY" -> grY["gridInt"],
    "gridInt" -> grX["gridInt"], "weights" -> mass,
    "hMin" -> grX["hMin"], "hMax" -> grX["hMax"], "uniform" -> grX["uniform"],
    "Nint" -> nx, "n" -> 3 n, "fbox" -> grX["fbox"],
    "bcType" -> "Dirichlet", "u" -> {ux, uy},
    "lNoise" -> lNoiseV, "pedestalFactor" -> pedFac,
    "pedestalContinuum" -> zet/(2 Pi Bv rho0v lNoiseV^2),
    "eigenvalues" -> eigvals, "maxReLambda" -> maxReLam,
    "minEigQ" -> minEigQ, "residual" -> resid,
    "lapIdentity" -> lapId, "antisymmetry" -> antisym,
    "derived" -> derived, "modelParams" -> modelParams,
    "scheme" -> "graded tensor grid; metric-aware adjoint pair Lhat = -Ghat^T Ghat \
with Ghat = Sf^(1/2) G M^(-1/2); divergence-form A21/A31; collocated cross terms; \
exact von Mises closure; Dirichlet reservoir wall; cutoff by GAUSSIAN QUADRATURE \
smoother normalised per node against the continuum coincidence limit (eq:lattice-bz \
does not apply on a graded grid)",
    "timings" -> timings|>];

nuExport[res_Association, outputDir_String, tag_String] :=
Module[{dir, binPath, metaPath},
  dir = If[DirectoryQ[outputDir], outputDir,
           Quiet@CreateDirectory[outputDir, CreateIntermediateDirectories -> True];
           outputDir];
  binPath = FileNameJoin[{dir, "covariance_" <> tag <> ".bin"}];
  metaPath = FileNameJoin[{dir, "covariance_" <> tag <> "_meta.m"}];
  Export[binPath, res["C"], {"Binary", "Real64"}];
  Export[metaPath,
     <|"n" -> res["n"], "Nint" -> res["Nint"], "fbox" -> res["fbox"],
       (* h is NOT a scalar here; gridInt/gridX/gridY carry the real geometry and
          binaryReadTools.m's covFn handles a non-uniform ListInterpolation grid
          unchanged.  hMin/hMax are reported so a reader can see the grading. *)
       "h" -> Missing["graded grid: see hMin, hMax, gridX, gridY"],
       "hMin" -> res["hMin"], "hMax" -> res["hMax"], "uniform" -> res["uniform"],
       "gridInt" -> res["gridInt"], "gridX" -> res["gridX"],
       "gridY" -> res["gridY"], "weights" -> res["weights"],
       "blocks" -> {"rho", "Q1", "Q2"},
       "ordering" -> "row-major Real64; grid index k = i + (j-1) Nint, x fastest",
       (* the health record travels WITH the data, not in a side file *)
       "maxReLambda" -> res["maxReLambda"], "minEigQ" -> res["minEigQ"],
       "residual" -> res["residual"], "lapIdentity" -> res["lapIdentity"],
       "antisymmetry" -> res["antisymmetry"],
       "lNoise" -> res["lNoise"], "bcType" -> res["bcType"], "u" -> res["u"],
       "pedestalContinuum" -> res["pedestalContinuum"],
       "derived" -> res["derived"], "modelParams" -> res["modelParams"],
       "scheme" -> res["scheme"], "timings" -> res["timings"]|>];
  {binPath, metaPath}];

End[];
EndPackage[];
