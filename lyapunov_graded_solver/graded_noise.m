(* ::Package:: *)

(* ============================================================================
   graded_noise.m -- the noise matrix and the physical UV cutoff on a graded
   grid.

   ------------------------------------------------------------------ Q ------
   In the symmetrised variables of graded_ops.m the congruence is EXACTLY the
   uniform one with G -> Ghat, and the 1/h^2 prefactor disappears.  Derivation:
   the flux lives on faces with <F F^T> = diag(Mf) Sf^-1 (the Sf^-1 being the
   grid form of delta(r-r') at a face, dean.tex:1561), the force is its
   divergence f = -M^-1 G^T Sf F, so

       Q    = M^-1 G^T Sf diag(Mf) G M^-1
       Qhat = M^(1/2) Q M^(1/2) = Ghat^T diag(Mf) Ghat.

   Uniform check: M = h^2, Sf = h^2, Ghat = G, so Qhat = G^T diag G while
   Q = (1/h^2) G^T diag G -- and M^(1/2) = h on each side reproduces exactly the
   pf = kap/h^2 of lyapunov_solver.m:1036.  The NON-divergence (rotational
   Lambda) noise is even simpler: Qloc = M^-1 diag(lam) so Qhat_loc = diag(lam),
   metric-free.  Both 1/h^2 factors were always 1/(cell area).

   --------------------------------------------------------- the cutoff ------
   eq:lattice-bz (dean.tex:1201) is a BRILLOUIN-ZONE average.  It assumes
   translation invariance and has no meaning on a graded grid, so lEff cannot be
   calibrated the way lyapunov_solver.m:449 does it.  It does not need to be.

   The bare discrete pedestal is Sigma_ii = T/(k_i w_i) -- mesh-dependent by
   construction, because delta(r-r') -> delta_ij / w_i.  Build the smoother from
   the ANALYTIC Gaussian sampled on the grid and weighted by the quadrature,

       W = Ghauss . M,   Ghauss[i,j] = g(x_i - x_j),   g of width lNoise/Sqrt[2]

   (the /Sqrt[2] because W carries symbol Exp[-k^2 l^2/4] and is applied twice).
   Then

       Sigma_ii = T Sum_j g(x_i - x_j)^2 w_j / k_j  ->  (T/k_i) / (2 Pi l^2),

   i.e. eq:iso-gauss: THE 1/w_i OF THE BARE PEDESTAL IS CANCELLED EXACTLY BY THE
   w_j INSIDE W.  Mesh-independent, flat, on any grid, with no calibration step.
   Trapezoid quadrature of a Gaussian converges SPECTRALLY (error ~
   Exp[-2 Pi^2 sigma^2/h^2], 2.6e-9 at h = sigma) where the MatrixExp form of the
   discrete Laplacian carries the algebraic 1 + 1/(4 t^2) of eq:lattice-bz.  That
   is why this is not merely a convenience: it is what makes a varying h legal.

   In hatted variables What = M^(1/2) Ghauss M^(1/2) is symmetric, stays
   separable as Ey (x) Ex, and costs what Emat (x) Emat costs now.  It is still a
   congruence, so Qhat -> What Qhat What preserves PSD.

   Dirichlet is imposed by IMAGES -- the current sine-basis W is precisely the
   Dirichlet heat kernel, and images reproduce it while keeping the spectral
   accuracy that a discrete-Laplacian exponential would lose.

   ------------------------------------------- the per-node normalisation ----
   ONE CORRECTION IS STILL NEEDED, and it is easy to miss.  Trapezoid quadrature
   of a Gaussian is spectrally accurate only at EQUAL spacing (Poisson summation);
   with h varying it reverts to O(h^2).  Measured on a 6:1 graded grid the raw
   coincidence limit is off by 3.6% -- against the 0.2% flatness the uniform N131
   reference achieves, that is not good enough.

   So each 1-D factor is row-scaled by

       c_i = Sqrt[ Pcont_i / Pdisc_i ],
       Pdisc_i = Sum_j gD(x_i,x_j)^2 w_j,   Pcont_i = Int gD(x_i,y)^2 dy,

   the continuum integral being evaluated once on a fine UNIFORM subgrid where
   the trapezoid rule IS spectral.  This makes the coincidence limit exact at
   every node by construction, at no cost and without breaking separability or
   the congruence (W -> diag(c) W is still a congruence; W need not be symmetric,
   only Q -> W Q W^T).

   Normalising against the CONTINUUM rather than against a constant is what keeps
   the physics: near the wall the image cancellation genuinely suppresses gD, and
   that suppression appears in Pcont and Pdisc alike, so it survives the
   correction.  Only the quadrature error is removed.  This is the graded-grid
   replacement for the scalar lEff calibration of lyapunov_solver.m:449, and like
   it, it is a PEDESTAL-CHANNEL statement: it makes the rho-sector coincidence
   limit exact and does not claim to remove lattice corrections from the dQ
   sector.
   ============================================================================ *)

BeginPackage["GradedNoise`", {"GradedGrid`", "GradedOps`"}];

nuKappaOfS::usage = "nuKappaOfS[s] inverts I1(k)/I0(k) = s for the von Mises concentration.";
nuG2::usage = "nuG2[s] = I0 I2 / I1^2, the exact-closure U coefficient.";
nuG3::usage = "nuG3[s] = I0^2 I3 / I1^3, the exact-closure W coefficient.";
nuSmoother1D::usage =
  "nuSmoother1D[gridAssoc, lNoise] returns the symmetric hatted 1-D smoothing \
factor Ehat[i,j] = Sqrt[w_i] gD(x_i,x_j) Sqrt[w_j], gD the Dirichlet heat kernel \
of width lNoise/Sqrt[2] by images.  lNoise = 0 gives the identity.";
nuPedestalFactor::usage =
  "nuPedestalFactor[gridX, gridY, lNoise] returns the LOCAL lattice factor, the \
graded-grid generalisation of eq:lattice-bz: 2 Pi lNoise^2 times (W*W)(r,r).  \
With the per-node normalisation in force it is 1 in the interior by \
construction and droops at the wall where the Dirichlet images physically \
suppress it.  Option \"Raw\" -> True reports the uncorrected quadrature instead, \
which on a uniform grid is the analogue of 2 Pi t^2 F(t)^2.  Returns an nx x ny \
array.";
nuAssembleQ::usage =
  "nuAssembleQ[ops, flds, mp] builds the hatted 3n x 3n noise matrix from the \
node- and face-sampled steady state in flds and the model parameters mp.";
nuSmoothQ::usage =
  "nuSmoothQ[Q, Ex, Ey, nx, ny] applies the congruence Q -> What Q What blockwise \
by contracting Ex, Ey onto the four grid slots of each n x n block.";

Begin["`Private`"];

(* ---- exact von Mises closure: verbatim from lyapunov_solver.m:388-394 ---- *)
nuKappaOfS[s_?NumericQ] := If[s < 1.*^-6, 0.,
   kk /. Quiet@FindRoot[BesselI[1, kk]/BesselI[0, kk] == s,
                        {kk, 2 s/(1 - s^2)}, MaxIterations -> 200]];
nuG2[s_?NumericQ] := If[s < 0.02, 1/2 + s^2/6 + 11 s^4/180,
   With[{k = nuKappaOfS[s]}, BesselI[0, k] BesselI[2, k]/BesselI[1, k]^2]];
nuG3[s_?NumericQ] := If[s < 0.02, 1/6 + s^2/8 + 23 s^4/480,
   With[{k = nuKappaOfS[s]}, BesselI[0, k]^2 BesselI[3, k]/BesselI[1, k]^3]];

(* ---- the smoother ------------------------------------------------------- *)
(* Dirichlet heat kernel on [-L, L] by images.  n = 0 already suffices at
   sigma << L; +-1 is carried because it costs nothing and makes the wall exact. *)
(* VECTORISED over the second argument -- Clip keeps the far tail from
   underflowing (which otherwise floods the log with General::munfl) and
   UnitStep zeroes it exactly. *)
gKerV[v_, sig_] := Module[{z = v^2/(2. sig^2)},
   Exp[-Clip[z, {0., 700.}]] UnitStep[700. - z]/(Sqrt[2. Pi] sig)];
gDirV[x_?NumericQ, ys_List, sig_, L_] :=
   Sum[gKerV[x - ys + 4. L nn, sig] - gKerV[x + ys + 2. L + 4. L nn, sig],
       {nn, -1, 1}];

(* Sum_j gD(x_i,x_j)^2 w_j -- the discrete coincidence limit, per axis *)
perAxisDisc[gr_Association, lNoise_?NumericQ] :=
Module[{g = gr["gridInt"], w = gr["w"], L = gr["fbox"], sig},
  sig = N[lNoise]/Sqrt[2.];
  Table[gDirV[g[[i]], g, sig, L]^2 . w, {i, Length[g]}]];

(* Int_-L^L gD(x_i,y)^2 dy on a fine UNIFORM subgrid, where trapezoid IS
   spectral.  nsub = 4000 gives dy/sigma ~ 0.006 at the parameters in use, so
   the quadrature error here is far below round-off. *)
perAxisCont[gr_Association, lNoise_?NumericQ, nsub_: 4000] :=
Module[{g = gr["gridInt"], L = gr["fbox"], sig, ys, dy, wq},
  sig = N[lNoise]/Sqrt[2.];
  ys = N@Subdivide[-L, L, nsub];
  dy = 2. L/nsub;
  wq = dy ConstantArray[1., nsub + 1];
  wq[[1]] /= 2.; wq[[-1]] /= 2.;
  Table[gDirV[g[[i]], ys, sig, L]^2 . wq, {i, Length[g]}]];

(* the row scaling that makes the coincidence limit exact at every node *)
nuSmootherScale[gr_Association, lNoise_?NumericQ] :=
   Sqrt[perAxisCont[gr, lNoise]/perAxisDisc[gr, lNoise]];

nuSmoother1D[gr_Association, lNoise_?NumericQ] :=
Module[{g = gr["gridInt"], w = gr["w"], L = gr["fbox"], sig, ms, c},
  If[lNoise == 0., Return[IdentityMatrix[gr["Nint"]]]];
  sig = N[lNoise]/Sqrt[2.];
  ms = Sqrt[w];
  c = nuSmootherScale[gr, lNoise];
  (* Ehat[i,j] = c_i Sqrt[w_i] gD(x_i,x_j) Sqrt[w_j] *)
  c * Transpose[ms * Transpose[ms * Table[gDirV[g[[i]], g, sig, L], {i, Length[g]}]]]];

(* The local lattice factor: the graded generalisation of eq:lattice-bz.  With
   "Raw" -> True it reports the UNCORRECTED quadrature, so the size of the
   normalisation is visible; by default it reports what the solver actually uses,
   which is 1 in the interior by construction and retains the physical Dirichlet
   droop at the wall. *)
Options[nuPedestalFactor] = {"Raw" -> False};
nuPedestalFactor[grX_Association, grY_Association, lNoise_?NumericQ,
                 OptionsPattern[]] :=
Module[{px, py},
  {px, py} = If[TrueQ@OptionValue["Raw"],
                {perAxisDisc[grX, lNoise], perAxisDisc[grY, lNoise]},
                {perAxisCont[grX, lNoise], perAxisCont[grY, lNoise]}];
  2. Pi N[lNoise]^2 Outer[Times, px, py]];

(* ---- the noise matrix --------------------------------------------------- *)
(* Hatted divergence block.  Identical in form to lyapunov_solver.m:1040-1042;
   the metric lives inside Ghat/Dhat and there is no 1/h^2 prefactor. *)
divBlk[op_, cxxF_, cyyF_, cxy_, cyx_] :=
   Transpose[op["Gx"]] . (cxxF * op["Gx"]) + Transpose[op["Gy"]] . (cyyF * op["Gy"]) +
   Transpose[op["Dx"]] . (cxy * op["Dy"]) + Transpose[op["Dy"]] . (cyx * op["Dx"]);

nuAssembleQ[op_Association, f_Association, mp_Association] :=
Module[{zet, xi0, rho0v, Lam, kap, locPf, n, dg,
        rhoV, Q1V, Q2V, rhoXE, Q1XE, Q2XE, rhoYE, Q1YE, Q2YE,
        g2V, g3V, g2XE, g3XE, g2YE, g3YE,
        h2c, h2s, h3c, h3s, e2, e3, h2cXE, h2sXE, h3cXE, h3sXE, e2XE, e3XE,
        h2cYE, h2sYE, h3cYE, h3sYE, e2YE, e3YE,
        Qrr, QQ1Q1, QQ2Q2, QQ1Q2, QQ1r, QQ2r, Qt},
  zet = mp["zeta"]; xi0 = mp["xi0"]; rho0v = mp["rho0"]; Lam = mp["Lambda"];
  kap = 2 zet/(xi0 rho0v);
  locPf = 2 Lam/rho0v;                  (* metric-free: see the header *)
  n = op["n"];
  dg[c_] := DiagonalMatrix[SparseArray@c];

  {rhoV, Q1V, Q2V} = f /@ {"rhoV", "Q1V", "Q2V"};
  {rhoXE, Q1XE, Q2XE} = f /@ {"rhoXE", "Q1XE", "Q2XE"};
  {rhoYE, Q1YE, Q2YE} = f /@ {"rhoYE", "Q1YE", "Q2YE"};
  {g2V, g3V, g2XE, g3XE, g2YE, g3YE} =
     f /@ {"g2V", "g3V", "g2XE", "g3XE", "g2YE", "g3YE"};

  h2c = Q1V^2 - Q2V^2;        h2s = 2 Q1V Q2V;
  h3c = Q1V^3 - 3 Q1V Q2V^2;  h3s = 3 Q1V^2 Q2V - Q2V^3;
  e2 = g2V/(2 rhoV);          e3 = g3V/(4 rhoV^2);
  h2cXE = Q1XE^2 - Q2XE^2;        h2sXE = 2 Q1XE Q2XE;
  h3cXE = Q1XE^3 - 3 Q1XE Q2XE^2; h3sXE = 3 Q1XE^2 Q2XE - Q2XE^3;
  e2XE = g2XE/(2 rhoXE);          e3XE = g3XE/(4 rhoXE^2);
  h2cYE = Q1YE^2 - Q2YE^2;        h2sYE = 2 Q1YE Q2YE;
  h3cYE = Q1YE^3 - 3 Q1YE Q2YE^2; h3sYE = 3 Q1YE^2 Q2YE - Q2YE^3;
  e2YE = g2YE/(2 rhoYE);          e3YE = g3YE/(4 rhoYE^2);

  Qrr = kap divBlk[op, Q1XE + rhoXE, -Q1YE + rhoYE, Q2V, Q2V];

  QQ1Q1 = kap divBlk[op, 3 Q1XE/4 + rhoXE/2 + e2XE h2cXE + e3XE h3cXE,
                         -3 Q1YE/4 + rhoYE/2 + e2YE h2cYE - e3YE h3cYE,
                         Q2V/4 + e3 h3s, Q2V/4 + e3 h3s] +
          locPf dg[2 rhoV - 4 e2 h2c];

  QQ2Q2 = kap divBlk[op, Q1XE/4 + rhoXE/2 - e2XE h2cXE - e3XE h3cXE,
                         -Q1YE/4 + rhoYE/2 - e2YE h2cYE + e3YE h3cYE,
                         3 Q2V/4 - e3 h3s, 3 Q2V/4 - e3 h3s] +
          locPf dg[2 rhoV + 4 e2 h2c];

  QQ1Q2 = kap divBlk[op, Q2XE/4 + e2XE h2sXE + e3XE h3sXE,
                         -Q2YE/4 + e2YE h2sYE - e3YE h3sYE,
                         Q1V/4 - e3 h3c, Q1V/4 - e3 h3c] +
          locPf dg[-4 e2 h2s];

  QQ1r = kap divBlk[op, Q1XE + rhoXE/2 + e2XE h2cXE,
                        Q1YE - rhoYE/2 - e2YE h2cYE, e2 h2s, e2 h2s];
  QQ2r = kap divBlk[op, Q2XE + e2XE h2sXE, Q2YE - e2YE h2sYE,
                        rhoV/2 - e2 h2c, rhoV/2 - e2 h2c];

  Qt = ArrayFlatten[{{Qrr, Transpose[QQ1r], Transpose[QQ2r]},
                     {QQ1r, QQ1Q1, QQ1Q2},
                     {QQ2r, Transpose[QQ1Q2], QQ2Q2}}];
  (Qt + Transpose[Qt])/2];

(* ---- the congruence Q -> What Q What, blockwise ------------------------- *)
(* k = i + (j-1) nx, so an n x n block reshapes to X[[j,i,j',i']] and the
   congruence is Ex, Ey contracted onto the four slots -- lyapunov_solver.m:1098. *)
smoothBlock[X_, Ex_, Ey_, nx_, ny_] :=
Module[{T = ArrayReshape[X, {ny, nx, ny, nx}]},
  T = Ey . T;
  T = Transpose[Ex . Transpose[T, {2, 1, 3, 4}], {2, 1, 3, 4}];
  T = Transpose[Ey . Transpose[T, {3, 2, 1, 4}], {3, 2, 1, 4}];
  T = Transpose[Ex . Transpose[T, {4, 2, 3, 1}], {4, 2, 3, 1}];
  ArrayReshape[T, {nx ny, nx ny}]];

nuSmoothQ[Q_, Ex_, Ey_, nx_, ny_] :=
Module[{n = nx ny, R, rng, S},
  R = Normal[Q];
  rng[a_] := (a - 1) n + 1 ;; a n;
  Do[Do[S = smoothBlock[R[[rng[a], rng[b]]], Ex, Ey, nx, ny];
        R[[rng[a], rng[b]]] = S, {b, 3}], {a, 3}];
  (R + Transpose[R])/2];

End[];
EndPackage[];
