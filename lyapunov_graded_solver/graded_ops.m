(* ::Package:: *)

(* ============================================================================
   graded_ops.m -- the adjoint pair on a graded tensor grid.

   THE PROBLEM.  lyapunov_solver.m rests on ONE identity (its :745 gate,
   dean.tex eq:G-lap-identity):

       Lap + (Gx^T Gx + Gy^T Gy) == 0

   with a PLAIN transpose.  That is the adjoint in the unweighted l2 inner
   product, and it is the right adjoint only because a uniform mesh has mass
   matrix h^2 I -- a scalar multiple of the identity, which cancels from both
   sides.  On a graded grid the natural inner product is <u,v> = Sum_k w_k u_k v_k
   with w varying, the adjoint of G becomes M^-1 G^T S_f, and a plain transpose
   silently breaks the discrete fluctuation-dissipation balance.  dean.tex:1526
   prices that failure at "exactly one half, at every h, for a scheme that is
   second-order accurate and passes every classical consistency check", and
   :1530 notes refinement cannot repair it.  So it must not be got wrong.

   THE FIX -- one similarity transform, not a rewrite.  With

       M   = diag(w)           node dual cells   (= the trapezoid weights)
       S_f = diag(d)           face dual cells   (= the node separations)
       G   : (G u)_m = (u_m - u_{m-1}) / d_m

   define the SYMMETRISED gradient and the scaled state

       Ghat = S_f^(1/2) G M^(-1/2),        psiHat = M^(1/2) psi.

   Then Ghat^T Ghat = -M^(1/2) L M^(-1/2), i.e.

       Lhat = -Ghat^T Ghat                with a PLAIN transpose,

   the uniform-grid identity verbatim.  Every one of the 13 identities
   lyapunov_solver.m relies on carries over under G -> Ghat, D -> Dhat.  When w
   is constant Ghat == G exactly, so the uniform solver is the special case and
   the regression test in verify.wls is meaningful.

   WHY L IS DEFINED, NOT DIFFERENCED.  L := -M^-1 G^T S_f G is the second-order
   finite-volume Laplacian.  A fourth-order FiniteDifferenceDerivative Laplacian
   -- what comoving_core.m:139 builds for the deterministic steady solve -- does
   NOT satisfy the Gram identity, so it would break R == 1 and with it the flat
   pedestal.  Second order here is a deliberate trade, and it is the same one the
   uniform Neumann branch already makes (lyapunov_solver.m:697).

   THE COLLOCATED DERIVATIVE.  D := M^-1 P^T S_f G, the graded analogue of the
   uniform D1 = Pav^T G1.  Two things make this the right choice and not merely
   a plausible one:
     - M D is the antisymmetric shift matrix (delta_{j,i+1} - delta_{j,i-1})/2
       INDEPENDENT of the spacing, so Dhat^T == -Dhat holds exactly.  That is
       identity #4, on which divBlk's transposed-gradient form depends
       (lyapunov_solver.m:1038-1039).
     - it evaluates to (u_{i+1} - u_{i-1})/(x_{i+1} - x_{i-1}), the natural
       graded central difference.

   ORDERING.  k = i + (j-1) Nint, x fastest -- the repo convention, matching
   lyapunov_solver.m:725 and comoving_core.m:151.  Hence Gx = I (x) G1x and
   Gy = G1y (x) I.  The y-metric factors cancel out of Ghat_x, so each hatted
   2-D operator is a Kronecker lift of a purely 1-D hatted factor.
   ============================================================================ *)

BeginPackage["GradedOps`", {"GradedGrid`"}];

nuOps1D::usage =
  "nuOps1D[gridAssoc] returns the 1-D operators for one axis: \"G\", \"P\", \
\"L\", \"D\" (unweighted, acting on plain nodal values) and \"Ghat\", \"Phat\", \
\"Lhat\", \"Dhat\" (symmetrised, acting on M^(1/2)-scaled values), plus \"M\", \
\"Sf\", \"msqrt\", \"misqrt\".";

nuOps2D::usage =
  "nuOps2D[gridX, gridY] lifts the 1-D operators to the 2-D tensor grid in \
x-fastest order and returns the hatted operators the solver uses: \"Gx\",\"Gy\",\
\"Px\",\"Py\",\"Dx\",\"Dy\",\"Lx\",\"Ly\",\"Lap\",\"Dxy\", together with \
\"mass\" (the nodal cell areas w^x_i w^y_j), \"msqrt\", \"misqrt\", the node and \
face coordinate vectors, and the grid associations.";

nuLapIdentity::usage =
  "nuLapIdentity[ops] returns the relative residual of Lap + (Gx^T Gx + Gy^T Gy) \
on the hatted operators.  Should be at round-off; it is the gate that says the \
scheme is an adjoint pair.";

nuAntisymmetry::usage =
  "nuAntisymmetry[ops] returns the relative residual of Dhat^T + Dhat for both \
axes -- identity #4, which divBlk's transposed-gradient form requires.";

Begin["`Private`"];

nuOps1D[gr_Association] :=
Module[{Nint = gr["Nint"], d = gr["d"], w = gr["w"], nF, G, P, Sf, M,
        msqrt, misqrt, sfsqrt, L, D, Ghat, Phat, Lhat, Dhat},
  nF = Nint + 1;
  (* (G u)_m = (u_m - u_{m-1})/d_m, with the Dirichlet ghosts u_0 = u_{Nint+1} = 0
     already eliminated.  Same sparsity as lyapunov_solver.m:719. *)
  G = SparseArray[Join[Table[{m, m} -> 1/d[[m]], {m, Nint}],
                       Table[{m + 1, m} -> -1/d[[m + 1]], {m, Nint}]], {nF, Nint}];
  P = SparseArray[Join[Table[{m, m} -> 1/2, {m, Nint}],
                       Table[{m + 1, m} -> 1/2, {m, Nint}]], {nF, Nint}];
  Sf = d;  M = w;
  msqrt = Sqrt[M];  misqrt = 1./msqrt;  sfsqrt = Sqrt[Sf];

  (* unweighted forms, for reference and for transforming results back *)
  L = SparseArray[-(1/M) (Transpose[G] . (Sf * G))];
  D = SparseArray[(1/M) (Transpose[P] . (Sf * G))];

  (* symmetrised: Ghat = Sf^(1/2) G M^(-1/2) *)
  Ghat = SparseArray[sfsqrt * (G . DiagonalMatrix[SparseArray@misqrt])];
  Phat = SparseArray[sfsqrt * (P . DiagonalMatrix[SparseArray@misqrt])];
  Lhat = SparseArray[-Transpose[Ghat] . Ghat];
  Dhat = SparseArray[Transpose[Phat] . Ghat];

  <|"Nint" -> Nint, "nFace" -> nF, "G" -> G, "P" -> P, "L" -> L, "D" -> D,
    "Ghat" -> Ghat, "Phat" -> Phat, "Lhat" -> Lhat, "Dhat" -> Dhat,
    "M" -> M, "Sf" -> Sf, "msqrt" -> msqrt, "misqrt" -> misqrt,
    "grid" -> gr|>];

nuOps2D[grX_Association, grY_Association] :=
Module[{ox, oy, nx, ny, Ix, Iy, Gx, Gy, Px, Py, Dx, Dy, Lx, Ly, Dxy,
        mass, xflat, yflat, xeX, xeY, yeX, yeY, gX, gY, fX, fY},
  ox = nuOps1D[grX];  oy = nuOps1D[grY];
  nx = ox["Nint"];    ny = oy["Nint"];
  Ix = IdentityMatrix[nx, SparseArray];
  Iy = IdentityMatrix[ny, SparseArray];

  (* x fastest => KroneckerProduct[y-factor, x-factor] *)
  Gx = KroneckerProduct[Iy, ox["Ghat"]];
  Gy = KroneckerProduct[oy["Ghat"], Ix];
  Px = KroneckerProduct[Iy, ox["Phat"]];
  Py = KroneckerProduct[oy["Phat"], Ix];
  Dx = KroneckerProduct[Iy, ox["Dhat"]];
  Dy = KroneckerProduct[oy["Dhat"], Ix];
  Lx = KroneckerProduct[Iy, ox["Lhat"]];
  Ly = KroneckerProduct[oy["Lhat"], Ix];
  Dxy = KroneckerProduct[oy["Dhat"], ox["Dhat"]];

  (* nodal cell areas, k-ordered -- the 2-D mass matrix *)
  mass = Flatten@Outer[Times, oy["M"], ox["M"]];

  gX = grX["gridInt"];  gY = grY["gridInt"];
  (* faces: face m of the x-axis sits between gridFull nodes m and m+1 *)
  fX = MovingAverage[grX["gridFull"], 2];
  fY = MovingAverage[grY["gridFull"], 2];

  (* node and face sample coordinates, in the orders Gx / Gy expect *)
  xflat = Flatten@ConstantArray[gX, ny];
  yflat = Flatten@Transpose@ConstantArray[gY, nx];
  xeX = Flatten@ConstantArray[fX, ny];
  xeY = Flatten@Transpose@ConstantArray[gY, Length[fX]];
  yeX = Flatten@ConstantArray[gX, Length[fY]];
  yeY = Flatten@Transpose@ConstantArray[fY, nx];

  <|"nx" -> nx, "ny" -> ny, "n" -> nx ny,
    "Gx" -> Gx, "Gy" -> Gy, "Px" -> Px, "Py" -> Py,
    "Dx" -> Dx, "Dy" -> Dy, "Lx" -> Lx, "Ly" -> Ly,
    "Lap" -> Lx + Ly, "Dxy" -> Dxy,
    "mass" -> mass, "msqrt" -> Sqrt[mass], "misqrt" -> 1./Sqrt[mass],
    "xflat" -> xflat, "yflat" -> yflat,
    "xeX" -> xeX, "xeY" -> xeY, "yeX" -> yeX, "yeY" -> yeY,
    "gridX" -> grX, "gridY" -> grY, "ops1Dx" -> ox, "ops1Dy" -> oy,
    "faceX" -> fX, "faceY" -> fY|>];

nuLapIdentity[op_Association] :=
Module[{r, s},
  r = op["Lap"] + (Transpose[op["Gx"]] . op["Gx"] + Transpose[op["Gy"]] . op["Gy"]);
  s = Max@Abs@op["Lap"]["NonzeroValues"];
  Max@Abs@Flatten@Normal[r]/s];

nuAntisymmetry[op_Association] :=
  Max @ Table[
     Module[{d = op[k], s},
        s = Max@Abs@d["NonzeroValues"];
        Max@Abs@Flatten@Normal[Transpose[d] + d]/s],
     {k, {"Dx", "Dy"}}];

End[];
EndPackage[];
