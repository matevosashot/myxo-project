(* ::Package:: *)

(* ============================================================================
   graded_grid.m -- 1-D tensor grids with a varying cell size, for the
   fluctuation solver.

   WHY A GRADED GRID.  Two constraints fight on a uniform mesh:

     stability   maxReLambda collapses onto h alone and changes sign between
                 h = 0.294 and h = 0.229 (measured, box15 series).  The unstable
                 mode is 85% Q2, 59% of its power inside r < 3, 90% inside r < 8
                 (claude_experiments/step0_unstable_mode.wls), so it is a CORE
                 object: only the core needs h <~ 0.23.

     pedestal    the far field must reproduce eq:iso-gauss.  With the Gaussian
                 quadrature smoother of graded_noise.m that needs only
                 h <~ lNoise/Sqrt[2] for spectral accuracy, i.e. hMax ~ 1.2 at
                 lNoise = 1.73 -- five times coarser than the core.

   Holding the core spacing everywhere costs N ~ 348 at box40 (n ~ 363k, a 1.1 TB
   dump, Bartels-Stewart impossible).  Grading turns that into N ~ 150.

   THE LAW.  h(x) = Min[hCore + growth Max[0, |x| - rFine], hMax, wall term],
   equidistributed
   exactly as comoving_core.m:87 cdGridPade does it: xi(x) = Int_0^x dx'/h(x') is
   the cell index, so sampling x at uniform xi gives spacing h(x) by construction.
   The point count is DERIVED from the law, not dialled in.

   REFINING AT THE WALL TOO.  Coarsening monotonically outward is wrong for this
   problem.  dean.tex:1670: "The dQ sector does [have a boundary layer]: its drift
   carries c_L L together with a mass term, so Var(dQ) has a genuine boundary
   layer whose width is the nematic correlation length and which survives h->0,
   while Var(drho) cannot have one at all."  A law that reaches hMax at the wall
   therefore puts its COARSEST cells exactly where a layer of width elld = 0.8
   lives, and under-resolves Var(dQ) there -- while Var(drho), which has no such
   layer, is unaffected.  So the law carries an optional third branch,

       hWall + growthWall Max[0, fbox - wallBand - |x|],

   which pulls the spacing back down to hWall within wallBand of the wall.
   hWall -> hMax (the default) switches it off and recovers the two-branch law.

   CONVENTION.  Returns the FULL grid including the two walls at +-fbox.  The
   Lyapunov state lives on the interior, gridFull[[2 ;; -2]], exactly as the
   uniform Dirichlet branch of lyapunov_solver.m:669-673.  Nint = n - 2 is odd
   and 0 is an interior node, which the pedestal diagnostics require.
   ============================================================================ *)

BeginPackage["GradedGrid`"];

nuGrid::usage =
  "nuGrid[fbox, hCore, hMax, rFine, growth, refine:1] returns an association \
describing a symmetric graded 1-D grid on [-fbox, fbox] whose spacing follows \
h(x) = Min[hCore + growth Max[0, |x| - rFine], hMax, wall term].  Options \
\"hWall\" (Automatic = hMax, i.e. no wall refinement), \"wallBand\" and \
\"growthWall\" refine the grid near the wall as well, which the dQ boundary \
layer of width elld needs (dean.tex:1670).  Keys: \"gridFull\" (all n = Nint+2 \
nodes, walls included), \"gridInt\" (the Nint interior nodes), \"Nint\", \
\"d\" (the Nint+1 face spacings), \"w\" (the Nint interior dual-cell widths), \
\"hMin\", \"hMax\", \"cells\", \"hOf\".";

nuGridUniform::usage =
  "nuGridUniform[fbox, Nint] returns the same association for the UNIFORM \
node-centred Dirichlet grid h = 2 fbox/(Nint+1), so the graded machinery can be \
regression-tested against lyapunov_solver.m.";

nuGridFromNodes::usage =
  "nuGridFromNodes[gridFull] builds the association from an explicit symmetric \
node list whose first and last entries are the walls.";

Begin["`Private`"];

(* d and w from the node list.  w_i is the dual cell of interior node i, i.e.
   the trapezoid weight (x_{i+1} - x_{i-1})/2 -- the same quantity cdGridPade
   returns, and exactly the mass matrix the adjoint pair needs. *)
nuGridFromNodes[gridFull_List] :=
Module[{g = N[gridFull], n, Nint, d, w},
  n = Length[g];
  If[n < 5, Message[nuGridFromNodes::short, n]; Return[$Failed]];
  Nint = n - 2;
  d = Differences[g];                                  (* Nint+1 faces *)
  w = Table[(g[[i + 2]] - g[[i]])/2, {i, Nint}];       (* interior dual cells *)
  <|"gridFull" -> g, "gridInt" -> g[[2 ;; -2]], "Nint" -> Nint,
    "d" -> d, "w" -> w,
    "hMin" -> Min[d], "hMax" -> Max[d], "fbox" -> Last[g],
    "uniform" -> (Max[d] - Min[d] <= 1.*^-12 Max[d])|>];

nuGridFromNodes::short = "a grid needs at least 5 nodes; got `1`.";

nuGridUniform[fboxV_?NumericQ, nInt_Integer] :=
Module[{h = (2 N[fboxV])/(nInt + 1)},
  Join[nuGridFromNodes[Table[-N[fboxV] + k h, {k, 0, nInt + 1}]],
       <|"cells" -> N[nInt + 1], "hOf" -> (h &)|>]];

Options[nuGrid] = {"hWall" -> Automatic, "wallBand" -> 2., "growthWall" -> 0.35};

nuGrid[fboxV_?NumericQ, hCore_?NumericQ, hMaxv_?NumericQ, rFine_?NumericQ,
       growth_?NumericQ, refine_: 1, OptionsPattern[]] :=
Module[{hOf, nsub = 20000, xs, dx, invh, cum, tot, m, xOf, gpos, g,
        hW, band, gW},
  hW = OptionValue["hWall"] /. Automatic -> hMaxv;
  band = OptionValue["wallBand"];
  gW = OptionValue["growthWall"];
  hOf[x_] := Min[hCore + growth Max[0., Abs[x] - rFine], hMaxv,
                 hW + gW Max[0., N[fboxV] - band - Abs[x]]];

  (* xi(x) = Int_0^x dx'/h on a fine subgrid, by trapezoid -- cdGridPade:101 *)
  xs = N@Subdivide[0., N[fboxV], nsub];
  dx = N[fboxV]/nsub;
  invh = 1./(hOf /@ xs);
  cum = Prepend[Accumulate[(Most[invh] + Rest[invh]) dx/2], 0.];
  tot = Last[cum];
  (* m+1 nodes on [0, fbox]; mirroring gives 2m+1 total, so Nint = 2m-1 is ODD
     and 0 is an interior node. *)
  m = Max[3, Ceiling[refine tot]];

  xOf = Interpolation[Transpose@{cum, xs}, InterpolationOrder -> 3];
  gpos = xOf /@ (Range[0, m] tot/m);
  gpos[[1]] = 0.; gpos[[-1]] = N[fboxV];
  g = Join[-Reverse[Rest[gpos]], gpos];

  Join[nuGridFromNodes[g], <|"cells" -> tot, "hOf" -> hOf|>]];

End[];
EndPackage[];
