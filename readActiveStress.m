(* ::Package:: *)

(* ============================================================================
   Active-force statistics from a dumped covariance matrix.

       Q = {{Q1, Q2}, {Q2, -Q1}},   f_a = d_b Q_ab
       f_x = dx Q1 + dy Q2,   f_y = dx Q2 - dy Q1

   The target is <|f|^2> = <f_x^2> + <f_y^2> for the FLUCTUATION field dQ, i.e.
   the variance of the active force, on the grid.

   Why this does not need readCorrQ.  <|f|^2>(r) is a COINCIDENT-POINT quantity:
   every term is a mixed second derivative of the two-point correlator evaluated
   on the diagonal r' = r,

       G^ab_(al)(be)(r) = <d_al dQ_a(r) d_be dQ_b(r)>
                        = d_al d_be' C_ab(r, r')|_{r'=r},     a,b in {1,2}

   and the discrete d_al is the solver's own 2-point central difference (Dx =
   KroneckerProduct[Idn, D1], Dy = KroneckerProduct[D1, Idn], see
   lyapunov_solver_module.m:230-262).  So only NINE shifted diagonals of each
   Sigma block pair are ever touched,

       Delta in {0, +-2, +-2 Nint, Nint +- 1, -Nint +- 1},

   read as Sigma[[(a-1) ng + k, (b-1) ng + k + Delta]] -- the pattern covDiag
   already supports.  That is 27 vectors of length ng = 4.3 MB at Nint = 141,
   versus 2.9 GiB PER BLOCK for readCorrQ (which returned $Aborted at Nint=141).
   Block pair (3,2) is not read: Sigma is symmetric, so G^21_(al)(be) =
   G^12_(be)(al) comes for free.

   Boundary.  Reading Sigma[[k, k + p + q Nint]] is geometrically correct only
   when i+p and j+q stay in [1, Nint]; otherwise the flat index k = i + (j-1) Nint
   wraps into the neighbouring y-row.  The valid set is i,j in [2, Nint-1], so the
   returned fields live on gridF = gridInt[[2 ;; -2]] -- one cell narrower each
   side than "grid".  USE "gridF", NOT "grid", TO CLAMP PLOT RANGES for these
   fields.  Dropping that ring also removes the only place where the Dirichlet and
   Neumann discretisations differ: in the interior both reduce to the same
   (f[i+1] - f[i-1])/(2h), so no BC branching is needed here.

   UV caveat.  <|f|^2> at coincident points is UV-divergent -- it is the Q-sector
   analogue of the rho pedestal Sigma^rr = zeta rho_ss/(B rho_0 h^2) of section5.tex
   -- so its value is set by the lattice cutoff.  The stencil choice therefore
   matters at LEADING order.  That is why the 2nd-order central difference of the
   solver is reproduced exactly rather than "improved", and why numbers are only
   comparable between runs at matched h.
   ============================================================================ *)

(* covStem / readCovMeta / covBase / covMat / covFn.  $InputFileName is set while
   this file is being Get; when it is empty (pasted into a notebook) DirectoryName
   gives "" and the load falls back to $Directory. *)
If[DownValues[covStem] === {},
   Get[FileNameJoin[{DirectoryName[$InputFileName], "binaryReadTools.m"}]]];

(* ---------------------------------------------------------------------------
   covDiags: several shifted diagonals in one pass.

   Returns <| Delta -> vector of length k2-k1+1 |>, entry k being
   Sigma[[k, k + Delta]] -- covDiag generalised from one shift to a list.  The
   shifts are sorted and Split into runs whose neighbours are at most "MaxGap"
   apart, and each run costs ONE seek plus one bulk read per row instead of one
   seek per shift.  For the nine force shifts this gives five runs per row:
   {-2N} {-N-1,-N+1} {-2,0,2} {N-1,N+1} {2N}.

   "MaxGap" is the read-strategy knob.  16 means five small reads per row;
   2 Nint + 2 collapses the nine shifts into ONE contiguous 4 Nint + 3 read per
   row -- five times fewer seeks for ~40x the bytes.  Seeks dominate by a wide
   margin, so readActiveStress defaults to the latter; measured cold-cache on
   /data (Nint ~ 140, ~26 GiB dumps, different files so neither was cached):

       "MaxGap" -> 16          139 s,   81 MiB peak,    6 MiB payload
       "MaxGap" -> 2 Nint + 2  9.4 s,  152 MiB peak,  253 MiB payload

   Both give bit-identical answers (checked in claude_experiments/active_stress).

   Columns outside [1, n] are zero-filled PER ELEMENT, not per run: a run is a
   read-batching device, so clipping a whole run because one of its members left
   the matrix would silently zero valid neighbours (Delta = Nint-1 in the Q2Q2
   block is used by the interior cut at rows where Delta = Nint+1 is not).
   --------------------------------------------------------------------------- *)
Options[covDiags] = {"MaxGap" -> 16};

covDiags[stem_String, n_Integer, k1_Integer, k2_Integer, shifts_List,
         OptionsPattern[]] :=
Module[{srt, runs, bounds, str, raw},
  srt = Sort @ DeleteDuplicates @ shifts;
  runs = Split[srt, #2 - #1 <= OptionValue["MaxGap"] &];
  bounds = {First[#], Last[#]} & /@ runs;
  str = OpenRead[stem <> ".bin", BinaryFormat -> True];
  raw = WithCleanup[
     Table[
        Table[
           Module[{cLo = k + rb[[1]], cHi = k + rb[[2]], a0, b0},
              If[cHi < 1 || cLo > n,
                 ConstantArray[0., cHi - cLo + 1],
                 a0 = Max[cLo, 1]; b0 = Min[cHi, n];
                 SetStreamPosition[str, 8 ((k - 1) n + (a0 - 1))];
                 If[a0 == cLo && b0 == cHi,
                    BinaryReadList[str, "Real64", b0 - a0 + 1],
                    Join[ConstantArray[0., a0 - cLo],
                         BinaryReadList[str, "Real64", b0 - a0 + 1],
                         ConstantArray[0., cHi - b0]]]]],
           {rb, bounds}],
        {k, k1, k2}],
     Close[str]];
  KeyTake[
     Association @ Flatten @ Table[
        With[{run = runs[[r]]},
           Table[s -> Developer`ToPackedArray[raw[[All, r, s - First[run] + 1]]],
                 {s, run}]],
        {r, Length[runs]}],
     shifts]];

(* Block pair (a,b) of Sigma, banded: <| Delta -> Nint x Nint array |> with
   B[[i,j]] = Sigma[[(a-1) ng + k, (b-1) ng + k + Delta]], k = i + (j-1) Nint.
   The 4-D analogue is covBlock; this is the diagonal-band analogue. *)
covBand[stem_String, n_Integer, ng_Integer, Nint_Integer, a_Integer, b_Integer,
        shifts_List, opts : OptionsPattern[covDiags]] :=
Module[{d},
  d = covDiags[stem, n, (a - 1) ng + 1, a ng, (b - a) ng + # & /@ shifts, opts];
  Association @ Table[s -> covMat[d[(b - a) ng + s], Nint], {s, shifts}]];

(* ---------------------------------------------------------------------------
   Stencils.  A stencil is a list of {{p, q}, c}: the grid-point offset and its
   coefficient, so (d_al u)[[i,j]] = Sum c u[[i+p, j+q]].  Second-order central
   difference, matching the solver.
   --------------------------------------------------------------------------- *)
stencilPts = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};

gradStencil[h_] := <|
   "x" -> {{{ 1, 0},  1/(2 h)}, {{-1, 0}, -1/(2 h)}},
   "y" -> {{{ 0, 1},  1/(2 h)}, {{ 0,-1}, -1/(2 h)}}|>;

(* The Delta's reachable by one left and one right stencil offset:
   {0, +-2, +-2 Nint, Nint +- 1, -Nint +- 1}. *)
forceShifts[Nint_Integer] := Union @ Flatten @
   Table[(u[[1]] - v[[1]]) + (u[[2]] - v[[2]]) Nint,
         {v, stencilPts}, {u, stencilPts}];

(* G(r) = <(d_L dQ_a)(r) (d_R dQ_b)(r)> on the interior i,j in [2, Nint-1].
   The term for output point (i,j) is B^Delta[[i+pL, j+qL]] with
   Delta = (pR - pL) + (qR - qL) Nint; since pL,qL in {-1,0,1} every index stays
   inside [1, Nint].  Result is (Nint-2) x (Nint-2). *)
gradCorr[band_Association, Nint_Integer, sL_List, sR_List] :=
   Total @ Flatten[
      Table[
         With[{pL = l[[1, 1]], qL = l[[1, 2]], cL = l[[2]],
               pR = r[[1, 1]], qR = r[[1, 2]], cR = r[[2]]},
            cL cR band[(pR - pL) + (qR - qL) Nint][[
                  2 + pL ;; Nint - 1 + pL, 2 + qL ;; Nint - 1 + qL]]],
         {l, sL}, {r, sR}],
      1];

(* ---------------------------------------------------------------------------
   readActiveStress[file] -- the whole gradient-correlator set plus the force.

   "G"        <| {a, al, b, be} -> (Nint-2)x(Nint-2) array |>, all 16 keys,
              G[{a,al,b,be}] = <d_al dQ_a  d_be dQ_b>; symmetric under
              (a,al) <-> (b,be), so 10 are independent.
   "fxfx" "fyfy" "fxfy" "f2"    interpolants on gridF, with *Values arrays.
   "gridF"    gridInt[[2 ;; -2]], the domain of all of the above.

       <f_x^2>  = G11xx + 2 G12xy + G22yy
       <f_y^2>  = G22xx - 2 G21xy + G11yy
       <f_xf_y> = G12xx - G11xy + G22yx - G21yy
       <|f|^2>  = <f_x^2> + <f_y^2>
   --------------------------------------------------------------------------- *)
readActiveStress[file_String, opts : OptionsPattern[covDiags]] :=
Module[{stem = covStem[file], meta, Nint, ng, n, h, grid, gridF, shifts, gap,
        band, st, gDir, g, fxfx, fyfy, fxfy, f2},
  meta = readCovMeta[stem];
  Nint = meta["Nint"];
  ng   = Nint^2;
  n    = meta["n"];
  h    = N[meta["h"]];
  grid = N[meta["gridInt"]];
  gridF  = grid[[2 ;; -2]];
  shifts = forceShifts[Nint];
  (* One contiguous read per row unless the caller says otherwise -- see the
     covDiags header for the measurement behind this default. *)
  gap = Lookup[Association @ Flatten @ {opts}, "MaxGap", 2 Nint + 2];

  (* Sigma blocks: 1 = rho, 2 = Q1, 3 = Q2.  (2,1) is reconstructed below. *)
  band = <|{1, 1} -> covBand[stem, n, ng, Nint, 2, 2, shifts, "MaxGap" -> gap],
           {2, 2} -> covBand[stem, n, ng, Nint, 3, 3, shifts, "MaxGap" -> gap],
           {1, 2} -> covBand[stem, n, ng, Nint, 2, 3, shifts, "MaxGap" -> gap]|>;

  st = gradStencil[h];
  gDir = Association @ Flatten @ Table[
     With[{a = ab[[1]], b = ab[[2]]},
        Table[{a, al, b, be} -> gradCorr[band[ab], Nint, st[al], st[be]],
              {al, {"x", "y"}}, {be, {"x", "y"}}]],
     {ab, {{1, 1}, {2, 2}, {1, 2}}}];
  (* (a,al) <-> (b,be) supplies the missing {2, al, 1, be} keys.  gDir last so a
     directly computed value always wins over its transposed twin. *)
  g = Join[
     Association @ KeyValueMap[{#1[[3]], #1[[4]], #1[[1]], #1[[2]]} -> #2 &, gDir],
     gDir];

  fxfx = g[{1, "x", 1, "x"}] + 2 g[{1, "x", 2, "y"}] + g[{2, "y", 2, "y"}];
  fyfy = g[{2, "x", 2, "x"}] - 2 g[{2, "x", 1, "y"}] + g[{1, "y", 1, "y"}];
  fxfy = g[{1, "x", 2, "x"}] - g[{1, "x", 1, "y"}]
       + g[{2, "y", 2, "x"}] - g[{2, "y", 1, "y"}];
  f2   = fxfx + fyfy;

  Join[covBase[meta],
     <|"gridF" -> gridF, "G" -> g,
       "fxfxValues" -> fxfx, "fxfx" -> covFn[fxfx, gridF],
       "fyfyValues" -> fyfy, "fyfy" -> covFn[fyfy, gridF],
       "fxfyValues" -> fxfy, "fxfy" -> covFn[fxfy, gridF],
       "f2Values"   -> f2,   "f2"   -> covFn[f2, gridF]|>]];
