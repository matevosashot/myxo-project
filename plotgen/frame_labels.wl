(* ::Package:: *)

(* frame_labels.wl -- axis labels at a set distance from the frame (2026-11-01).

   FrameLabel puts its labels about 8-11 pt (design) beyond the tick labels, and no option
   changes that gap; with a smaller ImagePadding it only clips them. frameLabels returns
   Insets for an Epilog instead, centred on the frame:
     frameLabels[{xlab, ylab}, {dx, dy}]
       the top of the x label's box dx pt below the frame, the right edge of the rotated
       y label's box dy pt left of it;
     frameLabels[{xlab, ylab, toplab}, {dx, dy, dt}]
       also a label whose bottom edge sits dt pt above the frame (a top axis).
   None skips a label. Distances are in design pt, before the 0.72 of the pdf export.
   The figure needs PlotRangeClipping -> False (the Epilog is clipped otherwise) and an
   explicit ImagePadding with room for the labels.
   Usage:  Get[FileNameJoin[{<repoRoot>, "plotgen", "frame_labels.wl"}]] *)

frameLabels[{xlab_, ylab_}, {dx_, dy_}] := frameLabels[{xlab, ylab, None}, {dx, dy, 0}];
frameLabels[{xlab_, ylab_, toplab_}, {dx_, dy_, dt_}] := {
   If[xlab === None, {}, Inset[xlab, Offset[{0, -dx}, Scaled[{0.5, 0}]], {Center, Top}]],
   (* rotated by 90 deg: the label's bottom edge faces the frame *)
   If[ylab === None, {},
    Inset[ylab, Offset[{-dy, 0}, Scaled[{0, 0.5}]], {Center, Bottom}, Automatic, {0, 1}]],
   If[toplab === None, {}, Inset[toplab, Offset[{0, dt}, Scaled[{0.5, 1}]], {Center, Bottom}]]};

(* MaTeX tick labels (2026-11-01). A MaTeX box carries about 1 pt of white margin on each
   side, which pushes tick labels away from the frame; trimSides[g, d] cuts d pt off its left
   and right (as varrho.wls does), so the labels sit where plain-text ones would.
   texTicks[vals, texfn] gives {v, label} pairs, the label texfn[str] with str the value
   printed as an integer or as given by fmt; minorTicks[vals, len] unlabelled ones. *)
trimSides[g_Graphics, d_] := Module[{pr = PlotRange /. Options[g],
    is = First@Cases[Options[g], (ImageSize -> v_) :> v]},
   Graphics[First[g], PlotRange -> {pr[[1]] + {d, -d}, pr[[2]]}, ImageSize -> is - {2 d, 0},
    AspectRatio -> Automatic, BaselinePosition -> (BaselinePosition /. Options[g])]];
tickString[v_] := If[Abs[v - Round[v]] < 10^-9, ToString[Round[v]],
   ToString[NumberForm[v, {Infinity, 1}]]];
texTicks[vals_List, texfn_, fmt_ : tickString] := {#, trimSides[texfn[fmt[#]], 1.]} & /@ vals;
minorTicks[vals_List, len_ : {0.004, 0}] := {#, "", len} & /@ vals;

(* Line widths (2026-11-02, the author's rule): no line prints thinner than about 1 px, taken
   as 0.75 pt; a line meant to look thinner than the frame (grid, axes, reference lines) is
   drawn in gray instead. f is the figure's design pt per printed pt, ImageSize / (width in
   manuscript.tex); frame and ticks take frameStyle[f], reference lines refStyle[f]. *)
thinLine[f_] := AbsoluteThickness[0.75 f];
frameStyle[f_] := Directive[Black, thinLine[f]];
refStyle[f_] := Directive[GrayLevel[0.6], thinLine[f]];
(* grid lines, including lines at zero drawn as grid lines so that they lie under the curves
   (2026-11-02, the author's later request): 0.6 pt, a darker gray than refStyle *)
gridStyle[f_] := Directive[GrayLevel[0.5], AbsoluteThickness[0.6 f]];
