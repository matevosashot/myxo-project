(* Compare MaTeX preambles against the revtex body font for \mathcal P *)
SetEnvironment["QT_QPA_PLATFORM" -> "offscreen"];
<< MaTeX`
s = "\\text{Layering probability}~~\\mathcal P(r) / \\mathcal P(\\infty),~ \\ell_{\\rm n}~(\\upmu\\rm m)";
preambles = {
  "newtx (current)" -> {"\\usepackage{amsmath,siunitx,newtxtext,newtxmath}"},
  "CM default (revtex)" -> {"\\usepackage{amsmath,siunitx,upgreek}"}};
rows = Table[{Style[p[[1]], 10], MaTeX[s, FontSize -> 11, "Preamble" -> p[[2]]]}, {p, preambles}];
Export["pick_font_mathematica.pdf", Grid[rows, Alignment -> Left, Spacings -> {2, 2}]];
