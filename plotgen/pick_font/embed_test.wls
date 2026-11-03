SetEnvironment["QT_QPA_PLATFORM" -> "offscreen"];
g[f_] := Plot[Sin[x], {x, 0, 3}, Frame -> True, ImageSize -> 250,
   FrameLabel -> {"Layering probability", "time t (min)"}, PlotLabel -> f,
   FrameStyle -> Directive[Black, FontFamily -> f, FontSize -> 11]];
Export["embed_cmu.pdf", g["CMU Serif"]];
Export["embed_times.pdf", g["Times"]];
Print[Options[$ExportFormats; "PDF"]];
