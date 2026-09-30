#!/usr/bin/env wolframscript

covFile = If[Length[$ScriptCommandLine] > 1, $ScriptCommandLine[[2]],
   "./data/fluctuations-set-2/covariance_box20_N31.bin"];
outputDir = If[Length[$ScriptCommandLine] > 2, $ScriptCommandLine[[3]], "./output"];
plotBox = If[Length[$ScriptCommandLine] > 3, ToExpression[$ScriptCommandLine[[4]]], 9];

meta = Import[StringReplace[covFile, ".bin" -> "_meta.m"]];
{n, Nint, fbox, h, gridInt} = meta /@ {"n", "Nint", "fbox", "h", "gridInt"};
gridInt = N[gridInt];
tag = "box" <> ToString[fbox] <> "_N" <> ToString[Nint];

Print["Reading ", covFile, "  n = ", n, ", Nint = ", Nint, ", fbox = ", fbox];
Cmat = ArrayReshape[BinaryReadList[covFile, "Real64"], {n, n}];

ng = Nint^2;
blockRange[b_] := (b - 1) ng + 1 ;; b ng;
block[a_, b_] := Cmat[[blockRange[a], blockRange[b]]];
{Crr, CrQ1, CrQ2} = block[1, #] & /@ {1, 2, 3};
{CQ1r, CQ1Q1, CQ1Q2} = block[2, #] & /@ {1, 2, 3};
{CQ2r, CQ2Q1, CQ2Q2} = block[3, #] & /@ {1, 2, 3};

sigRhoMat = Partition[Diagonal[Crr], Nint];
sigQ1Mat = Partition[Diagonal[CQ1Q1], Nint];
sigQ2Mat = Partition[Diagonal[CQ2Q2], Nint];
sigRhoFn = ListInterpolation[Transpose[sigRhoMat], {gridInt, gridInt},
   InterpolationOrder -> 2];

Print["Var(d rho) min = ", Min@sigRhoMat, "  max = ", Max@sigRhoMat,
      "  mean = ", Mean@Flatten@sigRhoMat];
Print["Var(d Q1) MinMax = ", MinMax@sigQ1Mat, "   Var(d Q2) MinMax = ", MinMax@sigQ2Mat];

{xLo, xHi} = {Max[-plotBox, First@gridInt], Min[plotBox, Last@gridInt]};
plotPts = Min[Nint, 200];
lbl = "Var(\[Delta]\[Rho])(r)  [Dirichlet BC, N=" <> ToString[Nint] <>
      ", h=" <> ToString@NumberForm[h, 3] <> "]";

plot2D = DensityPlot[sigRhoFn[x, y], {x, xLo, xHi}, {y, xLo, xHi},
   PlotPoints -> plotPts, MaxRecursion -> 0,
   PlotRange -> {{xLo, xHi}, {xLo, xHi}}, PlotLegends -> Automatic,
   FrameLabel -> {"x", "y"}, PlotLabel -> lbl, ColorFunction -> "SunsetColors"];

plot3D = Plot3D[sigRhoFn[x, y], {x, xLo, xHi}, {y, xLo, xHi},
   PlotPoints -> plotPts, MaxRecursion -> 0, Mesh -> None,
   AxesLabel -> {"x", "y", "Var(\[Delta]\[Rho])"}, PlotLabel -> lbl,
   ColorFunction -> "SunsetColors"];

If[!DirectoryQ[outputDir], CreateDirectory[outputDir, CreateIntermediateDirectories -> True]];
out2D = FileNameJoin[{outputDir, "sigma_rho_dirichlet_" <> tag <> ".png"}];
out3D = FileNameJoin[{outputDir, "sigma_rho_surface_" <> tag <> ".png"}];
Export[out2D, plot2D, ImageResolution -> 150];
Export[out3D, plot3D, ImageResolution -> 150];
Print["Saved: ", out2D];
Print["Saved: ", out3D];
