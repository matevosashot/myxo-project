params = {S0 -> 1, l -> 1};
W = 10;
S[r_] := S0 f[r/l]
f[x_] := x Sqrt[(0.34 + 0.07 x^2)/(1 + 0.41 x^2 + 0.07 x^4)]
\[Theta][\[Phi]_] := \[Phi]/2
Q[x_, y_] := {{S[r] Cos[2 \[Theta][\[Phi]]], 
    S[r] Sin[2 \[Theta][\[Phi]]]}, {S[r] Sin[
      2 \[Theta][\[Phi]]], -S[r] Cos[2 \[Theta][\[Phi]]]}} /. {r -> 
    Sqrt[x^2 + y^2], \[Phi] -> ArcTan[x, y]}

A[x_, y_] := 3 IdentityMatrix[2] + Q[x, y]
\[CapitalLambda][x_, y_] := 3 IdentityMatrix[2] + 0 Q[x, y]


G[x1_, y1_, x2_, 
   y2_] = (A[x1, 
        y1] Exp[-{x1 - x2, y1 - y2} . \[CapitalLambda][x1, 
            y1] . {x1 - x2, y1 - y2}/2]
      + A[x2, 
        y2] Exp[-{x1 - x2, y1 - y2} . \[CapitalLambda][x2, 
            y2] . {x1 - x2, y1 - y2}/2]) /. params // Simplify;