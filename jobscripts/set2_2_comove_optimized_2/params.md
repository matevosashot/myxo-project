paramRelations = {ell == Sqrt[L1/(-a1)], 
   S0 == Sqrt[(-a1)/(2 b)], \[Gamma] == -8 a1/\[Xi]r, \[Xi]r == 
    3  \[Xi]0,
   	a1 == a + \[CapitalLambda]  \[Xi]r, 
   L1 == L + (\[Zeta]  \[Xi]r)/(4  \[Xi]0), \[Gamma]0 == 
    B/\[Xi]0, \[Zeta] == 2.5  \[Xi]0};
knownParams = {ell -> 1, 
   S0 -> 0.9, \[Xi]0 -> 5, \[Rho]0 -> 0.25, \[CapitalLambda] -> 0.13};

guessedParams = {B -> 140000., \[Gamma] -> 10};

inferedParams = 
  First@Solve[
    paramRelations /. Join[knownParams, guessedParams], {b, a, L, a1, 
     L1, \[Gamma]0, \[Xi]r, \[Zeta]}];
derivedNames = {B, a, b, 
   L, \[Zeta], \[Xi]0, \[Xi]r, \[Rho]0, \[CapitalLambda], a1, 
   L1, \[Gamma]0};

modelParams = 
  Append[Thread[
    derivedNames -> (derivedNames /. 
       Join[knownParams, guessedParams, inferedParams])], 
   lNoise -> 0.5];
modelParams


Nint in [140 150 180]

B in [110000, 100000, 90000,120000, 130000]



for the steady state use
Import[FileNameJoin[repoRoot, "comoving_defect", 
   "comoving_steady_Qfree_solver.m"]];
solSteadyRhoNeumann = 
  SolveActiveNematicComovingSteadyQfree[modelParams, solverParams, 
   "RhoNeumann"];

Use Dirichlet for the Lyapunof solver.

