load "modular_param.m";
load "divisor_checker.m";

//Compute space of cusp forms, and select the 2-dimensional space of rational + cusp forms over Q(sqrt(2))
M := CuspForms(41^2);
B := Basis(M);
T := HeckeOperator(M, 2);
killer := T - 1;
ker := Kernel(killer);
ks := Basis(ker);
fs := [&+[Eltseq(k)[i]*B[i] : i in [1 .. #B]] : k in ks];

//The shadow is 1716D4 + 1464D8 + 12D16 - D', where D' is cuspidal of degree 6384
print(T2PointSq(41));

//The points in the upper half plane corresponding to CM by discriminants: 4, 8, 16. Only one point is needed for each discriminant due to the action of Atkin--Lehner
tau1 := (Sqrt(-1) - 378)/41^2;
tau2 := (Sqrt(-2) - 71)/41^2;
tau3 := (2*Sqrt(-1) - 756)/41^2;

//The integrals of each element of the eigenspace from infinity to the marked points
ints1 := modularIntegral(fs, tau1 : tol := 1e-7);
ints2 := modularIntegral(fs, tau2 : tol := 1e-7);
ints3 := modularIntegral(fs, tau3 : tol := 1e-7);

//Torsion bound for this isogeny factor is 4
print(torsBound(B, 200));

pt := [4*2*(143 * ints1[i] + 122 * ints2[i] + ints3[i]) : i in [1, 2]];

A := ModularAbelianVariety("1681A");
Ps := Periods(A, 6610);

//Check if the integrals lie in the period lattice
print(latticeCheck(pt, Ps : tol := 1e-3));
