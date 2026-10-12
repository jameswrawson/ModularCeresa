load "modular_param.m";
load "divisor_checker.m";

//Generate space of cusp forms, and select the 2 dimensional rational subspace of + forms over Q(sqrt(5))
M := CuspForms(29^2);
B := Basis(M);
T := HeckeOperator(M, 11);
killer := T^2 + 5*T + 5;
ker := Kernel(killer);
ks := Basis(ker);
fs := [&+[Eltseq(k)[i]*B[i] : i in [1 .. #B]] : k in ks];

//Check the shadow is 870D4 + 1368D7 + 12D16 - D', D' a cuspidal divisor of degree 4500
print(T2PointSq(29));

//Points in the upper half plane with CM discriminants 4, 7 and 16
tau1 := (Sqrt(-1) - 41)/29^2;
tau2 := (((-1 + Sqrt(-7))/2) + 109)/29^2;
tau3 := (2*Sqrt(-1) + 1600)/29^2;
ints1 := modularIntegral(fs, tau1);
ints2 := modularIntegral(fs, tau2);
ints3 := modularIntegral(fs, tau3);

//The torsion bound for this isogeny factor is 1 
print(torsBound(B, 10));

//The above divisor can be rescaled by a factor of 6. It is then multiplied by 2 (to account for the Atkin--Lehner action)
pt := [2 * (145 * ints1[i] + 228 * ints2[i] + 2 * ints3[i]) : i in [1, 2]];
print(pt);
A := ModularAbelianVariety("841A");
Ps := Periods(A, 4000);

print(latticeCheck(pt, Ps));
