load "modular_param.m";
load "divisor_checker.m";

//Compute the 2-dimensional rational space of + cusp forms over Q(sqrt(2))
M := CuspForms(31^2);
B := Basis(M);
T := HeckeOperator(M, 3);
killer := T^2 + 2*T - 1;
ker := Kernel(killer);
ks := Basis(ker);
fs := [];
fs := [&+[Eltseq(k)[i]*B[i] : i in [1 .. #B]] : k in ks];

//The shadow point is 232D3 + 24D12 - D', where D' is a cuspidal divisor of degree 512
print(T2PointSq(31));

//Points in the upper half plane of discriminant 3 and 12
tau1 := (((-1 + Sqrt(-3))/2) - 521)/31^2;
tau2 := (Sqrt(-3) - 82)/31^2;
ints1 := modularIntegral(fs, tau1 : tol := 1e-8);
ints2 := modularIntegral(fs, tau2 : tol := 1e-8);

//The rational torsion in this isogeny factor is trivial
print(torsBound(B, 10));

//The above point can be rescaled by a factor of 8
pt := [2*(29 * ints1[i] + 3 * ints2[i]) : i in [1, 2]];
A := ModularAbelianVariety("961A");
Ps := Periods(A, 3070);

print(latticeCheck(pt, Ps));
