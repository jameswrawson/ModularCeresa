load "modular_param.m";

//Compute a basis of cusp forms of level 7^3, and isolate the 3-dimensional rational subspace of the Atkin-Lehner + subspace
M := CuspForms(3^5);
B := Basis(M);
T := HeckeOperator(M, 7);
killer := T + 4;
ker := Kernel(killer);
ks := Basis(ker);
fs := [&+[Eltseq(k)[i]*B[i] : i in [1 .. #B]] : k in ks];
print(fs);

//The 6*shadow of T2 is 216*(D8 + 2D0) - 38*D', for a cuspidal divisor D' of degree 36
//Up to a cuspidal divisor (hence torsion), the shadow is 216*D8 - 216c, for a cusp c
//We can check D8 - 2c

//The points in the upper half plane corresponding to D8
tau1 := (22 + Sqrt(-2))/243;
ints1 := modularIntegral(fs, tau1);

//Torsion bound for the isogeny factor of J0 corresponding to the above subspace
print(torsBound(fs, 20));

//Rescale the point above to clear the torsion. The extra factor of 2 is from having 2 points in each D_i, which are swapped by Atkin--Lehner (and so map to the same point in the + quotient)
pt := [3*2*(ints1[1])];
print(pt);
A := ModularAbelianVariety("243A");
Ps := Periods(A, 740);

print(latticeCheck(pt, Ps));
