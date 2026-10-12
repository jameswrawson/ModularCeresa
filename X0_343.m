load "modular_param.m";

//Compute a basis of cusp forms of level 7^3, and isolate the 3-dimensional rational subspace of the Atkin-Lehner + subspace
M := CuspForms(7^3);
B := Basis(M);
T := HeckeOperator(M, 2);
killer := T^3  + 4*T^2 + 3*T - 1;
ker := Kernel(killer);
ks := Basis(ker);
fs := [&+[Eltseq(k)[i]*B[i] : i in [1 .. #B]] : k in ks];

//The 6*shadow of T3 is 300(2D0 + D3 + D12) - 32(D' - 4D3) - 2(4D' - 4D3 - 12D27) + 8(D' - 4D3), for a cuspidal divisor D'
//Up to a cuspidal divisor (hence torsion), the shadow is -1456c + 404 D3 + 300D12 + 24D27, for a cusp c
//We can check 101D3 + 75D12 + 6D27 - 364c

//The points in the upper half plane corresponding to D3, D12 and D27
tau1 := (((-1 + Sqrt(-3))/2) - 18)/343;
tau2 := (Sqrt(-3) - 37)/343;
tau3 := (((-3 + 3*Sqrt(-3))/2) - 54)/343;
ints1 := modularIntegral(fs, tau1);
ints2 := modularIntegral(fs, tau2);
ints3 := modularIntegral(fs, tau3);

//Torsion bound for the isogeny factor of J0 corresponding to the above subspace
print(torsBound(fs, 20));

//Rescale the point above to clear the torsion. The extra factor of 2 is from having 2 points in each D_i, which are swapped by Atkin--Lehner (and so map to the same point in the + quotient)
pt := [1*2*(101 * ints1[i] + 75 * ints2[i] + 6 * ints3[i]) : i in [1, 2, 3]];
print(pt);
A := ModularAbelianVariety("343A");
Ps := Periods(A, 1757);

print(latticeCheck(pt, Ps));
