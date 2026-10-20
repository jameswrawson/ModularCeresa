load "modular_param.m";

//Compute a basis of cusp forms of level 17^2, and isolate the 3-dimensional rational subspace of the Atkin-Lehner + subspace
M := CuspForms(17^2);
B := Basis(M);
T := HeckeOperator(M, 5);
killer := T - 2;
ker := Kernel(killer);
ks := Basis(ker);
fs := [&+[Eltseq(k)[i]*B[i] : i in [1 .. #B]] : k in ks];
print(fs);

//The 6*shadow of T2 is 6*32*(D4 + D8 + 2D0) - 40*(D' - 3D4) - 2(3D' - 3D4 - 6D16) + 6*(D' - 3D4), for a cuspidal divisor D' of degree 198 
//Up to a cuspidal divisor (hence torsion), the shadow is 300D4 + 192D8 + 12D16 - 1008c, for a cusp c
//We can check 25D4 + 16D8 + D16 - 84c

//The points in the upper half plane corresponding to D4, D8 and D16
tau1 := (38 + Sqrt(-1))/289;
tau2 := (24 + Sqrt(-2))/289;
tau3 := (76 + 2*Sqrt(-1))/289;
ints1 := modularIntegral(fs, tau1);
ints2 := modularIntegral(fs, tau2);
ints3 := modularIntegral(fs, tau3);

//Torsion bound for the isogeny factor of J0 corresponding to the above subspace
print(torsBound(fs, 20));

//Rescale the point above to clear the torsion. The extra factor of 2 is from having 2 points in each D_i, which are swapped by Atkin--Lehner (and so map to the same point in the + quotient)
pt := [4*2*(25*ints1[1] + 16*ints2[1] + ints3[1])];
print(pt);
A := ModularAbelianVariety("289A");
Ps := Periods(A, 1268);

print(latticeCheck(pt, Ps));
