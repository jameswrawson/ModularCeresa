load "modular_param.m";

//Compute a basis of cusp forms of level 17^2, and isolate the 3-dimensional rational subspace of the Atkin-Lehner + subspace
M := CuspForms(19^2);
B := Basis(M);
T := HeckeOperator(M, 3);
killer := T;
ker := Kernel(killer);
ks := Basis(ker);
fs := [&+[Eltseq(k)[i]*B[i] : i in [1 .. #B]] : k in ks];
print(fs);

//The 6*shadow of T2 is 6*42*(D8 + 2D0) - 42*(D' - 4D3) - 2(3D' - 12D12) + 6*(D' - 4D3), for a cuspidal divisor D' of degree 260 
//Up to a cuspidal divisor (hence torsion), the shadow is 144D3 + 252D8 + 24D12 - 840c, for a cusp c
//We can check 12D3 + 21D8 + 2D12 - 70c

//The points in the upper half plane corresponding to D3, D8 and D12
tau1 := (68 + ((-1 + Sqrt(-3))/2))/361;
tau2 := (63 + Sqrt(-2))/361;
tau3 := (137 + Sqrt(-3))/361;
ints1 := modularIntegral(fs, tau1);
ints2 := modularIntegral(fs, tau2);
ints3 := modularIntegral(fs, tau3);

//Torsion bound for the isogeny factor of J0 corresponding to the above subspace
print(torsBound(fs, 20));

//Rescale the point above to clear the torsion. The extra factor of 2 is from having 2 points in each D_i, which are swapped by Atkin--Lehner (and so map to the same point in the + quotient)
pt := [1*2*(12*ints1[1] + 21*ints2[1] + 2*ints3[1])];
print(pt);
A := ModularAbelianVariety("361A");
Ps := Periods(A, 1853);

print(latticeCheck(pt, Ps));
