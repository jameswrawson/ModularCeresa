load "modular_param.m";

//Compute a basis of cusp forms of level 37^2, and isolate a 1-dimensional rational subspace of the Atkin-Lehner + subspace
M := CuspForms(37^2);
B := Basis(M);
T := HeckeOperator(M, 5);
killer := T - 1;
ker := Kernel(killer);
ks := Basis(ker);
fs := [&+[Eltseq(k)[i]*B[i] : i in [1 .. #B]] : k in ks];
print(fs);

//The 6*shadow of T2 is 6*194*(2D0 + D4 + 2D7) - 82*(D' - 4D3 - 3D4) - 2(3D' - 12D12 - 3D4 - 6D16) + 6*(D' - 4D3 - 3D4), for a cuspidal divisor D' of degree 1178
//Up to a cuspidal divisor (hence torsion), the shadow is 304D3 + 1398D4 + 2328D7 + 24D12 + 12D16 - 4066c, for a cusp c
//We can check 152D3 + 699D4 + 1164D7 + 12D12 + 6D16 - 2033c

//The points in the upper half plane corresponding to D3, D4, D7, D12, D16
tau1 := (581 + ((-1 + Sqrt(-3))/2))/1369;
tau2 := (117 + Sqrt(-1))/1369;
tau3 := (398 + ((-1 + Sqrt(-7))/2))/1369;
tau4 := (206 + Sqrt(-3))/1369;
tau5 := (334 + 2*Sqrt(-1))/1369;
ints1 := modularIntegral(fs, tau1 : tol := 1e-5);
ints2 := modularIntegral(fs, tau2 : tol := 1e-5);
ints3 := modularIntegral(fs, tau3 : tol := 1e-5);
ints4 := modularIntegral(fs, tau4 : tol := 1e-5);
ints5 := modularIntegral(fs, tau5 : tol := 1e-5);

//Torsion bound for the isogeny factor of J0 corresponding to the above subspace
print(torsBound(fs, 20));

//Rescale the point above to clear the torsion. The extra factor of 2 is from having 2 points in each D_i, which are swapped by Atkin--Lehner (and so map to the same point in the + quotient)
pt := [1*2*(152*ints1[1] + 699*ints2[1] + 1164*ints3[1] + 12*ints4[1] + 6*ints5[1])];
print(pt);
A := ModularAbelianVariety("1369B");
Ps := Periods(A, 4463);

print(latticeCheck(pt, Ps : tol := 1e-1));
