load "modular_param.m";

M := CuspForms(41^2);
B := Basis(M);
T := HeckeOperator(M, 2);
killer := T - 1;
ker := Kernel(killer);
ks := Basis(ker);
fs := [];
for k in ks do
	f := 0 * B[1];
	elts := Eltseq(k);
	for i := 1 to #elts do
		f := f + elts[i]*B[i];
	end for;
	Append(~fs, f);
end for;
print("Basis computed");

tau1 := (Sqrt(-1) - 378)/41^2;
tau2 := (Sqrt(-2) - 71)/41^2;
tau3 := (2*Sqrt(-1) - 756)/41^2;
ints1 := modularIntegral(fs, tau1 : tol := 1e-8);
ints2 := modularIntegral(fs, tau2 : tol := 1e-8);
ints3 := modularIntegral(fs, tau3 : tol := 1e-8);

pt := [(1716 * ints1[i] + 1464 * ints2[i] + 12 * ints3[i]) : i in [1, 2]];
print(pt);
A := ModularAbelianVariety("1681A");
Ps := Periods(A, 6610);

print(latticeCheck(pt, Ps));
