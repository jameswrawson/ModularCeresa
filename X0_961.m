load "modular_param.m";

M := CuspForms(31^2);
B := Basis(M);
T := HeckeOperator(M, 3);
killer := T^2 + 2*T - 1;
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

tau1 := (((-1 + Sqrt(-3))/2) - 521)/31^2;
tau2 := (Sqrt(-3) - 82)/31^2;
ints1 := modularIntegral(fs, tau1 : tol := 1e-5);
ints2 := modularIntegral(fs, tau2 : tol := 1e-5);

pt := [(29 * ints1[i] + 3 * ints2[i]) : i in [1, 2]];
print(pt);
A := ModularAbelianVariety("961A");
Ps := Periods(A, 3070);

print(latticeCheck(pt, Ps));
