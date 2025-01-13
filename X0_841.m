load "modular_param.m";

M := CuspForms(29^2);
B := Basis(M);
T := HeckeOperator(M, 11);
killer := T^2 + 5*T + 5;
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

tau1 := (Sqrt(-1) - 41)/29^2;
tau2 := (((-1 + Sqrt(-7))/2) + 109)/29^2;
tau3 := (2*Sqrt(-1) + 1600)/29^2;
ints1 := modularIntegral(fs, tau1);
ints2 := modularIntegral(fs, tau2);
ints3 := modularIntegral(fs, tau3);

pt := [245 * 2 * (145 * ints1[i] + 228 * ints2[i] + 2 * ints3[i]) : i in [1, 2]];
print(pt);
A := ModularAbelianVariety("841A");
Ps := Periods(A, 4000);

print(latticeCheck(pt, Ps));
