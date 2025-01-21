load "modular_param.m";

M := CuspForms(7^3);
B := Basis(M);
T := HeckeOperator(M, 2);
killer := T^3  + 4*T^2 + 3*T - 1;
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

tau1 := (((-1 + Sqrt(-3))/2) - 18)/343;
tau2 := (Sqrt(-3) - 37)/343;
tau3 := (((-3 + 3*Sqrt(-3))/2) - 54)/343;
ints1 := modularIntegral(fs, tau1);
ints2 := modularIntegral(fs, tau2);
ints3 := modularIntegral(fs, tau3);

pt := [(35 * ints1[i] + 25 * ints2[i] + 2 * ints3[i]) : i in [1, 2, 3]];
print(pt);
A := ModularAbelianVariety("343A");
Ps := Periods(A, 1757);

print(latticeCheck(pt, Ps));
