//Given a list of cusp forms, f, and a point in the upper halfplane tau, compute the integral of f from infinity to tau, using the first n terms of the q-expansion
function modularIntegral(f, tau : tol := 1e-10)
	CC := ComplexField(30);
	q := Exp(2*Pi(CC)*CC.1 * tau);
	abs_q := Abs(q);
	n := Ceiling(Log((1 - abs_q) * tol / 2) / Log(abs_q));
	print(n);

	ints := [CC ! 0 : i in [1..#f]];
	for i := 1 to n do
		for j := 1 to #f do
			ints[j] := ints[j] + (Coefficient(f[j], i) / i) * q^i;
		end for;
	end for;
	
	return ints;
end function;

//return true if point is in the period lattice
function latticeCheck(point, periods : tol := 1e-4)
	latmat := [];
	N := #periods;
	print(N);
	for i := 1 to N do
		latent := [];
		for j := 1 to Floor(N / 2) do
			latent := latent cat [Re(periods[i][j]), Im(periods[i][j])];
		end for;
		Append(~latmat, latent);
	end for;
	
	latmat := Matrix(RealField(30), N, N, latmat);
	
	test_pt := [];
	for j := 1 to Floor(N / 2) do
		test_pt := test_pt cat [Re(point[j]), Im(point[j])];
	end for;
	
	sol, spc := Solution(latmat, Vector(RealField(30), N, test_pt));
	print(sol);
	good := true;
	for i := 1 to N do
		if Abs(sol[i] - Round(sol[i])) gt tol then
			good := false;
		end if;
	end for;
	
	return good;
end function;

function torsBound(basis, primbd)
	ps := PrimesInInterval(3, primbd);
	level := Level(basis[1]);
	S := CuspForms(level, 2);
	Sv, phi := VectorSpace(S);
	psi := Inverse(phi);
	V := sub<Sv | [psi(b) : b in basis]>;
	print("Converter computed");
	
	order := -1;
	for p in ps do
		print(p);
		if (level mod p) eq 0 then
			continue;
		end if;
		
		img_basis := [Coordinates(V, psi(HeckeOperator(p, b))) : b in basis];
		hecke_mat := Matrix(Integers(), #basis, #basis, img_basis);
		f := CharacteristicPolynomial(hecke_mat);
		bd := Evaluate(f, p + 1);
		
		if order eq -1 then
			order := bd;
		else
			order := GCD(order, bd);
		end if;
		
		if order eq 1 then
			break;
		end if;

	end for;
	return order;
end function;
