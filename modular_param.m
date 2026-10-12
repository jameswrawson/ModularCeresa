//Given a list of cusp forms, f, and a point in the upper halfplane tau, compute the integral of f from infinity to tau, using the first n terms of the q-expansion so that the value is correct up to the value of tol
//NB tol cannot be smaller than 1e-30, and probably should be at least a little larger to allow for any rounding errors
function modularIntegral(f, tau : tol := 1e-10)
	CC := ComplexField(30); //instantiate a complex field with 30 decimal digits of accuracy
	q := Exp(2*Pi(CC)*CC.1 * tau);
	abs_q := Abs(q);
	n := Ceiling(Log((1 - abs_q) * tol / 2) / Log(abs_q)); //The number of terms required to compute the integral sufficiently accurately, based on a comparison with geometric series (Section 6.3)

	ints := [CC ! 0 : i in [1..#f]];
	for i := 1 to n do
		for j := 1 to #f do
			ints[j] := ints[j] + (Coefficient(f[j], i) / i) * q^i;
		end for;
	end for;
	
	return ints;
end function;

//return true if point is in the period lattice, within the stated tolerance
function latticeCheck(point, periods : tol := 1e-4)
	latmat := [];
	N := #periods;

	//Build a matrix whose entries are the basis for the period lattice
	for i := 1 to N do
		latent := [];
		for j := 1 to Floor(N / 2) do
			latent := latent cat [Re(periods[i][j]), Im(periods[i][j])];
		end for;
		Append(~latmat, latent);
	end for;
	
	latmat := Matrix(RealField(30), N, N, latmat);
	
	//Determine the coordinates of the input point in terms of the lattice basis
	test_pt := [];
	for j := 1 to Floor(N / 2) do
		test_pt := test_pt cat [Re(point[j]), Im(point[j])];
	end for;
	
	sol, spc := Solution(latmat, Vector(RealField(30), N, test_pt));

	//Check if this coordinate vector is (approximately) integral, i.e. the vector is (approximately) in the lattice
	good := true;
	for i := 1 to N do
		if Abs(sol[i] - Round(sol[i])) gt tol then
			good := false;
		end if;
	end for;
	
	return good;
end function;

//Compute an upper bound for the order of rational torsion on X_0(N) using primes between 3 and the prime bound given, along with a basis for the cusp forms on X0(N)
//This functions by computing the order of the torsion modulo p (via characteristic polynomials of T_p), and taking the GCD of these values
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
		if (level mod p) eq 0 then
			continue;
		end if;
		
		//Compute the matrix of the Hecke operators in this basis
		img_basis := [Coordinates(V, psi(HeckeOperator(p, b))) : b in basis];
		hecke_mat := Matrix(Integers(), #basis, #basis, img_basis);
		f := CharacteristicPolynomial(hecke_mat);

		//Evaluate the characteristic polynomial of the Hecke operator at p + 1 (this is the torsion order at p by Eichler--Shimura)
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
