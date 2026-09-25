//Compute a basis for the space of level N weight 2 cusp forms to precision prec. 
//The basis is returned as two lists: the first is + forms for the Atkin-Lehner operator w_p, and the second is - forms
function ALBasis(N, p, prec)
	S := CuspForms(N);
	B := Basis(S, prec);
	M := AtkinLehnerOperator(S, p);
	
	//Compute spaces with given Atkin--Lehner sign
	spc1 := Kernel(M - ScalarMatrix(#B, 1));
	spc2 := Kernel(M + ScalarMatrix(#B, 1));
	
	b1s := Basis(spc1);
	b2s := Basis(spc2);

	//Convert the basis vectors to q-expansions
	basis1 := [&+[b[i] * B[i] : i in [1 .. #B]] : b in b1s];
	basis2 := [&+[b[i] * B[i] : i in [1 .. #B]] : b in b2s];
	
	return basis1, basis2;
end function;

//Compute a canonical model for the modular curve X_0(N), assuming it is only cut out by quadrics.
//By Petri's theorem, this is a correct model unless the curve is: hyperelliptic, trigonal or a smooth plane quintic.
//This returns 3 values: a list of the quadric equations; the coordinates fixed under the Atkin--Lehner involution w_p; the coordinates where w_p acts by -1
//By default, the precision used is an approximation to the Sturm bound for weight 4 forms of level N
function QuadricModel(N, p : prec := 1 + Integers() ! Ceiling((N + 1)/3))
	//Compute the Atkin--Lehner basis of forms
	Qq<q> := PowerSeriesRing(Rationals(), prec);
	basis1, basis2 := ALBasis(N, p, prec);
	basis := basis1 cat basis2;
	n := #basis;

	Qx<[x]> := PolynomialRing(Rationals(), n);

	//Construct a basis for homogenous degree 2 polynomials
	mons := MonomialsOfDegree(Qx, 2);
	mons := SetToSequence(mons);
	
	//Evaluate the basis on the basis of modular forms to deduce the relations
	monImages := [Evaluate(mon, basis) : mon in mons];
	monImages := [[Coefficient(m, i) : i in [1 .. prec]] : m in monImages];
	monImages := Matrix(monImages);
	BB := Basis(Kernel(monImages));
	
	//Form the relations and then clear denominators
	qds := [&+[Eltseq(b)[i]*mons[i] : i in [1 .. #mons]] : b in BB];

	for i := 1 to #qds do
		f := qds[i];
		dnm := LCM([Denominator(cf) : cf in Coefficients(f)]);
		qds[i] := qds[i]*dnm;
	end for;
	
	return qds, [x[i] : i in [1 .. #basis1]], [x[i] : i in [1 + #basis1 .. n]];
end function;

//Given the equations, the fixed variables (under Atkin--Lehner involution), the anti-invariant variables and an auxiliary prime q, compute the order of the Atkin--Lehner shadow point modulo q.
function ALShadowOrder(eqs, fix, neg, q)
	R := Parent(eqs[1]);
		
	//Construct the reduction of the modular curve
	Fpx<[x]> := PolynomialRing(GF(q), Rank(R));
	red_eq := [Fpx ! eqn : eqn in eqs];
	
	C := Curve(ProjectiveSpace(Fpx), red_eq);
	assert Genus(C) eq Rank(R);

	//We compute the divisors explicitly. To be invariant under Atkin--Lehner, either the coordinates which are fixed by Atkin--Lehner are all zero, or the coordinates which are anti-invariant are all zero.
	I := Ideal([Fpx ! f : f in fix]);
	D1 := Divisor(C,I);

	I2 := Ideal([Fpx ! f : f in neg]);
	D2 := Divisor(C,I2);

	F := D1 + D2; 

	Cl, phi, psi := ClassGroup(C);

	K := CanonicalDivisor(C);

	Sh := Degree(K)*F - Degree(F)*K; //Form the shadow divisor from the fixed point divisor and the canonical
	return Order(psi(Sh));
end function;
