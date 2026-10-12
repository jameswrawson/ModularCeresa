//Compute the shadow of T2 for X_0(p)
//The output is a list of tuples. The first value of each tuple is the CM discriminant (0 represents cusps, rather than infinity in the paper), and the second is the multiplicity
//NB: The degree calculations all make use of the fact that D_n has degree 2 in prime level, for all n where it is non-zero
function T2Point(p)
	//Compute the fixed point divisor. For prime level this is 2D_infinity + D_4 + D_8 + 2D_7
	F := [<0, 2>];
	if LegendreSymbol(-1, p) eq 1 then
		Append(~F, <4, 1>);
	end if;
	
	if LegendreSymbol(-2, p) eq 1 then
		Append(~F, <8, 1>);
	end if;

	if LegendreSymbol(-7, p) eq 1 then
		Append(~F, <7, 2>);
	end if;
	
	//Set up lists for the shadow divisor, the canonical, and the image of the canonical divisor under T2
	D := [];
	KC := [];
	T2KC := [];
	
	//for each residue class, use the formula for the canonical divisor (or its multiple), and thus find its image under T2. Then, scale the fixed point divisor accordingly
	case (p mod 12):
		when 11:
			KC := [<0, (p - 11)/12>];
			T2KC := [<0, (p - 11)/4>];
			D := [<t[1], (p - 11) * t[2]/6> : t in F];
		when 7:
			KC := [<0, (p - 11)/4>, <3, -2>];
			T2KC := [<0, 3*(p - 11)/4>, <12, -6>];
			D := [<t[1], (p - 19) * t[2]/2> : t in F];
		when 5:
			KC := [<0, (p - 11)/6>, <4, -1>];
			T2KC := [<0, (p - 11)/2>, <4, -1>, <16, -2>];
			D := [<t[1], (p - 17) * t[2]/3> : t in F];
		when 1:
			KC := [<0, (p - 11)/2>, <3, -4>, <4, -3>];
			T2KC := [<0, 3*(p - 11)/2>, <12, -12>, <4, -3>, <16, -6>];
			D := [<t[1], (p - 25) * t[2]> : t in F];
	end case;
	
	//total degree of the fixed point divisor
	tot := 2*&+[t[2] : t in F];
	//Calculate D as (2g - 2)F - (deg F - 6)K_C - 2(T_2 K_C)
	D := D cat [<t[1], -(tot - 6) * t[2]> : t in KC] cat [<t[1], -2*t[2]> : t in T2KC];

	//Simplify D by combining any duplicate terms
	keys := [];
	flatD := [];
	for tup in D do
		if tup[1] in keys then
			flatD[Index(keys, tup[1])] := flatD[Index(keys, tup[1])] + tup[2];
		else
			Append(~keys, tup[1]);
			Append(~flatD, tup[2]);
		end if;
	end for;
	
	outD := [];
	for i := 1 to #keys do
		if flatD[i] ne 0 then
			Append(~outD, <keys[i], flatD[i]>);
		end if;
	end for;
	
	return outD;	
end function;

//Compute the shadow of T2 for level p^2
//The multiple on infinity (0) is half its degree, rather than its multiplicity. See comments above T2Point for some insight on this
function T2PointSq(p)
	//Compute the fixed point divisor
	F := [<0, (p + 1)>];
	if LegendreSymbol(-1, p) eq 1 then
		Append(~F, <4, 1>);
	end if;
	
	if LegendreSymbol(-2, p) eq 1 then
		Append(~F, <8, 1>);
	end if;

	if LegendreSymbol(-7, p) eq 1 then
		Append(~F, <7, 2>);
	end if;
	
	D := [];
	KC := [];
	T2KC := [];
	
	//Compute six times the canonical (and its image under T2)
	case (p mod 12):
		when 11:
			KC := [<0, (p - 6)*(p + 1)/2>];
			T2KC := [<0, 3*(p - 6)*(p + 1)/2>];
			D := [<t[1], (p - 6)*(p + 1) * t[2]> : t in F];
		when 7:
			KC := [<0, (p - 6)*(p + 1)/2>, <3, -4>];
			T2KC := [<0, 3*(p - 6)*(p + 1)/2>, <12, -12>];
			D := [<t[1], ((p - 6)*(p + 1) - 8) * t[2]> : t in F];
		when 5:
			KC := [<0, (p - 6)*(p + 1)/2>, <4, -3>];
			T2KC := [<0, 3*(p - 6)*(p + 1)/2>, <4, -3>, <16, -6>];
			D := [<t[1], ((p - 6)*(p + 1) - 6)* t[2]> : t in F];
		when 1:
			KC := [<0, (p - 6)*(p + 1) / 2>, <3, -4>, <4, -3>];
			T2KC := [<0, 3*(p - 6)*(p + 1) / 2>, <12, -12>, <4, -3>, <16, -6>];
			D := [<t[1], ((p - 6)*(p + 1) - 14) * t[2]> : t in F];
	end case;
	
	//total degree of the fixed point divisor, and then use it to compute the shadow point
	tot := 2*&+[t[2] : t in F];
	D := D cat [<t[1], -(tot - 6) * t[2]> : t in KC] cat [<t[1], -2*t[2]> : t in T2KC];

	//Flatten D (i.e. remove any redundant terms)
	keys := [];
	flatD := [];
	for tup in D do
		if tup[1] in keys then
			flatD[Index(keys, tup[1])] := flatD[Index(keys, tup[1])] + tup[2];
		else
			Append(~keys, tup[1]);
			Append(~flatD, tup[2]);
		end if;
	end for;

	outD := [];
	for i := 1 to #keys do
		if flatD[i] ne 0 then
			Append(~outD, <keys[i], flatD[i]>);
		end if;
	end for;
	
	return outD;	
end function;

//Compute whether the reduction of the divisor modulo a prime is trivial. It returns the primes where it is trivial as a list, with multiplicity. The parameter ssj is the (supersingular) j-invariant to be used in the test (cf Proposition 5.3, 5.4).
//A second j-invariant can also be specified, which can be advantageous for large degree divisors: instead of factorising the resulting value, it factorises the GCD of two values.
function jTest(divisor, ssj : ssj2 := 0)
	//Create an array whose values are the j-invariants with CM discriminant stored in the key
	js := AssociativeArray();
	js[3] := 0;
	js[4] := 1728;
	js[7] := -3375;
	js[8] := 8000;
	js[12] := 54000;
	js[16] := 287496;

	tot1 := 1;
	for D in divisor do
		tot1 := tot1 * (ssj - js[D[1]])^(D[2]);
	end for;
	K<a> := Parent(ssj);
	
	//Only the numerator is necessary
	tot1 := Numerator(tot1);
	//Isolate the Galois anti-invariant component of the numerator (i.e. what cannot be corrected by subtracting integers)
	tot1 := (tot1 - Conjugate(tot1)) / (a - Conjugate(a));
	tot1 := Integers() ! Numerator(tot1);
	out := tot1;
	
	//If a second j-invariant is provided, repeat the above, and take the GCD
	if ssj2 ne 0 then
		L<b> := Parent(ssj2);
		tot2 := 1;
		for D in divisor do
			tot2 := tot2 * (ssj2 - js[D[1]])^(D[2]);
		end for;
		tot2 := Numerator(tot2);
		tot2 := Numerator((tot2 - Conjugate(tot2)) / (b - Conjugate(b)));
		tot2 := Integers() ! tot2;

		out := GCD(tot1, tot2);
	end if;
	
	return Factorisation(out);
end function;

//This carries out the same calculations as jTest (namely those of Proposition 5.3) but is tailored to the calculation for a single level/prime. If only a single level is under consideration, this is the faster calculation
//It returns true if the divisor is non-zero, false otherwise
function jTest_modp(divisor, ssj)
	js := AssociativeArray();
	js[3] := 0;
	js[4] := 1728;
	js[7] := -3375;
	js[8] := 8000;
	js[12] := 54000;
	js[16] := 287496;

	tot1 := Parent(ssj) ! 1;
	for D in divisor do
		if D[1] ne 0 then
			tot1 := tot1 * (ssj - js[D[1]])^(Integers() ! D[2]);
		end if;
	end for;
	
	err := tot1 - Frobenius(tot1);
	
	return (err ne 0);	
end function;

//Use T2 to check if a X0(p) has non-trivial Ceresa cycle.
//It returns two values: the first is whether there are any supersingular j-invariants that are not Fp rational, and the second is whether the divisor is non-zero.
function checkPrime(p)
	//Compute the shadow
	D := T2Point(p);
	if D eq [] then
		return true, false;
	end if;
	
	//Check if D is as reduced as possible, if not, then rescale
	new_D := [];
	scale := Integers() ! D[1][2];
	for pt in D do
		scale := GCD(scale, Integers() ! pt[2]);
		if pt[1] ne 0 then
			Append(~new_D, pt);
		end if;
	end for;
	
	D := new_D;
	for i := 1 to #D do
		D[i][2] := Integers() ! (D[i][2] / scale);
	end for;
	
	//Generate the supersingular j-invariants
	ssp := SupersingularPolynomial(p);
	K := GF(p^2);
	rts := Roots(ssp, K);
	has_quad_ss := false;
	
	for rt in rts do
		//test if the j-invariant is defined over Fp
		err := rt[1] - Frobenius(rt[1]);
		if err ne 0 then
			has_quad_ss := true;
			//Check if the T2 shadow is trivial mod p
			ctpr := jTest_modp(D, rt[1]);
			if ctpr then 
				return true, true;
			end if;
		end if;
	end for;

	return has_quad_ss, false;
end function;
