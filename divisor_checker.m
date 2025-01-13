function T2Point(p)
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
	
	D := [];
	KC := [];
	T2KC := [];
	
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
	
	tot := 2*&+[t[2] : t in F];
	D := D cat [<t[1], -(tot - 6) * t[2]> : t in KC] cat [<t[1], -2*t[2]> : t in T2KC];

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

function T2PointSq(p)
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
	
	print(F);
	
	D := [];
	KC := [];
	T2KC := [];
	
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

	tot := 2*&+[t[2] : t in F];
	D := D cat [<t[1], -(tot - 6) * t[2]> : t in KC] cat [<t[1], -2*t[2]> : t in T2KC];

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

function jTest(divisor, ssj : ssj2 := 0)
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
	print(a^2);
	
	tot1 := Numerator(tot1);
	tot1 := (tot1 - Conjugate(tot1)) / (a - Conjugate(a));
	print(tot1);
	tot1 := Integers() ! Numerator(tot1);
	out := tot1;

	if ssj2 ne 0 then
		L<b> := Parent(ssj2);
		print(b^2);
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

function checkPrime(p)
	D := T2Point(p);
	if D eq [] then
		return true, false;
	end if;
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
	ssp := SupersingularPolynomial(p);
	K := GF(p^2);
	rts := Roots(ssp, K);
	has_quad_ss := false;
	for rt in rts do
		err := rt[1] - Frobenius(rt[1]);
		if err ne 0 then
			has_quad_ss := true;
			ctpr := jTest_modp(D, rt[1]);
			if ctpr then 
				return true, true;
			end if;
			print("Bad good root");
		end if;
	end for;

	return has_quad_ss, false;
end function;
