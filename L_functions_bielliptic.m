//Compute L'(f x f x g, 2), where g is assumed to be a rank 1 cusp form, and both f and g are of level p
function tripleProd(f, g, p : prec := 5)
	L1 := LSeries(f : Precision := prec + 3);
	L2 := LSeries(g : Precision := prec + 3);

	Lsym :=	SymmetricPower(L1, 2);
	Ltrip := TensorProduct(Lsym, L2);
	
	R<x> := PolynomialRing(Rationals());
	loc := (1 + x)*(1 + p * x); //Local factor for the L-function at p, computed from e.g. Gross--Kudla
	Ltrip_g := ChangeLocalInformation(Ltrip, p, 4, loc); //Correct the local factor AND conductor
	LSetPrecision(Ltrip_g, prec); //Truncate the precision desired
	print(LCfRequired(Ltrip_g)); //Number of terms required to verify the functional equation holds
	print(CFENew(Ltrip_g)); //Verify the functional equation holds i.e. the local factor was correct
	LSetPrecision(L2, prec); 

	centmain := Evaluate(Ltrip_g, 2); //Evaluate L(g x Sym^2 f, 2)
	centsec := Evaluate(L2, 1 : Derivative := 1); //L'(g, 1)
	return centsec * centmain;
end function;

//The bi-elliptic levels of interest
ps := [43, 53, 61, 79, 83, 101, 131];

for p in ps do
	print(p);
	
	M := CuspForms(p, 2);
	N := Newforms(M); //List of new forms
	g := N[1][1]; //The unique + form
	for i := 2 to #N do
		fs := N[i][1];
		embeds := ComplexEmbeddings(fs)[1]; //Run over all embeddings of the cusp form
		for embed in embeds do
			f := embed;
			Lval := tripleProd(f, g, p : prec := 2); //Compute the central derivative to at least 2 decimal digits precision
			print("L-value is ");
			print(Lval);
		end for;
	end for;
	print("-----------");
end for;

