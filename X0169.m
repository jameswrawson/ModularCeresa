
// verify X0(169)

// we compute an AL invariant/anti-invariant basis for X0(169)

function ALBasis(N, p, prec)
	S := CuspForms(N);
	B := Basis(S, prec);
	M := AtkinLehnerOperator(S, p);
	spc1 := Kernel(M - ScalarMatrix(#B, 1));
	spc2 := Kernel(M + ScalarMatrix(#B, 1));
	b1s := Basis(spc1);
	b2s := Basis(spc2);

	basis1 := [];
	for b in b1s do
		f := 0;
		for i := 1 to #B do
			f := f + b[i] * B[i];
		end for;
		Append(~basis1, f);
	end for;
	
	basis2 := [];
	for b in b2s do
		f := 0;
		for i := 1 to #B do
			f := f + b[i] * B[i];
		end for;
		Append(~basis2, f);
	end for;

	return basis1, basis2;
end function;


// we use the above basis to compute a canonical model of X0(169)

Qq<q>:=PowerSeriesRing(Rationals(),100);

basis1, basis2 := ALBasis(169,169,100) ;
basis := basis1  cat basis2 ;
n := #basis;


Qx<[x]>:=PolynomialRing(Rationals(),n);

mons:=MonomialsOfDegree(Qx,2);
V:=VectorSpace(Rationals(),#mons);
monImages:=[Evaluate(mon,basis) : mon in mons];

W:=VectorSpace(Rationals(),40);
monImages:=[W![Coefficient(m,i) : i in [0..39] ] : m in monImages];

h:=hom<V->W | monImages>;

K:=Kernel(h);
K := Basis(K);
v:= [ Eltseq(V!(K[i])): i in [1..#K]];

Q:= [ &+[a[i]*mons[i] : i in [1..#mons]] : a in v ];    // qudratic cutting out our curve in P^7;
Q := [ 2*q : q in Q] ;
P := ProjectiveSpace(Qx);
C := Curve(P, Q) ;

Zx<[x]> := PolynomialRing(Integers(), #basis);
Q := [ Zx ! f : f in Q];

// mod 5
Zx<[x]> := PolynomialRing(GF(5), #basis);          
Q1 := [ Zx ! f : f in Q];
P := ProjectiveSpace(Zx);
C := Curve(P, Q) ;

I := Ideal([x[1], x[2], x[3]]) ;
D1 := Divisor(C,I) ;

I2 := Ideal([ x[4], x[5], x[6], x[7], x[8]]) ;
D2 := Divisor(C,I2);

D := D1 + D2 ;
Cl, phi, psi := ClassGroup(C) ;

K := CanonicalDivisor(C) ;
T := 7*D -3*K;

print "the order of the shadow mod 5 is:";
Order(psi(T));


// mod 7 
Zx<[x]> := PolynomialRing(GF(7), #basis);          
Q1 := [ Zx ! f : f in Q];
P := ProjectiveSpace(Zx);
C := Curve(P, Q) ;

I := Ideal([x[1], x[2], x[3]]) ;
D1 := Divisor(C,I) ;

I2 := Ideal([x[4], x[5], x[6], x[7], x[8]]) ;
D2 := Divisor(C,I2);

D := D1 + D2 ;
Cl, phi, psi := ClassGroup(C) ;

K := CanonicalDivisor(C) ;
T := 7*D -3*K;

print "the order of the shadow mod 7 is:";
Order(psi(T));

