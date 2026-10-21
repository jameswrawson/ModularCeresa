// verify X0(121)

load "modular_models.m";
// we use our basis to compute a canonical model of the curve 
pols, fix, neg := QuadricModel(121, 121);
R := Parent(pols[1]);

// base change to GF(5)
F5x<[x]> := PolynomialRing(GF(5), Rank(R));
C := Curve(ProjectiveSpace(F5x), [F5x ! pol : pol in pols]);

//We compute the divisors explicitly. To be invariant under Atkin--Lehner, either the coordinates which are fixed by Atkin--Lehner are all zero, or the coordinates which are anti-invariant are all zero.

I := Ideal([F5x ! f : f in fix]);
D1 := Divisor(C,I);

I2 := Ideal([F5x ! f : f in neg]);
D2 := Divisor(C,I2);

F := D1 + D2; 

Cl, phi, psi := ClassGroup(C);

K := CanonicalDivisor(C);

scle := GCD(Degree(F), Degree(K));
d1 := Integers() ! (Degree(K) / scle);
d2 := Integers() ! (Degree(F) / scle);
T := d1*F - d2*K; //Form the (scaled) shadow divisor from the fixed point divisor and the canonical
print "The order of the shadow mod 5 is:";
Order(psi(T));

// base change to GF(13)
Zx<[x]> := PolynomialRing(GF(13), #basis);
Q1 := [ Zx ! f : f in Q];
P := ProjectiveSpace(Zx);
C := Curve(P, Q1) ;


// we compute  the divisors explicitly 

I := Ideal([x[1], x[2]]) ;
D1 := Divisor(C,I) ;

I2 := Ideal([x[3], x[4], x[5], x[6]]) ;
D2 := Divisor(C,I2);

D := D1 + D2 ;

Cl, phi, psi := ClassGroup(C);

K := CanonicalDivisor(C) ;

T := 5*D -3*K;
print "The order of the shadow mod 13 is:";
Order(psi(T));

