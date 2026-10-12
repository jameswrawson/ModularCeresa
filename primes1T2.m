load "divisor_checker.m"

//Compute the primes for which the methods of Section 5.3 do not succeed, which depends on whether p is a quadratic residue mod 11 or not
P1s := PrimesInInterval(168, 150889);
P1s := [p : p in P1s | p mod 24 eq 1];
P1s := [p : p in P1s | LegendreSymbol(p, 7) eq 1];
P1s := [p : p in P1s | LegendreSymbol(p, 11) eq -1];

P2s := PrimesInInterval(168, 244897);
P2s := [p : p in P2s | p mod 24 eq 1];
P2s := [p : p in P2s | LegendreSymbol(p, 7) eq 1];
P2s := [p : p in P2s | LegendreSymbol(p, 11) eq 1];

Ps := P1s cat P2s;

no_quad := [];
triv := []; 

//For each prime, check the non-vanishing of the shadow associated to T2. For failures, we separate these by failure mode: insufficient supersingular j-invariants or the shadow is zero
for p in Ps do
	quad, tr := checkPrime(p);
	if not quad then
		Append(~no_quad, p);
	end if;

	if not tr then
		Append(~triv, p);
	end if;
end for;
