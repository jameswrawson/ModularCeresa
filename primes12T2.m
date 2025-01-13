load "divisor_checker.m"

P1s := PrimesInInterval(193, 150889);
P1s := [p : p in P1s | p mod 24 eq 1];
P1s := [p : p in P1s | LegendreSymbol(p, 7) eq 1];
P1s := [p : p in P1s | LegendreSymbol(p, 11) eq -1];

P2s := PrimesInInterval(193, 244897);
P2s := [p : p in P2s | p mod 24 eq 1];
P2s := [p : p in P2s | LegendreSymbol(p, 7) eq 1];
P2s := [p : p in P2s | LegendreSymbol(p, 11) eq 1];

Ps := P1s cat P2s;

for p in Ps do
	print(p);
	quad, tr := checkPrime(p);
	if not quad then
		Append(~no_quad, p);
	end if;

	if not tr then
		Append(~triv, p);
	end if;
end for;
