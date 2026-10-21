//Given an elliptic curve, EC, and points in the upper half plane, pts, compute whether the image of these points is trivial
function checkAnalPoint(EC, pts, mults : tol := 1e-5);
	Ps := Periods(EC);
	lat_pts := ModularParametrisation(EC, pts);
	tot := 0;
	for i := 1 to #pts do
		tot := tot + 2 * mults[i] * lat_pts[i]; //multiple of 2 from CM points coming in pairs
	end for;

	n2 := Im(tot) / Im(Ps[2]);
	if Abs(n2 - Round(n2)) gt tol then
		return true;
	end if;
	
	tot := tot - Round(n2) * Ps[2];
	n1 := RealField() ! (Re(tot) / Ps[1]);

	if Abs(n1 - Round(n1)) gt tol then
		return true;
	else
		return false;
	end if;
end function;
