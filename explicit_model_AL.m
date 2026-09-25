load "modular_models.m";

//Show the Atkin--Lehner shadow for X0(74) is non-zero
N := 74;
eqs, fix, neg := QuadricModel(N, N);
print("Order of the Atkin--Lehner shadow of X0(74) modulo 3 is");
m3 := ALShadowOrder(eqs, fix, neg, 3);
print(m3);
print("Order of the Atkin--Lehner shadow of X0(74) modulo 5 is");
m5 := ALShadowOrder(eqs, fix, neg, 5);
print(m5);
assert(m3 ne m5);
print("--------------");

//Show the Atkin--Lehner shadow for X0(121) is non-zero
N := 121;
eqs, fix, neg := QuadricModel(N, N);
print("Order of the Atkin--Lehner shadow of X0(121) modulo 3 is");
m3 := ALShadowOrder(eqs, fix, neg, 3);
print(m3);
print("Order of the Atkin--Lehner shadow of X0(121) modulo 5 is");
m5 := ALShadowOrder(eqs, fix, neg, 5);
print(m5);
assert(m3 ne m5);
print("--------------");

//Show the Atkin--Lehner shadow for X0(125) is non-zero
N := 125;
eqs, fix, neg := QuadricModel(N, N);
print("Order of the Atkin--Lehner shadow of X0(125) modulo 3 is");
m3 := ALShadowOrder(eqs, fix, neg, 3);
print(m3);
print("Order of the Atkin--Lehner shadow of X0(125) modulo 7 is");
m7 := ALShadowOrder(eqs, fix, neg, 7);
print(m7);
assert(m3 ne m7);
print("--------------");

//Show the Atkin--Lehner shadow for X0(169) is non-zero
N := 169;
eqs, fix, neg := QuadricModel(N, N);
print("Order of the Atkin--Lehner shadow of X0(169) modulo 3 is");
m3 := ALShadowOrder(eqs, fix, neg, 3);
print(m3);
print("Order of the Atkin--Lehner shadow of X0(169) modulo 5 is");
m5 := ALShadowOrder(eqs, fix, neg, 5);
print(m5);
assert(m3 ne m5);
print("--------------");
