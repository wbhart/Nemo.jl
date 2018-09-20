function test_gfp_fmpz_poly_constructors()
   print("gfp_fmpz_poly.constructors...")
 
   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   @test elem_type(S) == gfp_fmpz_poly
   @test elem_type(GFPFmpzPolyRing) == gfp_fmpz_poly
   @test parent_type(gfp_fmpz_poly) == GFPFmpzPolyRing

   @test typeof(S) <: GFPFmpzPolyRing

   @test isa(x, PolyElem)

   f = x^3 + 2x^2 + x + 1

   @test isa(f, PolyElem)

   g = S(2)

   @test isa(g, PolyElem)

   h = S(x^2 + 2x + 1)

   @test isa(h, PolyElem)

   k = S([R(1), R(0), R(3)])

   @test isa(k, PolyElem)

   l = S()

   @test isa(l, PolyElem)

   m = S(fmpz(123))

   @test isa(m, PolyElem)

   n = S([fmpz(1), fmpz(0), fmpz(3)])

   @test isa(n, PolyElem)

   T, y = PolynomialRing(ZZ, "y")

   p = 3y^3 + 2y - 1
   q = S(p)

   @test isa(q, PolyElem)

   r = S([1, 2, 3])

   @test isa(r, PolyElem)

   println("PASS")
end

function test_gfp_fmpz_poly_printing()
   print("gfp_fmpz_poly.printing...")
 
   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")
   f = x^3 + 2x^2 + x + 1

   @test string(f) == "x^3+2*x^2+x+1"

   println("PASS")
end

function test_gfp_fmpz_poly_manipulation()
   print("gfp_fmpz_poly.manipulation...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   @test iszero(zero(S))
   
   @test isone(one(S))

   @test isgen(gen(S))
   
   @test isunit(one(S))

   f = x^2 + 2x + 1

   @test lead(f) == 1

   @test degree(f) == 2

   @test length(f) == 3

   @test coeff(f, 1) == 2

   @test_throws DomainError coeff(f, -1)

   @test canonical_unit(-x + 1) == -1

   @test deepcopy(f) == f

   setcoeff!(f, 1, UInt(2))

   @test coeff(f, 1) == 2

   setcoeff!(f, 1, 3)

   @test coeff(f, 1) == 3

   setcoeff!(f, 1, fmpz(2)^100)

   @test coeff(f, 1) == 32146634986640907030

   @test modulus(x) == 123456789012345678949

   @test modulus(R) == 123456789012345678949

   println("PASS")
end

function test_gfp_fmpz_poly_binary_ops()
   print("gfp_fmpz_poly.binary_ops...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 2

   @test f + g == x^3+x^2+5*x+3

   @test f*g == x^5+2*x^4+4*x^3+8*x^2+7*x+2

   @test f - g == 123456789012345678948*x^3+x^2+123456789012345678948*x+123456789012345678948

   println("PASS")
end

function test_gfp_fmpz_poly_adhoc_binary()
   print("gfp_fmpz_poly.adhoc_binary...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 2

   @test f*12 == 12*x^2+24*x+12

   @test 7*g == 7*x^3+21*x+14

   @test fmpz(3)*g == 3*x^3+9*x+6

   @test f*fmpz(2) == 2*x^2+4*x+2

   @test f + 12 == x^2+2*x+13

   @test f + fmpz(12) == x^2+2*x+13

   @test f - 12 == x^2+2*x+123456789012345678938

   @test f - fmpz(12) == x^2+2*x+123456789012345678938

   @test 12 + g == x^3+3*x+14

   @test fmpz(12) + g == x^3+3*x+14

   @test 12 - g == 123456789012345678948*x^3+123456789012345678946*x+10

   @test fmpz(12) - g == 123456789012345678948*x^3+123456789012345678946*x+10

   @test f + R(12) == x^2+2*x+13

   @test R(12) + g == x^3+3*x+14

   @test f - R(12) == x^2+2*x+123456789012345678938

   @test R(12) - g == 123456789012345678948*x^3+123456789012345678946*x+10

   @test R(7)*g == 7*x^3+21*x+14

   @test f*R(12) == 12*x^2+24*x+12

   println("PASS")
end

function test_gfp_fmpz_poly_comparison()
   print("gfp_fmpz_poly.comparison...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 2

   @test f != g

   @test isequal(f, deepcopy(f))

   println("PASS")
end

function test_gfp_fmpz_poly_adhoc_comparison()
   print("gfp_fmpz_poly.adhoc_comparison...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1

   @test f != 1 

   @test 1 != f 

   @test S(7) == fmpz(7)

   @test fmpz(7) != f

   @test S(7) == R(7)

   @test R(7) != x + 1

   println("PASS")
end

function test_gfp_fmpz_poly_unary_ops()
   print("gfp_fmpz_poly.unary_ops...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1

   @test -f == 123456789012345678948*x^2+123456789012345678947*x+123456789012345678948

   println("PASS")
end

function test_gfp_fmpz_poly_truncation()
   print("gfp_fmpz_poly.truncation...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 1

   @test truncate(f, 2) == 2*x+1

   @test_throws DomainError truncate(f, -1)

   @test mullow(f, g, 3) == 7*x^2+5*x+1

   @test_throws DomainError mullow(f, g, -1)

   println("PASS")
end

function test_gfp_fmpz_poly_reverse()
   print("gfp_fmpz_poly.reverse...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 3

   @test reverse(f) == 3x^2 + 2x + 1

   println("PASS")
end

function test_gfp_fmpz_poly_shift()
   print("gfp_fmpz_poly.shift...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1

   @test shift_left(f, 3) == x^5 + 2x^4 + x^3

   @test_throws DomainError shift_left(f, -1)

   @test shift_right(f, 1) == x + 2

   @test_throws DomainError shift_right(f, -1)

   println("PASS")
end

function test_gfp_fmpz_poly_powering()
   print("gfp_fmpz_poly.powering...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1

   @test f^6 == x^12+12*x^11+66*x^10+220*x^9+495*x^8+792*x^7+924*x^6+792*x^5+495*x^4+220*x^3+66*x^2+12*x+1

   @test_throws DomainError f^-1

   println("PASS")
end

function test_gfp_fmpz_poly_exact_division()
   print("gfp_fmpz_poly.exact_division...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 1

   @test divexact(f*g, f) == g

   println("PASS")
end

function test_gfp_fmpz_poly_adhoc_exact_division()
   print("gfp_fmpz_poly.adhoc_exact_division...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1
   
   @test divexact(3*f, fmpz(3)) == f

   @test divexact(3*f, 3) == f

   @test divexact(R(7)*f, R(7)) == f
   
   println("PASS")
end

function test_gfp_fmpz_poly_modular_arithmetic()
   print("gfp_fmpz_poly.modular_arithmetic...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = 3*x^2 + x + 2
   g = 5*x^2 + 2*x + 1
   h = 3*x^3 + 2*x^2 + x + 7

   @test invmod(f, h) == 112883663504991137175*x^2+86761824016232857498*x+48511987621979662257

   @test mulmod(f, g, h) == 82304526008230452642*x^2+41152263004115226286*x+41152263004115226316

   @test powmod(f, 10, h) == 118470346535924950143*x^2+97790722831392543222*x+115967716915690326718

   @test powmod(f, fmpz(10), h) == 118470346535924950143*x^2+97790722831392543222*x+115967716915690326718

   @test powmod(f, -10, g) == 78305338116088931412*x+91239060941924718463

   @test powmod(f, -fmpz(10), g) == 78305338116088931412*x+91239060941924718463

   println("PASS")
end

function test_gfp_fmpz_poly_euclidean_division()
   print("gfp_fmpz_poly.euclidean_division...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 1

   @test mod(g, f) == 6*x+3

   @test divrem(g, f) == (x+123456789012345678947, 6*x+3)
 
   println("PASS")
end

function test_gfp_fmpz_poly_gcd()
   print("gfp_fmpz_poly.gcd...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 1
   h = x^2 + 1

   @test gcd(f*h, g*h) == x^2+1
 
   @test gcdx(f*h, g*h) == (x^2+1, 41152263004115226317*x^2+41152263004115226316*x+2,82304526008230452632*x+123456789012345678948)
   println("PASS")
end

function test_gfp_fmpz_poly_gcdinv()
   print("gfp_fmpz_poly.gcdinv...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 1

   @test gcdinv(f, g) == (1, 41152263004115226317*x^2+41152263004115226316*x+2)
 
   println("PASS")
end

function test_gfp_fmpz_poly_evaluation()
   print("gfp_fmpz_poly.evaluation...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1

   @test evaluate(f, 3) == 16
   
   @test evaluate(f, fmpz(10)) == 121

   @test evaluate(f, R(10)) == 121

if VERSION >= v"0.5.0-dev+3171"

   @test f(3) == 16

   @test f(fmpz(10)) == 121

   @test f(R(10)) == 121

end

   println("PASS")
end

function test_gfp_fmpz_poly_composition()
   print("gfp_fmpz_poly.composition...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 1

   @test compose(f, g) == x^6+6*x^4+4*x^3+9*x^2+12*x+4

   println("PASS")
end

function test_gfp_fmpz_poly_derivative()
   print("gfp_fmpz_poly.derivative...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1

   @test derivative(f) == 2x + 2

   println("PASS")
end

function test_gfp_fmpz_poly_integral()
   print("gfp_fmpz_poly.integral...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1
   
   @test integral(f) == 82304526008230452633*x^3+x^2+x

   println("PASS")
end

function test_gfp_fmpz_poly_resultant()
   print("gfp_fmpz_poly.resultant...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = 5x^2 + 2x + 1
   g = 2x^3 + 3x + 1

   @test resultant(f, g) == 212

   println("PASS")
end

function test_gfp_fmpz_poly_discriminant()
   print("gfp_fmpz_poly.discriminant...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1

   @test discriminant(f) == 0

   println("PASS")
end

function test_gfp_fmpz_poly_lift()
   print("gfp_fmpz_poly.lift...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   T, y = PolynomialRing(ZZ, "y")

   f = x^2 + 2x + 1

   @test lift(T, f) == y^2 + 2y + 1

   println("PASS")
end

function test_gfp_fmpz_poly_isirreducible()
   print("gfp_fmpz_poly.isirreducible...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1

   @test isirreducible(f) == false

   println("PASS")
end

function test_gfp_fmpz_poly_issquarefree()
   print("gfp_fmpz_poly.issquarefree...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1

   @test issquarefree(f) == false

   println("PASS")
end

function test_gfp_fmpz_poly_factor()
   print("gfp_fmpz_poly.factor...")

   R = GF(fmpz(123456789012345678949))
   S, x = PolynomialRing(R, "x")

   f = 3*(x^2 + 2x + 1)
   g = x^3 + 3x + 1

   R = factor(f*g)

   @test f*g == unit(R) * prod([ p^e for (p, e) in R])

   R = factor_squarefree(f*g)

   @test f*g == unit(R) * prod([ p^e for (p, e) in R])

   R = factor_distinct_deg((x + 1)*g*(x^5+x+1))

   @test length(R) == 2
   @test R == Dict(1 => x^3+2*x^2+2*x+1,
                3 => x^6+123456789012345678948*x^5+3*x^4+123456789012345678948*x^3+123456789012345678948*x^2+3*x+1)

   R = factor_shape(f*g)
   
   @test length(R) == 2
   @test R == Dict(3=>1, 1=>2)

   println("PASS")
end

function test_gfp_fmpz_poly_remove_valuation()
   print("gfp_fmpz_poly.remove_valuation...")

   R = GF(fmpz(123456789012345678949))
   S, y = PolynomialRing(R, "y")

   f = 7y^2 + 3y + 2
   g = f^5*(11y^3 - 2y^2 + 5)

   v, h = remove(g, f)

   @test valuation(g, f) == 5
   @test v == 5
   @test h == (11y^3 - 2y^2 + 5)

   v, q = divides(f*g, f)

   @test v
   @test q == g

   v, q = divides(f*g + 1, f)

   @test !v

   println("PASS")
end



function test_gfp_fmpz_poly()
   test_gfp_fmpz_poly_constructors()
   test_gfp_fmpz_poly_printing()
   test_gfp_fmpz_poly_manipulation()
   test_gfp_fmpz_poly_binary_ops()
   test_gfp_fmpz_poly_adhoc_binary()
   test_gfp_fmpz_poly_comparison()
   test_gfp_fmpz_poly_adhoc_comparison()
   test_gfp_fmpz_poly_unary_ops()
   test_gfp_fmpz_poly_truncation()
   test_gfp_fmpz_poly_reverse()
   test_gfp_fmpz_poly_shift()
   test_gfp_fmpz_poly_powering()
   test_gfp_fmpz_poly_exact_division()
   test_gfp_fmpz_poly_adhoc_exact_division()
   test_gfp_fmpz_poly_modular_arithmetic()
   test_gfp_fmpz_poly_euclidean_division()
   test_gfp_fmpz_poly_gcdinv()
   test_gfp_fmpz_poly_evaluation()
   test_gfp_fmpz_poly_composition()
   test_gfp_fmpz_poly_derivative()
   test_gfp_fmpz_poly_integral()
   test_gfp_fmpz_poly_resultant()
   test_gfp_fmpz_poly_discriminant()
   test_gfp_fmpz_poly_lift()
   test_gfp_fmpz_poly_isirreducible()
   test_gfp_fmpz_poly_issquarefree()
   test_gfp_fmpz_poly_factor()
   test_gfp_fmpz_poly_remove_valuation()

   println("")
end

