function test_qadic_constructors()
   print("qadic.constructors...")

   R = QadicField(7, 1, 30)

   @test elem_type(R) == qadic
   @test elem_type(FlintQadicField) == qadic
   @test parent_type(qadic) == FlintQadicField

   @test isa(R, FlintQadicField)

   S = QadicField(fmpz(1009), 1, 30)

   @test isa(S, FlintQadicField)

   @test isa(R(), qadic)

   @test isa(R(1), qadic)

   @test isa(R(ZZ(123)), qadic)

   @test isa(R(ZZ(1)//7^2), qadic)

   @test isa(1 + 2*7 + 4*7^2 + O(R, 7^3), qadic)

   @test isa(13 + 357*fmpz(1009) + O(S, fmpz(1009)^12), qadic)

   @test isa(fmpz(1)//7^2 + fmpz(2)//7 + 3 + 4*7 + O(R, 7^2), qadic)

   s = R()

   t = deepcopy(s)

   @test isa(t, qadic)

   @test parent(t) === R

   println("PASS")
end

function test_qadic_printing()
   print("qadic.printing...")

   R = QadicField(7, 1, 30)

   a = 1 + 2*7 + 4*7^2 + O(R, 7^3)

   b = string(a)

   @test b isa String

   println("PASS")
end

function test_qadic_manipulation()
   print("qadic.manipulation...")

   R = QadicField(7, 1, 30)

   a = 1 + 2*7 + 4*7^2 + O(R, 7^3)
   b = 7^2 + 3*7^3 + O(R, 7^5)
   c = R(2)

   @test isone(one(R))

   @test iszero(zero(R))

   @test precision(a) == 3

   @test prime(R) == 7

   @test valuation(b) == 2

   println("PASS")
end

function test_qadic_unary_ops()
   print("qadic.unary_ops...")

   R = QadicField(7, 1, 30)

   a = 1 + 2*7 + 4*7^2 + O(R, 7^3)
   b = R(0)

   @test -a == 6 + 4*7^1 + 2*7^2 + O(R, 7^3)

   @test iszero(-b)

   println("PASS")
end

function test_qadic_binary_ops()
   print("qadic.binary_ops...")

   R = QadicField(7, 1, 30)

   a = 1 + 2*7 + 4*7^2 + O(R, 7^3)
   b = 7^2 + 3*7^3 + O(R, 7^5)
   c = O(R, 7^3)
   d = R(2)

   @test a + b == 1 + 2*7^1 + 5*7^2 + O(R, 7^3)

   @test a - b == 1 + 2*7^1 + 3*7^2 + O(R, 7^3)

   @test a*b == 1*7^2 + 5*7^3 + 3*7^4 + O(R, 7^5)

   @test b*c == O(R, 7^5)

   @test a*d == 2 + 4*7^1 + 1*7^2 + O(R, 7^3)

   println("PASS")
end

function test_qadic_adhoc_binary()
   print("qadic.adhoc_binary...")

   R = QadicField(7, 1, 30)

   a = 1 + 2*7 + 4*7^2 + O(R, 7^3)
   b = 7^2 + 3*7^3 + O(R, 7^5)
   c = O(R, 7^3)
   d = R(2)

   @test a + 2 == 3 + 2*7^1 + 4*7^2 + O(R, 7^3)

   @test 3 - b == 3 + 6*7^2 + 3*7^3 + 6*7^4 + O(R, 7^5)

   @test a*fmpz(5) == 5 + 3*7^1 + O(R, 7^3)

   @test fmpz(3)*c == O(R, 7^3)

   @test 2*d == 4

   @test 2 + d == 4

   @test iszero(d - fmpz(2))

   @test a + fmpz(1)//7^2 == fmpz(1)//7^2 + 1 + 2*7^1 + 4*7^2 + O(R, 7^3)

   @test (fmpz(12)//11)*b == 3*7^2 + 3*7^3 + O(R, 7^5)

   @test c*(fmpz(1)//7) == O(R, 7^2)

   println("PASS")
end

function test_qadic_comparison()
   print("qadic.comparison...")

   R = QadicField(7, 1, 30)

   a = 1 + 2*7 + 4*7^2 + O(R, 7^3)
   b = 3*7^3 + O(R, 7^5)
   c = O(R, 7^3)
   d = R(2)

   @test a == 1 + 2*7 + O(R, 7^2)

   @test b == c

   @test c == R(0)

   @test d == R(2)

   println("PASS")
end

function test_qadic_adhoc_comparison()
   print("qadic.adhoc_comparison...")

   R = QadicField(7, 1, 30)

   a = 1 + O(R, 7^3)
   b = O(R, 7^5)
   c = R(2)

   @test a == 1

   @test b == ZZ(0)

   @test c == 2

   @test fmpz(2) == c

   @test a == fmpz(344)//1

   println("PASS")
end

function test_qadic_powering()
   print("qadic.powering...")

   R = QadicField(7, 1, 30)

   a = 1 + 7 + 2*7^2 + O(R, 7^3)
   b = O(R, 7^5)
   c = R(2)

   @test a^5 == 1 + 5*7^1 + 6*7^2 + O(R, 7^3)

   @test b^3 == O(R, 7^5)

   @test c^7 == 2 + 4*7^1 + 2*7^2

   println("PASS")
end

function test_qadic_inversion()
   print("qadic.inversion...")

   R = QadicField(7, 1, 30)

   a = 1 + 7 + 2*7^2 + O(R, 7^3)
   b = 2 + 3*7 + O(R, 7^5)
   c = 7^2 + 2*7^3 + O(R, 7^4)
   d = 7 + 2*7^2 + O(R, 7^5)

   @test inv(a) == 1 + 6*7^1 + 5*7^2 + O(R, 7^3)

   @test inv(b) == 4 + 4*7^1 + 3*7^2 + 1*7^3 + 1*7^4 + O(R, 7^5)

   @test inv(c) == fmpz(1)//7^2 + fmpz(5)//7 + O(R, 7^0)

   @test inv(d) == fmpz(1)//7 + 5 + 3*7^1 + 6*7^2 + O(R, 7^3)

   @test inv(R(1)) == 1

   println("PASS")
end

function test_qadic_exact_division()
   print("qadic.exact_division...")

   R = QadicField(7, 1, 30)

   a = 1 + 7 + 2*7^2 + O(R, 7^3)
   b = 2 + 3*7 + O(R, 7^5)
   c = 7^2 + 2*7^3 + O(R, 7^4)
   d = 7 + 2*7^2 + O(R, 7^5)

   @test divexact(a, b) == 4 + 1*7^1 + 2*7^2 + O(R, 7^3)

   @test divexact(c, d) == 1*7^1 + O(R, 7^3)

   @test divexact(d, R(7^3)) == fmpz(1)//7^2 + fmpz(2)//7 + O(R, 7^2)

   @test divexact(R(34), R(17)) == 2

   println("PASS")
end

function test_qadic_adhoc_exact_division()
   print("qadic.adhoc_exact_division...")

   R = QadicField(7, 1, 30)

   a = 1 + 7 + 2*7^2 + O(R, 7^3)
   b = 2 + 3*7 + O(R, 7^5)
   c = 7^2 + 2*7^3 + O(R, 7^4)
   d = 7 + 2*7^2 + O(R, 7^5)

   @test divexact(a, 2) == 4 + 1*7^2 + O(R, 7^3)

   @test divexact(b, fmpz(7)) == fmpz(2)//7 + 3 + O(R, 7^4)

   @test divexact(c, fmpz(12)//7^2) == 3*7^4 + 5*7^5 + O(R, 7^6)

   @test divexact(2, d) == fmpz(2)//7 + 3 + 6*7^2 + O(R, 7^3)

   @test divexact(R(3), 3) == 1

   @test divexact(fmpz(5)//7, R(5)) == fmpz(1)//7

   println("PASS")
end

function test_qadic_divides()
   print("qadic.divides...")

   R = QadicField(7, 1, 30)

   a = 1 + 7 + 2*7^2 + O(R, 7^3)
   b = 2 + 3*7 + O(R, 7^5)

   flag, q = divides(a, b)

   @test flag
   @test q == divexact(a, b)

   println("PASS")
end

function test_qadic_gcd()
   print("qadic.adhoc_gcd...")

   R = QadicField(7, 1, 30)

   a = 1 + 7 + 2*7^2 + O(R, 7^3)
   b = 2 + 3*7 + O(R, 7^5)

   @test gcd(a, b) == 1

   @test gcd(zero(R), zero(R)) == 0

   println("PASS")
end

function test_qadic_square_root()
   print("qadic.square_root...")

   R = QadicField(7, 1, 30)

   a = 1 + 7 + 2*7^2 + O(R, 7^3)
   b = 2 + 3*7 + O(R, 7^5)
   c = 7^2 + 2*7^3 + O(R, 7^4)

   @test sqrt(a)^2 == a

   @test sqrt(b)^2 == b

   @test sqrt(c)^2 == c

   @test sqrt(R(121))^2 == R(121)

   println("PASS")
end

function test_qadic_special_functions()
   print("qadic.special_functions...")

   R = QadicField(7, 1, 30)

   a = 1 + 7 + 2*7^2 + O(R, 7^3)
   b = 2 + 5*7 + 3*7^2 + O(R, 7^3)
   c = 3*7 + 2*7^2 + O(R, 7^5)

   @test exp(c) == 1 + 3*7^1 + 3*7^2 + 4*7^3 + 4*7^4 + O(R, 7^5)

   @test log(a) == 1*7^1 + 5*7^2 + O(R, 7^3)

   @test exp(R(0)) == 1

   @test log(R(1)) == 0

   @test teichmuller(b) == 2 + 4*7^1 + 6*7^2 + O(R, 7^3)

   println("PASS")
end

function test_qadic()
   test_qadic_constructors()
   test_qadic_printing()
   test_qadic_manipulation()
   test_qadic_unary_ops()
   test_qadic_binary_ops()
   test_qadic_adhoc_binary()
   test_qadic_comparison()
   test_qadic_adhoc_comparison()
   test_qadic_powering()
   test_qadic_inversion()
   test_qadic_exact_division()
   test_qadic_adhoc_exact_division()
   test_qadic_divides()
   test_qadic_gcd()
   test_qadic_square_root()
   test_qadic_special_functions()

   println("")
end
