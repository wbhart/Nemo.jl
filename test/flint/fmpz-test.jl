@testset "fmpz.abstract_types..." begin
   @test fmpz <: RingElem

   @test FlintIntegerRing <: Nemo.Ring

   @test elem_type(FlintIntegerRing()) == fmpz
   @test elem_type(FlintIntegerRing) == fmpz
   @test parent_type(fmpz) == FlintIntegerRing
end

@testset "fmpz.constructors..." begin
   a = fmpz(-123)
   @test isa(a, RingElem)

   b = fmpz(12.0)
   @test isa(b, RingElem)

   c = fmpz("-1234567876545678376545678900000000000000000000000000")
   @test isa(c, RingElem)

   d = fmpz(c)
   @test isa(d, RingElem)

   e = deepcopy(c)
   @test isa(e, RingElem)

   f = fmpz(BigFloat(10)^100)
   @test isa(f, RingElem)

   g = fmpz()
   @test isa(f, RingElem)
end

@testset "fmpz.rand..." begin
   @test rand(FlintZZ, 1:9) isa fmpz
   Random.seed!(rng, 0)
   t = rand(rng, FlintZZ, 1:9)
   @test t isa fmpz
   Random.seed!(rng, 0)
   @test t == rand(rng, FlintZZ, 1:9)
end

@testset "fmpz.printing..." begin
   a = fmpz(-123)

   @test string(a) == "-123"
end

@testset "fmpz.convert..." begin
   a = fmpz(-123)
   b = fmpz(12)

   @test Int(a) == -123
   @test UInt(b) == UInt(12)
   @test BigInt(a) == BigInt(-123)
   @test Float64(a) == Float64(-123)
   @test Float32(a) == Float32(-123)
   @test Float16(a) == Float16(-123)
   @test BigFloat(a) == BigFloat(-123)

   @test_throws InexactError Int(fmpz(1234484735687346876324432764872))
   @test_throws InexactError UInt(fmpz(typemin(Int)))
end

@testset "fmpz.manipulation..." begin
   a = one(FlintIntegerRing())
   b = zero(FlintIntegerRing())
   c = zero(fmpz)

   @test isa(a, RingElem)

   @test isa(b, RingElem)

   @test isa(c, RingElem)

   @test sign(a) == 1

   @test fits(Int, a)

   @test fits(UInt, a)

   @test size(a) == 1

   @test canonical_unit(fmpz(-12)) == -1

   @test isunit(fmpz(-1))

   @test iszero(b)

   @test isone(a)

   @test numerator(fmpz(12)) == fmpz(12)

   @test denominator(fmpz(12)) == fmpz(1)

   @test floor(fmpz(12)) == fmpz(12)

   @test ceil(fmpz(12)) == fmpz(12)

   @test iseven(fmpz(12))
   @test isodd(fmpz(13))
   b = big(2)
   x = rand(-b^rand(1:1000):b^rand(1:1000))
   y = fmpz(x)
   @test iseven(x) == iseven(y)
   @test isodd(x) == isodd(y)

   @test characteristic(ZZ) == 0
end

@testset "fmpz.binary_ops..." begin
   a = fmpz(12)
   b = fmpz(26)

   @test a + b == 38

   @test a - b == -14

   @test a*b == 312

   @test b%a == 2

   @test b&a == 8

   @test b|a == 30

   @test xor(b, a) == 22
end

@testset "fmpz.division..." begin
   a = fmpz(12)
   b = fmpz(26)

   @test fdiv(b, a) == 2

   @test cdiv(b, a) == 3

   @test tdiv(b, a) == 2

   @test div(b, a) == 2
end

@testset "fmpz.remainder..." begin
   a = fmpz(12)
   b = fmpz(26)

   @test mod(b, a) == 2

   @test mod(fmpz(3), fmpz(-2)) == fmpz(-1)

   @test rem(b, a) == 2

   @test mod(b, 12) == 2

   @test rem(b, 12) == 2
end

@testset "fmpz.exact_division..." begin
   @test divexact(fmpz(24), fmpz(12)) == 2
   @test_throws ArgumentError divexact(fmpz(24), fmpz(11))
end

@testset "fmpz.divides..." begin
   flag, q = divides(fmpz(12), fmpz(0))
   @test flag == false
   @test divides(fmpz(12), fmpz(6)) == (true, fmpz(2))
   @test divides(fmpz(0), fmpz(0)) == (true, fmpz(0))
end

@testset "fmpz.gcd_lcm..." begin
   a = fmpz(12)
   b = fmpz(26)

   @test gcd(a, b) == 2
   @test gcd(a, 26) == 2
   @test gcd(12, b) == 2

   @test_throws ErrorException gcd(fmpz[])

   @test gcd(fmpz[8]) == 8

   @test gcd([fmpz(10), fmpz(2)]) == 2

   @test gcd([fmpz(1), fmpz(2), fmpz(3)]) == 1

   @test gcd([fmpz(9), fmpz(27), fmpz(3)]) == 3

   @test lcm(a, b) == 156
   @test lcm(12, b) == 156
   @test lcm(a, 26) == 156

   @test_throws ErrorException lcm(fmpz[])

   @test lcm(fmpz[2]) == 2

   @test lcm(fmpz[2, 3]) == 6

   @test lcm(fmpz[2, 2, 2]) == 2

   @test lcm(fmpz[2, 3, 2]) == 6
end

@testset "fmpz.logarithm..." begin
   a = fmpz(12)
   b = fmpz(26)

   @test flog(b, a) == 1

   @test_throws DomainError flog(b, -a)

   @test flog(b, 12) == 1

   @test_throws DomainError flog(b, -12)

   @test clog(b, a) == 2

   @test_throws DomainError clog(b, -a)

   @test clog(b, 12) == 2

   @test_throws DomainError clog(b, -12)
end

@testset "fmpz.adhoc_binary..." begin
   a = fmpz(-12)

   @test 3 + a == -9

   @test a + 3 == -9

   @test a - 3 == -15

   @test 5 - a == 17

   @test a*5 == -60

   @test 5*a == -60

   @test a%5 == -2
end

@testset "fmpz.adhoc_division..." begin
   a = fmpz(-12)

   @test fdiv(a, 5) == -3

   @test tdiv(a, 7) == -1

   @test cdiv(a, 7) == -1

   @test div(a, 3) == -4

   @test div(-12, fmpz(3)) == -4

   @test mod(-12, fmpz(3)) == 0

   @test isa(mod(fmpz(2), -3), fmpz)

   @test mod(fmpz(2), -3) == -1

   @test rem(-12, fmpz(3)) == 0

   @test_throws ArgumentError divexact(ZZ(2), 3)
end

@testset "fmpz.shift.." begin
   a = fmpz(-12)

   @test a >> 3 == -2

   @test fdivpow2(a, 2) == -3

   @test_throws DomainError fdivpow2(a, -1)

   @test cdivpow2(a, 2) == -3

   @test_throws DomainError cdivpow2(a, -1)

   @test tdivpow2(a, 2) == -3

   @test_throws DomainError tdivpow2(a, -1)

   @test a << 4 == -192
end

@testset "fmpz.powering..." begin
   a = fmpz(-12)

   @test a^5 == -248832

   @test a^1 == a
   @test a^1 !== a

   @test isone(a^0)

   for a in fmpz.(-5:5)
      for e = -5:-1
         if a != 1 && a != -1
            @test_throws DomainError a^e
         end
      end
      @test a^1 == a
      @test a^1 !== a
   end

   a = fmpz(1)
   for e = -2:2
      @test isone(a^e)
      @test a^e !== a
   end

   a = fmpz(-1)
   for e = [-3, -1, 1, 3, 5]
      @test a^e == a
      @test a^e !== a
   end
   for e = [-2, 0, 2, 4]
      @test isone(a^e )
   end
end

@testset "fmpz.comparison..." begin
   a = fmpz(-12)
   b = fmpz(5)

   @test a < b

   @test b > a

   @test b >= a

   @test a <= b

   @test a == fmpz(-12)

   @test a != b

   @test isequal(a, fmpz(-12))

   @test cmpabs(a, b) == 1

   @test cmp(a, b) == -1

   @test fmpz(2) < 47632748687326487326487326487326

   @test fmpz(2) < 476327486873264873264873264873264837624982
end

@testset "fmpz.adhoc_comparison..." begin
   a = fmpz(-12)

   @test a < 7

   @test a > -40

   @test 7 > a

   @test -40 < a

   @test a <= 7

   @test a >= -40

   @test 7 >= a

   @test -40 <= a

   @test a == -12

   @test a != 4

   @test -12 == a

   @test 4 != a

   a = fmpz(2)

   @test a < UInt(7)

   @test a > UInt(1)

   @test UInt(7) > a

   @test UInt(1) < a

   @test a <= UInt(7)

   @test a >= UInt(2)

   @test UInt(7) >= a

   @test UInt(1) <= a

   @test a == UInt(2)

   @test a != UInt(4)

   @test UInt(2) == a

   @test UInt(4) != a
end

@testset "fmpz.unary_ops..." begin
   @test -fmpz(12) == -12

   @test ~fmpz(-5) == 4
end

@testset "fmpz.abs..." begin
   @test abs(fmpz(-12)) == 12
end

@testset "fmpz.divrem..." begin
   @test fdivrem(fmpz(12), fmpz(5)) == (fmpz(2), fmpz(2))

   @test tdivrem(fmpz(12), fmpz(5)) == (fmpz(2), fmpz(2))

   @test divrem(fmpz(12), fmpz(5)) == (fmpz(2), fmpz(2))
end

@testset "fmpz.roots..." begin
   @test isqrt(fmpz(12)) == 3

   @test_throws DomainError isqrt(-fmpz(12))

   @test isqrtrem(fmpz(12)) == (3, 3)

   @test_throws DomainError isqrtrem(-fmpz(12))

   @test root(fmpz(1000), 3) == 10

   @test_throws DomainError root(-fmpz(1000), 4)

   @test_throws DomainError root(fmpz(1000), -3)
end

@testset "fmpz.extended_gcd..." begin
   @test gcdx(fmpz(12), fmpz(5)) == (1, -2, 5)
   @test gcdx(fmpz(12), 5) == (1, -2, 5)
   @test gcdx(12, fmpz(5)) == (1, -2, 5)

   @test gcdinv(fmpz(5), fmpz(12)) == (1, 5)
   @test gcdinv(fmpz(5), 12) == (1, 5)
   @test gcdinv(5, fmpz(12)) == (1, 5)

   @test_throws DomainError gcdinv(-fmpz(5), fmpz(12))

   @test_throws DomainError gcdinv(fmpz(13), fmpz(12))

   for i = -10:10
      for j = -10:10
         @test gcdx(fmpz(i), fmpz(j)) == gcdx(i, j)
      end
   end
end

@testset "fmpz.bit_twiddling..." begin
   a = fmpz(12)

   @test popcount(a) == 2

   @test nextpow2(a) == 16

   @test prevpow2(a) == 8

   @test trailing_zeros(a) == 2

   combit!(a, 2)

   @test a == 8

   @test_throws DomainError combit!(a, -1)

   setbit!(a, 0)

   @test a == 9

   @test_throws DomainError setbit!(a, -1)

   clrbit!(a, 0)

   @test a == 8

   @test_throws DomainError clrbit!(a, -1)
end

@testset "fmpz.bases..." begin
   a = fmpz(12)

   @test bin(a) == "1100"

   @test oct(a) == "14"

   @test dec(a) == "12"

   @test hex(a) == "c"

   @test base(a, 13) == "c"

   @test nbits(a) == 4

   @test ndigits(a, 3) == 3
end

@testset "fmpz.string_io..." begin
   a = fmpz(12)

   @test string(a) == "12"
end

@testset "fmpz.modular_arithmetic..." begin
   @test powmod(fmpz(12), fmpz(110), fmpz(13)) == 1

   @test_throws DomainError powmod(fmpz(12), fmpz(110), fmpz(-1))

   @test powmod(fmpz(12), 110, fmpz(13)) == 1

   @test_throws DomainError powmod(fmpz(12), 110, fmpz(-1))

   @test invmod(fmpz(12), fmpz(13)) == 12

   @test_throws DomainError invmod(fmpz(12), fmpz(-13))

   @test sqrtmod(fmpz(12), fmpz(13)) == 5

   @test_throws DomainError sqrtmod(fmpz(12), fmpz(-13))

   @test crt(fmpz(5), fmpz(13), fmpz(7), fmpz(37), true) == 44

   @test crt(fmpz(5), fmpz(13), 7, 37, false) == 44

   @test_throws DomainError crt(fmpz(5), fmpz(13), -7, 37, true)

   @test_throws DomainError crt(fmpz(5), fmpz(13), 7, -37, true)

   @test_throws DomainError crt(fmpz(5), fmpz(13), -7, -37, true)
end

@testset "fmpz.factor..." begin
   a = fmpz(-3*5*7*11*13^10)

   fact = factor(a)

   b = unit(fact)

   for (p, e) in fact
      b = b*p^e
   end

   @test b == a

   @test fact[fmpz(3)] == 1
   @test fact[fmpz(5)] == 1
   @test fact[fmpz(7)] == 1
   @test fact[fmpz(11)] == 1
   @test fact[fmpz(13)] == 10
   @test 3 in fact
   @test !(2 in fact)

   fact = factor(fmpz(-1))

   @test fact.fac == Dict{fmpz, Int}()

   fact = factor(fmpz(-2))

   @test fact.fac == Dict(fmpz(2) => 1)
   @test unit(fact) == -1

   @test_throws ArgumentError factor(fmpz(0))

   n = fmpz(2 * 1125899906842679)
   b, f = Nemo.ecm(n)
   @test mod(n, f) == 0

   n = factorial(ZZ(50))
   d, u = Nemo._factor_trial_range(n, 0, 50)
   @test isone(u)
   @test prod(p^e for (p, e) in d) == n
end

@testset "fmpz.number_theoretic..." begin
   @test isprime(fmpz(13))

   @test isprime(13)

   @test isprobable_prime(fmpz(13))

   @test divisible(fmpz(12), fmpz(6))

   @test issquare(fmpz(36))

   @test factorial(ZZ(100)) == fmpz("93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000")

   @test divisor_sigma(fmpz(128), 10) == fmpz("1181745669222511412225")

   @test_throws DomainError divisor_sigma(fmpz(1), -1)

   @test euler_phi(fmpz(12480)) == 3072

   @test fibonacci(2) == 1

   @test fibonacci(0) == 0

   @test fibonacci(-2) == -1

   @test fibonacci(fmpz(2)) == 1

   @test fibonacci(fmpz(-2)) == -1

   @test_throws DomainError  euler_phi(-fmpz(12480))

   @test remove(fmpz(12), fmpz(2)) == (2, 3)

   @test valuation(fmpz(12), fmpz(2)) == 2

   @test valuation(fmpz(12), 2) == 2

   @test valuation(12, 2) == 2

   @test divisor_lenstra(fmpz(12), fmpz(4), fmpz(5)) == 4

   @test_throws DomainError divisor_lenstra(fmpz(12), -fmpz(4), fmpz(5))
   @test_throws DomainError divisor_lenstra(fmpz(1), fmpz(4), fmpz(5))
   @test_throws DomainError divisor_lenstra(fmpz(10), fmpz(4), fmpz(3))

   @test rising_factorial(fmpz(12), 5) == 524160

   @test_throws DomainError rising_factorial(fmpz(12), -1)

   @test rising_factorial(12, 5) == 524160

   @test_throws DomainError rising_factorial(12, -1)

   @test primorial(7) == 210

   @test_throws DomainError primorial(-7)

   @test binomial(ZZ(12), ZZ(5)) == 792

   @test bell(12) == 4213597

   @test_throws DomainError bell(-1)

   @test moebius_mu(fmpz(13)) == -1

   @test_throws DomainError moebius_mu(-fmpz(1))

   @test jacobi_symbol(fmpz(2), fmpz(5)) == -1

   @test_throws DomainError jacobi_symbol(fmpz(5), fmpz(-2))

   @test_throws DomainError jacobi_symbol(fmpz(5), fmpz(2))

   @test jacobi_symbol(2, 3) == -1

   @test_throws DomainError jacobi_symbol(2, 0)

   @test_throws DomainError jacobi_symbol(-5, 4)

   if !Nemo.iswindows64()

      @test number_of_partitions(10) == 42

      @test number_of_partitions(fmpz(1000)) == fmpz("24061467864032622473692149727991")

      @test number_of_partitions(0) == 1

      @test number_of_partitions(-1) == 0

      @test number_of_partitions(fmpz(-2)) == 0
   end
end

@testset "fmpz.square_root..." begin
   @test sqrt(fmpz(4)) == 2

   @test_throws DomainError sqrt(-fmpz(4))

   @test sqrt(fmpz(0)) == 0
end
