function randelem(R::Nemo.GaloisField, n)
   return rand(R, -n:n)
end

function randmat(R::GFPMatSpace, d::Int)
   m = R.rows
   n = R.cols
   r = R()
   for i = 1:m
      for j = 1:n
         r[i, j] = randelem(base_ring(R), d)
      end
   end
   return r
end

function test_gfp_mat_constructors()
  print("gfp_mat.constructors...")

  Z2 = GF(2)
  Z3 = GF(3)

  R = GFPMatSpace(Z2, 2, 2)

  @test elem_type(R) == gfp_mat
  @test elem_type(GFPMatSpace) == gfp_mat
  @test parent_type(gfp_mat) == GFPMatSpace

  @test isa(R, GFPMatSpace)

  @test base_ring(R) == Z2

  S = GFPMatSpace(Z3, 2, 2)

  @test isa(S, GFPMatSpace)

  RR = GFPMatSpace(Z2, 2, 2)

  @test isa(RR, GFPMatSpace)

  @test R == RR

  @test_throws ErrorException GFPMatSpace(Z2, 2, -1)
  @test_throws ErrorException GFPMatSpace(Z2, -1, 2)
  @test_throws ErrorException GFPMatSpace(Z2, -1, -1)

  a = R()

  @test isa(a, gfp_mat)
  @test parent(a) == R

  ar = [ BigInt(1) BigInt(1); BigInt(1) BigInt(1) ]

  b = R(ar)

  @test isa(b, gfp_mat)
  @test parent(b) == R
  @test rows(b) == 2 && cols(b) == 2
  @test_throws ErrorConstrDimMismatch R(reshape(ar,1,4))
  @test b == R([BigInt(1), BigInt(1), BigInt(1), BigInt(1)])
  @test_throws ErrorConstrDimMismatch R([BigInt(1) BigInt(1)])
  @test_throws ErrorConstrDimMismatch R([BigInt(1) BigInt(1) ; BigInt(1) BigInt(1) ;
                                 BigInt(1) BigInt(1)])
  @test_throws ErrorConstrDimMismatch R([BigInt(1), BigInt(1), BigInt(1)])
  @test_throws ErrorConstrDimMismatch R([BigInt(1), BigInt(1),
                                  BigInt(1), BigInt(1), BigInt(1)])

  ar = [ ZZ(1) ZZ(1); ZZ(1) ZZ(1) ]

  c = R(ar)
  @test isa(c, gfp_mat)
  @test parent(c) == R
  @test rows(c) == 2 && cols(c) == 2
  @test_throws ErrorConstrDimMismatch R(reshape(ar,4,1))
  @test c == R([ ZZ(1), ZZ(1), ZZ(1), ZZ(1)])
  @test_throws ErrorConstrDimMismatch R([ZZ(1) ZZ(1)])
  @test_throws ErrorConstrDimMismatch R([ZZ(1) ZZ(1) ; ZZ(1) ZZ(1) ; ZZ(1) ZZ(1)])
  @test_throws ErrorConstrDimMismatch R([ZZ(1), ZZ(1), ZZ(1)])
  @test_throws ErrorConstrDimMismatch R([ZZ(1), ZZ(1), ZZ(1), ZZ(1), ZZ(1)])

  ar = [ 1 1; 1 1]

  d = R(ar)

  @test isa(d, gfp_mat)
  @test parent(d) == R
  @test rows(d) == 2 && cols(d) == 2
  @test_throws ErrorConstrDimMismatch R(reshape(ar,1,4))
  @test d == R([1,1,1,1])
  @test_throws ErrorConstrDimMismatch R([1 1 ])
  @test_throws ErrorConstrDimMismatch R([1 1 ; 1 1 ; 1 1 ])
  @test_throws ErrorConstrDimMismatch R([1, 1, 1])
  @test_throws ErrorConstrDimMismatch R([1, 1, 1, 1, 1])

  ar = MatrixSpace(ZZ, 2, 2)([ 1 1; 1 1])

  e = R(ar)

  @test isa(e, gfp_mat)
  @test parent(e) == R
  @test rows(e) == 2 && cols(e) == 2

  ar = MatrixSpace(ZZ, 3, 3)([ 1 1 1 ; 1 1 1; 1 1 1])

  @test_throws ErrorException R(ar)

  ar = [ Z2(1) Z2(1); Z2(1) Z2(1) ]

  f = R(ar)

  @test isa(f, gfp_mat)
  @test parent(f) == R
  @test rows(f) == 2 && cols(f) == 2
  @test_throws ErrorConstrDimMismatch R(reshape(ar,4,1))
  @test f == R([Z2(1), Z2(1), Z2(1), Z2(1)])
  @test_throws ErrorConstrDimMismatch R([Z2(1) Z2(1) ])
  @test_throws ErrorConstrDimMismatch R([Z2(1) Z2(1) ; Z2(1) Z2(1) ; Z2(1) Z2(1) ])
  @test_throws ErrorConstrDimMismatch R([Z2(1), Z2(1), Z2(1)])
  @test_throws ErrorConstrDimMismatch R([Z2(1), Z2(1), Z2(1), Z2(1), Z2(1)])

  @test isa(S(1), gfp_mat)

  @test isa(S(fmpz(1)), gfp_mat)

  @test isa(S(Z3(1)), gfp_mat)

  g = deepcopy(e)

  @test b == c
  @test c == d
  @test d == e
  @test e == f
  @test g == e

   arr = [1 2; 3 4]
   arr2 = [1, 2, 3, 4, 5, 6]

   for T in [Z3, fmpz, Int, BigInt]
      M = matrix(Z3, map(T, arr))
      @test isa(M, gfp_mat)
      @test M.base_ring == Z3

      M2 = matrix(Z3, 2, 3, map(T, arr2))
      @test isa(M2, gfp_mat)
      @test M2.base_ring == Z3
      @test rows(M2) == 2
      @test cols(M2) == 3
      @test_throws ErrorConstrDimMismatch matrix(Z3, 2, 2, map(T, arr2))
      @test_throws ErrorConstrDimMismatch matrix(Z3, 2, 4, map(T, arr2))
   end

   M3 = zero_matrix(Z3, 2, 3)

   @test isa(M3, gfp_mat)
   @test M3.base_ring == Z3

   M4 = identity_matrix(Z3, 3)

   @test isa(M4, gfp_mat)
   @test M4.base_ring == Z3

  println("PASS")
end

function test_gfp_mat_printing()
  print("gfp_mat.printing...")

  Z2 = GF(2)
  R = GFPMatSpace(Z2, 2, 2)

  a = R(1)

  @test string(a) == "[1 0]\n[0 1]"

  println("PASS")
end

function test_gfp_mat_manipulation()
  print("gfp_mat.manipulation...")
  Z11 = GF(11)
  R = GFPMatSpace(Z11, 2, 2)
  Z23 = GF(23)
  S = GFPMatSpace(Z23, 2, 2)

  ar = [ 1 2; 3 4]

  a = R(ar)
  aa = S(ar)

  @test rows(a) == 2
  @test cols(a) == 2

  b = deepcopy(a)

  c = R([ 1 3; 2 4])

  @test a[1,1] == Z11(1)

  a[1,1] = UInt(2)

  @test a[1,1] == Z11(2)
  @test_throws BoundsError a[0,-1] = Z11(2)

  a[2,1] = ZZ(3)

  @test a[2,1] == Z11(ZZ(3))
  @test_throws BoundsError a[-10,-10] = ZZ(3)

  a[2,2] = Z11(4)

  @test a[2,2] == Z11(4)
  @test_throws BoundsError a[-2,2] = Z11(4)

  a[1,2] = 5

  @test a[1,2] == Z11(5)
  @test_throws BoundsError a[-2,2] = 5

  @test a != b

  d = one(R)

  @test isa(d, gfp_mat)

  e = zero(R)

  @test isa(e, gfp_mat)

  @test iszero(e)

  @test_throws ErrorException one(MatrixSpace(GF(2), 1, 2))

  @test issquare(a)

  @test a == a
  @test a == deepcopy(a)
  @test a != aa

  @test transpose(b) == c

  transpose!(b)

  @test b == c

  @test transpose(MatrixSpace(Z11,1,2)([ 1 2; ])) ==
          MatrixSpace(Z11,2,1)(reshape([ 1 ; 2],2,1))

  @test_throws ErrorConstrDimMismatch transpose!(R([ 1 2 ;]))

  println("PASS")
end

function test_gfp_mat_unary_ops()
  print("gfp_mat.unary_ops...")

  Z17 = GF(17)

  R = MatrixSpace(Z17, 3, 4)
  RR = MatrixSpace(Z17, 4, 3)
  Z2 = GF(2)
  S = MatrixSpace(Z2, 3, 4)

  a = R([ 1 2 3 1; 3 2 1 2; 1 3 2 0])

  b = R([ 2 1 0 1; 0 0 0 0; 0 1 2 0 ])

  c = R()

  d = -b

  @test d == R([ 15 16 0 16; 0 0 0 0; 0 16 15 0])

  println("PASS")
end

function test_gfp_mat_binary_ops()
  print("gfp_mat.binary_ops...")

  Z17 = GF(17)

  R = MatrixSpace(Z17, 3, 4)

  a = R([ 1 2 3 1; 3 2 1 2; 1 3 2 0])

  b = R([ 2 1 0 1; 0 0 0 0; 0 1 2 0 ])

  c = R()

  d = a + b

  @test d == R([3 3 3 2; 3 2 1 2; 1 4 4 0])

  d = a - b

  @test d == R([ 16 1 3 0; 3 2 1 2; 1 2 0 0 ])

  d = a*transpose(a)

  @test d == MatrixSpace(Z17, 3, 3)([15 12 13; 12 1 11; 13 11 14])

  d = transpose(a)*a

  @test d == MatrixSpace(Z17, 4, 4)([11 11 8 7; 11 0 14 6; 8 14 14 5; 7 6 5 5])

  println("PASS")
end

function test_gfp_mat_adhoc_binary()
  print("gfp_mat.adhoc_binary...")

  Z17 = GF(17)

  R = MatrixSpace(Z17, 3, 4)
  Z2 = GF(2)

  a = R([ 1 2 3 1; 3 2 1 2; 1 3 2 0])

  b = R([ 2 1 0 1; 0 0 0 0; 0 1 2 0 ])

  c = R()
  d = UInt(2)*a
  dd = a*UInt(2)

  @test d == R([ 2 4 6 2; 6 4 2 4; 2 6 4 0])
  @test dd == R([ 2 4 6 2; 6 4 2 4; 2 6 4 0])


  d = 2*a
  dd = a*2

  @test d == R([ 2 4 6 2; 6 4 2 4; 2 6 4 0])
  @test dd == R([ 2 4 6 2; 6 4 2 4; 2 6 4 0])

  d = ZZ(2)*a
  dd = a*ZZ(2)

  @test d == R([ 2 4 6 2; 6 4 2 4; 2 6 4 0])
  @test dd == R([ 2 4 6 2; 6 4 2 4; 2 6 4 0])

  d = Z17(2)*a
  dd = a*Z17(2)

  @test d == R([ 2 4 6 2; 6 4 2 4; 2 6 4 0])
  @test dd == R([ 2 4 6 2; 6 4 2 4; 2 6 4 0])

  @test_throws ErrorException Z2(1)*a

  println("PASS")
end

function test_gfp_mat_comparison()
  print("gfp_mat.comparison...")

  Z17 = GF(17)

  R = MatrixSpace(Z17, 3, 4)

  a = R([ 1 2 3 1; 3 2 1 2; 1 3 2 0])

  @test a == a

  @test deepcopy(a) == a

  @test a != R([0 1 3 1; 2 1 4 2; 1 1 1 1])

  println("PASS")
end

function test_gfp_mat_adhoc_comparison()
  print("gfp_mat.comparison...")

  Z17 = GF(17)

  R = MatrixSpace(Z17, 3, 4)

  @test R(5) == 5
  @test R(5) == fmpz(5)
  @test R(5) == Z17(5)

  @test 5 == R(5)
  @test fmpz(5) == R(5)
  @test Z17(5) == R(5)

  println("PASS")
end

function test_gfp_mat_powering()
  print("gfp_mat.powering...")

  Z17 = GF(17)

  R = MatrixSpace(Z17, 3, 4)

  a = R([ 1 2 3 1; 3 2 1 2; 1 3 2 0])

  f = a*transpose(a)

  g = f^1000

  @test g == MatrixSpace(Z17, 3, 3)([1 2 2; 2 13 12; 2 12 15])

  g = f^ZZ(1000)

  @test g == MatrixSpace(Z17, 3, 3)([1 2 2; 2 13 12; 2 12 15])

  @test_throws ErrorException f^(ZZ(2)^1000)

  println("PASS")
end

function test_gfp_mat_row_echelon_form()
  print("gfp_mat.row_echelon_form...")

  Z17 = GF(17)
  R = MatrixSpace(Z17, 3, 4)
  RR = MatrixSpace(Z17, 4, 3)
  Z2 = GF(2)
  S = MatrixSpace(Z2, 3, 4)

  a = R([ 1 2 3 1; 3 2 1 2; 1 3 2 0])

  b = R([ 2 1 0 1; 0 0 0 0; 0 1 2 0 ])

  b = transpose(b)

  c = a*transpose(a)

  r, d = rref(a)

  @test d == R([ 1 0 0 8; 0 1 0 15; 0 0 1 16])
  @test r == 3

  r = rref!(a)

  @test a == R([ 1 0 0 8; 0 1 0 15; 0 0 1 16])
  @test r == 3

  r, d = rref(b)

  @test d == parent(b)([ 1 0 0 ; 0 0 1; 0 0 0; 0 0 0])
  @test r == 2

  println("PASS")
end

function test_gfp_mat_howell_form()
  print("gfp_mat.howell_form...")

  Z17 = GF(17)
  R = MatrixSpace(Z17, 3, 3)

  a = R([4 1 0; 0 0 5; 0 0 0 ])

  b = R([8 5 5; 0 9 8; 0 0 10])

  c = R([1 13 0; 0 0 1; 0 0 0])

  d = R([4 0 0; 0 0 1; 0 0 0])

  @test howell_form(a) == c
  @test howell_form(b) == 1
  @test strong_echelon_form(d) == R([1 0 0; 0 0 0; 0 0 1])

  a = matrix(Z17, [0 1 0 0; 0 0 1 0; 0 0 0 1; 0 0 0 0])
  @test strong_echelon_form(a) == matrix(Z17, [0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1])

  println("PASS")
end

function test_gfp_mat_trace_det()
  print("gfp_mat.trace_det...")

  Z17 = GF(17)
  R = MatrixSpace(Z17, 3, 4)
  RR = MatrixSpace(Z17, 4, 3)
  Z2 = GF(2)
  S = MatrixSpace(Z2, 3, 4)

  a = R([ 1 2 3 1; 3 2 1 2; 1 3 2 0])

  aa = MatrixSpace(Z17,3,3)([ 1 2 3; 3 2 1; 1 1 2])

  b = R([ 2 1 0 1; 0 0 0 0; 0 1 2 0 ])

  a = transpose(a)*a

  c = tr(a)

  @test c == Z17(13)

  @test_throws ErrorException tr(b)

  c = det(a)

  @test c == zero(Z17)

  @test_throws ErrorException det(b)

  c = det(aa)

  @test c == Z17(13)

  a = R([ 1 2 3 1; 3 2 1 2; 1 3 2 0])

  println("PASS")
end

function test_gfp_mat_rank()
  print("gfp_mat.rank...")

  Z17 = GF(17)
  R = MatrixSpace(Z17, 3, 4)
  RR = MatrixSpace(Z17, 4, 3)
  Z2 = GF(2)
  S = MatrixSpace(Z2, 3, 4)

  a = R([ 1 2 3 1; 3 2 1 2; 1 3 2 0])

  aa = MatrixSpace(Z17,3,3)([ 1 2 3; 3 2 1; 1 1 2])

  b = R([ 2 1 0 1; 0 0 0 0; 0 1 2 0 ])

  c = rank(a)

  @test c == 3

  c = rank(aa)

  @test c == 3

  c = rank(b)

  @test c == 2

  println("PASS")
end

function test_gfp_mat_inv()
  print("gfp_mat.inv...")

  Z17 = GF(17)
  R = MatrixSpace(Z17, 3, 4)
  RR = MatrixSpace(Z17, 4, 3)
  Z2 = GF(2)
  S = MatrixSpace(Z2, 3, 4)

  a = R([ 1 2 3 1; 3 2 1 2; 1 3 2 0])

  aa = MatrixSpace(Z17,3,3)([ 1 2 3; 3 2 1; 1 1 2])

  b = R([ 2 1 0 1; 0 0 0 0; 0 1 2 0 ])

  c = inv(aa)

  @test c == parent(aa)([12 13 1; 14 13 15; 4 4 1])

  @test_throws ErrorException inv(a)

  @test_throws ErrorException inv(transpose(a)*a)

  println("PASS")
end

function test_gfp_mat_solve()
  print("gfp_mat.solve...")

  Z17 = GF(17)
  R = MatrixSpace(Z17, 3, 3)
  S = MatrixSpace(Z17, 3, 4)

  a = R([ 1 2 3 ; 3 2 1 ; 0 0 2 ])

  b = S([ 2 1 0 1; 0 0 0 0; 0 1 2 0 ])

  c = a*b

  d = solve(a,c)

  @test d == b

  a = zero(R)

  @test_throws ErrorException  solve(a,c)

  println("PASS")
end

function test_gfp_mat_lu()
  print("gfp_mat.lu...")


  Z17 = GF(17)
  R = MatrixSpace(Z17, 3, 3)
  S = MatrixSpace(Z17, 3, 4)

  a = R([ 1 2 3 ; 3 2 1 ; 0 0 2 ])

  b = S([ 2 1 0 1; 0 0 0 0; 0 1 2 0 ])

  r, P, l, u = lu(a)

  @test l*u == P*a

  r, P, l, u = lu(b)

  @test l*u == S([ 2 1 0 1; 0 1 2 0; 0 0 0 0])

  @test l*u == P*b

  c = matrix(Z17, 6, 3, [0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1])

  r, P, l, u = lu(c)

  @test r == 3
  @test l*u == P*c

  println("PASS")
end

function test_gfp_mat_swap_rows()
  print("gfp_mat.swap_rows...")

  Z17 = GF(17)

  A = matrix(Z17, 5, 1, [1, 2, 3, 4, 5])

  B = swap_rows(A, 3, 4)
  @test B == matrix(Z17, 5, 1, [1, 2, 4, 3, 5])

  swap_rows!(A, 3, 4)
  @test A == matrix(Z17, 5, 1, [1, 2, 4, 3, 5])

  @test_throws ErrorException swap_rows(A, 0, 5)
  @test_throws ErrorException swap_rows(A, 4, 6)

  println("PASS")
end
function test_gfp_mat_view()
  print("gfp_mat.view...")

  Z17 = GF(17)
  R = MatrixSpace(Z17, 3, 3)
  S = MatrixSpace(Z17, 3, 4)

  a = R([ 1 2 3 ; 3 2 1 ; 0 0 2 ])

  b = S([ 2 1 0 1; 0 0 0 0; 0 1 2 0 ])

  t = view(a, 1, 1, 3, 3)

  @test t == a

  @test view(a, 1, 1, 3, 3) == view(a, 1:3, 1:3)
  @test view(a, 1, 1, 3, 3) == sub(a, 1, 1, 3, 3)
  @test view(a, 1, 1, 3, 3) == sub(a, 1:3, 1:3)

  t = view(a, 1, 1, 2, 2)

  @test t == Z17[1 2; 3 2]

  t = view(a, 2, 2, 3, 2)

  @test t == transpose(Z17[2 0])

  @test view(a, 2, 2, 3, 2) == view(a, 2:3,  2:2)
  @test view(a, 2, 2, 3, 2) == sub(a, 2, 2, 3, 2)
  @test view(a, 2, 2, 3, 2) == sub(a, 2:3, 2:2)

  @test_throws BoundsError view(a, 4, 4, 1, 1)

  @test_throws ErrorException view(a, 2, 2, 1, 1)

  S = MatrixSpace(Z17, 3, 3)

  A = S([1 2 3; 4 5 6; 7 8 9])

  B = @inferred view(A, 1, 1, 2, 2)

  @test typeof(B) == gfp_mat
  @test B == MatrixSpace(Z17, 2, 2)([1 2; 4 5])

  B[1, 1] = 10
  @test A[1, 1] == 10

  C = @inferred view(B, 1:2, 1:2)

  @test typeof(C) == gfp_mat
  @test C == MatrixSpace(Z17, 2, 2)([10 2; 4 5])

  C[1, 1] = 20
  @test B[1, 1] == 20
  @test A[1, 1] == 20

  A = 0
  GC.gc()
  @test B[1, 1] == 20

  println("PASS")
end

function test_gfp_mat_sub()
   print("gfp_mat.sub...")

   Z17 = GF(17)
   S = MatrixSpace(Z17, 3, 3)

   A = S([1 2 3; 4 5 6; 7 8 9])

   B = @inferred sub(A, 1, 1, 2, 2)

   @test typeof(B) == gfp_mat
   @test B == MatrixSpace(Z17, 2, 2)([1 2; 4 5])

   B[1, 1] = 10
   @test A == S([1 2 3; 4 5 6; 7 8 9])

   C = @inferred sub(B, 1:2, 1:2)

   @test typeof(C) == gfp_mat
   @test C == MatrixSpace(Z17, 2, 2)([10 2; 4 5])

   C[1, 1] = 20
   @test B == MatrixSpace(Z17, 2, 2)([10 2; 4 5])
   @test A == S([1 2 3; 4 5 6; 7 8 9])

   println("PASS")
end

function test_gfp_mat_concatenation()
  print("gfp_mat.concatenation...")

  Z17 = GF(17)
  R = MatrixSpace(Z17, 3, 3)
  S = MatrixSpace(Z17, 3, 4)

  a = R([ 1 2 3 ; 3 2 1 ; 0 0 2 ])

  b = S([ 2 1 0 1; 0 0 0 0; 0 1 2 0 ])

  c = hcat(a,a)

  @test c == MatrixSpace(Z17, 3, 6)([1, 2, 3, 1, 2, 3,
                                     3, 2, 1, 3, 2, 1,
                                     0, 0, 2, 0, 0, 2])

  c = hcat(a,b)

  @test c == MatrixSpace(Z17, 3, 7)([1, 2, 3, 2, 1, 0, 1,
                                     3, 2, 1, 0, 0, 0, 0,
                                     0, 0, 2, 0, 1, 2, 0])

  @test_throws ErrorException c = hcat(a,transpose(b))

  c = vcat(a,transpose(b))

  @test c == MatrixSpace(Z17, 7, 3)([1, 2, 3,
                                     3, 2, 1,
                                     0, 0, 2,
                                     2, 0, 0,
                                     1, 0, 1,
                                     0, 0, 2,
                                     1, 0, 0])

  @test_throws ErrorException vcat(a,b)

  println("PASS")
end

function test_gfp_mat_conversion()
  print("gfp_mat.conversion...")

  Z17 = GF(17)
  R = MatrixSpace(Z17, 3, 3)
  S = MatrixSpace(ZZ, 3, 3)

  c = S()

  a = R([ 1 2 3 ; 3 2 1 ; 0 0 2 ])

  @test Array(a) == [Z17(1) Z17(2) Z17(3);
                     Z17(3) Z17(2) Z17(1);
                     Z17(0) Z17(0) Z17(2) ]

  b = lift(a)

  @test b == S([ 1 2 3; 3 2 1; 0 0 2])
  @test parent(b) == S

  lift!(c,a)

  @test c == S([ 1 2 3; 3 2 1; 0 0 2])

  println("PASS")
end

function test_gfp_mat_charpoly()
   print("gfp_mat.charpoly...")

   R = GF(17)

   for dim = 0:5
      S = MatrixSpace(R, dim, dim)
      U, x = PolynomialRing(R, "x")

      for i = 1:10
         M = randmat(S, 5)

         p1 = charpoly(U, M)
         p2 = charpoly_danilevsky!(U, M)

         @test p1 == p2
      end
   end

   println("PASS")
end

function test_gfp_mat()
  test_gfp_mat_constructors()
  test_gfp_mat_printing()
  test_gfp_mat_manipulation()
  test_gfp_mat_unary_ops()
  test_gfp_mat_binary_ops()
  test_gfp_mat_adhoc_binary()
  test_gfp_mat_comparison()
  test_gfp_mat_adhoc_comparison()
  test_gfp_mat_powering()
  test_gfp_mat_row_echelon_form()
  test_gfp_mat_howell_form()
  test_gfp_mat_trace_det()
  test_gfp_mat_rank()
  test_gfp_mat_inv()
  test_gfp_mat_lu()
  test_gfp_mat_view()
  test_gfp_mat_swap_rows()
  test_gfp_mat_sub()
  test_gfp_mat_concatenation()
  test_gfp_mat_conversion()
  test_gfp_mat_charpoly()

  println("")
end
