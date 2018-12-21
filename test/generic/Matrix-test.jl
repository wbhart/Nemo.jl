function test_Matrix_binary_ops_delayed_reduction()
   print("Matrix.binary_ops_delayed_reduction...")

   R, t = PolynomialRing(QQ, "t")
   K, a = NumberField(t^3 + 3t + 1, "a")
   S = MatrixSpace(K, 5, 5)

   for iter = 1:10
      f = rand(S, -10:10)
      g = rand(S, -10:10)
      h = rand(S, -10:10)

      @test f*(g + h) == f*g + f*h
      @test f*(g - h) == f*g - f*h
   end

   println("PASS")
end

function test_Matrix_lu_delayed_reduction()
   print("Matrix.lu_delayed_reduction...")

   R, t = PolynomialRing(QQ, "t")
   K, a = NumberField(t^3 + 3t + 1, "a")
   S = MatrixSpace(K, 5, 5)

   for iter = 1:10
      rk = rand(0:5)
      A = randmat_with_rank(S, rk, -10:10)

      r, P, L, U = lu(A)

      @test r == rk
      @test P*A == L*U
   end

   println("PASS")
end

function test_Matrix_fflu_delayed_reduction()
   print("Matrix.fflu_delayed_reduction...")

   R, t = PolynomialRing(QQ, "t")
   K, a = NumberField(t^3 + 3t + 1, "a")
   S = MatrixSpace(K, 5, 5)

   for iter = 1:20
      rk = rand(0:5)
      A = randmat_with_rank(S, rk, -10:10)

      r, d, P, L, U = fflu(A)

      if r == 5
         D = S()
         D[1, 1] = inv(U[1, 1])
         D[2, 2] = inv(U[1, 1]*U[2, 2])
         D[3, 3] = inv(U[2, 2]*U[3, 3])
         D[4, 4] = inv(U[3, 3]*U[4, 4])
         D[5, 5] = inv(U[4, 4])
      end

      @test r == rk && (r < 5 || P*A == L*D*U)
   end

   println("PASS")
end

function test_Matrix_minpoly_delayed_reduction()
   print("Matrix.minpoly_delayed_reduction...")

   # Tests reduce_row!

   R, t = PolynomialRing(QQ, "t")
   K, a = NumberField(t^3 + 3t + 1, "a")
   S = MatrixSpace(K, 6, 6)
   U, z = PolynomialRing(K, "z")

   M = S()
   for i = 1:3
      for j = 1:3
         M[i, j] = rand(K, -10:10)
         M[i + 3, j + 3] = deepcopy(M[i, j])
      end
   end

   p1 = minpoly(U, M)

   for i = 1:10
      similarity!(M, rand(1:6), rand(K, -3:3))
   end

   p2 = minpoly(U, M)

   @test p1 == p2

   println("PASS")
end

function test_Matrix_solve_fflu_delayed_reduction()
   print("Matrix.solve_fflu_delayed_reduction...")

   R, t = PolynomialRing(QQ, "t")
   K, a = NumberField(t^3 + 3t + 1, "a")

   for dim = 0:5
      S = MatrixSpace(K, dim, dim)
      T = MatrixSpace(K, dim, rand(1:5))

      M = randmat_with_rank(S, dim, -100:100)
      b = rand(T, -100:100)

      x, d = Generic.solve_fflu(M, b)

      @test divexact(M, d)*x == b
   end

   println("PASS")
end

function test_Matrix_solve_lu_delayed_reduction()
   print("Matrix.solve_lu_delayed_reduction...")

   R, t = PolynomialRing(QQ, "t")
   K, a = NumberField(t^3 + 3t + 1, "a")

   for dim = 0:5
      S = MatrixSpace(K, dim, dim)
      T = MatrixSpace(K, dim, rand(1:5))

      M = randmat_with_rank(S, dim, -100:100)
      b = rand(T, -100:100)

      x = Generic.solve_lu(M, b)

      @test M*x == b
   end

   println("PASS")
end

#=
   TODO: Add tests for the following when there are rings that are not fields
         that have delayed reduction
     - fflu over a ring
     - rref over a ring
     - minpoly over an integrally closed domain

   Note: backsolve! is also not tested as it appears to be unused.
=#

function test_Matrix()
   test_Matrix_binary_ops_delayed_reduction()
   test_Matrix_lu_delayed_reduction()
   test_Matrix_fflu_delayed_reduction()
   test_Matrix_minpoly_delayed_reduction()
   test_Matrix_solve_fflu_delayed_reduction()
   test_Matrix_solve_lu_delayed_reduction()

   println("")
end

