function benchmark_minpoly_finite_field()
   print("benchmark_minpoly_finite_field ... ")

   F, s = FlintFiniteField(103, 2, "s")

   M = MatrixSpace(F, 80, 80)()

   for i in 1:40
      for j in 1:40
         r = F(rand(1:103)) + s*(rand(1:103))
         M[i, j] = r
         M[40 + i, 40 + j] = r
      end
   end

   for i in 1:10
      similarity!(M, 80, F(rand(-3:3)))
   end

   tt = @elapsed minpoly(PolynomialRing(F, "x")[1], M)
   println("$tt")
end

