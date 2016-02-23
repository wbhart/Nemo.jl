function benchmark_det_poly_ring()
   print("benchmark_det_poly_ring ... ")
   ZZx, x = PolynomialRing(FlintZZ, "x")
   M = MatrixSpace(ZZx, 40, 40)()

   for i in 1:40
      for j in 1:40
         M[i, j] = ZZx(map(fmpz, rand(-20:20, 3)))
      end
   end

   tt = @elapsed determinant(M)
   println("$tt")
end
