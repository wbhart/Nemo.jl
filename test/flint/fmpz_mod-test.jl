function test_fmpz_mod_gcdx()
   print("fmpz_mod.gcdx...")

   for i = 1:100
      n = rand(1:24)
      R = ResidueRing(ZZ, ZZ(n))

      for iter = 1:100
         a = rand(R, 0:n-1)
         b = rand(R, 0:n-1)

         g, s, t = gcdx(a, b)

         @test g == s*a + t*b
      end
   end

   println("PASS")
end

function test_fmpz_mod()
   test_fmpz_mod_gcdx()

   println("")
end
