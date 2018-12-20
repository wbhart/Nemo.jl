include("generic/Poly-test.jl")
include("generic/MPoly-test.jl")
include("generic/Matrix-test.jl")

function test_generic()
   test_Poly()
   test_MPoly()
   test_Matrix()
end
