include("Rings-test.jl")
include("Fields-test.jl")
include("Generic-test.jl")
include("Benchmark-test.jl")

function test_all()
   test_rings()
   test_fields()
   test_generic()
   test_benchmarks()
end
