@static if VERSION < v"1.3.0"

  using BinaryProvider # requires BinaryProvider 0.3.0 or later

  # Parse some basic command-line arguments
  const verbose = "--verbose" in ARGS
  const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))

  products = [ ]

  dependencies = [
  "build_GMP.v6.1.2.jl",
  "build_MPFR.v4.0.2.jl",
  "build_FLINT.v0.0.1.jl",
  "build_Arb.v2.17.0.jl",
  "build_Antic.v0.1.0.jl",
   ]

  for file in dependencies
       build_file = joinpath(@__DIR__, file)
       m = @eval module $(gensym()); include($build_file); end
       if !occursin("FLINT", file) && !occursin("GMP", file) && !occursin("MPFR", file)
         append!(products, m.products)
       end
  end

  # Finally, write out a deps.jl file
  write_deps_file(joinpath(@__DIR__, "deps.jl"), Array{Product,1}(products), verbose=verbose)

end # VERSION
