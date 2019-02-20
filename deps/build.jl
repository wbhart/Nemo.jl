using Libdl

using BinaryProvider

# Parse some basic command-line arguments
const verbose = "--verbose" in ARGS

# Dependencies that must be installed before this package can be built
dependencies = [
  "https://github.com/JuliaMath/GMPBuilder/releases/download/v6.1.2-2/build_GMP.v6.1.2.jl",
  "https://github.com/JuliaMath/MPFRBuilder/releases/download/v4.0.1-3/build_MPFR.v4.0.1.jl",
  "https://github.com/thofma/Flint2Builder/releases/download/2b8f8acb/build_libflint.v2.0.0-b8f8acb317c265db99f828e7baf3266f07f92a7.jl",
  "https://github.com/thofma/ArbBuilder/releases/download/v2.16.0/build_libarb.v2.16.0.jl",
  "https://github.com/thofma/AnticBuilder/releases/download/v0.2.0/build_libantic.v0.0.1.jl"
 ]

const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))

const prefixpath = joinpath(@__DIR__, "usr")

products = []

for url in dependencies
		build_file = joinpath(@__DIR__, basename(url))
		if !isfile(build_file)
				download(url, build_file)
		end
end

# Execute the build scripts for the dependencies in an isolated module to avoid overwriting
# any variables/constants here
for url in dependencies
		build_file = joinpath(@__DIR__, basename(url))
		m = @eval module $(gensym()); include($build_file); end
		append!(products, m.products)
end

push!(Libdl.DL_LOAD_PATH, joinpath(prefixpath, "lib"), joinpath(prefixpath, "bin"))
