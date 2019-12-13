using Documenter, Nemo, AbstractAlgebra

makedocs(
         format   = :html,
         sitename = "Nemo.jl",
         modules = [Nemo, AbstractAlgebra],
         clean = true,
         checkdocs = :none,
         doctest = false,
         pages    = [
             "index.md",
             "about.md",
             "types.md",
             "constructors.md",
             "Rings" => [ "integer.md",
                          "polynomial.md",
                          "mpolynomial.md",
                          "series.md",
                          "puiseux.md",
                          "residue.md"],
             "Fields" => [ "fraction.md",
                           "rational.md",
                           "arb.md",
                           "acb.md",
                           "gfp.md",
                           "finitefield.md",
			   "ff_embedding.md",
                           "numberfield.md",
                           "padic.md",
                           "qadic.md"],
             "matrix.md",
             "factor.md"
         ]
)

deploydocs(
   julia = "1.0",
   repo   = "github.com/Nemocas/Nemo.jl.git",
   target = "build",
   deps = nothing,
   make   = nothing,
   osname = "linux"
)
