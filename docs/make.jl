using Documenter, Nemo

makedocs(
         format   = :html,
         sitename = "Nemo.jl",
         pages    = [
             "index.md",
             "about.md",
             "types.md",
             "constructors.md",
             "Rings" => [ "integer.md",
                          "polynomial.md",
                          "series.md",
                          "residue.md"],
             "Fields" => [ "fraction.md",
                           "rational.md",
                           "arb.md",
                           "acb.md",
                           "finitefield.md",
                           "numberfield.md",
                           "padic.md"],
             "matrix.md",
             "factor.md"
         ]
)

deploydocs(
   repo   = "github.com/Nemocas/Nemo.jl.git",
   target = "build",
   deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
   make   = nothing
)

