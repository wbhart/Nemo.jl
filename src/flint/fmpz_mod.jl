###############################################################################
#
#   fmpz_mod.jl : Nemo fmpz_mod (integers modulo large n)
#
###############################################################################

###############################################################################
#
#   GCD
#
###############################################################################

@doc Markdown.doc"""
    gcdx(a::ResElem{fmpz}, b::ResElem{fmpz})
> Compute the extended gcd with the Euclidean structure inherited from
> $\mathbb{Z}$.
"""
function gcdx(a::ResElem{fmpz}, b::ResElem{fmpz})
  m = modulus(a)
  R = parent(a)
  g, u, v = gcdx(fmpz(a.data), fmpz(b.data))
  G, U, V = gcdx(g, fmpz(m))
  return R(G), R(U)*R(u), R(U)*R(v)
end
