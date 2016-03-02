#################################################################################
#
#   Factorisation
#
#################################################################################

<<<<<<< HEAD
function _Factor(f::PariFactor, R::Ring)
=======
function pari_factor_to_dict(f::PariFactor, R::Ring)
>>>>>>> 508a753c6b74f2b54f565639dbfcd6fd476618ab
  D = Dict{typeof(zero(R)), Int}()
  for i=1:f.len
    p, n = f[i]
    D[R(p)] = n
  end
  return D
end

function factor(n::fmpz)
   f = factor(pari(n))
<<<<<<< HEAD
   return _Factor(f, FlintZZ)
=======
   return pari_factor_to_dict(f, FlintZZ)
>>>>>>> 508a753c6b74f2b54f565639dbfcd6fd476618ab
end

function factor(g::fmpz_poly)
   h = pari(g)
   f = factor(h)
<<<<<<< HEAD
   return _Factor(f, g.parent)
=======
   return pari_factor_to_dict(f, g.parent)
>>>>>>> 508a753c6b74f2b54f565639dbfcd6fd476618ab
end

function factor(g::fmpq_poly)
   h = pari(g)
   f = factor(h)
<<<<<<< HEAD
   return _Factor(f, g.parent)
=======
   return pari_factor_to_dict(f, g.parent)
>>>>>>> 508a753c6b74f2b54f565639dbfcd6fd476618ab
end

