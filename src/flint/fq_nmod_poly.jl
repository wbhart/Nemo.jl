################################################################################
#
#  fq_nmod_poly.jl: Flint fq_mod_poly (Polynomials over FqNmodFiniteField)
#
################################################################################

export fq_nmod_poly, FqNmodPolyRing

################################################################################
#
#  Type and parent object methods
#
################################################################################

parent_type(::Type{fq_nmod_poly}) = FqNmodPolyRing

elem_type(::FqNmodPolyRing) = fq_nmod_poly

base_ring(a::FqNmodPolyRing) = a.base_ring

parent(a::fq_nmod_poly) = a.parent

var(a::FqNmodPolyRing) = a.S

function check_parent(a::fq_nmod_poly, b::fq_nmod_poly) 
   a.parent != b.parent &&
         error("Operations on distinct polynomial rings not supported")
end

################################################################################
#
#   Basic manipulation
#
################################################################################
   
length(x::fq_nmod_poly) = ccall((:fq_nmod_poly_length, :libflint), Int,
                                (Ptr{fq_nmod_poly},), &x)

set_length!(x::fq_nmod_poly, n::Int) = ccall((:_fq_nmod_poly_set_length, :libflint), Void,
                              (Ptr{fq_nmod_poly}, Int), &x, n)

function coeff(x::fq_nmod_poly, n::Int)
   n < 0 && throw(DomainError())
   F = (x.parent).base_ring
   temp = F(1)
   ccall((:fq_nmod_poly_get_coeff, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod_poly}, Int, Ptr{FqNmodFiniteField}),
         &temp, &x, n, &F)
   return temp
end

zero(a::FqNmodPolyRing) = a(zero(base_ring(a)))

one(a::FqNmodPolyRing) = a(one(base_ring(a)))

gen(a::FqNmodPolyRing) = a([zero(base_ring(a)), one(base_ring(a))])

iszero(x::fq_nmod_poly) = ccall((:fq_nmod_poly_is_zero, :libflint), Bool,
                              (Ptr{fq_nmod_poly}, Ptr{FqNmodFiniteField}),
                              &x, &base_ring(x.parent))

isone(x::fq_nmod_poly) = ccall((:fq_nmod_poly_is_one, :libflint), Bool,
                              (Ptr{fq_nmod_poly}, Ptr{FqNmodFiniteField}),
                              &x, &base_ring(x.parent))

isgen(x::fq_nmod_poly) = ccall((:fq_nmod_poly_is_gen, :libflint), Bool,
                              (Ptr{fq_nmod_poly}, Ptr{FqNmodFiniteField}),
                              &x, &base_ring(x.parent))

degree(f::fq_nmod_poly) = f.length - 1

function deepcopy(a::fq_nmod_poly)
   z = fq_nmod_poly(a)
   z.parent = a.parent
   return z
end

################################################################################
#
#   Canonicalisation
#
################################################################################

canonical_unit(a::fq_nmod_poly) = canonical_unit(lead(a))
  
################################################################################
#
#  AbstractString{} I/O
#
################################################################################

function show(io::IO, x::fq_nmod_poly)
   if length(x) == 0
      print(io, "0")
   else
      cstr = ccall((:fq_nmod_poly_get_str_pretty, :libflint), Ptr{UInt8}, 
                  (Ptr{fq_nmod_poly}, Ptr{UInt8}, Ptr{FqNmodFiniteField}),
                  &x, bytestring(string(var(parent(x)))),
                  &((x.parent).base_ring))
      print(io, bytestring(cstr))
      ccall((:flint_free, :libflint), Void, (Ptr{UInt8},), cstr)
   end
end

function show(io::IO, R::FqNmodPolyRing)
  print(io, "Univariate Polynomial Ring in ")
  print(io, string(var(R)))
  print(io, " over ")
  show(io, base_ring(R))
end

show_minus_one(::Type{fq_nmod_poly}) = show_minus_one(fq_nmod)

################################################################################
#
#  Unary operations
#
################################################################################

function -(x::fq_nmod_poly)
   z = parent(x)()
   ccall((:fq_nmod_poly_neg, :libflint), Void,
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Ptr{FqNmodFiniteField}),
         &z, &x, &base_ring(parent(x)))
   return z
end

################################################################################
#
#  Binary operations
#
################################################################################

function +(x::fq_nmod_poly, y::fq_nmod_poly)
   check_parent(x,y)
   z = parent(x)()
   ccall((:fq_nmod_poly_add, :libflint), Void, 
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly},
         Ptr{fq_nmod_poly}, Ptr{FqNmodFiniteField}),
         &z, &x, &y, &base_ring(parent(x)))
   return z
end

function -(x::fq_nmod_poly, y::fq_nmod_poly)
   check_parent(x,y)
   z = parent(x)()
   ccall((:fq_nmod_poly_sub, :libflint), Void,
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly},
         Ptr{fq_nmod_poly}, Ptr{FqNmodFiniteField}),
         &z, &x, &y, &base_ring(parent(x)))
   return z
end

function *(x::fq_nmod_poly, y::fq_nmod_poly)
   check_parent(x,y)
   z = parent(x)()
   ccall((:fq_nmod_poly_mul, :libflint), Void,
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly},
         Ptr{fq_nmod_poly}, Ptr{FqNmodFiniteField}),
         &z, &x, &y, &base_ring(parent(x)))
   return z
end

################################################################################
#
#   Ad hoc binary operators
#
################################################################################

function *(x::fq_nmod, y::fq_nmod_poly)
   parent(x) != base_ring(parent(y)) &&
         error("Coefficient rings must be equal")
   z = parent(y)()
   ccall((:fq_nmod_poly_scalar_mul_fq_nmod, :libflint), Void,
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly},
         Ptr{fq_nmod}, Ptr{FqNmodFiniteField}),
         &z, &y, &x, &parent(x))
  return z
end

*(x::fq_nmod_poly, y::fq_nmod) = y*x

*(x::fmpz, y::fq_nmod_poly) = base_ring(parent(y))(x) * y

*(x::fq_nmod_poly, y::fmpz) = y*x

*(x::Integer, y::fq_nmod_poly) = fmpz(x)*y

*(x::fq_nmod_poly, y::Integer) = y*x

+(x::fq_nmod, y::fq_nmod_poly) = parent(y)(x) + y

+(x::fq_nmod_poly, y::fq_nmod) = y + x

+(x::fmpz, y::fq_nmod_poly) = base_ring(parent(y))(x) + y

+(x::fq_nmod_poly, y::fmpz) = y + x

+(x::fq_nmod_poly, y::Integer) = x + fmpz(y)

+(x::Integer, y::fq_nmod_poly) = y + x

-(x::fq_nmod, y::fq_nmod_poly) = parent(y)(x) - y

-(x::fq_nmod_poly, y::fq_nmod) = x - parent(x)(y)

-(x::fmpz, y::fq_nmod_poly) = base_ring(parent(y))(x) - y

-(x::fq_nmod_poly, y::fmpz) = x - base_ring(parent(x))(y)

-(x::fq_nmod_poly, y::Integer) = x - fmpz(y)

-(x::Integer, y::fq_nmod_poly) = fmpz(x) - y

################################################################################
#
#   Powering
#
################################################################################

function ^(x::fq_nmod_poly, y::Int)
   y < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fq_nmod_poly_pow, :libflint), Void,
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Int, Ptr{FqNmodFiniteField}), 
         &z, &x, y, &base_ring(parent(x)))
   return z
end

################################################################################
#
#   Comparisons
#
################################################################################

function ==(x::fq_nmod_poly, y::fq_nmod_poly)
   check_parent(x,y)
   r = ccall((:fq_nmod_poly_equal, :libflint), Cint,
             (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Ptr{FqNmodFiniteField}),
             &x, &y, &base_ring(parent(x)))
   return Bool(r)
end

################################################################################
#
#   Ad hoc comparisons
#
################################################################################

function ==(x::fq_nmod_poly, y::fq_nmod) 
   base_ring(parent(x)) != parent(y) && return false
   if length(x) > 1
      return false
   elseif length(x) == 1 
      r = ccall((:fq_nmod_poly_equal_fq_nmod, :libflint), Cint, 
                (Ptr{fq_nmod_poly}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}),
                &x, &y, &base_ring(parent(x)))
      return Bool(r)
   else
      return y == 0
  end 
end

==(x::fq_nmod, y::fq_nmod_poly) = y == x

==(x::fq_nmod_poly, y::fmpz) = x == base_ring(parent(x))(y)

==(x::fmpz, y::fq_nmod_poly) = y == x

==(x::fq_nmod_poly, y::Integer) = x == fmpz(y)

==(x::Integer, y::fq_nmod_poly) = y == x

################################################################################
#
#   Truncation
#
################################################################################

function truncate(x::fq_nmod_poly, n::Int)
   n < 0 && throw(DomainError())
   if length(x) <= n
      return x
   end
   z = parent(x)()
   ccall((:fq_nmod_poly_set_trunc, :libflint), Void,
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Int, Ptr{FqNmodFiniteField}),
         &z, &x, n, &base_ring(parent(x)))
   return z
end

function mullow(x::fq_nmod_poly, y::fq_nmod_poly, n::Int)
   check_parent(x,y)
   n < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fq_nmod_poly_mullow, :libflint), Void,
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly},
         Int, Ptr{FqNmodFiniteField}),
         &z, &x, &y, n, &base_ring(parent(x)))
   return z
end

################################################################################
#
#   Reversal
#
################################################################################

function reverse(x::fq_nmod_poly, len::Int)
   len < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fq_nmod_poly_reverse, :libflint), Void,
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Int, Ptr{FqNmodFiniteField}),
         &z, &x, len, &base_ring(parent(x)))
   return z
end

################################################################################
#
#   Shifting
#
################################################################################

function shift_left(x::fq_nmod_poly, len::Int)
   len < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fq_nmod_poly_shift_left, :libflint), Void,
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Int, Ptr{FqNmodFiniteField}),
         &z, &x, len, &base_ring(parent(x)))
   return z
end

function shift_right(x::fq_nmod_poly, len::Int)
   len < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fq_nmod_poly_shift_right, :libflint), Void,
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Int, Ptr{FqNmodFiniteField}),
         &z, &x, len, &base_ring(parent(x)))
   return z
end

################################################################################
#
#   Euclidean division
#
################################################################################

function div(x::fq_nmod_poly, y::fq_nmod_poly)
   check_parent(x,y)
   z = parent(x)()
   ccall((:fq_nmod_poly_div_basecase, :libflint), Void,
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly},
         Ptr{FqNmodFiniteField}), &z, &x, &y, &base_ring(parent(x)))
  return z
end

function rem(x::fq_nmod_poly, y::fq_nmod_poly)
   check_parent(x,y)
   z = parent(x)()
   ccall((:fq_nmod_poly_rem, :libflint), Void,
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly},
         Ptr{FqNmodFiniteField}), &z, &x, &y, &base_ring(parent(x)))
  return z
end

mod(x::fq_nmod_poly, y::fq_nmod_poly) = rem(x, y)

function divrem(x::fq_nmod_poly, y::fq_nmod_poly)
   check_parent(x,y)
   z = parent(x)()
   r = parent(x)()
   ccall((:fq_nmod_poly_divrem, :libflint), Void, (Ptr{fq_nmod_poly},
         Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly},
         Ptr{FqNmodFiniteField}), &z, &r, &x, &y, &base_ring(parent(x)))
   return z,r
end

################################################################################
#
#   Modular arithmetic
#
################################################################################

function powmod(x::fq_nmod_poly, n::Int, y::fq_nmod_poly)
   check_parent(x,y)
   z = parent(x)()
   ccall((:fq_nmod_poly_powmod_ui_binexp, :libflint), Void,
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Int, Ptr{fq_nmod_poly},
         Ptr{FqNmodFiniteField}), &z, &x, n, &y, &base_ring(parent(x)))
  return z
end

################################################################################
#
#   GCD
#
################################################################################

function gcd(x::fq_nmod_poly, y::fq_nmod_poly)
   check_parent(x,y)
   z = parent(x)()
   ccall((:fq_nmod_poly_gcd, :libflint), Void,
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly},
         Ptr{FqNmodFiniteField}), &z, &x, &y, &base_ring(parent(x)))
   return z
end

function gcdinv(x::fq_nmod_poly, y::fq_nmod_poly)
   check_parent(x,y)
   z = parent(x)()
   s = parent(x)()
   t = parent(x)()
   ccall((:fq_nmod_poly_xgcd, :libflint), Void,
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, 
          Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly},
           Ptr{FqNmodFiniteField}), &z, &s, &t, &x, &y, &base_ring(parent(x)))
   return z, s
end

function gcdx(x::fq_nmod_poly, y::fq_nmod_poly)
   check_parent(x,y)
   z = parent(x)()
   s = parent(x)()
   t = parent(x)()
   ccall((:fq_nmod_poly_xgcd, :libflint), Void,
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, 
          Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly},
           Ptr{FqNmodFiniteField}), &z, &s, &t, &x, &y, &base_ring(parent(x)))
   return z, s, t
end

################################################################################
#
#   Evaluation
#
################################################################################

function evaluate(x::fq_nmod_poly, y::fq_nmod)
   base_ring(parent(x)) != parent(y) && error("Incompatible coefficient rings")
   z = parent(y)()
   ccall((:fq_nmod_poly_evaluate_fq_nmod, :libflint), Void,
         (Ptr{fq_nmod}, Ptr{fq_nmod_poly}, Ptr{fq_nmod},
         Ptr{FqNmodFiniteField}), &z, &x, &y, &base_ring(parent(x)))
   return z
end

################################################################################
#
#   Composition
#
################################################################################

function compose(x::fq_nmod_poly, y::fq_nmod_poly)
   check_parent(x,y)
   z = parent(x)()
   ccall((:fq_nmod_poly_compose, :libflint), Void, 
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly},
         Ptr{FqNmodFiniteField}), &z, &x, &y, &base_ring(parent(x)))
   return z
end

################################################################################
#
#   Derivative
#
################################################################################

function derivative(x::fq_nmod_poly)
   z = parent(x)()
   ccall((:fq_nmod_poly_derivative, :libflint), Void, 
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Ptr{FqNmodFiniteField}),
         &z, &x, &base_ring(parent(x)))
   return z
end

################################################################################
#
#  Inflation and deflation
#
################################################################################

function inflate(x::fq_nmod_poly, n::Int)
   z = parent(x)()
   ccall((:fq_nmod_poly_inflate, :libflint), Void, (Ptr{fq_nmod_poly},
         Ptr{fq_nmod_poly}, Culong, Ptr{FqNmodFiniteField}),
         &z, &x, UInt(n), &base_ring(parent(x)))
   return z
end

function deflate(x::fq_nmod_poly, n::Int)
   z = parent(x)()
   ccall((:fq_nmod_poly_deflate, :libflint), Void,
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Culong, Ptr{FqNmodFiniteField}),
         &z, &x, UInt(n), &base_ring(parent(x)))
  return z
end

################################################################################
#
#  Squarefree testing
#
################################################################################

function issquarefree(x::fq_nmod_poly)
   return Bool(ccall((:fq_nmod_poly_is_squarefree, :libflint), Int32, 
       (Ptr{fq_nmod_poly}, Ptr{FqNmodFiniteField}), &x, &base_ring(parent(x))))
end

################################################################################
#
#  Factorization
#
################################################################################

function factor(x::fq_nmod_poly)
   R = parent(x)
   F = base_ring(R)
   a = F()
   fac = fq_nmod_poly_factor(F)
   ccall((:fq_nmod_poly_factor, :libflint), Void, (Ptr{fq_nmod_poly_factor},
         Ptr{fq_nmod}, Ptr{fq_nmod_poly}, Ptr{FqNmodFiniteField}),
         &fac, &a, &x, &F)
   res = Dict{fq_nmod_poly,Int}()
   for i in 1:fac.num
      f = R()
      ccall((:fq_nmod_poly_factor_get_poly, :libflint), Void,
            (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly_factor}, Int,
            Ptr{FqNmodFiniteField}), &f, &fac, i-1, &F)
      e = unsafe_load(fac.exp,i)
      res[f] = e
   end
   return res 
end  

function factor_distinct_deg(x::fq_nmod_poly)
   R = parent(x)
   F = base_ring(R)
   fac = fq_nmod_poly_factor(F)
   tmp = fq_nmod_poly_factor(F)
   ccall((:fq_nmod_poly_factor_fit_length, :libflint), Void,
         (Ptr{fq_nmod_poly_factor}, Int, Ptr{FqNmodFiniteField}),
         &tmp, degree(x), &F)
   ccall((:fq_nmod_poly_factor_distinct_deg, :libflint), Void, 
         (Ptr{fq_nmod_poly_factor}, Ptr{fq_nmod_poly}, Ptr{Int},
         Ptr{FqNmodFiniteField}), &fac, &x, &tmp.exp, &F)
   res = Dict{fq_nmod_poly, Int}()
   for i in 1:fac.num
      f = R()
      ccall((:fq_nmod_poly_factor_get_poly, :libflint), Void,
            (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly_factor}, Int,
            Ptr{FqNmodFiniteField}), &f, &fac, i-1, &F)
      d = unsafe_load(tmp.exp,i)
      res[f] = d
   end
   return res
end

################################################################################
#
#   Unsafe functions
#
################################################################################

function fit!(z::fq_nmod_poly, n::Int)
   ccall((:fq_nmod_poly_fit_length, :libflint), Void, 
         (Ptr{fq_nmod_poly}, Int, Ptr{FqNmodFiniteField}),
         &z, n, &base_ring(parent(z)))
end

function setcoeff!(z::fq_nmod_poly, n::Int, x::fq_nmod)
   ccall((:fq_nmod_poly_set_coeff, :libflint), Void, 
         (Ptr{fq_nmod_poly}, Int, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}),
         &z, n, &x, &base_ring(parent(z)))
end

function mul!(z::fq_nmod_poly, x::fq_nmod_poly, y::fq_nmod_poly)
   ccall((:fq_nmod_poly_mul, :libflint), Void, 
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly},
         Ptr{FqNmodFiniteField}), &z, &x, &y, &base_ring(parent(x)))
end

function add!(z::fq_nmod_poly, x::fq_nmod_poly, y::fq_nmod_poly)
   ccall((:fq_nmod_poly_add, :libflint), Void, 
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly},
         Ptr{FqNmodFiniteField}), &z, &x, &y, &base_ring(parent(x)))
end

function sub!(z::fq_nmod_poly, x::fq_nmod_poly, y::fq_nmod_poly)
   ccall((:fq_nmod_poly_sub, :libflint), Void, 
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly},
         Ptr{FqNmodFiniteField}), &z, &x, &y, &base_ring(parent(x)))
end


function addeq!(z::fq_nmod_poly, x::fq_nmod_poly)
   ccall((:fq_nmod_poly_add, :libflint), Void, 
         (Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly}, Ptr{fq_nmod_poly},
         Ptr{FqNmodFiniteField}), &z, &z, &x, &base_ring(parent(x)))
end

################################################################################
#
#  Promotion rules
#
################################################################################

Base.promote_rule{V <: Integer}(::Type{fq_nmod_poly}, ::Type{V}) = fq_nmod_poly

Base.promote_rule(::Type{fq_nmod_poly}, ::Type{fmpz}) = fq_nmod_poly

Base.promote_rule(::Type{fq_nmod_poly}, ::Type{fq_nmod}) = fq_nmod_poly

###############################################################################
#
#   Polynomial substitution
#
###############################################################################

function Base.call(f::fq_nmod_poly, a::fq_nmod)
   if parent(a) != base_ring(f)
      return subst(f, a)
   end
   return evaluate(f, a)
end

################################################################################
#
#   Parent object call overloads
#
################################################################################

function Base.call(R::FqNmodPolyRing)
   z = fq_nmod_poly()
   z.parent = R
   return z
end

function Base.call(R::FqNmodPolyRing, x::fq_nmod)
  z = fq_nmod_poly(x)
  z.parent = R
  return z
end

function Base.call(R::FqNmodPolyRing, x::fmpz)
   return R(base_ring(R)(x))
end

function Base.call(R::FqNmodPolyRing, x::Integer)
   return R(fmpz(x))
end

function Base.call(R::FqNmodPolyRing, x::Array{fq_nmod, 1})
   length(x) == 0 && error("Array must be non-empty")
   base_ring(R) != parent(x[1]) && error("Coefficient rings must coincide")
   z = fq_nmod_poly(x)
   z.parent = R
   return z
end

function Base.call(R::FqNmodPolyRing, x::Array{fmpz, 1})
   length(x) == 0 && error("Array must be non-empty")
   z = fq_nmod_poly(x, base_ring(R))
   z.parent = R
   return z
end

function Base.call{T <: Integer}(R::FqNmodPolyRing, x::Array{T, 1})
   length(x) == 0 && error("Array must be non-empty")
   return R(map(fmpz, x))
end

function Base.call(R::FqNmodPolyRing, x::fmpz_poly)
   z = fq_nmod_poly(x, base_ring(R))
   z.parent = R
   return z
end

function Base.call(R::FqNmodPolyRing, x::fq_nmod_poly)
  parent(x) != R && error("Unable to coerce to polynomial")
  return x
end

################################################################################
#
#   PolynomialRing constructor
#
################################################################################

function PolynomialRing(R::FqNmodFiniteField, s::AbstractString{})
   S = symbol(s)
   parent_obj = FqNmodPolyRing(R, S)
   return parent_obj, parent_obj([R(0), R(1)])
end

