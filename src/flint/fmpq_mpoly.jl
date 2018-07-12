###############################################################################
#
#   fmpq_mpoly.jl : Flint multivariate polynomials over fmpz
#
###############################################################################

export FmpqMPolyRing, fmpq_mpoly, degrees, symbols

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{fmpq_mpoly}) = FmpqMPolyRing

elem_type(::Type{FmpqMPolyRing}) = fmpq_mpoly

elem_type(::FmpqMPolyRing) = fmpq_mpoly

symbols(a::FmpqMPolyRing) = a.S

parent(a::fmpq_mpoly) = a.parent

function check_parent(a::fmpq_mpoly, b::fmpq_mpoly)
   parent(a) != parent(b) &&
      error("Incompatible polynomial rings in polynomial operation")
end

nvars(a::FmpqMPolyRing) = ccall((:fmpq_mpoly_ctx_nvars, :libflint), Int,
                                (Ref{FmpqMPolyRing}, ), a)

base_ring(a::FmpqMPolyRing) = a.base_ring

function ordering(a::FmpqMPolyRing)
   b = ccall((:fmpq_mpoly_ctx_ord, :libflint), Cint, (Ref{FmpqMPolyRing}, ), a)
   return flint_orderings[b + 1]
end

function gens(R::FmpqMPolyRing)
   A = Array{fmpq_mpoly}(R.nvars)
   for i = 1:R.nvars
      z = R()
      ccall((:fmpq_mpoly_gen, :libflint), Void,
            (Ref{fmpq_mpoly}, Int, Ref{FmpqMPolyRing}), z, i - 1, R)
      A[i] = z
   end
   return A
end

function gen(R::FmpqMPolyRing, i::Int)
   n = nvars(R)
   (i <= 0 || i > n) && error("Index must be between 1 and $n")
   z = R()
   ccall((:fmpq_mpoly_gen, :libflint), Void,
         (Ref{fmpq_mpoly}, Int, Ref{FmpqMPolyRing}), z, i - 1, R)
   return z
end

function isgen(a::fmpq_mpoly, i::Int)
   n = nvars(parent(a))
   (i <= 0 || i > n) && error("Index must be between 1 and $n")
   R = parent(a)
   return Bool(ccall((:fmpq_mpoly_is_gen, :libflint), Cint,
                     (Ref{fmpq_mpoly}, Int, Ref{FmpqMPolyRing}),
                     a, i - 1, a.parent))
end

function isgen(a::fmpq_mpoly)
   n = nvars(parent(a))
   for i in 1:n
      isgen(a, i) && return true
   end
   return false
end

function deepcopy_internal(a::fmpq_mpoly, dict::ObjectIdDict)
   z = parent(a)()
   ccall((:fmpq_mpoly_set, :libflint), Void,
         (Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}),
         z, a, a.parent)
   return z
end

function length(a::fmpq_mpoly)
   n = ccall((:fmpq_mpoly_length, :libflint), Int, (Ref{fmpq_mpoly}, ), a)
   return n
end

function one(R::FmpqMPolyRing)
   z = R()
   ccall((:fmpq_mpoly_one, :libflint), Void,
         (Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}), z, R)
   return z
end

function zero(R::FmpqMPolyRing)
   z = R()
   ccall((:fmpq_mpoly_zero, :libflint), Void,
         (Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}), z, R)
   return z
end

function isone(a::fmpq_mpoly)
   b = ccall((:fmpq_mpoly_is_one, :libflint), Cint,
             (Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}), a, a.parent)
   return Bool(b)
end

function iszero(a::fmpq_mpoly)
   b = ccall((:fmpq_mpoly_is_zero, :libflint), Cint,
             (Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}), a, a.parent)
   return Bool(b)
end

function ismonomial(a::fmpq_mpoly)
   return length(a) == 1
end

function isconstant(a::fmpq_mpoly)
   b = ccall((:fmpq_mpoly_is_fmpq, :libflint), Cint,
             (Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}), a, parent(a))
   return Bool(b)
end

################################################################################
#
#  Getting coefficients
#
################################################################################

function coeff(a::fmpq_mpoly, i::Int)
   z = fmpq()
   n = length(a)
   (i < 1 || i > n) && error("Index must be between 1 and $(length(a))")
   ccall((:fmpq_mpoly_get_termcoeff_fmpq, :libflint), Void,
         (Ref{fmpq}, Ref{fmpq_mpoly}, Int, Ref{FmpqMPolyRing}),
         z, a, i - 1, a.parent)
   return z
end

function coeff(a::fmpq_mpoly, b::fmpq_mpoly)
   check_parent(a, b)
   !isone(length(b)) && error("Second argument must be a monomial")
   z = fmpq()
   ccall((:fmpq_mpoly_get_coeff_fmpq_monomial, :libflint), Void,
         (Ref{fmpq}, Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}),
         z, a, b, parent(a))
   return z
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function degree_int(a::fmpq_mpoly, i::Int)
   n = nvars(parent(a))
   (i <= 0 || i > n) && error("Index must be between 1 and $n")
   d = ccall((:fmpq_mpoly_degree_si, :libflint), Int,
             (Ref{fmpq_mpoly}, Int, Ref{FmpqMPolyRing}), a, i - 1, a.parent)
   return d
end

function degree(a::fmpq_mpoly, i::Int)
   n = nvars(parent(a))
   (i <= 0 || i > n) && error("Index must be between 1 and $n")
   d = fmpz()
   ccall((:fmpq_mpoly_degree_fmpz, :libflint), Void,
         (Ref{fmpz}, Ref{fmpq_mpoly}, Int, Ref{FmpqMPolyRing}),
         d, a, i - 1, a.parent)
   return d
end

function degrees_fit_int(a::fmpq_mpoly)
   b = ccall((:fmpq_mpoly_degrees_fit_si, :libflint), Cint,
             (Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}), a, a.parent)
   return Bool(b)
end

function degrees_int(a::fmpq_mpoly)
   degs = Vector{Int}(nvars(parent(a)))
   ccall((:fmpq_mpoly_degrees_si, :libflint), Void,
         (Ptr{Int}, Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}),
         degs, a, a.parent)
   return degs
end

function degrees(a::fmpq_mpoly)
   n = nvars(parent(a))
   degs = Vector{fmpz}(n)
   for i in 1:n
      degs[i] = fmpz()
   end
   ccall((:fmpq_mpoly_degrees_fmpz, :libflint), Void,
         (Ptr{Ref{fmpz}}, Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}),
         degs, a, a.parent)
   return degs
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, x::fmpq_mpoly)
   if length(x) == 0
      print(io, "0")
   else
      cstr = ccall((:fmpq_mpoly_get_str_pretty, :libflint), Ptr{UInt8},
          (Ref{fmpq_mpoly}, Ptr{Ptr{UInt8}}, Ref{FmpqMPolyRing}),
          x, [string(s) for s in symbols(parent(x))], x.parent)
      print(io, unsafe_string(cstr))

      ccall((:flint_free, :libflint), Void, (Ptr{UInt8},), cstr)
   end
end

function show(io::IO, p::FmpqMPolyRing)
   const max_vars = 5 # largest number of variables to print
   n = nvars(p)
   print(io, "Multivariate Polynomial Ring in ")
   if n > max_vars
      print(io, nvars(p))
      print(io, " variables ")
   end
   for i = 1:min(n - 1, max_vars - 1)
      print(io, string(p.S[i]), ", ")
   end
   if n > max_vars
      print(io, "..., ")
   end
   print(io, string(p.S[n]))
   print(io, " over ")
   show(io, base_ring(p))
end

###############################################################################
#
#   Basic arithmetic
#
###############################################################################

function -(a::fmpq_mpoly)
   z = parent(a)()
   ccall((:fmpq_mpoly_neg, :libflint), Void,
       (Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}),
       z, a, a.parent)
   return z
end

function +(a::fmpq_mpoly, b::fmpq_mpoly)
   check_parent(a, b)
   z = parent(a)()
   ccall((:fmpq_mpoly_add, :libflint), Void,
       (Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}),
       z, a, b, a.parent)
   return z
end

function -(a::fmpq_mpoly, b::fmpq_mpoly)
   check_parent(a, b)
   z = parent(a)()
   ccall((:fmpq_mpoly_sub, :libflint), Void,
       (Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}),
       z, a, b, a.parent)
   return z
end

function *(a::fmpq_mpoly, b::fmpq_mpoly)
   check_parent(a, b)
   z = parent(a)()
   ccall((:fmpq_mpoly_mul, :libflint), Void,
       (Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}),
       z, a, b, a.parent)
   return z
end

###############################################################################
#
#   Ad hoc arithmetic
#
###############################################################################

for (jT, cN, cT) in ((fmpq, :fmpq, Ref{fmpq}), (fmpz, :fmpz, Ref{fmpz}),
                     (Int, :si, Int))
   @eval begin
      function +(a::fmpq_mpoly, b::($jT))
         z = parent(a)()
         ccall(($(string(:fmpq_mpoly_add_, cN)), :libflint), Void,
               (Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, ($cT), Ref{FmpqMPolyRing}),
               z, a, b, parent(a))
         return z
      end

      +(a::($jT), b::fmpq_mpoly) = b + a

      function -(a::fmpq_mpoly, b::($jT))
         z = parent(a)()
         ccall(($(string(:fmpq_mpoly_sub_, cN)), :libflint), Void,
               (Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, ($cT), Ref{FmpqMPolyRing}),
               z, a, b, parent(a))
         return z
      end

      -(a::($jT), b::fmpq_mpoly) = - (b - a)

      function *(a::fmpq_mpoly, b::($jT))
         z = parent(a)()
         ccall(($(string(:fmpq_mpoly_scalar_mul_, cN)), :libflint), Void,
               (Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, ($cT), Ref{FmpqMPolyRing}),
               z, a, b, parent(a))
         return z
      end

      *(a::($jT), b::fmpq_mpoly) = b * a

      function divexact(a::fmpq_mpoly, b::($jT))
         z = parent(a)()
         ccall(($(string(:fmpq_mpoly_scalar_div_, cN)), :libflint), Void,
               (Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, ($cT), Ref{FmpqMPolyRing}),
               z, a, b, parent(a))
         return z
      end

      //(a::fmpq_mpoly, b::($jT)) = divexact(a, b)
   end
end

+(a::fmpq_mpoly, b::Integer) = a + fmpz(b)

+(a::Integer, b::fmpq_mpoly) = b + a

+(a::fmpq_mpoly, b::Rational{<:Integer}) = a + fmpq(b)

+(a::Rational{<:Integer}, b::fmpq_mpoly) = b + a

*(a::fmpq_mpoly, b::Integer) = a * fmpz(b)

*(a::Integer, b::fmpq_mpoly) = b * a

*(a::fmpq_mpoly, b::Rational{<:Integer}) = a * fmpq(b)

*(a::Rational{<:Integer}, b::fmpq_mpoly) = b * a

divexact(a::fmpq_mpoly, b::Integer) = divexact(a, fmpz(b))

divexact(a::fmpq_mpoly, b::Rational{<:Integer}) = divexact(a, fmpq(b))

//(a::fmpq_mpoly, b::Integer) = //(a, fmpz(b))

//(a::fmpq_mpoly, b::Rational{<:Integer}) = //(a, fmpq(b))

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::fmpq_mpoly, b::Int)
   b < 0 && throw(DomainError())
   z = parent(a)()
   ccall((:fmpq_mpoly_pow_si, :libflint), Void,
         (Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, Int, Ref{FmpqMPolyRing}),
         z, a, b, parent(a))
   return z
end

function ^(a::fmpq_mpoly, b::fmpz)
   b < 0 && throw(DomainError())
   z = parent(a)()
   ccall((:fmpq_mpoly_pow_fmpz, :libflint), Void,
         (Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, Ref{fmpz}, Ref{FmpqMPolyRing}),
         z, a, b, parent(a))
   return z
end

################################################################################
#
#   GCD
#
################################################################################

function gcd(a::fmpq_mpoly, b::fmpq_mpoly)
   check_parent(a, b)
   z = parent(a)()
   ccall((:fmpq_mpoly_gcd, :libflint), Cint,
         (Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}),
         z, a, b, a.parent)
   return z
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::fmpq_mpoly, b::fmpq_mpoly)
   check_parent(a, b)
   return Bool(ccall((:fmpq_mpoly_equal, :libflint), Cint,
               (Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}),
               a, b, a.parent))
end

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

function ==(a::fmpq_mpoly, b::fmpq)
   return Bool(ccall((:fmpq_mpoly_equal_fmpq, :libflint), Cint,
                     (Ref{fmpq_mpoly}, Ref{fmpq}, Ref{FmpqMPolyRing}),
                     a, b, a.parent))
end

==(a::fmpq, b::fmpq_mpoly) = b == a

function ==(a::fmpq_mpoly, b::fmpz)
   return Bool(ccall((:fmpq_mpoly_equal_fmpz, :libflint), Cint,
                     (Ref{fmpq_mpoly}, Ref{fmpz}, Ref{FmpqMPolyRing}),
                     a, b, a.parent))
end

function ==(a::fmpq_mpoly, b::Int)
   return Bool(ccall((:fmpq_mpoly_equal_si, :libflint), Cint,
               (Ref{fmpq_mpoly}, Int, Ref{FmpqMPolyRing}),
               a, b, a.parent))
end

==(a::Int, b::fmpq_mpoly) = b == a

==(a::fmpq_mpoly, b::Integer) = a == fmpz(b)

==(a::Integer, b::fmpq_mpoly) = b == a

==(a::fmpq_mpoly, b::Rational{<:Integer}) = a == fmpq(b)

==(a::Rational{<:Integer}, b::fmpq_mpoly) = b == a

###############################################################################
#
#   Divisibility
#
###############################################################################

function divides(a::fmpq_mpoly, b::fmpq_mpoly)
   check_parent(a, b)
   iszero(b) && error("Cannot divide by zero")
   z = parent(a)()
   d = ccall((:fmpq_mpoly_divides, :libflint), Cint,
       (Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}),
       z, a, b, a.parent)
   return isone(d), z
end

###############################################################################
#
#   Division with remainder
#
###############################################################################

function div(a::fmpq_mpoly, b::fmpq_mpoly)
   check_parent(a, b)
   q = parent(a)()
   ccall((:fmpq_mpoly_div, :libflint), Void,
       (Ref{fmpq_mpoly}, Ref{fmpq_mpoly},
        Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}),
       q, a, b, a.parent)
   return q
end

function divrem(a::fmpq_mpoly, b::fmpq_mpoly)
   check_parent(a, b)
   q = parent(a)()
   r = parent(a)()
   ccall((:fmpq_mpoly_divrem, :libflint), Void,
       (Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, Ref{fmpq_mpoly},
        Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}),
       q, r, a, b, a.parent)
   return q, r
end

function divrem(a::fmpq_mpoly, b::Array{fmpq_mpoly, 1})
   len = length(b)
   q = [parent(a)() for i in 1:len]
   r = parent(a)()
   ccall((:fmpq_mpoly_divrem_ideal, :libflint), Void,
         (Ptr{Ref{fmpq_mpoly}}, Ref{fmpq_mpoly}, Ref{fmpq_mpoly},
          Ptr{Ref{fmpq_mpoly}}, Int, Ref{FmpqMPolyRing}),
       q, r, a, b, len, a.parent)
   return q, r
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::fmpq_mpoly, b::fmpq_mpoly)
   check_parent(a, b)
   b, q = divides(a, b)
   !b && error("Division is not exact in divexact")
   return q
end

###############################################################################
#
#   Calculus
#
###############################################################################

function derivative(a::fmpq_mpoly, i::Int)
   n = nvars(parent(a))
   (i <= 0 || i > n) && error("Index must be between 1 and $n")
   z = parent(a)()
   ccall((:fmpq_mpoly_derivative, :libflint), Void,
         (Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, Int, Ref{FmpqMPolyRing}),
         z, a, i - 1, parent(a))
   return z
end

function integral(a::fmpq_mpoly, i::Int)
   n = nvars(parent(a))
   (i <= 0 || i > n) && error("Index must be between 1 and $n")
   z = parent(a)()
   ccall((:fmpq_mpoly_integral, :libflint), Void,
         (Ref{fmpq_mpoly}, Ref{fmpq_mpoly}, Int, Ref{FmpqMPolyRing}),
         z, a, i - 1, parent(a))
   return z
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function addeq!(a::fmpq_mpoly, b::fmpq_mpoly)
   ccall((:fmpq_mpoly_add, :libflint), Void,
         (Ref{fmpq_mpoly}, Ref{fmpq_mpoly},
          Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}), a, a, b, a.parent)
   return a
end

function mul!(a::fmpq_mpoly, b::fmpq_mpoly, c::fmpq_mpoly)
   ccall((:fmpq_mpoly_mul, :libflint), Void,
         (Ref{fmpq_mpoly}, Ref{fmpq_mpoly},
          Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}), a, b, c, a.parent)
   return a
end

function setcoeff!(a::fmpq_mpoly, i::Int, c::fmpq)
   ccall((:fmpq_mpoly_set_coeff_fmpq, :libflint), Void,
         (Ref{fmpq_mpoly}, Int, Ref{fmpq}, Ref{FmpqMPolyRing}),
         a, i - 1, c, a.parent)
   return a
end

setcoeff!(a::fmpq_mpoly, i::Int, c::fmpz) = setcoeff!(a, i, fmpq(c))

setcoeff!(a::fmpq_mpoly, i::Int, c::Integer) = setcoeff!(a, i, fmpq(c))

setcoeff!(a::fmpq_mpoly, i::Int, c::Rational{<:Integer}) =
         setcoeff!(a, i, fmpq(c))

###############################################################################
#
#   Manipulating terms and monomials
#
###############################################################################

function _termexp_fits_ui(a::fmpq_mpoly, i::Int)
   b = ccall((:fmpq_mpoly_termexp_fits_ui, :libflint), Cint,
             (Ref{fmpq_mpoly}, Int, Ref{FmpqMPolyRing}), a, i - 1, a.parent)
   return Bool(b)
end

function _get_termexp_ui(a::fmpq_mpoly, i::Int)
   z = Vector{UInt}(nvars(parent(a)))
   ccall((:fmpq_mpoly_get_termexp_ui, :libflint), Void,
         (Ptr{UInt}, Ref{fmpq_mpoly}, Int, Ref{FmpqMPolyRing}),
         z, a, i - 1, parent(a))
   return z
end

function _get_termexp_fmpz(a::fmpq_mpoly, i::Int)
   n = nvars(parent(a))
   z = Vector{fmpz}(n)
   for j in 1:n
      z[j] = fmpz()
   end
   ccall((:fmpq_mpoly_get_termexp_fmpz, :libflint), Void,
         (Ptr{Ref{fmpz}}, Ref{fmpq_mpoly}, Int, Ref{FmpqMPolyRing}),
         z, a, i - 1, parent(a))
   return z
end

function _set_termexp_ui!(a::fmpq_mpoly, i::Int, exps::Vector{UInt})
   ccall((:fmpq_mpoly_set_termexp_ui, :libflint), Void,
         (Ref{fmpq_mpoly}, Int, Ptr{UInt}, Ref{FmpqMPolyRing}),
         a, i - 1, exps, parent(a))
   return a
end

function _set_termexp_fmpz!(a::fmpq_mpoly, i::Int, exps::Vector{fmpz})
   ccall((:fmpq_mpoly_set_termexp_fmpz, :libflint), Void,
         (Ref{fmpq_mpoly}, Int, Ptr{Ref{fmpz}}, Ref{FmpqMPolyRing}),
         a, i - 1, exps, parent(a))
   return a
end

function _get_term(a::fmpq_mpoly, exps::Vector{UInt})
   z = fmpq()
   ccall((:fmpq_mpoly_get_term_fmpq_ui, :libflint), Void,
         (Ref{fmpq}, Ref{fmpq_mpoly}, Ptr{UInt}, Ref{FmpqMPolyRing}),
         z, a, exps, parent(a))
   return z
end

function _set_term!(a::fmpq_mpoly, exps::Vector{UInt}, b::fmpq)
   ccall((:fmpq_mpoly_set_term_fmpq_ui, :libflint), Void,
         (Ref{fmpq_mpoly}, Ref{fmpq}, Ptr{UInt}, Ref{FmpqMPolyRing}),
         a, b, exps, parent(a))
   return a
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule{V <: Integer}(::Type{fmpq_mpoly}, ::Type{V}) = fmpq_mpoly

promote_rule{V <: Integer}(::Type{fmpq_mpoly}, ::Type{Rational{V}}) = fmpq_mpoly

promote_rule(::Type{fmpq_mpoly}, ::Type{fmpz}) = fmpq_mpoly

promote_rule(::Type{fmpq_mpoly}, ::Type{fmpq}) = fmpq_mpoly

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (R::FmpqMPolyRing)()
   z = fmpq_mpoly(R)
   return z
end

function (R::FmpqMPolyRing)(b::fmpq)
   z = fmpq_mpoly(R, b)
   return z
end

function (R::FmpqMPolyRing)(b::fmpz)
   z = fmpq_mpoly(R, b)
   return z
end

function (R::FmpqMPolyRing)(b::Int)
   z = fmpq_mpoly(R, b)
   return z
end

function (R::FmpqMPolyRing)(b::UInt)
   z = fmpq_mpoly(R, b)
   return z
end

function (R::FmpqMPolyRing)(b::Integer)
   return R(fmpz(b))
end

function (R::FmpqMPolyRing)(b::Rational{<:Integer})
   return R(fmpq(b))
end

function (R::FmpqMPolyRing)(a::fmpq_mpoly)
   parent(a) != R && error("Unable to coerce polynomial")
   return a
end

function (R::FmpqMPolyRing)(a::Vector{fmpq}, b::Vector{Vector{T}}) where {T <: Union{fmpz, UInt}}
   length(a) != length(b) && error("Coefficient and exponent vector must have the same length")

   for i in 1:length(b)
     length(b[i]) != nvars(R) && error("Exponent vector $i has length $(length(b[i])) (expected $(nvars(R))")
   end

   z = fmpq_mpoly(R, a, b)
   return z
end

function (R::FmpqMPolyRing)(a::Vector{Any}, b::Vector{Vector{T}}) where T
   n = nvars(R)
   length(a) != length(b) && error("Coefficient and exponent vector must have the same length")
   newa = map(FlintQQ, a)
   newb = map(x -> map(FlintZZ, x), b)
   newaa = convert(Vector{fmpq}, newa)
   newbb = convert(Vector{Vector{fmpz}}, newb)

   for i in 1:length(newbb)
      length(newbb[i]) != n && error("Exponent vector $i has length $(length(newbb[i])) (expected $(nvars(R)))")
   end

   return R(newaa, newbb)
end

###############################################################################
#
#   PolynomialRing constructor
#
###############################################################################

function PolynomialRing(R::FlintRationalField, s::Array{String, 1}; cached::Bool = true, ordering::Symbol = :lex)
   U = [Symbol(x) for x in s]
   parent_obj = FmpqMPolyRing(U, ordering, cached)
   return tuple(parent_obj, gens(parent_obj))
end
