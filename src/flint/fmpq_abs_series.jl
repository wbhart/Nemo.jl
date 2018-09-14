###############################################################################
#
#   fmpq_abs_series.jl : Power series over flint fmpq rationals (using fmpq_poly)
#
###############################################################################

export fmpq_abs_series, FmpqAbsSeriesRing, tan, tanh, sin, sinh, asin, asinh,
       atan, atanh, log

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

function O(a::fmpq_abs_series)
   prec = length(a) - 1
   prec < 0 && throw(DomainError("Precision must be non-negative: $prec"))
   z = parent(a)()
   z.prec = prec
   z.parent = parent(a)
   return z
end

elem_type(::Type{FmpqAbsSeriesRing}) = fmpq_abs_series

parent_type(::Type{fmpq_abs_series}) = FmpqAbsSeriesRing

base_ring(R::FmpqAbsSeriesRing) = R.base_ring

var(a::FmpqAbsSeriesRing) = a.S

###############################################################################
#
#   Basic manipulation
#
###############################################################################

max_precision(R::FmpqAbsSeriesRing) = R.prec_max

function normalise(a::fmpq_abs_series, len::Int)
   if len > 0
      c = fmpq()
      ccall((:fmpq_poly_get_coeff_fmpq, :libflint), Nothing,
         (Ref{fmpq}, Ref{fmpq_abs_series}, Int), c, a, len - 1)
   end
   while len > 0 && iszero(c)
      len -= 1
      if len > 0
         ccall((:fmpq_poly_get_coeff_fmpq, :libflint), Nothing,
            (Ref{fmpq}, Ref{fmpq_abs_series}, Int), c, a, len - 1)
      end
   end

   return len
end

function coeff(x::fmpq_abs_series, n::Int)
   if n < 0
      return fmpq(0)
   end
   z = fmpq()
   ccall((:fmpq_poly_get_coeff_fmpq, :libflint), Nothing,
         (Ref{fmpq}, Ref{fmpq_abs_series}, Int), z, x, n)
   return z
end

function length(x::fmpq_abs_series)
   return ccall((:fmpq_poly_length, :libflint), Int, (Ref{fmpq_abs_series},), x)
end

precision(x::fmpq_abs_series) = x.prec

zero(R::FmpqAbsSeriesRing) = R(0)

one(R::FmpqAbsSeriesRing) = R(1)

function gen(R::FmpqAbsSeriesRing)
   z = fmpq_abs_series([fmpq(0), fmpq(1)], 2, max_precision(R))
   z.parent = R
   return z
end

function deepcopy_internal(a::fmpq_abs_series, dict::IdDict)
   z = fmpq_abs_series(a)
   z.prec = a.prec
   z.parent = parent(a)
   return z
end

function isgen(a::fmpq_abs_series)
   return precision(a) == 0 || ccall((:fmpq_poly_is_x, :libflint), Bool,
                            (Ref{fmpq_abs_series},), a)
end

iszero(a::fmpq_abs_series) = length(a) == 0

isunit(a::fmpq_abs_series) = valuation(a) == 0 && isunit(coeff(a, 0))

function isone(a::fmpq_abs_series)
   return precision(a) == 0 || ccall((:fmpq_poly_is_one, :libflint), Bool,
                                (Ref{fmpq_abs_series},), a)
end

# todo: write an fmpq_poly_valuation
function valuation(a::fmpq_abs_series)
   for i = 1:length(a)
      if !iszero(coeff(a, i - 1))
         return i - 1
      end
   end
   return precision(a)
end

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function show(io::IO, a::FmpqAbsSeriesRing)
   print(io, "Univariate power series ring in ", var(a), " over ")
   show(io, base_ring(a))
end

show_minus_one(::Type{fmpq_abs_series}) = show_minus_one(fmpq)

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(x::fmpq_abs_series)
   z = parent(x)()
   ccall((:fmpq_poly_neg, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}),
               z, x)
   z.prec = x.prec
   return z
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::fmpq_abs_series, b::fmpq_abs_series)
   check_parent(a, b)
   lena = length(a)
   lenb = length(b)

   prec = min(a.prec, b.prec)

   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   z = parent(a)()
   z.prec = prec
   ccall((:fmpq_poly_add_series, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, a, b, lenz)
   return z
end

function -(a::fmpq_abs_series, b::fmpq_abs_series)
   check_parent(a, b)
   lena = length(a)
   lenb = length(b)

   prec = min(a.prec, b.prec)

   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   z = parent(a)()
   z.prec = prec
   ccall((:fmpq_poly_sub_series, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, a, b, lenz)
   return z
end

function *(a::fmpq_abs_series, b::fmpq_abs_series)
   check_parent(a, b)
   lena = length(a)
   lenb = length(b)

   aval = valuation(a)
   bval = valuation(b)

   prec = min(a.prec + bval, b.prec + aval)
   prec = min(prec, max_precision(parent(a)))

   lena = min(lena, prec)
   lenb = min(lenb, prec)

   z = parent(a)()
   z.prec = prec

   if lena == 0 || lenb == 0
      return z
   end

   lenz = min(lena + lenb - 1, prec)

   ccall((:fmpq_poly_mullow, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, a, b, lenz)
   return z
end


###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(x::Int, y::fmpq_abs_series)
   z = parent(y)()
   z.prec = y.prec
   ccall((:fmpq_poly_scalar_mul_si, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, y, x)
   return z
end

function *(x::fmpz, y::fmpq_abs_series)
   z = parent(y)()
   z.prec = y.prec
   ccall((:fmpq_poly_scalar_mul_fmpz, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Ref{fmpz}),
               z, y, x)
   return z
end

function *(x::fmpq, y::fmpq_abs_series)
   z = parent(y)()
   z.prec = y.prec
   ccall((:fmpq_poly_scalar_mul_fmpq, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Ref{fmpq}),
               z, y, x)
   return z
end

*(x::fmpq_abs_series, y::Int) = y*x

*(x::fmpq_abs_series, y::fmpz) = y*x

*(x::fmpq_abs_series, y::fmpq) = y*x

*(x::fmpq_abs_series, y::Integer) = x*fmpz(y)

*(x::Integer, y::fmpq_abs_series) = fmpz(x)*y

*(x::fmpq_abs_series, y::Rational) = x*fmpq(y)

*(x::Rational, y::fmpq_abs_series) = fmpq(x)*y

+(x::fmpq_abs_series, y::Rational) = x + fmpq(y)

+(x::Rational, y::fmpq_abs_series) = fmpq(x) + y

-(x::fmpq_abs_series, y::Rational) = x - fmpq(y)

-(x::Rational, y::fmpq_abs_series) = fmpq(x) - y

###############################################################################
#
#   Shifting
#
###############################################################################

function shift_left(x::fmpq_abs_series, len::Int)
   len < 0 && throw(DomainError("Shift must be non-negative: $len"))
   xlen = length(x)
   z = parent(x)()
   z.prec = x.prec + len
   ccall((:fmpq_poly_shift_left, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, x, len)
   return z
end

function shift_right(x::fmpq_abs_series, len::Int)
   len < 0 && throw(DomainError("Shift must be non-negative: $len"))
   xlen = length(x)
   z = parent(x)()
   if len >= xlen
      z.prec = max(0, x.prec - len)
   else
      z.prec = x.prec - len
      ccall((:fmpq_poly_shift_right, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, x, len)
   end
   return z
end

###############################################################################
#
#   Truncation
#
###############################################################################

function truncate(x::fmpq_abs_series, prec::Int)
   prec < 0 && throw(DomainError("Precision must be non-negative: $prec"))
   if x.prec <= prec
      return x
   end
   z = parent(x)()
   z.prec = prec
   ccall((:fmpq_poly_set_trunc, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, x, prec)
   return z
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::fmpq_abs_series, b::Int)
   b < 0 && throw(DomainError("Exponent must be non-negative: $b"))
   # special case powers of x for constructing power series efficiently
   if precision(a) > 0 && isgen(a) && b > 0
      return shift_left(a, b - 1)
   elseif length(a) == 1
      z = parent(a)(coeff(a, 0)^b)
      set_prec!(z, precision(a))
      return z
   elseif b == 0
      z = one(parent(a))
      set_prec!(z, precision(a))
      return z
   else
      bit = ~((~UInt(0)) >> 1)
      while (UInt(bit) & b) == 0
         bit >>= 1
      end
      z = a
      bit >>= 1
      while bit !=0
         z = z*z
         if (UInt(bit) & b) != 0
            z *= a
         end
         bit >>= 1
      end
      return z
   end
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(x::fmpq_abs_series, y::fmpq_abs_series)
   check_parent(x, y)
   prec = min(x.prec, y.prec)

   n = max(length(x), length(y))
   n = min(n, prec)

   return Bool(ccall((:fmpq_poly_equal_trunc, :libflint), Cint,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               x, y, n))
end

function isequal(x::fmpq_abs_series, y::fmpq_abs_series)
   if parent(x) != parent(y)
      return false
   end
   if x.prec != y.prec || length(x) != length(y)
      return false
   end
   return Bool(ccall((:fmpq_poly_equal, :libflint), Cint,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               x, y, length(x)))
end

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

==(x::fmpq_abs_series, y::Rational{T}) where T <: Union{Int, BigInt} = x == fmpq(y)

==(x::fmpq_abs_series, y::Integer) = x == fmpz(y)

==(x::Rational{T}, y::fmpq_abs_series) where T <: Union{Int, BigInt} = y == x

==(x::Integer, y::fmpq_abs_series) = y == x

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(x::fmpq_abs_series, y::fmpq_abs_series)
   check_parent(x, y)
   iszero(y) && throw(DivideError())
   v2 = valuation(y)
   v1 = valuation(x)
   if v2 != 0
      if v1 >= v2
         x = shift_right(x, v2)
         y = shift_right(y, v2)
      end
   end
   !isunit(y) && error("Unable to invert power series")
   prec = min(x.prec, y.prec - v2 + v1)
   z = parent(x)()
   z.prec = prec
   ccall((:fmpq_poly_div_series, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, x, y, prec)
   return z
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(x::fmpq_abs_series, y::Int)
   y == 0 && throw(DivideError())
   z = parent(x)()
   z.prec = x.prec
   ccall((:fmpq_poly_scalar_div_si, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, x, y)
   return z
end

function divexact(x::fmpq_abs_series, y::fmpz)
   iszero(y) && throw(DivideError())
   z = parent(x)()
   z.prec = x.prec
   ccall((:fmpq_poly_scalar_div_fmpz, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Ref{fmpz}),
               z, x, y)
   return z
end

function divexact(x::fmpq_abs_series, y::fmpq)
   iszero(y) && throw(DivideError())
   z = parent(x)()
   z.prec = x.prec
   ccall((:fmpq_poly_scalar_div_fmpq, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Ref{fmpq}),
               z, x, y)
   return z
end

divexact(x::fmpq_abs_series, y::Integer) = divexact(x, fmpz(y))

divexact(x::fmpq_abs_series, y::Rational{T}) where T <: Union{Int, BigInt} = divexact(x, fmpq(y))

###############################################################################
#
#   Inversion
#
###############################################################################

function inv(a::fmpq_abs_series)
  iszero(a) && throw(DivideError())
   !isunit(a) && error("Unable to invert power series")
   ainv = parent(a)()
   ainv.prec = a.prec
   ccall((:fmpq_poly_inv_series, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               ainv, a, a.prec)
   return ainv
end

###############################################################################
#
#   Special functions
#
###############################################################################

function Base.exp(a::fmpq_abs_series)
   !iszero(coeff(a, 0)) && error("Constant term not zero in exp")
   if length(a) == 0 || a.prec == 1
      return parent(a)([fmpq(1)], 1, a.prec)
   end
   z = parent(a)()
   z.prec = a.prec
   ccall((:fmpq_poly_exp_series, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, a, a.prec)
   return z
end

@doc Markdown.doc"""
   log(a::fmpq_abs_series)
> Return log$(a)$. Requires the constant term to be one.
"""
function log(a::fmpq_abs_series)
   !isone(coeff(a, 0)) && error("Constant term not one in log")
   if length(a) == 1 || a.prec < 2
      return parent(a)()
   end
   z = parent(a)()
   z.prec = a.prec
   ccall((:fmpq_poly_log_series, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, a, a.prec)
   return z
end

@doc Markdown.doc"""
   tan(a::fmpq_abs_series)
> Return tan$(a)$. Requires a zero constant term.
"""
function tan(a::fmpq_abs_series)
   !iszero(coeff(a, 0)) && error("Constant term not zero in tan")
   if iszero(a) || a.prec < 2
      return parent(a)()
   end
   z = parent(a)()
   z.prec = a.prec
   ccall((:fmpq_poly_tan_series, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, a, a.prec)
   return z
end

@doc Markdown.doc"""
   tanh(a::fmpq_abs_series)
> Return tanh$(a)$. Requires a zero constant term.
"""
function tanh(a::fmpq_abs_series)
   !iszero(coeff(a, 0)) && error("Constant term not zero in tanh")
   if iszero(a) || a.prec < 2
      return parent(a)()
   end
   z = parent(a)()
   z.prec = a.prec
   ccall((:fmpq_poly_tanh_series, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, a, a.prec)
   return z
end

@doc Markdown.doc"""
   sin(a::fmpq_abs_series)
> Return sin$(a)$. Requires a zero constant term.
"""
function sin(a::fmpq_abs_series)
   !iszero(coeff(a, 0)) && error("Constant term not zero in sin")
   if iszero(a) || a.prec < 2
      return parent(a)()
   end
   z = parent(a)()
   z.prec = a.prec
   ccall((:fmpq_poly_sin_series, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, a, a.prec)
   return z
end

@doc Markdown.doc"""
   sinh(a::fmpq_abs_series)
> Return sinh$(a)$. Requires a zero constant term.
"""
function sinh(a::fmpq_abs_series)
   !iszero(coeff(a, 0)) && error("Constant term not zero in sinh")
   if iszero(a) || a.prec < 2
      return parent(a)()
   end
   z = parent(a)()
   z.prec = a.prec
   ccall((:fmpq_poly_sinh_series, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, a, a.prec)
   return z
end

@doc Markdown.doc"""
   cos(a::fmpq_abs_series)
> Return cos$(a)$. Requires a zero constant term.
"""
function cos(a::fmpq_abs_series)
   !iszero(coeff(a, 0)) && error("Constant term not zero in cos")
   if length(a) == 0 || a.prec == 1
      return one(parent(a))
   end
   z = parent(a)()
   z.prec = a.prec
   ccall((:fmpq_poly_cos_series, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, a, a.prec)
   return z
end

@doc Markdown.doc"""
   cosh(a::fmpq_abs_series)
> Return cosh$(a)$. Requires a zero constant term.
"""
function cosh(a::fmpq_abs_series)
   !iszero(coeff(a, 0)) && error("Constant term not zero in cosh")
   if length(a) == 0 || a.prec == 1
      return one(parent(a))
   end
   z = parent(a)()
   z.prec = a.prec
   ccall((:fmpq_poly_cosh_series, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, a, a.prec)
   return z
end

@doc Markdown.doc"""
   asin(a::fmpq_abs_series)
> Return asin$(a)$. Requires a zero constant term.
"""
function asin(a::fmpq_abs_series)
   !iszero(coeff(a, 0)) && error("Constant term not zero in asin")
   if iszero(a) || a.prec < 2
      return parent(a)()
   end
   z = parent(a)()
   z.prec = a.prec
   ccall((:fmpq_poly_asin_series, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, a, a.prec)
   return z
end

@doc Markdown.doc"""
   asinh(a::fmpq_abs_series)
> Return asinh$(a)$. Requires a zero constant term.
"""
function asinh(a::fmpq_abs_series)
   !iszero(coeff(a, 0)) && error("Constant term not zero in asinh")
   if iszero(a) || a.prec < 2
      return parent(a)()
   end
   z = parent(a)()
   z.prec = a.prec
   ccall((:fmpq_poly_asinh_series, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, a, a.prec)
   return z
end

@doc Markdown.doc"""
   atan(a::fmpq_abs_series)
> Return atan$(a)$. Requires a zero constant term.
"""
function atan(a::fmpq_abs_series)
   !iszero(coeff(a, 0)) && error("Constant term not zero in atan")
   if iszero(a) || a.prec < 2
      return parent(a)()
   end
   z = parent(a)()
   z.prec = a.prec
   ccall((:fmpq_poly_atan_series, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, a, a.prec)
   return z
end

@doc Markdown.doc"""
   atanh(a::fmpq_abs_series)
> Return atanh$(a)$. Requires a zero constant term.
"""
function atanh(a::fmpq_abs_series)
   !iszero(coeff(a, 0)) && error("Constant term not zero in atanh")
   if iszero(a) || a.prec < 2
      return parent(a)()
   end
   z = parent(a)()
   z.prec = a.prec
   ccall((:fmpq_poly_atanh_series, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, a, a.prec)
   return z
end

@doc Markdown.doc"""
   sqrt(a::fmpq_abs_series)
> Return the power series square root of $a$. Requires a constant term equal to
> one.
"""
function Base.sqrt(a::fmpq_abs_series)
   !isone(coeff(a, 0)) && error("Constant term not one in sqrt")
   z = parent(a)()
   z.prec = a.prec
   ccall((:fmpq_poly_sqrt_series, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, a, a.prec)
   return z
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function setcoeff!(z::fmpq_abs_series, n::Int, x::fmpq)
   ccall((:fmpq_poly_set_coeff_fmpq, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Int, Ref{fmpq}),
               z, n, x)
   return z
end

function mul!(z::fmpq_abs_series, a::fmpq_abs_series, b::fmpq_abs_series)
   lena = length(a)
   lenb = length(b)

   aval = valuation(a)
   bval = valuation(b)

   prec = min(a.prec + bval, b.prec + aval)
   prec = min(prec, max_precision(parent(a)))

   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = min(lena + lenb - 1, prec)
   if lenz < 0
      lenz = 0
   end

   z.prec = prec
   ccall((:fmpq_poly_mullow, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, a, b, lenz)
   return z
end

function addeq!(a::fmpq_abs_series, b::fmpq_abs_series)
   lena = length(a)
   lenb = length(b)

   prec = min(a.prec, b.prec)

   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   a.prec = prec
   ccall((:fmpq_poly_add_series, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               a, a, b, lenz)
   return a
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{fmpq_abs_series}, ::Type{T}) where {T <: Integer} = fmpq_abs_series

promote_rule(::Type{fmpq_abs_series}, ::Type{Rational{T}}) where T <: Union{Int, BigInt} = fmpq_abs_series

promote_rule(::Type{fmpq_abs_series}, ::Type{fmpz}) = fmpq_abs_series

promote_rule(::Type{fmpq_abs_series}, ::Type{fmpq}) = fmpq_abs_series

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::FmpqAbsSeriesRing)()
   z = fmpq_abs_series()
   z.prec = a.prec_max
   z.parent = a
   return z
end

function (a::FmpqAbsSeriesRing)(b::Integer)
   if b == 0
      z = fmpq_abs_series()
      z.prec = a.prec_max
   else
      z = fmpq_abs_series([fmpq(b)], 1, a.prec_max)
   end
   z.parent = a
   return z
end

function (a::FmpqAbsSeriesRing)(b::fmpz)
   if iszero(b)
      z = fmpq_abs_series()
      z.prec = a.prec_max
   else
      z = fmpq_abs_series([fmpq(b)], 1, a.prec_max)
   end
   z.parent = a
   return z
end

function (a::FmpqAbsSeriesRing)(b::fmpq)
   if iszero(b)
      z = fmpq_abs_series()
      z.prec = a.prec_max
   else
      z = fmpq_abs_series([b], 1, a.prec_max)
   end
   z.parent = a
   return z
end

(a::FmpqAbsSeriesRing)(b::Rational{T}) where T <: Union{Int, BigInt} = a(fmpq(b))

function (a::FmpqAbsSeriesRing)(b::fmpq_abs_series)
   parent(b) != a && error("Unable to coerce power series")
   return b
end

function (a::FmpqAbsSeriesRing)(b::Array{fmpq, 1}, len::Int, prec::Int)
   z = fmpq_abs_series(b, len, prec)
   z.parent = a
   return z
end

###############################################################################
#
#   PowerSeriesRing constructor
#
###############################################################################

function PowerSeriesRing(R::FlintRationalField, prec::Int, s::AbstractString; model=:capped_relative, cached = true)
   S = Symbol(s)

   if model == :capped_relative
      parent_obj = FmpqRelSeriesRing(prec, S, cached)
   elseif model == :capped_absolute
      parent_obj = FmpqAbsSeriesRing(prec, S, cached)
   else
      error("Unknown model")
   end

   return parent_obj, gen(parent_obj)
end
