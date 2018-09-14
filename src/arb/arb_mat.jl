###############################################################################
#
#   arb_mat.jl : Arb matrices over arb
#
###############################################################################

export rows, cols, zero, one, deepcopy, -, transpose, +, *, &, ==, !=,
       strongequal, overlaps, contains, inv, divexact, charpoly, det,
       lu, lu!, solve, solve!, solve_lu_precomp, solve_lu_precomp!,
       swap_rows, swap_rows!, bound_inf_norm

###############################################################################
#
#   Similar
#
###############################################################################

function similar(x::arb_mat)
   z = arb_mat(rows(x), cols(x))
   z.base_ring = x.base_ring
   return z
end

function similar(x::arb_mat, r::Int, c::Int)
   z = arb_mat(r, c)
   z.base_ring = x.base_ring
   return z
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{arb_mat}) = ArbMatSpace

elem_type(::Type{ArbMatSpace}) = arb_mat

base_ring(a::ArbMatSpace) = a.base_ring

base_ring(a::arb_mat) = a.base_ring

parent(x::arb_mat, cached::Bool = true) =
      MatrixSpace(base_ring(x), rows(x), cols(x))

prec(x::ArbMatSpace) = prec(x.base_ring)

function getindex!(z::arb, x::arb_mat, r::Int, c::Int)
  GC.@preserve x begin
     v = ccall((:arb_mat_entry_ptr, :libarb), Ptr{arb},
                 (Ref{arb_mat}, Int, Int), x, r - 1, c - 1)
     ccall((:arb_set, :libarb), Nothing, (Ref{arb}, Ptr{arb}), z, v)
  end
  return z
end

@inline function getindex(x::arb_mat, r::Int, c::Int)
  @boundscheck Generic._checkbounds(x, r, c)

  z = base_ring(x)()
  GC.@preserve x begin
     v = ccall((:arb_mat_entry_ptr, :libarb), Ptr{arb},
                 (Ref{arb_mat}, Int, Int), x, r - 1, c - 1)
     ccall((:arb_set, :libarb), Nothing, (Ref{arb}, Ptr{arb}), z, v)
  end
  return z
end

for T in [Int, UInt, fmpz, fmpq, Float64, BigFloat, arb, AbstractString]
   @eval begin
      @inline function setindex!(x::arb_mat, y::$T, r::Int, c::Int)
         @boundscheck Generic._checkbounds(x, r, c)

         GC.@preserve x begin
            z = ccall((:arb_mat_entry_ptr, :libarb), Ptr{arb},
                      (Ref{arb_mat}, Int, Int), x, r - 1, c - 1)
            Nemo._arb_set(z, y, prec(base_ring(x)))
         end
      end
   end
end

Base.@propagate_inbounds setindex!(x::arb_mat, y::Integer,
                                 r::Int, c::Int) =
         setindex!(x, fmpz(y), r, c)

Base.@propagate_inbounds setindex!(x::arb_mat, y::Rational{T},
                                 r::Int, c::Int) where {T <: Integer} =
         setindex!(x, fmpz(y), r, c)

zero(a::ArbMatSpace) = a()

function one(x::ArbMatSpace)
  z = x()
  ccall((:arb_mat_one, :libarb), Nothing, (Ref{arb_mat}, ), z)
  return z
end

rows(a::arb_mat) = a.r

cols(a::arb_mat) = a.c

function deepcopy_internal(x::arb_mat, dict::IdDict)
  z = arb_mat(rows(x), cols(x))
  ccall((:arb_mat_set, :libarb), Nothing, (Ref{arb_mat}, Ref{arb_mat}), z, x)
  z.base_ring = x.base_ring
  return z
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, a::ArbMatSpace)
   print(io, "Matrix Space of ")
   print(io, a.rows, " rows and ", a.cols, " columns over ")
   print(io, base_ring(a))
end

function show(io::IO, a::arb_mat)
   r = rows(a)
   c = cols(a)
   if r*c == 0
      print(io, "$r by $c matrix")
   end
   for i = 1:r
      print(io, "[")
      for j = 1:c
         print(io, a[i, j])
         if j != c
            print(io, " ")
         end
      end
      print(io, "]")
      if i != r
         println(io, "")
      end
   end
end

################################################################################
#
#  Unary operations
#
################################################################################

function -(x::arb_mat)
  z = similar(x)
  ccall((:arb_mat_neg, :libarb), Nothing, (Ref{arb_mat}, Ref{arb_mat}), z, x)
  return z
end

################################################################################
#
#  Transpose
#
################################################################################

function transpose(x::arb_mat)
  z = similar(x, cols(x), rows(x))
  ccall((:arb_mat_transpose, :libarb), Nothing,
              (Ref{arb_mat}, Ref{arb_mat}), z, x)
  return z
end

################################################################################
#
#  Binary operations
#
################################################################################

function +(x::arb_mat, y::arb_mat)
  check_parent(x, y)
  z = similar(x)
  ccall((:arb_mat_add, :libarb), Nothing,
              (Ref{arb_mat}, Ref{arb_mat}, Ref{arb_mat}, Int),
              z, x, y, prec(parent(x)))
  return z
end

function -(x::arb_mat, y::arb_mat)
  check_parent(x, y)
  z = similar(x)
  ccall((:arb_mat_sub, :libarb), Nothing,
              (Ref{arb_mat}, Ref{arb_mat}, Ref{arb_mat}, Int),
              z, x, y, prec(parent(x)))
  return z
end

function *(x::arb_mat, y::arb_mat)
  cols(x) != rows(y) && error("Matrices have wrong dimensions")
  z = similar(x, rows(x), cols(y))
  ccall((:arb_mat_mul, :libarb), Nothing,
              (Ref{arb_mat}, Ref{arb_mat}, Ref{arb_mat}, Int),
              z, x, y, prec(base_ring(x)))
  return z
end

################################################################################
#
#   Ad hoc binary operators
#
################################################################################

function ^(x::arb_mat, y::UInt)
  rows(x) != cols(x) && error("Matrix must be square")
  z = similar(x)
  ccall((:arb_mat_pow_ui, :libarb), Nothing,
              (Ref{arb_mat}, Ref{arb_mat}, UInt, Int),
              z, x, y, prec(base_ring(x)))
  return z
end

function *(x::arb_mat, y::Int)
  z = similar(x)
  ccall((:arb_mat_scalar_mul_si, :libarb), Nothing,
              (Ref{arb_mat}, Ref{arb_mat}, Int, Int),
              z, x, y, prec(base_ring(x)))
  return z
end

*(x::Int, y::arb_mat) = y*x

function *(x::arb_mat, y::fmpz)
  z = similar(x)
  ccall((:arb_mat_scalar_mul_fmpz, :libarb), Nothing,
              (Ref{arb_mat}, Ref{arb_mat}, Ref{fmpz}, Int),
              z, x, y, prec(base_ring(x)))
  return z
end

*(x::fmpz, y::arb_mat) = y*x

function *(x::arb_mat, y::arb)
  z = similar(x)
  ccall((:arb_mat_scalar_mul_arb, :libarb), Nothing,
              (Ref{arb_mat}, Ref{arb_mat}, Ref{arb}, Int),
              z, x, y, prec(base_ring(x)))
  return z
end

*(x::arb, y::arb_mat) = y*x

for T in [Integer, fmpz, fmpq, arb]
   @eval begin
      function +(x::arb_mat, y::$T)
         z = deepcopy(x)
         for i = 1:min(rows(x), cols(x))
            z[i, i] += y
         end
         return z
      end

      +(x::$T, y::arb_mat) = y + x

      function -(x::arb_mat, y::$T)
         z = deepcopy(x)
         for i = 1:min(rows(x), cols(x))
            z[i, i] -= y
         end
         return z
      end

      function -(x::$T, y::arb_mat)
         z = -y
         for i = 1:min(rows(y), cols(y))
            z[i, i] += x
         end
         return z
      end
   end
end

function +(x::arb_mat, y::Rational{T}) where T <: Union{Int, BigInt}
   z = deepcopy(x)
   for i = 1:min(rows(x), cols(x))
      z[i, i] += y
   end
   return z
end

+(x::Rational{T}, y::arb_mat) where T <: Union{Int, BigInt} = y + x

function -(x::arb_mat, y::Rational{T}) where T <: Union{Int, BigInt}
   z = deepcopy(x)
   for i = 1:min(rows(x), cols(x))
      z[i, i] -= y
   end
   return z
end

function -(x::Rational{T}, y::arb_mat) where T <: Union{Int, BigInt}
   z = -y
   for i = 1:min(rows(y), cols(y))
      z[i, i] += x
   end
   return z
end

###############################################################################
#
#   Shifting
#
###############################################################################

@doc Markdown.doc"""
    ldexp(x::acb_mat, y::Int)
> Return $2^yx$. Note that $y$ can be positive, zero or negative.
"""
function ldexp(x::arb_mat, y::Int)
  z = similar(x)
  ccall((:arb_mat_scalar_mul_2exp_si, :libarb), Nothing,
              (Ref{arb_mat}, Ref{arb_mat}, Int), z, x, y)
  return z
end

###############################################################################
#
#   Comparisons
#
###############################################################################

@doc Markdown.doc"""
    isequal(x::arb_mat, y::arb_mat)
> Return `true` if the matrices of balls $x$ and $y$ are precisely equal,
> i.e. if all matrix entries have the same midpoints and radii.
"""
function isequal(x::arb_mat, y::arb_mat)
  r = ccall((:arb_mat_equal, :libarb), Cint,
              (Ref{arb_mat}, Ref{arb_mat}), x, y)
  return Bool(r)
end

function ==(x::arb_mat, y::arb_mat)
  check_parent(x, y)
  r = ccall((:arb_mat_eq, :libarb), Cint, (Ref{arb_mat}, Ref{arb_mat}), x, y)
  return Bool(r)
end

function !=(x::arb_mat, y::arb_mat)
  r = ccall((:arb_mat_ne, :libarb), Cint, (Ref{arb_mat}, Ref{arb_mat}), x, y)
  return Bool(r)
end

@doc Markdown.doc"""
    overlaps(x::arb_mat, y::arb_mat)
> Returns `true` if all entries of $x$ overlap with the corresponding entry of
> $y$, otherwise return `false`.
"""
function overlaps(x::arb_mat, y::arb_mat)
  r = ccall((:arb_mat_overlaps, :libarb), Cint,
              (Ref{arb_mat}, Ref{arb_mat}), x, y)
  return Bool(r)
end

@doc Markdown.doc"""
    contains(x::arb_mat, y::arb_mat)
> Returns `true` if all entries of $x$ contain the corresponding entry of
> $y$, otherwise return `false`.
"""
function contains(x::arb_mat, y::arb_mat)
  r = ccall((:arb_mat_contains, :libarb), Cint,
              (Ref{arb_mat}, Ref{arb_mat}), x, y)
  return Bool(r)
end

###############################################################################
#
#   Ad hoc comparisons
#
###############################################################################

@doc Markdown.doc"""
    contains(x::arb_mat, y::fmpz_mat)
> Returns `true` if all entries of $x$ contain the corresponding entry of
> $y$, otherwise return `false`.
"""
function contains(x::arb_mat, y::fmpz_mat)
  r = ccall((:arb_mat_contains_fmpz_mat, :libarb), Cint,
              (Ref{arb_mat}, Ref{fmpz_mat}), x, y)
  return Bool(r)
end


@doc Markdown.doc"""
    contains(x::arb_mat, y::fmpq_mat)
> Returns `true` if all entries of $x$ contain the corresponding entry of
> $y$, otherwise return `false`.
"""
function contains(x::arb_mat, y::fmpq_mat)
  r = ccall((:arb_mat_contains_fmpq_mat, :libarb), Cint,
              (Ref{arb_mat}, Ref{fmpq_mat}), x, y)
  return Bool(r)
end

==(x::arb_mat, y::Integer) = x == parent(x)(y)

==(x::Integer, y::arb_mat) = y == x

==(x::arb_mat, y::fmpz) = x == parent(x)(y)

==(x::fmpz, y::arb_mat) = y == x

==(x::arb_mat, y::fmpz_mat) = x == parent(x)(y)

==(x::fmpz_mat, y::arb_mat) = y == x

###############################################################################
#
#   Inversion
#
###############################################################################

@doc Markdown.doc"""
    inv(M::arb_mat)
> Given a  $n\times n$ matrix of type `arb_mat`, return an
> $n\times n$ matrix $X$ such that $AX$ contains the
> identity matrix. If $A$ cannot be inverted numerically an exception is raised.
"""
function inv(x::arb_mat)
  cols(x) != rows(x) && error("Matrix must be square")
  z = similar(x)
  r = ccall((:arb_mat_inv, :libarb), Cint,
              (Ref{arb_mat}, Ref{arb_mat}, Int), z, x, prec(base_ring(x)))
  Bool(r) ? (return z) : error("Matrix cannot be inverted numerically")
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(x::arb_mat, y::arb_mat)
   cols(x) != cols(y) && error("Incompatible matrix dimensions")
   x*inv(y)
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(x::arb_mat, y::Int)
  y == 0 && throw(DivideError())
  z = similar(x)
  ccall((:arb_mat_scalar_div_si, :libarb), Nothing,
              (Ref{arb_mat}, Ref{arb_mat}, Int, Int),
              z, x, y, prec(base_ring(x)))
  return z
end

function divexact(x::arb_mat, y::fmpz)
  z = similar(x)
  ccall((:arb_mat_scalar_div_fmpz, :libarb), Nothing,
              (Ref{arb_mat}, Ref{arb_mat}, Ref{fmpz}, Int),
              z, x, y, prec(base_ring(x)))
  return z
end

function divexact(x::arb_mat, y::arb)
  z = similar(x)
  ccall((:arb_mat_scalar_div_arb, :libarb), Nothing,
              (Ref{arb_mat}, Ref{arb_mat}, Ref{arb}, Int),
              z, x, y, prec(base_ring(x)))
  return z
end

################################################################################
#
#  Characteristic polynomial
#
################################################################################

function charpoly(x::ArbPolyRing, y::arb_mat)
  base_ring(y) != base_ring(x) && error("Base rings must coincide")
  z = x()
  ccall((:arb_mat_charpoly, :libarb), Nothing,
              (Ref{arb_poly}, Ref{arb_mat}, Int), z, y, prec(base_ring(y)))
  return z
end

###############################################################################
#
#   Determinant
#
###############################################################################

function det(x::arb_mat)
  cols(x) != rows(x) && error("Matrix must be square")
  z = base_ring(x)()
  ccall((:arb_mat_det, :libarb), Nothing,
              (Ref{arb}, Ref{arb_mat}, Int), z, x, prec(base_ring(x)))
  return z
end

################################################################################
#
#   Exponential function
#
################################################################################

@doc Markdown.doc"""
    exp(x::arb_mat)
> Returns the exponential of the matrix $x$.
"""
function Base.exp(x::arb_mat)
  cols(x) != rows(x) && error("Matrix must be square")
  z = similar(x)
  ccall((:arb_mat_exp, :libarb), Nothing,
              (Ref{arb_mat}, Ref{arb_mat}, Int), z, x, prec(base_ring(x)))
  return z
end

###############################################################################
#
#   Linear solving
#
###############################################################################

function lu!(P::Generic.perm, x::arb_mat)
  cols(x) != rows(x) && error("Matrix must be square")
  parent(P).n != rows(x) && error("Permutation does not match matrix")
  P.d .-= 1
  r = ccall((:arb_mat_lu, :libarb), Cint,
              (Ptr{Int}, Ref{arb_mat}, Ref{arb_mat}, Int),
              P.d, x, x, prec(base_ring(x)))
  r == 0 && error("Could not find $(rows(x)) invertible pivot elements")
  P.d .+= 1
  inv!(P)
  return rows(x)
end

function lu(x::arb_mat, P = PermGroup(rows(x)))
  p = P()
  R = base_ring(x)
  L = similar(x)
  U = deepcopy(x)
  n = cols(x)
  r = lu!(p, U)
  for i = 1:n
    for j = 1:n
      if i > j
        L[i, j] = U[i, j]
        U[i, j] = R()
      elseif i == j
        L[i, j] = one(R)
      else
        L[i, j] = R()
      end
    end
  end
  return r, p, L, U
end

function solve!(z::arb_mat, x::arb_mat, y::arb_mat)
  r = ccall((:arb_mat_solve, :libarb), Cint,
              (Ref{arb_mat}, Ref{arb_mat}, Ref{arb_mat}, Int),
              z, x, y, prec(base_ring(x)))
  r == 0 && error("Matrix cannot be inverted numerically")
  nothing
end

function solve(x::arb_mat, y::arb_mat)
  cols(x) != rows(x) && error("First argument must be square")
  cols(x) != rows(y) && error("Matrix dimensions are wrong")
  z = similar(y)
  solve!(z, x, y)
  return z
end

function solve_lu_precomp!(z::arb_mat, P::Generic.perm, LU::arb_mat, y::arb_mat)
  Q = inv(P)
  ccall((:arb_mat_solve_lu_precomp, :libarb), Nothing,
              (Ref{arb_mat}, Ptr{Int}, Ref{arb_mat}, Ref{arb_mat}, Int),
              z, Q.d .- 1, LU, y, prec(base_ring(LU)))
  nothing
end

function solve_lu_precomp(P::Generic.perm, LU::arb_mat, y::arb_mat)
  cols(LU) != rows(y) && error("Matrix dimensions are wrong")
  z = similar(y)
  solve_lu_precomp!(z, P, LU, y)
  return z
end

################################################################################
#
#   Row swapping
#
################################################################################

function swap_rows(x::arb_mat, i::Int, j::Int)
  Generic._checkbounds(rows(x), i) || throw(BoundsError())
  Generic._checkbounds(rows(x), j) || throw(BoundsError())
  z = deepcopy(x)
  swap_rows!(z, i, j)
  return z
end

function swap_rows!(x::arb_mat, i::Int, j::Int)
  ccall((:arb_mat_swap_rows, :libarb), Nothing,
              (Ref{arb_mat}, Ptr{Nothing}, Int, Int),
              x, C_NULL, i - 1, j - 1)
end

################################################################################
#
#   Norm
#
################################################################################

@doc Markdown.doc"""
    bound_inf_norm(x::arb_mat)
> Returns a nonnegative element $z$ of type `arb`, such that $z$ is an upper
> bound for the infinity norm for every matrix in $x$
"""
function bound_inf_norm(x::arb_mat)
  z = arb()
  GC.@preserve x z begin
     t = ccall((:arb_rad_ptr, :libarb), Ptr{mag_struct}, (Ref{arb}, ), z)
     ccall((:arb_mat_bound_inf_norm, :libarb), Nothing,
                 (Ptr{mag_struct}, Ref{arb_mat}), t, x)
     s = ccall((:arb_mid_ptr, :libarb), Ptr{arf_struct}, (Ref{arb}, ), z)
     ccall((:arf_set_mag, :libarb), Nothing,
                 (Ptr{arf_struct}, Ptr{mag_struct}), s, t)
     ccall((:mag_zero, :libarb), Nothing,
                 (Ptr{mag_struct},), t)
  end
  return base_ring(x)(z)
end

################################################################################
#
#   Unsafe functions
#
################################################################################

for (s,f) in (("add!","arb_mat_add"), ("mul!","arb_mat_mul"),
              ("sub!","arb_mat_sub"))
  @eval begin
    function ($(Symbol(s)))(z::arb_mat, x::arb_mat, y::arb_mat)
      ccall(($f, :libarb), Nothing,
                  (Ref{arb_mat}, Ref{arb_mat}, Ref{arb_mat}, Int),
                  z, x, y, prec(base_ring(x)))
      return z
    end
  end
end

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function (x::ArbMatSpace)()
  z = arb_mat(x.rows, x.cols)
  z.base_ring = x.base_ring
  return z
end

function (x::ArbMatSpace)(y::fmpz_mat)
  (x.cols != cols(y) || x.rows != rows(y)) &&
      error("Dimensions are wrong")
  z = arb_mat(y, prec(x))
  z.base_ring = x.base_ring
  return z
end

function (x::ArbMatSpace)(y::Array{T, 2}) where {T <: Union{Int, UInt, fmpz, fmpq, Float64, BigFloat, arb, AbstractString}}
  _check_dim(x.rows, x.cols, y)
  z = arb_mat(x.rows, x.cols, y, prec(x))
  z.base_ring = x.base_ring
  return z
end

function (x::ArbMatSpace)(y::Array{T, 1}) where {T <: Union{Int, UInt, fmpz, fmpq, Float64, BigFloat, arb, AbstractString}}
  _check_dim(x.rows, x.cols, y)
  z = arb_mat(x.rows, x.cols, y, prec(x))
  z.base_ring = x.base_ring
  return z
end

function (x::ArbMatSpace)(y::Union{Int, UInt, fmpz, fmpq, Float64,
                          BigFloat, arb, AbstractString})
  z = x()
  for i in 1:rows(z)
      for j = 1:cols(z)
         if i != j
            z[i, j] = zero(base_ring(x))
         else
            z[i, j] = y
         end
      end
   end
   return z
end

(x::ArbMatSpace)(y::arb_mat) = y

###############################################################################
#
#   Matrix constructor
#
###############################################################################

function matrix(R::ArbField, arr::Array{T, 2}) where {T <: Union{Int, UInt, fmpz, fmpq, Float64, BigFloat, arb, AbstractString}}
   z = arb_mat(size(arr, 1), size(arr, 2), arr, prec(R))
   z.base_ring = R
   return z
end

function matrix(R::ArbField, r::Int, c::Int, arr::Array{T, 1}) where {T <: Union{Int, UInt, fmpz, fmpq, Float64, BigFloat, arb, AbstractString}}
   _check_dim(r, c, arr)
   z = arb_mat(r, c, arr, prec(R))
   z.base_ring = R
   return z
end

function matrix(R::ArbField, arr::Array{<: Integer, 2})
   arr_fmpz = map(fmpz, arr)
   return matrix(R, arr_fmpz)
end

function matrix(R::ArbField, r::Int, c::Int, arr::Array{<: Integer, 1})
   arr_fmpz = map(fmpz, arr)
   return matrix(R, r, c, arr_fmpz)
end

function matrix(R::ArbField, arr::Array{Rational{T}, 2}) where {T <: Integer}
   arr_fmpz = map(fmpq, arr)
   return matrix(R, arr_fmpz)
end

function matrix(R::ArbField, r::Int, c::Int, arr::Array{Rational{T}, 1}) where {T <: Integer}
   arr_fmpz = map(fmpq, arr)
   return matrix(R, r, c, arr_fmpz)
end

###############################################################################
#
#  Zero matrix
#
###############################################################################

function zero_matrix(R::ArbField, r::Int, c::Int)
   if r < 0 || c < 0
     error("dimensions must not be negative")
   end
   z = arb_mat(r, c)
   z.base_ring = R
   return z
end

###############################################################################
#
#  Identity matrix
#
###############################################################################

function identity_matrix(R::ArbField, n::Int)
   if n < 0
     error("dimension must not be negative")
   end
   z = arb_mat(n, n)
   ccall((:arb_mat_one, :libarb), Nothing, (Ref{arb_mat}, ), z)
   z.base_ring = R
   return z
end

###############################################################################
#
#   Promotions
#
###############################################################################

promote_rule(::Type{arb_mat}, ::Type{T}) where {T <: Integer} = arb_mat

promote_rule(::Type{arb_mat}, ::Type{Rational{T}}) where T <: Union{Int, BigInt} = arb_mat

promote_rule(::Type{arb_mat}, ::Type{fmpz}) = arb_mat

promote_rule(::Type{arb_mat}, ::Type{fmpq}) = arb_mat

promote_rule(::Type{arb_mat}, ::Type{arb}) = arb_mat

promote_rule(::Type{arb_mat}, ::Type{Float64}) = arb_mat

promote_rule(::Type{arb_mat}, ::Type{BigFloat}) = arb_mat

promote_rule(::Type{arb_mat}, ::Type{fmpz_mat}) = arb_mat

promote_rule(::Type{arb_mat}, ::Type{fmpq_mat}) = arb_mat

###############################################################################
#
#   MatrixSpace constructor
#
###############################################################################

function MatrixSpace(R::ArbField, r::Int, c::Int; cached = true)
  (r <= 0 || c <= 0) && error("Dimensions must be positive")
  return ArbMatSpace(R, r, c, cached)
end
