###############################################################################
#
#   fmpq_mat.jl : Flint matrices over the rationals
#
###############################################################################

export fmpq_mat, FmpqMatSpace, gso, hilbert

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

elem_type(::Type{FmpqMatSpace}) = fmpq_mat

parent_type(::Type{fmpq_mat}) = FmpqMatSpace

base_ring(a::FmpqMatSpace) = a.base_ring

base_ring(a::fmpq_mat) = a.base_ring

parent(a::fmpq_mat, cached::Bool = true) =
      FmpqMatSpace(rows(a), cols(a), cached)

function check_parent(a::fmpq_mat, b::fmpq_mat)
   (rows(a) != rows(b) || cols(a) != cols(b) || base_ring(a) != base_ring(b)) &&
      error("Incompatible matrices")
end

###############################################################################
#
#   Similar
#
###############################################################################

function similar(x::fmpq_mat)
   z = fmpq_mat(rows(x), cols(x))
   z.base_ring = x.base_ring
   return z
end

function similar(x::fmpq_mat, r::Int, c::Int)
   z = fmpq_mat(r, c)
   z.base_ring = x.base_ring
   return z
end

###############################################################################
#
#   Windows - handle with care!!!
#
###############################################################################

function Base.view(x::fmpq_mat, r1::Int, c1::Int, r2::Int, c2::Int)
  Generic._checkbounds(x, r1, c1)
  Generic._checkbounds(x, r2, c2)
  (r1 > r2 || c1 > c2) && error("Invalid parameters")
  b = fmpq_mat()
  b.base_ring = x.base_ring
  ccall((:fmpq_mat_window_init, :libflint), Void,
        (Ref{fmpq_mat}, Ref{fmpq_mat}, Int, Int, Int, Int),
            b, x, r1 - 1, c1 - 1, r2, c2)
  finalizer(b, _fmpq_mat_window_clear_fn)
  return b
end

function Base.view(x::fmpq_mat, r::UnitRange{Int}, c::UnitRange{Int})
  return Base.view(x, r.start, c.start, r.stop, c.stop)
end

function _fmpq_mat_window_clear_fn(a::fmpq_mat)
   ccall((:fmpq_mat_window_clear, :libflint), Void, (Ref{fmpq_mat},), a)
end

function sub(x::fmpq_mat, r1::Int, c1::Int, r2::Int, c2::Int)
   return deepcopy(view(x, r1, c1, r2, c2))
end

function sub(x::fmpq_mat, r::UnitRange{Int}, c::UnitRange{Int})
   return deepcopy(view(x, r, c))
end

getindex(x::fmpq_mat, r::UnitRange{Int}, c::UnitRange{Int}) = sub(x, r, c)

################################################################################
#
#   Size
#
################################################################################

size(x::fmpq_mat) = tuple(rows(x), cols(x))

size(t::fmpq_mat, d) = d <= 2 ? size(t)[d] : 1

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function getindex!(v::fmpq, a::fmpq_mat, r::Int, c::Int)
   z = ccall((:fmpq_mat_entry, :libflint), Ptr{fmpq},
             (Ref{fmpq_mat}, Int, Int), a, r - 1, c - 1)
   ccall((:fmpq_set, :libflint), Void, (Ref{fmpq}, Ptr{fmpq}), v, z)
end

@inline function getindex(a::fmpq_mat, r::Int, c::Int)
   @boundscheck Generic._checkbounds(a, r, c)
   v = fmpq()
   z = ccall((:fmpq_mat_entry, :libflint), Ptr{fmpq},
             (Ref{fmpq_mat}, Int, Int), a, r - 1, c - 1)
   ccall((:fmpq_set, :libflint), Void, (Ref{fmpq}, Ptr{fmpq}), v, z)
   return v
end

@inline function setindex!(a::fmpq_mat, d::fmpz, r::Int, c::Int)
   @boundscheck Generic._checkbounds(a, r, c)
   z = ccall((:fmpq_mat_entry_num, :libflint), Ptr{fmpz},
             (Ref{fmpq_mat}, Int, Int), a, r - 1, c - 1)
   ccall((:fmpz_set, :libflint), Void, (Ptr{fmpz}, Ref{fmpz}), z, d)
   z = ccall((:fmpq_mat_entry_den, :libflint), Ptr{fmpz},
             (Ref{fmpq_mat}, Int, Int), a, r - 1, c - 1)
   ccall((:fmpz_set_si, :libflint), Void, (Ptr{fmpz}, Int), z, 1)
end

@inline function setindex!(a::fmpq_mat, d::fmpq, r::Int, c::Int)
   @boundscheck Generic._checkbounds(a, r, c)
   z = ccall((:fmpq_mat_entry, :libflint), Ptr{fmpq},
             (Ref{fmpq_mat}, Int, Int), a, r - 1, c - 1)
   ccall((:fmpq_set, :libflint), Void, (Ptr{fmpq}, Ref{fmpq}), z, d)
end

Base.@propagate_inbounds setindex!(a::fmpq_mat, d::Integer,
                                 r::Int, c::Int) =
         setindex!(a, fmpq(d), r, c)

@inline function setindex!(a::fmpq_mat, d::Int, r::Int, c::Int)
   @boundscheck Generic._checkbounds(a, r, c)
   z = ccall((:fmpq_mat_entry, :libflint), Ptr{fmpq},
             (Ref{fmpq_mat}, Int, Int), a, r - 1, c - 1)
   ccall((:fmpq_set_si, :libflint), Void, (Ptr{fmpq}, Int, Int), z, d, 1)
end

Base.@propagate_inbounds setindex!(a::fmpq_mat, d::Rational,
                                 r::Int, c::Int) =
         setindex!(a, fmpq(d), r, c)

rows(a::fmpq_mat) = a.r

cols(a::fmpq_mat) = a.c

zero(a::FmpqMatSpace) = a()

one(a::FmpqMatSpace) = a(1)

iszero(a::fmpq_mat) = ccall((:fmpq_mat_is_zero, :libflint), Bool,
                            (Ref{fmpq_mat},), a)

isone(a::fmpq_mat) = ccall((:fmpq_mat_is_one, :libflint), Bool,
                           (Ref{fmpq_mat},), a)

function deepcopy_internal(d::fmpq_mat, dict::ObjectIdDict)
   z = fmpq_mat(d)
   z.base_ring = d.base_ring
   return z
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(a::fmpq_mat) = canonical_unit(a[1, 1])

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function show(io::IO, a::FmpqMatSpace)
   print(io, "Matrix Space of ")
   print(io, a.rows, " rows and ", a.cols, " columns over ")
   print(io, "Rational Field")
end

function show(io::IO, a::fmpq_mat)
   r = rows(a)
   c = cols(a)
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

show_minus_one(::Type{fmpq_mat}) = show_minus_one(fmpq)

###############################################################################
#
#   Unary operations
#
###############################################################################

function -(x::fmpq_mat)
   z = similar(x)
   ccall((:fmpq_mat_neg, :libflint), Void,
         (Ref{fmpq_mat}, Ref{fmpq_mat}), z, x)
   return z
end

###############################################################################
#
#   transpose
#
###############################################################################

function transpose(x::fmpq_mat)
   z = similar(x, cols(x), rows(x))
   ccall((:fmpq_mat_transpose, :libflint), Void,
         (Ref{fmpq_mat}, Ref{fmpq_mat}), z, x)
   return z
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function +(x::fmpq_mat, y::fmpq_mat)
   check_parent(x, y)
   z = similar(x)
   ccall((:fmpq_mat_add, :libflint), Void,
                (Ref{fmpq_mat}, Ref{fmpq_mat},  Ref{fmpq_mat}),
               z, x, y)
   return z
end

function -(x::fmpq_mat, y::fmpq_mat)
   check_parent(x, y)
   z = similar(x)
   ccall((:fmpq_mat_sub, :libflint), Void,
                (Ref{fmpq_mat}, Ref{fmpq_mat},  Ref{fmpq_mat}),
               z, x, y)
   return z
end

function *(x::fmpq_mat, y::fmpq_mat)
   cols(x) != rows(y) && error("Incompatible matrix dimensions")
   z = similar(x, rows(x), cols(y))
   ccall((:fmpq_mat_mul, :libflint), Void,
                (Ref{fmpq_mat}, Ref{fmpq_mat},  Ref{fmpq_mat}),
               z, x, y)
   return z
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(x::fmpz, y::fmpq_mat)
   z = similar(y)
   ccall((:fmpq_mat_scalar_mul_fmpz, :libflint), Void,
                (Ref{fmpq_mat}, Ref{fmpq_mat}, Ref{fmpz}), z, y, x)
   return z
end

function *(x::fmpq, y::fmpq_mat)
   z = similar(y)
   ccall((:fmpq_mat_scalar_mul_fmpz, :libflint), Void,
                (Ref{fmpq_mat}, Ref{fmpq_mat}, Ref{fmpq}), z, y, numerator(x))
   ccall((:fmpq_mat_scalar_div_fmpz, :libflint), Void,
                (Ref{fmpq_mat}, Ref{fmpq_mat}, Ref{fmpq}), z, z, denominator(x))
   return z
end

*(x::fmpq_mat, y::fmpq) = y*x

*(x::fmpq_mat, y::fmpz) = y*x

*(x::Integer, y::fmpq_mat) = fmpz(x)*y

*(x::fmpq_mat, y::Integer) = fmpz(y)*x

*(x::Rational, y::fmpq_mat) = fmpq(x)*y

*(x::fmpq_mat, y::Rational) = fmpq(y)*x

for T in [Integer, fmpz, fmpq]
   @eval begin
      function +(x::fmpq_mat, y::$T)
         z = deepcopy(x)
         for i = 1:min(rows(x), cols(x))
            z[i, i] += y
         end
         return z
      end

      +(x::$T, y::fmpq_mat) = y + x

      function -(x::fmpq_mat, y::$T)
         z = deepcopy(x)
         for i = 1:min(rows(x), cols(x))
            z[i, i] -= y
         end
         return z
      end

      function -(x::$T, y::fmpq_mat)
         z = -y
         for i = 1:min(rows(y), cols(y))
            z[i, i] += x
         end
         return z
      end
   end
end

function +(x::fmpq_mat, y::Rational)
   z = deepcopy(x)
   for i = 1:min(rows(x), cols(x))
      z[i, i] += y
   end
   return z
end

+(x::Rational, y::fmpq_mat) = y + x

function -(x::fmpq_mat, y::Rational)
   z = deepcopy(x)
   for i = 1:min(rows(x), cols(x))
      z[i, i] -= y
   end
   return z
end

function -(x::Rational, y::fmpq_mat)
   z = -y
   for i = 1:min(rows(y), cols(y))
      z[i, i] += x
   end
   return z
end

###############################################################################
#
#   Comparisons
#
###############################################################################

function ==(x::fmpq_mat, y::fmpq_mat)
   check_parent(x, y)
   ccall((:fmpq_mat_equal, :libflint), Bool,
                                       (Ref{fmpq_mat}, Ref{fmpq_mat}), x, y)
end

###############################################################################
#
#   Ad hoc comparisons
#
###############################################################################

function ==(x::fmpq_mat, y::Integer)
   for i = 1:min(rows(x), cols(x))
      if x[i, i] != y
         return false
      end
   end
   for i = 1:rows(x)
      for j = 1:cols(x)
         if i != j && x[i, j] != 0
            return false
         end
      end
   end
   return true
end

==(x::Integer, y::fmpq_mat) = y == x

==(x::fmpq_mat, y::Rational{T}) where T <: Union{Int, BigInt} = x == fmpq(y)

==(x::Rational{T}, y::fmpq_mat) where T <: Union{Int, BigInt} = y == x

###############################################################################
#
#   Inversion
#
###############################################################################

function inv(x::fmpq_mat)
   z = similar(x)
   success = ccall((:fmpq_mat_inv, :libflint), Cint,
         (Ref{fmpq_mat}, Ref{fmpq_mat}), z, x)
   success == 0 && error("Matrix not invertible")
   return z
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(x::fmpq_mat, y::fmpq_mat)
   cols(x) != cols(y) && error("Incompatible matrix dimensions")
   x*inv(y)
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(x::fmpq_mat, y::fmpq)
   z = similar(x)
   ccall((:fmpq_mat_scalar_div_fmpz, :libflint), Void,
                (Ref{fmpq_mat}, Ref{fmpq_mat}, Ref{fmpz}), z, x, numerator(y))
   ccall((:fmpq_mat_scalar_mul_fmpz, :libflint), Void,
                (Ref{fmpq_mat}, Ref{fmpq_mat}, Ref{fmpz}), z, z, denominator(y))
   return z
end

function divexact(x::fmpq_mat, y::fmpz)
   z = similar(x)
   ccall((:fmpq_mat_scalar_div_fmpz, :libflint), Void,
                (Ref{fmpq_mat}, Ref{fmpq_mat}, Ref{fmpz}), z, x, y)
   return z
end

divexact(x::fmpq_mat, y::Integer) = divexact(x, fmpz(y))

divexact(x::fmpq_mat, y::Rational{T}) where T <: Union{Int, BigInt} = divexact(x, fmpq(y))

###############################################################################
#
#   Characteristic polynomial
#
###############################################################################

function charpoly(R::FmpqPolyRing, x::fmpq_mat)
   rows(x) != cols(x) && error("Non-square")
   z = R()
   ccall((:fmpq_mat_charpoly, :libflint), Void,
                (Ref{fmpq_poly}, Ref{fmpq_mat}), z, x)
   return z
end

###############################################################################
#
#   Minimal polynomial
#
###############################################################################

function minpoly(R::FmpqPolyRing, x::fmpq_mat)
   rows(x) != cols(x) && error("Non-square")
   z = R()
   ccall((:fmpq_mat_minpoly, :libflint), Void,
                (Ref{fmpq_poly}, Ref{fmpq_mat}), z, x)
   return z
end

###############################################################################
#
#   Determinant
#
###############################################################################

function det(x::fmpq_mat)
   rows(x) != cols(x) && error("Non-square matrix")
   z = fmpq()
   ccall((:fmpq_mat_det, :libflint), Void,
                (Ref{fmpq}, Ref{fmpq_mat}), z, x)
   return z
end

###############################################################################
#
#   Gram-Schmidt orthogonalisation
#
###############################################################################

doc"""
    gso(x::fmpq_mat)
> Return the Gram-Schmidt Orthogonalisation of the matrix $x$.
"""
function gso(x::fmpq_mat)
   z = similar(x)
   ccall((:fmpq_mat_gso, :libflint), Void,
                (Ref{fmpq_mat}, Ref{fmpq_mat}), z, x)
   return z
end

###############################################################################
#
#   Hilbert matrix
#
###############################################################################

doc"""
    hilbert(R::FmpqMatSpace)
> Return the Hilbert matrix in the given matrix space. This is the matrix with
> entries $H_{i,j} = 1/(i + j - 1)$.
"""
function hilbert(R::FmpqMatSpace)
   z = R()
   ccall((:fmpq_mat_hilbert_matrix, :libflint), Bool,
                   (Ref{fmpq_mat},), z)
   return z
end

###############################################################################
#
#   Rank
#
###############################################################################

function rank(x::fmpq_mat)
   z = similar(x)
   r = ccall((:fmpq_mat_rref, :libflint), Int,
         (Ref{fmpq_mat}, Ref{fmpq_mat}), z, x)
   return r
end

###############################################################################
#
#   Reduced row echelon form
#
###############################################################################

function rref(x::fmpq_mat)
   z = similar(x)
   r = ccall((:fmpq_mat_rref, :libflint), Int,
         (Ref{fmpq_mat}, Ref{fmpq_mat}), z, x)
   return r, z
end

###############################################################################
#
#   Linear solving
#
###############################################################################

function solve(a::fmpq_mat, b::fmpq_mat)
   rows(a) != cols(a) && error("Not a square matrix in solve")
   rows(b) != rows(a) && error("Incompatible dimensions in solve")
   z = similar(b)
   nonsing = ccall((:fmpq_mat_solve_fraction_free, :libflint), Bool,
      (Ref{fmpq_mat}, Ref{fmpq_mat}, Ref{fmpq_mat}), z, a, b)
   !nonsing && error("Singular matrix in solve")
   return z
end

doc"""
    solve_dixon(a::fmpq_mat, b::fmpq_mat)
> Solve $ax = b$ by clearing denominators and using Dixon's algorithm. This is
> usually faster for large systems.
"""
function solve_dixon(a::fmpq_mat, b::fmpq_mat)
   rows(a) != cols(a) && error("Not a square matrix in solve")
   rows(b) != rows(a) && error("Incompatible dimensions in solve")
   z = similar(b)
   nonsing = ccall((:fmpq_mat_solve_dixon, :libflint), Bool,
      (Ref{fmpq_mat}, Ref{fmpq_mat}, Ref{fmpq_mat}), z, a, b)
   !nonsing && error("Singular matrix in solve")
   return z
end

###############################################################################
#
#   Trace
#
###############################################################################

function trace(x::fmpq_mat)
   rows(x) != cols(x) && error("Not a square matrix in trace")
   d = fmpq()
   ccall((:fmpq_mat_trace, :libflint), Void,
                (Ref{fmpq}, Ref{fmpq_mat}), d, x)
   return d
end

###############################################################################
#
#   Concatenation
#
###############################################################################

function hcat(a::fmpq_mat, b::fmpq_mat)
  rows(a) != rows(b) && error("Incompatible number of rows in hcat")
  c = similar(a, rows(a), cols(a) + cols(b))
  ccall((:fmpq_mat_concat_horizontal, :libflint), Void,
        (Ref{fmpq_mat}, Ref{fmpq_mat}, Ref{fmpq_mat}), c, a, b)
  return c
end

function vcat(a::fmpq_mat, b::fmpq_mat)
  cols(a) != cols(b) && error("Incompatible number of columns in vcat")
  c = similar(a, rows(a) + rows(b), cols(a))
  ccall((:fmpq_mat_concat_vertical, :libflint), Void,
        (Ref{fmpq_mat}, Ref{fmpq_mat}, Ref{fmpq_mat}), c, a, b)
  return c
end

###############################################################################
#
#   Similarity
#
###############################################################################

function similarity!(z::fmpq_mat, r::Int, d::fmpq)
   ccall((:fmpq_mat_similarity, :libflint), Void, 
         (Ref{fmpq_mat}, Int, Ref{fmpq}), z, r - 1, d)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function mul!(z::fmpq_mat, x::fmpq_mat, y::fmpq_mat)
   ccall((:fmpq_mat_mul, :libflint), Void,
                (Ref{fmpq_mat}, Ref{fmpq_mat}, Ref{fmpq_mat}), z, x, y)
   return z
end

function mul!(y::fmpq_mat, x::Int)
   ccall((:fmpq_mat_scalar_mul_fmpz, :libflint), Void,
                (Ref{fmpq_mat}, Ref{fmpq_mat}, Ref{fmpq}), y, y, fmpz(x))
   return y
end

function mul!(y::fmpq_mat, x::fmpz)
   ccall((:fmpq_mat_scalar_mul_fmpz, :libflint), Void,
                (Ref{fmpq_mat}, Ref{fmpq_mat}, Ref{fmpz}), y, y, x)
   return y
end

function addeq!(z::fmpq_mat, x::fmpq_mat)
   ccall((:fmpq_mat_add, :libflint), Void,
                (Ref{fmpq_mat}, Ref{fmpq_mat}, Ref{fmpq_mat}), z, z, x)
   return z
end

function zero!(z::fmpq_mat)
   ccall((:fmpq_mat_zero, :libflint), Void,
                (Ref{fmpq_mat},), z)
   return z
end

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function (a::FmpqMatSpace)()
   z = fmpq_mat(a.rows, a.cols)
   z.base_ring = a.base_ring
   return z
end

function (a::FmpqMatSpace)(arr::Array{fmpq, 2})
   _check_dim(a.rows, a.cols, arr)
   z = fmpq_mat(a.rows, a.cols, arr)
   z.base_ring = a.base_ring
   return z
end

function (a::FmpqMatSpace)(arr::Array{fmpz, 2})
   _check_dim(a.rows, a.cols, arr)
   z = fmpq_mat(a.rows, a.cols, arr)
   z.base_ring = a.base_ring
   return z
end


function (a::FmpqMatSpace)(arr::Array{T, 2}) where {T <: Integer}
   _check_dim(a.rows, a.cols, arr)
   z = fmpq_mat(a.rows, a.cols, arr)
   z.base_ring = a.base_ring
   return z
end

function (a::FmpqMatSpace)(arr::Array{Rational{T}, 2}) where {T <: Integer}
   _check_dim(a.rows, a.cols, arr)
   z = fmpq_mat(a.rows, a.cols, map(fmpq, arr))
   z.base_ring = a.base_ring
   return z
end

function (a::FmpqMatSpace)(arr::Array{fmpq, 1})
   _check_dim(a.rows, a.cols, arr)
   z = fmpq_mat(a.rows, a.cols, arr)
   z.base_ring = a.base_ring
   return z
end

function (a::FmpqMatSpace)(arr::Array{fmpz, 1})
   _check_dim(a.rows, a.cols, arr)
   z = fmpq_mat(a.rows, a.cols, arr)
   z.base_ring = a.base_ring
   return z
end

function (a::FmpqMatSpace)(arr::Array{T, 1}) where {T <: Integer}
   _check_dim(a.rows, a.cols, arr)
   z = fmpq_mat(a.rows, a.cols, arr)
   z.base_ring = a.base_ring
   return z
end

function (a::FmpqMatSpace)(arr::Array{Rational{T}, 1}) where {T <: Integer}
   _check_dim(a.rows, a.cols, arr)
   z = fmpq_mat(a.rows, a.cols, map(fmpq, arr))
   z.base_ring = a.base_ring
   return z
end

function (a::FmpqMatSpace)(d::fmpq)
   z = fmpq_mat(a.rows, a.cols, d)
   z.base_ring = a.base_ring
   return z
end

function (a::FmpqMatSpace)(d::fmpz)
   z = fmpq_mat(a.rows, a.cols, fmpq(d))
   z.base_ring = a.base_ring
   return z
end

function (a::FmpqMatSpace)(d::Integer)
   z = fmpq_mat(a.rows, a.cols, fmpq(d))
   z.base_ring = a.base_ring
   return z
end

(a::FmpqMatSpace)(d::Rational) = a(fmpq(d))

function (a::FmpqMatSpace)(M::fmpz_mat)
   (a.cols == cols(M) && a.rows == rows(M)) || error("wrong matrix dimension")
   z = a()
   ccall((:fmpq_mat_set_fmpz_mat, :libflint), Void, (Ref{fmpq_mat}, Ref{fmpz_mat}), z, M)
   return z
end

(a::FmpqMatSpace)(d::fmpq_mat) = d

###############################################################################
#
#   Promotions
#
###############################################################################

promote_rule(::Type{fmpq_mat}, ::Type{T}) where {T <: Integer} = fmpq_mat

promote_rule(::Type{fmpq_mat}, ::Type{fmpq}) = fmpq_mat

promote_rule(::Type{fmpq_mat}, ::Type{fmpz}) = fmpq_mat

promote_rule(::Type{fmpq_mat}, ::Type{Rational{T}}) where T <: Union{Int, BigInt} = fmpq_mat

###############################################################################
#
#   Matrix constructor
#
###############################################################################

function matrix(R::FlintRationalField, arr::Array{fmpq, 2})
   z = fmpq_mat(size(arr, 1), size(arr, 2), arr)
   z.base_ring = FlintQQ
   return z
end

function matrix(R::FlintRationalField, arr::Array{<: Union{fmpz, Int, BigInt}, 2})
   z = fmpq_mat(size(arr, 1), size(arr, 2), arr)
   z.base_ring = FlintQQ
   return z
end

function matrix(R::FlintRationalField, arr::Array{Rational{T}, 2}) where {T <: Integer}
   z = fmpq_mat(size(arr, 1), size(arr, 2), map(fmpq, arr))
   z.base_ring = FlintQQ
   return z
end

function matrix(R::FlintRationalField, r::Int, c::Int, arr::Array{fmpq, 1})
   _check_dim(r, c, arr)
   z = fmpq_mat(r, c, arr)
   z.base_ring = FlintQQ
   return z
end

function matrix(R::FlintRationalField, r::Int, c::Int, arr::Array{<: Union{fmpz, Int, BigInt}, 1})
   _check_dim(r, c, arr)
   z = fmpq_mat(r, c, arr)
   z.base_ring = FlintQQ
   return z
end

function matrix(R::FlintRationalField, r::Int, c::Int, arr::Array{Rational{T}, 1}) where {T <: Union{fmpz, Int, BigInt}}
   _check_dim(r, c, arr)
   z = fmpq_mat(r, c, map(fmpq, arr))
   z.base_ring = FlintQQ
   return z
end

###############################################################################
#
#  Zero matrix
#
###############################################################################

function zero_matrix(R::FlintRationalField, r::Int, c::Int)
   z = fmpq_mat(r, c)
   z.base_ring = FlintQQ
   return z
end

###############################################################################
#
#  Identity matrix
#
###############################################################################

function identity_matrix(R::FlintRationalField, n::Int)
   z = fmpq_mat(n, n)
   ccall((:fmpq_mat_one, :libflint), Void, (Ref{fmpq_mat}, ), z)
   z.base_ring = FlintQQ
   return z
end

###############################################################################
#
#   MatrixSpace constructor
#
###############################################################################

function MatrixSpace(R::FlintRationalField, r::Int, c::Int; cached = true)
   return FmpqMatSpace(r, c, cached)
end
