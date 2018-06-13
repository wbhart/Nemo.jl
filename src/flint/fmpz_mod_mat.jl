################################################################################
#
#  fmpz_mod_mat.jl: flint fmpz_mod_mat types in julia
#
################################################################################

export fmpz_mod_mat, FmpzModMatSpace, getindex, setindex!, set_entry!, deepcopy, rows, 
       cols, parent, base_ring, zero, one, show, transpose,
       transpose!, rref, rref!, trace, det, rank, inv, solve, lufact,
       sub, hcat, vcat, Array, lift, lift!, MatrixSpace, check_parent,
       howell_form, howell_form!, strong_echelon_form, strong_echelon_form!

################################################################################
#
#  Data type and parent object methods
#
################################################################################

parent_type(::Type{fmpz_mod_mat}) = FmpzModMatSpace

elem_type(::Type{FmpzModMatSpace}) = fmpz_mod_mat

function check_parent(x::fmpz_mod_mat, y::fmpz_mod_mat)
  base_ring(x) != base_ring(y) && error("Residue rings must be equal")
  (cols(x) != cols(y)) && (rows(x) != rows(y)) &&
          error("Matrices have wrong dimensions")
  return nothing
end

###############################################################################
#
#   Similar
#
###############################################################################

function similar(x::fmpz_mod_mat)
   z = fmpz_mod_mat(rows(x), cols(x), x.n)
   z.base_ring = x.base_ring
   return z
end

function similar(x::fmpz_mod_mat, r::Int, c::Int)
   z = fmpz_mod_mat(r, c, x.n)
   z.base_ring = x.base_ring
   return z
end

################################################################################
#
#  Manipulation
#
################################################################################

@inline function getindex(a::fmpz_mod_mat, i::Int, j::Int)
  @boundscheck Generic._checkbounds(a, i, j)
  u = ccall((:fmpz_mod_mat_get_entry, :libflint), UInt,
              (Ref{fmpz_mod_mat}, Int, Int), a, i - 1 , j - 1)
  return base_ring(a)(u)
end

#as above, but as a plain UInt
function getindex_raw(a::fmpz_mod_mat, i::Int, j::Int)
  return ccall((:fmpz_mod_mat_get_entry, :libflint), UInt,
                 (Ref{fmpz_mod_mat}, Int, Int), a, i - 1, j - 1)
end

@inline function setindex!(a::fmpz_mod_mat, u::UInt, i::Int, j::Int)
  @boundscheck Generic._checkbounds(a, i, j)
  set_entry!(a, i, j, u)
end

@inline function setindex!(a::fmpz_mod_mat, u::fmpz, i::Int, j::Int)
  @boundscheck Generic._checkbounds(a, i, j)
  set_entry!(a, i, j, u)
end

@inline function setindex!(a::fmpz_mod_mat, u::nmod, i::Int, j::Int)
  @boundscheck Generic._checkbounds(a, i, j)
  (base_ring(a) != parent(u)) && error("Parent objects must coincide") 
  set_entry!(a, i, j, u.data)
end

setindex!(a::fmpz_mod_mat, u::Integer, i::Int, j::Int) =
        setindex!(a, fmpz(u), i, j)

setindex_t!(a::fmpz_mod_mat, u::T, i::Int, j::Int) where {T<:Union{RingElem, Integer}} =
  setindex!(a, u, j, i)

function set_entry!(a::fmpz_mod_mat, i::Int, j::Int, u::fmpz)
  ccall((:fmpz_mod_mat_set_entry, :libflint), Void,
        (Ref{fmpz_mod_mat}, Int, Int, Ref{fmpz}), a, i - 1, j - 1, u)
end

function set_entry!(a::fmpz_mod_mat, i::Int, j::Int, u::Integer)
  ccall((:fmpz_mod_mat_set_entry, :libflint), Void,
        (Ref{fmpz_mod_mat}, Int, Int, Ref{fmpz}), a, i - 1, j - 1, fmpz(u))
end

set_entry!(a::fmpz_mod_mat, i::Int, j::Int, u::nmod) =
        set_entry!(a, i, j, u.data)

set_entry_t!(a::fmpz_mod_mat, i::Int, j::Int, u::T) where {T<:Union{RingElem, Integer}} =
  set_entry!(a, j, i, u)
 
function deepcopy_internal(a::fmpz_mod_mat, dict::ObjectIdDict)
  z = fmpz_mod_mat(rows(a), cols(a), a.n)
  if isdefined(a, :base_ring)
    z.base_ring = a.base_ring
  end
  ccall((:fmpz_mod_mat_set, :libflint), Void,
          (Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}), z, a)
  return z
end

rows(a::fmpz_mod_mat) = a.r

cols(a::fmpz_mod_mat) = a.c

parent(a::fmpz_mod_mat, cached::Bool = true) = MatrixSpace(base_ring(a), rows(a), cols(a), cached)

base_ring(a::FmpzModMatSpace) = a.base_ring

base_ring(a::fmpz_mod_mat) = a.base_ring

zero(a::FmpzModMatSpace) = a()

function one(a::FmpzModMatSpace)
  (a.rows != a.cols) && error("Matrices must be quadratic")
  z = a()
  ccall((:fmpz_mod_mat_one, :libflint), Void, (Ref{fmpz_mod_mat}, ), z)
  return z
end

function iszero(a::fmpz_mod_mat)
  r = ccall((:fmpz_mod_mat_is_zero, :libflint), Cint, (Ref{fmpz_mod_mat}, ), a)
  return Bool(r)
end

################################################################################
#
#  AbstractString I/O
#
################################################################################

function show(io::IO, a::FmpzModMatSpace)
   print(io, "Matrix Space of ")
   print(io, a.rows, " rows and ", a.cols, " columns over ")
   print(io, a.base_ring)
end

function show(io::IO, a::fmpz_mod_mat)
   rows = a.r
   cols = a.c
   if rows * cols == 0
      print(io, "$rows by $cols matrix")
   end
   for i = 1:rows
      print(io, "[")
      for j = 1:cols
         print(io, a[i, j])
         if j != cols
            print(io, " ")
         end
      end
      print(io, "]")
      if i != rows
         println(io, "")
      end
   end
end

################################################################################
#
#  Comparison
#
################################################################################

==(a::fmpz_mod_mat, b::fmpz_mod_mat) = (a.base_ring == b.base_ring) &&
        Bool(ccall((:fmpz_mod_mat_equal, :libflint), Cint,
                (Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}), a, b))

################################################################################
#
#  Transpose
#
################################################################################

function transpose(a::fmpz_mod_mat)
  z = similar(a, cols(a), rows(a))
  ccall((:fmpz_mod_mat_transpose, :libflint), Void,
          (Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}), z, a)
  return z
end

function transpose!(a::fmpz_mod_mat)
  !issquare(a) && error("Matrix must be a square matrix")
  ccall((:fmpz_mod_mat_transpose, :libflint), Void,
          (Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}), a, a)
end

################################################################################
#
#  Unary operators
#
################################################################################

function -(x::fmpz_mod_mat)
  z = similar(x)
  ccall((:fmpz_mod_mat_neg, :libflint), Void,
          (Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}), z, x)
  return z
end

################################################################################
#
#  Binary operators
#
################################################################################

function +(x::fmpz_mod_mat, y::fmpz_mod_mat)
  check_parent(x,y)
  z = similar(x)
  ccall((:fmpz_mod_mat_add, :libflint), Void,
          (Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}), z, x, y)
  return z
end

function -(x::fmpz_mod_mat, y::fmpz_mod_mat)
  check_parent(x,y)
  z = similar(x)
  ccall((:fmpz_mod_mat_sub, :libflint), Void,
          (Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}), z, x, y)
  return z
end

function *(x::fmpz_mod_mat, y::fmpz_mod_mat)
  (base_ring(x) != base_ring(y)) && error("Base ring must be equal")
  (cols(x) != rows(y)) && error("Dimensions are wrong")
  z = similar(x, rows(x), cols(y))
  ccall((:fmpz_mod_mat_mul, :libflint), Void,
          (Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}), z, x, y)
  return z
end


################################################################################
#
#  Unsafe operations
#
################################################################################

function mul!(a::fmpz_mod_mat, b::fmpz_mod_mat, c::fmpz_mod_mat)
  ccall((:fmpz_mod_mat_mul, :libflint), Void, (Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}), a, b, c)
  return a
end

function add!(a::fmpz_mod_mat, b::fmpz_mod_mat, c::fmpz_mod_mat)
  ccall((:fmpz_mod_mat_add, :libflint), Void, (Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}), a, b, c)
  return a
end

function zero!(a::fmpz_mod_mat)
  ccall((:fmpz_mod_mat_zero, :libflint), Void, (Ref{fmpz_mod_mat}, ), a)
  return a
end

################################################################################
#
#  Ad hoc binary operators
#
################################################################################

function *(x::fmpz_mod_mat, y::UInt)
  z = similar(x)
  ccall((:fmpz_mod_mat_scalar_mul_ui, :libflint), Void,
          (Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}, UInt), z, x, y)
  return z
end

*(x::UInt, y::fmpz_mod_mat) = y*x

function *(x::fmpz_mod_mat, y::Int)
  z = similar(x)
  ccall((:fmpz_mod_mat_scalar_mul_si, :libflint), Void,
          (Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}, Int), z, x, y)
  return z
end

*(x::Int, y::fmpz_mod_mat) = y*x

function *(x::fmpz_mod_mat, y::fmpz)
  z = similar(x)
  ccall((:fmpz_mod_mat_scalar_mul_fmpz, :libflint), Void,
          (Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}, UInt), z, x, y)
  return x*tt
end

*(x::fmpz, y::fmpz_mod_mat) = y*x

function *(x::fmpz_mod_mat, y::Integer)
  return x*fmpz(y)
end

*(x::Integer, y::fmpz_mod_mat) = y*x

function *(x::fmpz_mod_mat, y::nmod)
  (base_ring(x) != parent(y)) && error("Parent objects must coincide")
  return x*y.data
end

*(x::Generic.Res{fmpz}, y::fmpz_mod_mat) = y*x

################################################################################
#
#  Powering
#
################################################################################

function ^(x::fmpz_mod_mat, y::UInt)
  z = similar(x)
  ccall((:fmpz_mod_mat_pow, :libflint), Void,
          (Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}, UInt), z, x, y)
  return z
end

function ^(x::fmpz_mod_mat, y::Int)
  ( y < 0 ) && error("Exponent must be positive")
  return x^UInt(y)
end

function ^(x::fmpz_mod_mat, y::fmpz)
  ( y < 0 ) && error("Exponent must be positive")
  ( y > fmpz(typemax(UInt))) &&
          error("Exponent must be smaller then ", fmpz(typemax(UInt)))
  return x^(UInt(y))
end

################################################################################
#
#  Row echelon form
#
################################################################################

function rref(a::fmpz_mod_mat)
  z = deepcopy(a)
  r = ccall((:fmpz_mod_mat_rref, :libflint), Int, (Ref{fmpz_mod_mat}, ), z)
  return r, z
end

function rref!(a::fmpz_mod_mat)
  r = ccall((:fmpz_mod_mat_rref, :libflint), Int, (Ref{fmpz_mod_mat}, ), a)
  return r
end

################################################################################
#
#  Strong echelon form and Howell form
#
################################################################################

function strong_echelon_form!(a::fmpz_mod_mat)
  ccall((:fmpz_mod_mat_strong_echelon_form, :libflint), Void, (Ref{fmpz_mod_mat}, ), a)
end

doc"""
    strong_echelon_form(a::fmpz_mod_mat)
> Return the strong echeleon form of $a$. The matrix $a$ must have at least as
> many rows as columns.
"""
function strong_echelon_form(a::fmpz_mod_mat)
  (rows(a) < cols(a)) &&
              error("Matrix must have at least as many rows as columns")
  z = deepcopy(a)
  strong_echelon_form!(z)
  return z
end

function howell_form!(a::fmpz_mod_mat)
  ccall((:fmpz_mod_mat_howell_form, :libflint), Void, (Ref{fmpz_mod_mat}, ), a)
end

doc"""
    howell_form(a::fmpz_mod_mat)
> Return the Howell normal form of $a$. The matrix $a$ must have at least as
> many rows as columns.
"""
function howell_form(a::fmpz_mod_mat)
  (rows(a) < cols(a)) &&
              error("Matrix must have at least as many rows as columns")

  z = deepcopy(a)
  howell_form!(z)
  return z
end

################################################################################
#
#  Trace
#
################################################################################

function trace(a::fmpz_mod_mat)
  !issquare(a) && error("Matrix must be a square matrix")
  t = fmpz()
  r = ccall((:fmpz_mod_mat_trace, :libflint), Void, (Ref{fmpz},
            Ref{fmpz_mod_mat}), t, a)
  return t
end

################################################################################
#
#  Determinant
#
################################################################################

function det(a::fmpz_mod_mat)
  !issquare(a) && error("Matrix must be a square matrix")
  if is_prime(a.n)
     r = ccall((:fmpz_mod_mat_det, :libflint), UInt, (Ref{fmpz_mod_mat}, ), a)
     return base_ring(a)(r)
  else
     try
        return det_fflu(a)
     catch
        return det_df(a)
     end
  end
end

################################################################################
#
#  Rank
#
################################################################################

function rank(a::fmpz_mod_mat)
  r = ccall((:fmpz_mod_mat_rank, :libflint), Int, (Ref{fmpz_mod_mat}, ), a)
  return r
end

################################################################################
#
#  Inverse
#
################################################################################

function inv(a::fmpz_mod_mat)
  !issquare(a) && error("Matrix must be a square matrix")
  z = similar(a)
  r = ccall((:fmpz_mod_mat_inv, :libflint), Int,
          (Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}), z, a)
  !Bool(r) && error("Matrix not invertible")
  return z
end

################################################################################
#
#  Linear solving
#
################################################################################

function solve(x::fmpz_mod_mat, y::fmpz_mod_mat)
  (base_ring(x) != base_ring(y)) && error("Matrices must have same base ring")
  !issquare(x)&& error("First argument not a square matrix in solve")
  (y.r != x.r) || y.c != 1 && ("Not a column vector in solve")
  z = similar(y)
  r = ccall((:fmpz_mod_mat_solve, :libflint), Int,
          (Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}), z, x, y)
  !Bool(r) && error("Singular matrix in solve")
  return z
end

################################################################################
#
#  LU decomposition
#
################################################################################

function lufact!(P::Generic.perm, x::fmpz_mod_mat)
  rank = ccall((:fmpz_mod_mat_lu, :libflint), Cint, (Ptr{Int}, Ref{fmpz_mod_mat}, Cint),
           P.d, x, 0)

  for i in 1:length(P.d)
    P.d[i] += 1
  end

  # flint does x == PLU instead of Px == LU (docs are wrong)
  inv!(P)

  return rank
end

function lufact(x::fmpz_mod_mat, P = PermGroup(rows(x)))
  m = rows(x)
  n = cols(x)
  P.n != m && error("Permutation does not match matrix")
  p = P()
  R = base_ring(x)
  U = deepcopy(x)

  L = similar(x, m, m)

  rank = lufact!(p, U)

  for i = 1:m
    for j = 1:n
      if i > j
        L[i, j] = U[i, j]
        U[i, j] = R()
      elseif i == j
        L[i, j] = R(1)
      elseif j <= m
        L[i, j] = R()
      end
    end
  end
  return rank, p, L, U
end

################################################################################
#
#  Windowing
#
################################################################################

function Base.view(x::fmpz_mod_mat, r1::Int, c1::Int, r2::Int, c2::Int)
  Generic._checkbounds(x, r1, c1)
  Generic._checkbounds(x, r2, c2)
  (r1 > r2 || c1 > c2) && error("Invalid parameters")
  z = fmpz_mod_mat()
  z.base_ring = x.base_ring
  ccall((:fmpz_mod_mat_window_init, :libflint), Void,
          (Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}, Int, Int, Int, Int),
          z, x, r1 - 1, c1 - 1, r2, c2)
  finalizer(z, _fmpz_mod_mat_window_clear_fn)
  return z
end

function Base.view(x::fmpz_mod_mat, r::UnitRange{Int}, c::UnitRange{Int})
  return Base.view(x, r.start, c.start, r.stop, c.stop)
end

function _fmpz_mod_mat_window_clear_fn(a::fmpz_mod_mat)
  ccall((:fmpz_mod_mat_window_clear, :libflint), Void, (Ref{fmpz_mod_mat}, ), a)
end

function sub(x::fmpz_mod_mat, r1::Int, c1::Int, r2::Int, c2::Int)
  return deepcopy(Base.view(x, r1, c1, r2, c2))
end

function sub(x::fmpz_mod_mat, r::UnitRange{Int}, c::UnitRange{Int})
  return deepcopy(Base.view(x, r, c))
end

getindex(x::fmpz_mod_mat, r::UnitRange{Int}, c::UnitRange{Int}) = sub(x, r, c)

################################################################################
#
#  Concatenation
#
################################################################################

function hcat(x::fmpz_mod_mat, y::fmpz_mod_mat)
  (base_ring(x) != base_ring(y)) && error("Matrices must have same base ring")
  (x.r != y.r) && error("Matrices must have same number of rows")
  z = similar(x, rows(x), cols(x) + cols(y))
  ccall((:fmpz_mod_mat_concat_horizontal, :libflint), Void,
          (Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}), z, x, y)
  return z
end

function vcat(x::fmpz_mod_mat, y::fmpz_mod_mat)
  (base_ring(x) != base_ring(y)) && error("Matrices must have same base ring")
  (x.c != y.c) && error("Matrices must have same number of columns")
  z = similar(x, rows(x) + rows(y), cols(x))
  ccall((:fmpz_mod_mat_concat_vertical, :libflint), Void,
          (Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}, Ref{fmpz_mod_mat}), z, x, y)
  return z
end

################################################################################
#
#  Conversion
#
################################################################################

function Array(b::fmpz_mod_mat)
  a = Array{Generic.Res{fmpz}}(b.r, b.c)
  for i = 1:b.r
    for j = 1:b.c
      a[i,j] = b[i,j]
    end
  end
  return a
end

################################################################################
#
#  Lifting
#
################################################################################

doc"""
    lift(a::fmpz_mod_mat)
> Return a lift of the matrix $a$ to a matrix over $\mathbb{Z}$, i.e. where the
> entries of the returned matrix are those of $a$ lifted to $\mathbb{Z}$.
"""
function lift(a::fmpz_mod_mat)
    return a.mat
end

function lift!(z::fmpz_mat, a::fmpz_mod_mat)
    z = a.mat
  return z 
end

################################################################################
#
#  Characteristic polynomial
#
################################################################################

function charpoly(R::FmpzModPolyRing, a::fmpz_mod_mat)
  m = deepcopy(a)
  p = R()
  ccall((:fmpz_mod_mat_charpoly, :libflint), Void,
          (Ref{nmod_poly}, Ref{fmpz_mod_mat}), p, m)
  return p
end

################################################################################
#
#  Minimal polynomial
#
################################################################################

function minpoly(R::FmpzModPolyRing, a::fmpz_mod_mat)
  p = R()
  ccall((:fmpz_mod_mat_minpoly, :libflint), Void,
          (Ref{nmod_poly}, Ref{fmpz_mod_mat}), p, a)
  return p
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{fmpz_mod_mat}, ::Type{V}) where {V <: Integer} = fmpz_mod_mat

promote_rule(::Type{fmpz_mod_mat}, ::Type{nmod}) = fmpz_mod_mat

promote_rule(::Type{fmpz_mod_mat}, ::Type{fmpz}) = fmpz_mod_mat

################################################################################
#
#  Parent object overloading
#
################################################################################

function (a::FmpzModMatSpace)()
  z = fmpz_mod_mat(a.rows, a.cols, a.n)
  z.base_ring = a.base_ring
  return z
end

function (a::FmpzModMatSpace)(b::Integer)
   M = a()
   for i = 1:a.rows
      for j = 1:a.cols
         if i != j
            M[i, j] = zero(base_ring(a))
         else
            M[i, j] = base_ring(a)(b)
         end
      end
   end
   return M
end

function (a::FmpzModMatSpace)(b::fmpz)
   M = a()
   for i = 1:a.rows
      for j = 1:a.cols
         if i != j
            M[i, j] = zero(base_ring(a))
         else
            M[i, j] = base_ring(a)(b)
         end
      end
   end
   return M
end

function (a::FmpzModMatSpace)(b::nmod)
   parent(b) != base_ring(a) && error("Unable to coerce to matrix")
   M = a()
   for i = 1:a.rows
      for j = 1:a.cols
         if i != j
            M[i, j] = zero(base_ring(a))
         else
            M[i, j] = fmpz(b.data)
         end
      end
   end
   return M
end

function (a::FmpzModMatSpace)(arr::Array{BigInt, 2}, transpose::Bool = false)
  _check_dim(a.rows, a.cols, arr, transpose)
  z = fmpz_mod_mat(a.rows, a.cols, a.n, arr, transpose)
  z.base_ring = a.base_ring
  return z
end

function (a::FmpzModMatSpace)(arr::Array{BigInt, 1}, transpose::Bool = false)
  _check_dim(a.rows, a.cols, arr)
  z = fmpz_mod_mat(a.rows, a.cols, a.n, arr, transpose)
  z.base_ring = a.base_ring
  return z
end

function (a::FmpzModMatSpace)(arr::Array{fmpz, 2}, transpose::Bool = false)
  _check_dim(a.rows, a.cols, arr, transpose)
  z = fmpz_mod_mat(a.rows, a.cols, a.n, arr, transpose)
  z.base_ring = a.base_ring
  return z
end

function (a::FmpzModMatSpace)(arr::Array{fmpz, 1}, transpose::Bool = false)
  _check_dim(a.rows, a.cols, arr)
  z = fmpz_mod_mat(a.rows, a.cols, a.n, arr, transpose)
  z.base_ring = a.base_ring
  return z
end

function (a::FmpzModMatSpace)(arr::Array{Int, 2}, transpose::Bool = false)
  _check_dim(a.rows, a.cols, arr, transpose)
  z = fmpz_mod_mat(a.rows, a.cols, a.n, arr, transpose)
  z.base_ring = a.base_ring
  return z
end

function (a::FmpzModMatSpace)(arr::Array{Int, 1}, transpose::Bool = false)
  _check_dim(a.rows, a.cols, arr)
  z = fmpz_mod_mat(a.rows, a.cols, a.n, arr, transpose)
  z.base_ring = a.base_ring
  return z
end

function (a::FmpzModMatSpace)(arr::Array{nmod, 2}, transpose::Bool = false)
  _check_dim(a.rows, a.cols, arr, transpose)
  (length(arr) > 0 && (base_ring(a) != parent(arr[1]))) && error("Elements must have same base ring")
  z = fmpz_mod_mat(a.rows, a.cols, a.n, arr, transpose)
  z.base_ring = a.base_ring
  return z
end

function (a::FmpzModMatSpace)(arr::Array{nmod, 1}, transpose::Bool = false)
  _check_dim(a.rows, a.cols, arr)
  (length(arr) > 0 && (base_ring(a) != parent(arr[1]))) && error("Elements must have same base ring")
  z = fmpz_mod_mat(a.rows, a.cols, a.n, arr, transpose)
  z.base_ring = a.base_ring
  return z
end

function (a::FmpzModMatSpace)(b::fmpz_mat)
  (a.cols != b.c || a.rows != b.r) && error("Dimensions do not fit")
  z = fmpz_mod_mat(a.n, b)
  z.base_ring = a.base_ring
  return z
end

###############################################################################
#
#   Matrix constructor
#
###############################################################################

function matrix(R::Generic.ResRing{fmpz}, arr::Array{<: Union{nmod, fmpz, Integer}, 2})
   z = fmpz_mod_mat(size(arr, 1), size(arr, 2), R.n, arr)
   z.base_ring = R
   return z
end

function matrix(R::Generic.ResRing{fmpz}, r::Int, c::Int, arr::Array{<: Union{nmod, fmpz, Integer}, 1})
   _check_dim(r, c, arr)
   z = fmpz_mod_mat(r, c, R.n, arr)
   z.base_ring = R
   return z
end

###############################################################################
#
#  Zero matrix
#
###############################################################################

function zero_matrix(R::Generic.ResRing{fmpz}, r::Int, c::Int)
   if r < 0 || c < 0
     error("dimensions must not be negative")
   end
   z = fmpz_mod_mat(r, c, R.n)
   z.base_ring = R
   return z
end

###############################################################################
#
#  Identity matrix
#
###############################################################################

function identity_matrix(R::Generic.ResRing{fmpz}, n::Int)
   z = zero_matrix(R, n, n)
   for i in 1:n
      z[i, i] = one(R)
   end
   z.base_ring = R
   return z
end

################################################################################
#
#  Matrix space constructor
#
################################################################################

function MatrixSpace(R::Generic.ResRing{fmpz}, r::Int, c::Int, cached::Bool = true)
  FmpzModMatSpace(R, r, c, cached)
end
