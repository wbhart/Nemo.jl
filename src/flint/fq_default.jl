###############################################################################
#
#   fq_default.jl : Flint finite fields
#
###############################################################################

export FqFiniteField, fq_default, FqDefaultFiniteField

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

parent_type(::Type{fq_default}) = FqDefaultFiniteField

elem_type(::Type{FqDefaultFiniteField}) = fq_default

@doc Markdown.doc"""
    base_ring(a::FqDefaultFiniteField)

Returns `Union{}` as this field is not dependent on another field.
"""
base_ring(a::FqDefaultFiniteField) = Union{}

@doc Markdown.doc"""
    base_ring(a::fq_default)

Returns `Union{}` as this field is not dependent on another field.
"""
base_ring(a::fq_default) = Union{}

@doc Markdown.doc"""
    parent(a::fq_default)

Returns the parent of the given finite field element.
"""
parent(a::fq_default) = a.parent

isdomain_type(::Type{fq_default}) = true

function check_parent(a::fq_default, b::fq_default)
   a.parent != b.parent && error("Operations on distinct finite fields not supported")
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::fq_default, h::UInt)
   b = 0xb310fb6ea97e1f1a%UInt
   for i in 0:degree(parent(a)) - 1
      b = xor(b, xor(hash(coeff(a, i), h), h))
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

@doc Markdown.doc"""
    coeff(x::fq_default, n::Int)

Return the degree $n$ coefficient of the polynomial representing the given
finite field element.
"""
function coeff(x::fq_default, n::Int)
   n < 0 && throw(DomainError(n, "Index must be non-negative"))
   z = fmpz()
   ccall((:fq_default_get_coeff_fmpz, libflint), Nothing,
               (Ref{fmpz}, Ref{fq_default}, Int, Ref{FqDefaultFiniteField}),
                                                           z, x, n, parent(x))
   return z
end

@doc Markdown.doc"""
    zero(a::FqDefaultFiniteField)

Return the additive identity, zero, in the given finite field.
"""
function zero(a::FqDefaultFiniteField)
   d = a()
   ccall((:fq_default_zero, libflint), Nothing, (Ref{fq_default}, Ref{FqDefaultFiniteField}), d, a)
   return d
end

@doc Markdown.doc"""
    one(a::FqDefaultFiniteField)

Return the multiplicative identity, one, in the given finite field.
"""
function one(a::FqDefaultFiniteField)
   d = a()
   ccall((:fq_default_one, libflint), Nothing, (Ref{fq_default}, Ref{FqDefaultFiniteField}), d, a)
   return d
end

@doc Markdown.doc"""
    gen(a::FqDefaultFiniteField)

Return the generator of the finite field. Note that this is only guaranteed
to be a multiplicative generator if the finite field is generated by a
Conway polynomial automatically.
"""
function gen(a::FqDefaultFiniteField)
   d = a()
   ccall((:fq_default_gen, libflint), Nothing, (Ref{fq_default}, Ref{FqDefaultFiniteField}), d, a)
   return d
end

@doc Markdown.doc"""
    iszero(a::fq_default)

Return `true` if the given finite field element is zero, otherwise return
`false`.
"""
iszero(a::fq_default) = ccall((:fq_default_is_zero, libflint), Bool,
                     (Ref{fq_default}, Ref{FqDefaultFiniteField}), a, a.parent)

@doc Markdown.doc"""
    isone(a::fq_default)

Return `true` if the given finite field element is one, otherwise return
`false`.
"""
isone(a::fq_default) = ccall((:fq_default_is_one, libflint), Bool,
                    (Ref{fq_default}, Ref{FqDefaultFiniteField}), a, a.parent)

@doc Markdown.doc"""
    isgen(a::fq_default)

Return `true` if the given finite field element is the generator of the
finite field, otherwise return `false`.
"""
isgen(a::fq_default) = a == gen(parent(a))

@doc Markdown.doc"""
    isunit(a::fq_default)

Return `true` if the given finite field element is invertible, i.e. nonzero,
otherwise return `false`.
"""
isunit(a::fq_default) = ccall((:fq_default_is_invertible, libflint), Bool,
                     (Ref{fq_default}, Ref{FqDefaultFiniteField}), a, a.parent)

@doc Markdown.doc"""
    characteristic(a::FqDefaultFiniteField)

Return the characteristic of the given finite field.
"""
function characteristic(a::FqDefaultFiniteField)
   d = fmpz()
   ccall((:fq_default_ctx_prime, libflint), Nothing,
         (Ref{fmpz}, Ref{FqDefaultFiniteField}), d, a)
   return d
end

@doc Markdown.doc"""
    order(a::FqDefaultFiniteField)

Return the order, i.e. the number of elements in, the given finite field.
"""
function order(a::FqDefaultFiniteField)
   d = fmpz()
   ccall((:fq_default_ctx_order, libflint), Nothing,
         (Ref{fmpz}, Ref{FqDefaultFiniteField}), d, a)
   return d
end

@doc Markdown.doc"""
    degree(a::FqDefaultFiniteField)

Return the degree of the given finite field.
"""
function degree(a::FqDefaultFiniteField)
   return ccall((:fq_default_ctx_degree, libflint), Int, (Ref{FqDefaultFiniteField},), a)
end

function deepcopy_internal(d::fq_default, dict::IdDict)
   z = fq_default(parent(d), d)
   return z
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(x::fq_default) = x

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function expressify(a::fq_default; context = nothing)
   x = a.parent.var
   d = degree(a.parent)

   sum = Expr(:call, :+)
   for k in (d - 1):-1:0
        c = coeff(a, k)
        if !iszero(c)
            xk = k < 1 ? 1 : k == 1 ? x : Expr(:call, :^, x, k)
            if isone(c)
                push!(sum.args, Expr(:call, :*, xk))
            else
                push!(sum.args, Expr(:call, :*, expressify(c, context = context), xk))
            end
        end
    end
    return sum
end

show(io::IO, a::fq_default) = print(io, AbstractAlgebra.obj_to_string(a, context = io))

function show(io::IO, a::FqDefaultFiniteField)
   print(io, "Finite field of degree ", degree(a))
   print(io, " over F_", characteristic(a))
end

###############################################################################
#
#   Unary operations
#
###############################################################################

function -(x::fq_default)
   z = parent(x)()
   ccall((:fq_default_neg, libflint), Nothing,
         (Ref{fq_default}, Ref{fq_default}, Ref{FqDefaultFiniteField}), z, x, x.parent)
   return z
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function +(x::fq_default, y::fq_default)
   check_parent(x, y)
   z = parent(y)()
   ccall((:fq_default_add, libflint), Nothing,
        (Ref{fq_default}, Ref{fq_default}, Ref{fq_default}, Ref{FqDefaultFiniteField}), z, x, y, y.parent)
   return z
end

function -(x::fq_default, y::fq_default)
   check_parent(x, y)
   z = parent(y)()
   ccall((:fq_default_sub, libflint), Nothing,
        (Ref{fq_default}, Ref{fq_default}, Ref{fq_default}, Ref{FqDefaultFiniteField}), z, x, y, y.parent)
   return z
end

function *(x::fq_default, y::fq_default)
   check_parent(x, y)
   z = parent(y)()
   ccall((:fq_default_mul, libflint), Nothing,
        (Ref{fq_default}, Ref{fq_default}, Ref{fq_default}, Ref{FqDefaultFiniteField}), z, x, y, y.parent)
   return z
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(x::Int, y::fq_default)
   z = parent(y)()
   ccall((:fq_default_mul_si, libflint), Nothing,
         (Ref{fq_default}, Ref{fq_default}, Int, Ref{FqDefaultFiniteField}), z, y, x, y.parent)
   return z
end

*(x::Integer, y::fq_default) = fmpz(x)*y

*(x::fq_default, y::Integer) = y*x

function *(x::fmpz, y::fq_default)
   z = parent(y)()
   ccall((:fq_default_mul_fmpz, libflint), Nothing,
         (Ref{fq_default}, Ref{fq_default}, Ref{fmpz}, Ref{FqDefaultFiniteField}),
                                            z, y, x, y.parent)
   return z
end

*(x::fq_default, y::fmpz) = y*x

+(x::fq_default, y::Integer) = x + parent(x)(y)

+(x::Integer, y::fq_default) = y + x

+(x::fq_default, y::fmpz) = x + parent(x)(y)

+(x::fmpz, y::fq_default) = y + x

-(x::fq_default, y::Integer) = x - parent(x)(y)

-(x::Integer, y::fq_default) = parent(y)(x) - y

-(x::fq_default, y::fmpz) = x - parent(x)(y)

-(x::fmpz, y::fq_default) = parent(y)(x) - y

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::fq_default, y::Int)
   if y < 0
      x = inv(x)
      y = -y
   end
   z = parent(x)()
   ccall((:fq_default_pow_ui, libflint), Nothing,
         (Ref{fq_default}, Ref{fq_default}, Int, Ref{FqDefaultFiniteField}), z, x, y, x.parent)
   return z
end

function ^(x::fq_default, y::fmpz)
   if y < 0
      x = inv(x)
      y = -y
   end
   z = parent(x)()
   ccall((:fq_default_pow, libflint), Nothing,
         (Ref{fq_default}, Ref{fq_default}, Ref{fmpz}, Ref{FqDefaultFiniteField}),
                                            z, x, y, x.parent)
   return z
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(x::fq_default, y::fq_default)
   check_parent(x, y)
   ccall((:fq_default_equal, libflint), Bool,
         (Ref{fq_default}, Ref{fq_default}, Ref{FqDefaultFiniteField}), x, y, y.parent)
end

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

==(x::fq_default, y::Integer) = x == parent(x)(y)

==(x::fq_default, y::fmpz) = x == parent(x)(y)

==(x::Integer, y::fq_default) = parent(y)(x) == y

==(x::fmpz, y::fq_default) = parent(y)(x) == y

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(x::fq_default, y::fq_default)
   check_parent(x, y)
   iszero(y) && throw(DivideError())
   z = parent(y)()
   ccall((:fq_default_div, libflint), Nothing,
        (Ref{fq_default}, Ref{fq_default}, Ref{fq_default}, Ref{FqDefaultFiniteField}), z, x, y, y.parent)
   return z
end

function divides(a::fq_default, b::fq_default)
   if iszero(a)
      return true, zero(parent(a))
   end
   if iszero(b)
      return false, zero(parent(a))
   end
   return true, divexact(a, b)
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

divexact(x::fq_default, y::Integer) = divexact(x, parent(x)(y))

divexact(x::fq_default, y::fmpz) = divexact(x, parent(x)(y))

divexact(x::Integer, y::fq_default) = divexact(parent(y)(x), y)

divexact(x::fmpz, y::fq_default) = divexact(parent(y)(x), y)

###############################################################################
#
#   Inversion
#
###############################################################################

@doc Markdown.doc"""
    inv(x::fq_default)

Return $x^{-1}$.
"""
function inv(x::fq_default)
   iszero(x) && throw(DivideError())
   z = parent(x)()
   ccall((:fq_default_inv, libflint), Nothing,
         (Ref{fq_default}, Ref{fq_default}, Ref{FqDefaultFiniteField}), z, x, x.parent)
   return z
end

###############################################################################
#
#   Special functions
#
###############################################################################

@doc Markdown.doc"""
    square_root(x::fq_default)

Return a square root of $x$ in the finite field. If $x$ is not a square
an exception is raised.
"""
function square_root(x::fq_default)
   z = parent(x)()
   res = Bool(ccall((:fq_default_sqrt, libflint), Cint,
                    (Ref{fq_default}, Ref{fq_default}, Ref{FqDefaultFiniteField}),
                    z, x, x.parent))
   res || error("Not a square in square_root")
   return z
end

@doc Markdown.doc"""
    issquare(x::fq_default)

Return `true` if $x$ is a square in the finite field (includes zero).
"""
function issquare(x::fq_default)
   return Bool(ccall((:fq_default_is_square, libflint), Cint,
                     (Ref{fq_default}, Ref{FqDefaultFiniteField}),
                     x, x.parent))
end

@doc Markdown.doc"""
    issquare_with_square_root(x::fq_default)

If $x$ is a square, return the tuple of `true` and a square root of $x$.
Otherwise, return the tuple of `false` and an unspecified element of the
finite field.
"""
function issquare_with_square_root(x::fq_default)
   z = parent(x)()
   flag = ccall((:fq_default_sqrt, libflint), Cint,
                (Ref{fq_default}, Ref{fq_default}, Ref{FqDefaultFiniteField}),
                z, x, x.parent)
   return (Bool(flag), z)
end

@doc Markdown.doc"""
    pth_root(x::fq_default)

Return the $p$-th root of $x$ in the finite field of characteristic $p$. This
is the inverse operation to the Frobenius map $\sigma_p$.
"""
function pth_root(x::fq_default)
   z = parent(x)()
   ccall((:fq_default_pth_root, libflint), Nothing,
         (Ref{fq_default}, Ref{fq_default}, Ref{FqDefaultFiniteField}), z, x, x.parent)
   return z
end

@doc Markdown.doc"""
    tr(x::fq_default)

Return the trace of $x$. This is an element of $\mathbb{F}_p$, but the value returned
is this value embedded in the original finite field.
"""
function tr(x::fq_default)
   z = fmpz()
   ccall((:fq_default_trace, libflint), Nothing,
         (Ref{fmpz}, Ref{fq_default}, Ref{FqDefaultFiniteField}), z, x, x.parent)
   return parent(x)(z)
end

@doc Markdown.doc"""
    norm(x::fq_default)

Return the norm of $x$. This is an element of $\mathbb{F}_p$, but the value returned
is this value embedded in the original finite field.
"""
function norm(x::fq_default)
   z = fmpz()
   ccall((:fq_default_norm, libflint), Nothing,
         (Ref{fmpz}, Ref{fq_default}, Ref{FqDefaultFiniteField}), z, x, x.parent)
   return parent(x)(z)
end

@doc Markdown.doc"""
    frobenius(x::fq_default, n = 1)

Return the iterated Frobenius $\sigma_p^n(x)$ where $\sigma_p$ is the
Frobenius map sending the element $a$ to $a^p$ in the finite field of
characteristic $p$. By default the Frobenius map is applied $n = 1$ times if
$n$ is not specified.
"""
function frobenius(x::fq_default, n = 1)
   z = parent(x)()
   ccall((:fq_default_frobenius, libflint), Nothing,
         (Ref{fq_default}, Ref{fq_default}, Int, Ref{FqDefaultFiniteField}), z, x, n, x.parent)
   return z
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(z::fq_default)
   ccall((:fq_default_zero, libflint), Nothing,
        (Ref{fq_default}, Ref{FqDefaultFiniteField}), z, z.parent)
   return z
end

function mul!(z::fq_default, x::fq_default, y::fq_default)
   ccall((:fq_default_mul, libflint), Nothing,
        (Ref{fq_default}, Ref{fq_default}, Ref{fq_default}, Ref{FqDefaultFiniteField}), z, x, y, y.parent)
   return z
end

function addeq!(z::fq_default, x::fq_default)
   ccall((:fq_default_add, libflint), Nothing,
        (Ref{fq_default}, Ref{fq_default}, Ref{fq_default}, Ref{FqDefaultFiniteField}), z, z, x, x.parent)
   return z
end

function add!(z::fq_default, x::fq_default, y::fq_default)
   ccall((:fq_default_add, libflint), Nothing,
        (Ref{fq_default}, Ref{fq_default}, Ref{fq_default}, Ref{FqDefaultFiniteField}), z, x, y, x.parent)
   return z
end

###############################################################################
#
#   Promotions
#
###############################################################################

promote_rule(::Type{fq_default}, ::Type{T}) where {T <: Integer} = fq_default

promote_rule(::Type{fq_default}, ::Type{fmpz}) = fq_default

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::FqDefaultFiniteField)()
   z = fq_default(a)
   return z
end

(a::FqDefaultFiniteField)(b::Integer) = a(fmpz(b))

function (a::FqDefaultFiniteField)(b::Int)
   z = fq_default(a, b)
   return z
end

function (a::FqDefaultFiniteField)(b::fmpz)
   z = fq_default(a, b)
   return z
end

function (a::FqDefaultFiniteField)(b::fq_default)
   parent(b) != a && error("Coercion between finite fields not implemented")
   return b
end

###############################################################################
#
#   FlintFiniteField constructor
#
###############################################################################

@doc Markdown.doc"""
    NGFiniteField(char::Union{fmpz, Integer}, deg::Int, s::AbstractString; cached = true)

Returns a tuple $S, x$ consisting of a finite field parent object $S$ and
generator $x$ for the finite field of the given characteristic and degree.
The string $s$ is used to designate how the finite field generator will be
printed. The characteristic must be prime. When a Conway polynomial is known,
the field is generated using the Conway polynomial. Otherwise a random
sparse, irreducible polynomial is used. The generator of the field is
guaranteed to be a multiplicative generator only if the field is generated by
a Conway polynomial. We require the degree to be positive.
"""
function NGFiniteField(char::Union{fmpz, Integer}, deg::Int, s::AbstractString; cached = true)
   char = fmpz(char)
   S = Symbol(s)
   parent_obj = FqDefaultFiniteField(char, deg, S, cached)

   return parent_obj, gen(parent_obj)
end

@doc Markdown.doc"""
    NGFiniteField(pol::Union{fmpz_mod_poly, gfp_fmpz_poly}, s::AbstractString; cached = true, check = true)

Returns a tuple $S, x$ consisting of a finite field parent object $S$ and
generator $x$ for the finite field over $F_p$ defined by the given
polynomial, i.e. $\mathbb{F}_p[t]/(pol)$. The characteristic is specified by
the modulus of `pol`. The polynomial is required to be irreducible, but this
is not checked. The base ring of the polynomial is required to be a field, which
is checked by default. Use `check = false` to disable the check.
The string $s$ is used to designate how the finite field
generator will be printed. The generator will not be multiplicative in
general.
"""
function NGFiniteField(pol::Union{fmpz_mod_poly, gfp_fmpz_poly},
                          s::AbstractString; cached = true, check::Bool=true)
   S = Symbol(s)
   parent_obj = FqDefaultFiniteField(pol, S, cached, check=check)

   return parent_obj, gen(parent_obj)
end
