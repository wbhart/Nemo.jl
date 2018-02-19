################################################################################
#
#   Fields.jl : generic fields
#
################################################################################

include("flint/fmpq.jl")

include("flint/fq.jl")

include("flint/fq_nmod.jl")

include("antic/nf_elem.jl")

include("arb/arb.jl")

include("arb/acb.jl")

//(a::T, b::T) where {T <: FieldElem} = divexact(a, b)

//(x::T, y::Union{Integer, Rational}) where {T <: RingElem} = x//parent(x)(y)
                                          
//(x::Union{Integer, Rational}, y::T) where {T <: RingElem} = parent(y)(x)//y

characteristic(R::Field) = 0


