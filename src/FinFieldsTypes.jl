################################################################################
#
#  FinFieldsTypes.jl : some types needed to play with embeddings
#
################################################################################

################################################################################
#
#  FinFieldsMorphism : embeddings
#
################################################################################

struct FinFieldMorphism{T}
    domain::T
    codomain::T
    f::Function
    inv::Function
end

domain(f::FinFieldMorphism) = f.domain

codomain(f::FinFieldMorphism) = f.codomain

function (f::FinFieldMorphism)(x)
    return f.f(x)::elem_type(codomain(f))
end

function Base.show(io::IO, f::FinFieldMorphism)
    print("Morphism from $(domain(f))\nto $(codomain(f))")
end

struct FinFieldSection{T}
    domain::T
    codomain::T
    f::Function
    inv::Function
end

domain(f::FinFieldSection) = f.domain

codomain(f::FinFieldSection) = f.codomain

function (f::FinFieldSection)(x)
    a = f.f(x)::elem_type(codomain(f))
    b = f.inv(a)
    if x == b
        return a
    else
        throw(ArgumentError(string("not an element in the subfield of degree ",
                                   degree(codomain(f)), " over F_",
                                   characteristic(codomain(f)))))
    end
end

function Base.show(io::IO, f::FinFieldSection)
    print("Section from $(domain(f))\nto $(codomain(f))")
end

section(f::FinFieldMorphism) = FinFieldSection(f.codomain, f.domain, f.inv, f.f)
