################################################################################
#
#  FinFieldsMorphism : some types needed to work with embeddings
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

struct FinFieldPreimage{T}
    domain::T
    codomain::T
    f::Function
    inv::Function
end

domain(f::FinFieldPreimage) = f.domain

codomain(f::FinFieldPreimage) = f.codomain

function (f::FinFieldPreimage)(x)
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

function Base.show(io::IO, f::FinFieldPreimage)
    print("Preimage from $(domain(f))\nto $(codomain(f))")
end

preimage(f::FinFieldMorphism) = FinFieldPreimage(f.codomain, f.domain, f.inv, f.f)
