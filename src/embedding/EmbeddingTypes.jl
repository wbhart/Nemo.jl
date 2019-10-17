################################################################################
#
#  FinFieldsMorphism : some types needed to work with embeddings
#
################################################################################

struct FinFieldMorphism{T} <: AbstractAlgebra.Map{T, T, AbstractAlgebra.SetMap,
                                                  FinFieldMorphism} 
    map::AbstractAlgebra.FunctionalMap
    preimage::AbstractAlgebra.FunctionalMap

    function FinFieldMorphism(domain::T, codomain::T, image_fn::Function,
                              inverse_fn::Function) where T
        map = AbstractAlgebra.FunctionalMap(domain, codomain, image_fn)
        preimage = AbstractAlgebra.FunctionalMap(codomain, domain, inverse_fn)
        return FinFieldMorphism(map, preimage)
    end
end

domain(f::FinFieldMorphism) = (f.map).domain
codomain(f::FinFieldMorphism) = (f.map).codomain
image_fn(f::FinFieldMorphism) = (f.map).image_fn
inverse_fn(f::FinFieldMorphism) = (f.preimage).image_fn

function (f::FinFieldMorphism)(x)
    return image_fn(f)(x)::elem_type(codomain(f))
end

function Base.show(io::IO, f::FinFieldMorphism)
    print("Morphism from $(domain(f))\nto $(codomain(f))")
end

struct FinFieldPreimage{T} <: AbstractAlgebra.Map{T, T, AbstractAlgebra.SetMap,
                                                  FinFieldMorphism} where T
    map::AbstractAlgebra.FunctionalMap
    preimage::AbstractAlgebra.FunctionalMap

    function FinFieldPreimage(domain::T, codomain::T, image_fn::Function,
                              inverse_fn::Function)
        map = AbstractAlgebra.FunctionalMap(domain, codomain, image_fn)
        preimage = AbstractAlgebra.FunctionalMap(codomain, domain, inverse_fn)
        return FinFieldPreimage(map, preimage)
    end
end

domain(f::FinFieldPreimage) = (f.map).domain
codomain(f::FinFieldPreimage) = (f.map).codomain
image_fn(f::FinFieldPreimage) = (f.map).image_fn
inverse_fn(f::FinFieldPreimage) = (f.preimage).image_fn

function (f::FinFieldPreimage)(x)
    a = inverse_fn(f)(x)::elem_type(domain(f))
    b = image_fn(f)(a)
    if x == b
        return a
    else
        throw(ArgumentError(string("not an element in the subfield of degree ",
                                   degree(domain(f)), " over F_",
                                   characteristic(domain(f)))))
    end
end

function Base.show(io::IO, f::FinFieldPreimage)
    print("Preimage of the morphism from $(domain(f))\nto $(codomain(f))")
end

preimage(f::FinFieldMorphism) = FinFieldPreimage(domain(f), codomain(f),
                                                 image_fn(f), inverse_fn(f))
