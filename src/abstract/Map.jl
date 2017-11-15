export Map, domain, codomain, endomorphism, GenMap, fun

################################################################################
#
#  Abstract map type
#
################################################################################

abstract type Map{D, C} end

################################################################################
#
#  Generic composition type
#
################################################################################

# If f is a map f : D -> C and g a map g : C -> B, then
# Then Comp(g, f) corresponds to the map D -> B, x -> g(f(x))
# So the first field corresponds to the first map you apply

mutable struct Comp{M, N, D, B} <: Map{D, B}
  second::N
  first::M
end

Comp{D, C, B}(f::Map{D, C}, g::Map{C, B}) =
    Comp{typeof(f), typeof(g), D, B}(g, f)

# implement the map interface

domain(f::Comp) = domain(f.first)

codomain(f::Comp) = codomain(f.second)

(f::Comp{M, N}){M, N}(x) = f.second(f.first(x))

# printing

function Base.show(io::IO, f::Comp)
  print(io, "Map from $(domain(f)) to $(codomain(f))\n")
  print(io, "Composition of $(f.second)\nand\n$(f.first)\n")
end

################################################################################
#
#  Homogenous composition type
#
################################################################################

mutable struct HomgComp{M, D} <: Map{D, D}
  maps::Array{M, 1}
end

function HomgComp{M}(f::Array{M, 1})
  length(f) == 0 && error("Number of maps must be positive")
  return HomgComp{M, typeof(domain(f))}(f)
end

HomgComp{M}(f::M, g::M) = HomgComp{M, typeof(domain(f))}([f, g])

map_type{M, D}(::HomgComp{M, D}) = M

Base.length(f::HomgComp) = length(f.maps)

# implement the map interface

domain(f::HomgComp) = domain(f.maps[end])

codomain(f::HomgComp) = codomain(f.maps[1])

function (f::HomgComp{M, D}){M, D}(x)
  for i in length(f.maps):-1:1
    x = f.maps[i](x)
  end
  return x
end

# printing

function Base.show(io::IO, f::HomgComp)
  print(io, "Map from $(domain(f)) to $(codomain(f))\n")
  print(io, "Composition of\n")
  for i in length(f):-1:1
    println(io, f.maps[i])
  end
end

################################################################################
#
#  Generic maps
#
################################################################################

# This is not optimal with the abstract types.
# If you parametrize, we cannot have partial intialization without the inverse.
# Another option would be too use something like FunctionsWrapper.jl

mutable struct GenMap{D, C} <: Map{D, C}
  domain::D
  codomain::C
  f::Function
  inv::Function

  GenMap{D, C}(R::D, S::C, f::Function, inv::Function) where {D, C} = new(R, S, f, inv)

  function GenMap{D, C}(R::D, S::C, f::Function) where {D, C}
    z = new()
    z.domain = R
    z.codomain = S
    z.f = f
    return z
  end
end

GenMap{D, C}(R::D, S::C, f::Function, inv::Function) = GenMap{D, C}(R, S, f, inv)

GenMap{D, C}(R::D, S::C, f::Function) = GenMap{D, C}(R, S, f)

fun{D, C}(R::D, S::C, f::Function) = GenMap(R, S, f)

# implement the map interface

domain(f::GenMap) = f.domain

codomain(f::GenMap) = f.codomain

function (f::GenMap)(x)
  return f.f(x)::elem_type(domain(f))
end

# printing

function Base.show(io::IO, f::GenMap)
  print(io, "Map from $(domain(f)) to $(codomain(f))\n")
  print(io, "given by julia function")
end

################################################################################
#
#  Composition
#
################################################################################

# Idea: Try to produce as many homogenous compositions as possible
# (But don't split compositions)

function Base.:*{D, C, B}(f::Map{D, C}, g::Map{C, B})
  if typeof(f) <: HomgComp
    if typeof(g) <: HomgComp
      if map_type(f) == map_type(g)
        return typeof(f)(vcat(f.maps, g.maps))
      else
        return Comp(f, g)
      end
    elseif typeof(g) == map_type(f)
      return typeof(f)(vcat(f.maps, g))
    else
      return Comp(f, g)
    end
  elseif typeof(g) <: HomgComp
    if typeof(f) == map_type(g)
      return typeof(g)(vcat(g.maps, f))
    else
      return Comp(f, g)
    end
  elseif typeof(f) == typeof(g)
    return HomgComp(f, g)
  else
    return Comp(f, g)
  end
end

################################################################################
#
#  Examples
#
################################################################################

# Example 1
# Endomorphisms of polynomial rings

mutable struct PolyRingToPolyRing{D, T} <: Map{D, D}
  domain::D
  codomain::D
  a::T # image of the generator of the domain
end

function endomorphism{T <: PolyRing}(R::T, a)
  return PolyRingToPolyRing{T, elem_type(T)}(R, R, R(a))
end

domain(f::PolyRingToPolyRing) = f.domain

codomain(f::PolyRingToPolyRing) = f.codomain

function Base.show(io::IO, f::PolyRingToPolyRing)
  print("Map from \n$(domain(f))\nto\n$(codomain(f))\n")
  print("mapping $(gen(domain(f))) to $(f.a)")
end

function (f::PolyRingToPolyRing)(x)
  parent(x) != domain(f) && error("Element not in domain")
  return subst(x, f.a)
end

#=
julia> Zx, x = PolynomialRing(ZZ,"x");

julia> f = endomorphism(Zx, x^2)
Map from 
Univariate Polynomial Ring in x over Integer Ring
to
Univariate Polynomial Ring in x over Integer Ring
mapping x to x^2

julia> domain(f) == Zx
true

julia> g = endomorphism(Zx, x + 1)
Map from 
Univariate Polynomial Ring in x over Integer Ring
to
Univariate Polynomial Ring in x over Integer Ring
mapping x to x+1

julia> h = f * g
Map from Univariate Polynomial Ring in x over Integer Ring to Univariate Polynomial Ring in x over Integer Ring
Composition of
Map from 
Univariate Polynomial Ring in x over Integer Ring
to
Univariate Polynomial Ring in x over Integer Ring
mapping x to x+1
Map from 
Univariate Polynomial Ring in x over Integer Ring
to
Univariate Polynomial Ring in x over Integer Ring
mapping x to x^2

julia> h(x)
x^2+1

=#

# Example 2
# Canonical injection of ZZ into any commutative unitary ring

mutable struct ZtoRing{D} <: Map{FlintIntegerRing, D}
  codomain::D
end

canonical_map{T}(::FlintIntegerRing, R::T) = ZtoRing{T}(R)

domain(::ZtoRing) = FlintZZ

codomain(f::ZtoRing) = f.codomain

(f::ZtoRing)(x::fmpz) = (codomain(f))(x)

(f::ZtoRing)(x::Integer) = (codomain(f))(x)

function Base.show(io::IO, f::ZtoRing)
  print("Canonical morphism from ZZ to\n$(codomain(f))")
end

################################################################################
#
#  Finite fields morphism
#
################################################################################

struct FinFieldMorphism{T} <: Map{T, T}
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

inv(f::FinFieldMorphism) = FinFieldMorphism(f.codomain, f.domain, f.inv, f.f)

function embed{T <: FinField}(k::T, K::T)
    M, N = embed_matrices(k, K)
    f(x) = embed_pre_mat(x, K, M)
    inv(y) = embed_pre_mat(y, k, N)
    return FinFieldMorphism(k, K, f, inv)
end

function project{T <: FinField}(K::T, k::T)
    M, N = embed_matrices(k, K)
    f(x) = embed_pre_mat(x, k, N)
    inv(y) = embed_pre_mat(y, K, M)
    return FinFieldMorphism(K, k, f, inv)
end
