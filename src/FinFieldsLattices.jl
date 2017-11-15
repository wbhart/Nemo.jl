################################################################################
#
#   FinFieldsLattices.jl : Finite Fields Lattices
#
################################################################################

export FieldNode, fieldNode

################################################################################
#
#   Base type : nodes
#
################################################################################

mutable struct FieldNode{T <: FinField}
    field::T
    overfields::Dict{Int, Array{FinFieldMorphism{T}, 1}}
    subfields::Dict{Int, Array{FinFieldMorphism{T}, 1}}
end

function fieldNode{T <: FinField}(k::T)
    overfields = Dict{Int, Array{FinFieldMorphism{T}, 1}}()
    subfields = Dict{Int, Array{FinFieldMorphism{T}, 1}}()

    return FieldNode(k, overfields, subfields)
end

field(f::FieldNode) = f.field
overfields(f::FieldNode) = f.overfields
subfields(f::FieldNode) = f.subfields

function addOverfield!{T <: FinField}(F::FieldNode{T},
                                   f::FinFieldMorphism{T})

    d = degree(codomain(f))
    over = overfields(F)

    if haskey(over, d)
        push!(over[d], f)
    else
        a = FinFieldMorphism{T}[f]
        over[d] = a
    end
end

function addSubfield!{T <: FinField}(F::FieldNode{T},
                                  f::FinFieldMorphism{T})

    d = degree(codomain(f))
    sub = subfields(F)

    if haskey(sub, d)
        push!(sub[d], f)
    else
        a = FinFieldMorphism{T}[f]
        sub[d] = a
    end
end

function Base.show(io::IO, k::FieldNode)  
    print(io, "Node composed of Finite Field of degree ", degree(k.field))
    print(io, " over F_", characteristic(k.field))
end

################################################################################
#
#  Modulus 
#
################################################################################

function modulus(k::FqNmodFiniteField)
    p::Int = characteristic(k)
    R = PolynomialRing(ResidueRing(ZZ, p), "T")[1]
    Q = R()
    P = ccall((:fq_nmod_ctx_modulus, :libflint), Ptr{nmod_poly},
                 (Ptr{FqNmodFiniteField},), &k)
    ccall((:nmod_poly_set, :libflint), Void,
          (Ptr{nmod_poly}, Ptr{nmod_poly}),
          &Q, P)

    return Q
end

################################################################################
#
#  Root Finding 
#
################################################################################

function anyRoot(Q::PolyElem)
    fact = factor(Q)
    for f in fact
        return -coeff(f[1], 0)
    end
end

################################################################################
#
#   Computing the defining polynomial
#
################################################################################

function coordinates_mat(f::FinFieldMorphism)
    E = domain(f)
    F = codomain(f)
    n = degree(E)
    m = degree(F)
    d::Int = m/n
    p::Int = characteristic(F)
    R = ResidueRing(ZZ, p)
    S = MatrixSpace(R, m, m)
    t = f(gen(E))
    u = gen(F)
    M = S()

    for i in 0:(d-1)
        for j in 0:(n-1)
            tmp = u^i*t^j
            for k in 0:(m-1)
                M[k+1, i*n+j+1] = coeff(tmp, k)
            end
        end
    end
    return M
end

function coordinates(x::fq_nmod, f::FinFieldMorphism)
    M = inv(coordinates_mat(f))
    R = base_ring(M)
    E = domain(f)
    F = codomain(f)
    n = degree(E)
    m = degree(F)
    d::Int = m/n
    S = MatrixSpace(R, m, 1)
    col = S()
    for i in 0:(m-1)
        col[i+1, 1] = coeff(x, i)
    end
    product = M*col
    res = Array{fq_nmod}(d)
    for i in 0:(d-1)
        tmp = E()
        for j in 0:(n-1)
            coeff!(tmp, j, Int(data(product[i*n+j+1, 1])))
        end
        res[i+1] = tmp
    end
    return res
end

function coordinates_pre(x::fq_nmod, f::FinFieldMorphism, M::nmod_mat)
    R = base_ring(M)
    E = domain(f)
    F = codomain(f)
    n = degree(E)
    m = degree(F)
    d::Int = m/n
    S = MatrixSpace(R, m, 1)
    col = S()
    for i in 0:(m-1)
        col[i+1, 1] = coeff(x, i)
    end
    product = M*col
    res = Array{fq_nmod}(d)
    for i in 0:(d-1)
        tmp = E()
        for j in 0:(n-1)
            coeff!(tmp, j, Int(data(product[i*n+j+1, 1])))
        end
        res[i+1] = tmp
    end
    return res
end

function defining_poly(f::FinFieldMorphism)

    E = domain(f)
    F = codomain(f)
    d::Int = degree(F)/degree(E)

    tmp = coordinates(gen(F)^d, f)

    R, T = PolynomialRing(E, "T")

    res = R(tmp)

    return T^d-res
end

function embed(k::FieldNode, K::FieldNode)

    E, F = field(k), field(K)
    S = PolynomialRing(F, "T")[1]
    Q = S()
    needy = false
    m, n = degree(E), degree(F)

    for l in keys(subfields(k))
        if haskey(subfields(K), l)
            f = inv(subfields(k)[l][1])
            g = inv(subfields(K)[l][1])
            P = embedPoly(defining_poly(f), g)
            if needy
                Q = gcd(Q, P)
            else
                Q = P
            end
            needy = true
        end
    end

    if needy
        defPol = modulus(E)
        t = anyRoot(Q)
        M, N = embed_matrices_pre(gen(E), t, defPol)
        f(x) = embed_pre_mat(x, F, M)
        g(y) = embed_pre_mat(y, E, N)
        morph = FinFieldMorphism(E, F, f, g)
        addOverfield!(k, morph)
        addSubfield!(K, inv(morph))
    else
        morph = embed(E, F)
        addOverfield!(k, morph)
        addSubfield!(K, inv(morph))
    end
    return morph
end
