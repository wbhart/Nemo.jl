################################################################################
#
#   FinFieldsLattices.jl : Finite Fields Lattices
#
################################################################################

export embed, section

################################################################################
#
#  Over/Sub fields
#
################################################################################

overfields(k::FqNmodFiniteField) = k.overfields
subfields(k::FqNmodFiniteField) = k.subfields

"""
    addOverfield!{T <: FinField}(F::T,
                                 f::FinFieldMorphism{T})

Add an overfield to `F`, represented by a morphism ``f: F->codomain(f)``.
"""
function addOverfield!{T <: FinField}(F::T,
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

"""
    addSubfield!{T <: FinField}(F::T,
                                f::FinFieldMorphism{T})

Add a subfield to `F`, represented by a morphism ``f: domain(f)->F``.
"""
function addSubfield!{T <: FinField}(F::T,
                                     f::FinFieldMorphism{T})

    d = degree(domain(f))
    sub = subfields(F)

    if haskey(sub, d)
        push!(sub[d], f)
    else
        a = FinFieldMorphism{T}[f]
        sub[d] = a
    end
end

################################################################################
#
#  Modulus 
#
################################################################################

"""
    modulus(k::FqNmodFiniteField)

Return the modulus defining the finite field `k`.
"""
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
#  Random monic and splitting polynomials
#
################################################################################

function rand(R::PolyRing, d::Int)
    x = gen(R)
    k = base_ring(R)
    r = one(R)
    for i in 1:d
        r *= x + rand(k)
    end
    return r
end

################################################################################
#
#  Root Finding 
#
################################################################################

anyRoot(x::PolyElem) = -coeff(linfactor(x), 0)

################################################################################
#
# Minimal polynomial of the generator
#
################################################################################

"""
    berlekampMassey{Y <: RingElem}(a::Array{Y,1}, n::Int)

Compute the minimal polynomial of a linear recurring sequence.
"""
function berlekampMassey{Y <: FieldElem}(a::Array{Y, 1}, n::Int)

  S, T = PolynomialRing(parent(a[1]), "T")
  m = 2*n-1
  R0 = T^(2*n)
  R1 = S(reverse(a))
  V0 = S(0)
  V1 = S(1)

  while n <= degree(R1)
    Q, R = divrem(R0,R1)
    V = V0-Q*V1
    V0 = V1
    V1 = V
    R0 = R1
    R1 = R
  end

  return V1*lead(V1)^(-1)
end

"""
    genminpoly(f::FinFieldMorphism)

Compute the minimal polynomial of the generator of the codomain 
of `f` over the domain of `f`.
"""
function genminpoly(f::FinFieldMorphism)

    E = domain(f)
    F = codomain(f)
    d::Int = degree(F)/degree(E)
    b = 2*d

    # We define `sec`, the section of the morphism `f`
    
    sec = f.inv
    x = gen(F)
    y = one(F)

    # We compute the recurring sequence sec(x^j) for j from 0 to 2d-1

    A = Array{typeof(x)}(b)

    for j in 1:(b-1)
        A[j] = sec(y)
        y *= x
    end
    A[b] = sec(y)

    # We apply Berlekamp Massey algorithm to the sequence

    return berlekampMassey(A, d)
end

################################################################################
#
#   Embedding a polynomial
#
################################################################################

function embedPoly(P::fq_nmod_poly, f::FinFieldMorphism)
    S = PolynomialRing(codomain(f), "T")[1]
    return S([f(coeff(P, j)) for j in 0:degree(P)])
end

################################################################################
#
#  Embedding 
#
################################################################################

"""
    is_embedded{T <: FinField}(k::T, K::T)

If `k` is embbeded in `K`, return the corresponding embedding.
"""
function is_embedded{T <: FinField}(k::T, K::T)

    d = degree(K)
    ov = overfields(k)

    # We look for an embedding that has k as domain and K as codomain

    if haskey(ov, d)
        for f in ov[d]
            if domain(f) == k && codomain(f) == K
                return f
            end
        end
    end
end

"""
    embed_no_cond{T <: FinField}(k::T, K::T)

Embed `k` in `K` without worrying about compatibility conditions.
"""
function embed_no_cond{T <: FinField}(k::T, K::T)

    # We call the Flint algorithms directly, currently this is based on
    # factorization

    M, N = embed_matrices(k, K)
    f(x) = embed_pre_mat(x, K, M)
    inv(y) = embed_pre_mat(y, k, N)

    return FinFieldMorphism(k, K, f, inv)
end

"""
    find_morph{T <: FinField}(k::T, K::T)

Return a compatible embedding from `k` to `K`.
"""
function find_morph{T <: FinField}(k::T, K::T)

    S = PolynomialRing(K, "T")[1]
    Q = S()
    needy = false
    m, n = degree(k), degree(K)

    # For each common subfield S of k and K, we compute the minimal polynomial
    # of the canonical generator of k over S, with coefficient seen in K and we
    # compute the gcd of all these polynomials 

    for l in keys(subfields(k))
        if haskey(subfields(K), l)
            f = subfields(k)[l][1]
            g = subfields(K)[l][1]
            P = embedPoly(genminpoly(f), g)
            if needy
                Q = gcd(Q, P)
            else
                Q = P
            end
            needy = true
        end
    end

    # If there is at least one common subfield, we define the embedding from k
    # to K by sending the canonical generator of k to a root of the gcd
    # computed above

    if needy
        defPol = modulus(k)
        t = anyRoot(Q)
        M, N = embed_matrices_pre(gen(k), t, defPol)
        f(x) = embed_pre_mat(x, K, M)
        g(y) = embed_pre_mat(y, k, N)
        morph = FinFieldMorphism(k, K, f, g)

    # If there is no common subfield, there is no compatibility condition to
    # fulfill
    else
        morph = embed_no_cond(k, K)
    end

    return morph
end

"""
    transitive_closure(f::FinFieldMorphism)

Compute the transitive closure.
"""
function transitive_closure(f::FinFieldMorphism)

    k = domain(f)
    K = codomain(f)

    # Subfields

    subk = subfields(k)
    subK = subfields(K)

    # We go through all subfields of k and check if they are also subfields of
    # K, we add them if they are not

    for d in keys(subk)
        if !haskey(subK, d)
            for g in subk[d]
                t(y) = f(g(y))
                tinv(x) = g.inv(f.inv(x))
                phi = FinFieldMorphism(domain(g), K, t, tinv)

                addSubfield!(K, phi)
                addOverfield!(domain(g), phi)
            end
        else
            val = FqNmodFiniteField[codomain(v) for v in subK[d]]
            
            for g in subk[d]
                if !(domain(g) in val)
                    t(y) = f(g(y))
                    tinv(x) = g.inv(f.inv(x))
                    phi = FinFieldMorphism(domain(g), K, t, tinv)

                    addSubfield!(K, phi)
                    addOverfield!(domain(g), phi)
               end
            end
        end
    end

    # Overfields

    ov = overfields(K)

    # We call the same procedure on the overfields

    for d in keys(ov)
        for g in ov[d]
            transitive_closure(g)
        end
    end
end

"""
    intersections{T <: FinField}(k::T, K::T)

For each subfield S of K, embed I in S and k, where I is the intersection
between S and k.
"""
function intersections{T <: FinField}(k::T, K::T)

    d = degree(k)
    subk = subfields(k)
    subK = subfields(K)
    needmore = true

    # We loop through the subfields of K and we have different cases
    for l in keys(subK)
        c = gcd(d, l)

        # The intersection may be trivial, I = F_p
        if c == 1

        # I = k, so k is a subfield of S and we embed k in S
        # In that case, we finally have the embeddings k in S and S in K, so by
        # transitive closure we have k in K and we do not need more work
        elseif c == d
            for g in subK[l]
                embed(k, domain(g))
            end
            needmore = false

        # I = S, so S is a subfield of k, and we embed S in k
        elseif c == l
            for g in subK[l]
                embed(domain(g), k)
            end

        # If I is a subfield of k, we embed it in S
        elseif haskey(subk, c)
            L = domain(subk[c][1])
            for h in subK[l]
                embed(L, domain(h))
            end

        # If I is a subfield of K, we embed it in k and S
        elseif haskey(subK, c)
            L = domain(subK[c][1])
            embed(L, k)
            for h in subK[l]
                embed(L, domain(h))
            end

        # Otherwise it means that there is no field I around so we create one
        # and we embed it in k and S
        else
            p::Int = characteristic(k)
            kc, xc = FiniteField(p, c, string("x", c))
            embed(kc, k)
            for g in subK[l]
                embed(kc, domain(g))
            end
        end
    end

    # We return a boolean to tell if some more work needs to be done
    return needmore
end

"""
    embed{T <: FinField}(k::T, K::T)

Embed `k` in `K`, with some additionnal computations in order to satisfy
compatibility conditions with previous and future embeddings.
"""
function embed{T <: FinField}(k::T, K::T)

    # If k == K we return the identity

    if k == K
        identity(x) = x
        morph = FinFieldMorphism(k, k, identity, identity)
        return morph
    end

    # If k is already embedded in K, we return the corresponding embedding

    tr = is_embedded(k, K)
    if tr != nothing
        return tr
    end

    # Prior to embed k in K, we compute the needed embeddings
    # The embedding might be computed during the process !

    needmore = intersections(k, K) # reccursive calls to embed

    # And, if the wanted embeddings has not been computed during the process, we
    # finally compute a compatible embedding

    if needmore 

        # We compute a compatible embedding
        morph = find_morph(k, K)

        # We had it to the over and sub fields of k and K
        addOverfield!(k, morph)
        addSubfield!(K, morph)

        # We compute the transitive closure induced by the new embedding
        transitive_closure(morph)

        # And return the embedding
        return morph
    else

        # If the embedding has already been computing, we return it
        return is_embedded(k, K)
    end
end

################################################################################
#
#   Sections
#
################################################################################

"""
    section{T <: FinField}(K::T, k::T)

Compute a section of the embedding of `k` into `K`.
"""
function section{T <: FinField}(K::T, k::T)
    f = embed(k, K)
    return section(f)
end
