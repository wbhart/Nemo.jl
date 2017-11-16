###############################################################################
#
#   fq_nmod_embed.jl : Flint finite fields embeddings
#
###############################################################################

export embed_matrices, embed

###############################################################################
#
#   Naive functions
#
###############################################################################

function embed_gens(k::FqNmodFiniteField, K::FqNmodFiniteField)
    a = k()
    b = K()
    p::Int = characteristic(k)
    R = ResidueRing(ZZ, p)
    PR = PolynomialRing(R, "T")[1]
    P = PR()

    ccall((:fq_nmod_embed_gens, :libflint), Void, (Ptr{fq_nmod}, Ptr{fq_nmod},
    Ptr{nmod_poly}, Ptr{FqNmodFiniteField}, Ptr{FqNmodFiniteField}), &a, &b,
    &P, &k, &K)

    return a, b, P
end

function embed_matrices(k::FqNmodFiniteField, K::FqNmodFiniteField)
    a, b, P = embed_gens(k, K)
    m, n = degree(k), degree(K)
    R = base_ring(P)
    S1 = MatrixSpace(R, n, m)
    S2 = MatrixSpace(R, m, n)
    s1 = S1()
    s2 = S2()

    ccall((:fq_nmod_embed_matrices, :libflint), Void, (Ptr{nmod_mat},
    Ptr{nmod_mat}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}, Ptr{fq_nmod},
    Ptr{FqNmodFiniteField}, Ptr{nmod_poly}), &s1, &s2, &a, &k, &b, &K, &P)
    return s1, s2
end

function embed_matrices_pre(a::fq_nmod, b::fq_nmod, P::nmod_poly)
    k = parent(a)
    K = parent(b)
    m, n = degree(k), degree(K)
    R = base_ring(P)
    S1 = MatrixSpace(R, n, m)
    S2 = MatrixSpace(R, m, n)
    s1 = S1()
    s2 = S2()

    ccall((:fq_nmod_embed_matrices, :libflint), Void, (Ptr{nmod_mat},
    Ptr{nmod_mat}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}, Ptr{fq_nmod},
    Ptr{FqNmodFiniteField}, Ptr{nmod_poly}), &s1, &s2, &a, &k, &b, &K, &P)
    return s1, s2
end

function coeff!(x::fq_nmod, j::Int, c::Int)
    ccall((:nmod_poly_set_coeff_ui, :libflint), Void,
          (Ptr{fq_nmod}, Int, UInt), &x, j, c)
end

function embed(x::fq_nmod, K::FqNmodFiniteField)

    k = parent(x)
    d = degree(k)
    M = embed_matrices(k, K)[1]
    S = MatrixSpace(base_ring(M), d, 1)
    col = S()

    for j in 0:(d-1)
        col[j+1, 1] = coeff(x, j)
    end

    product = M*col
    res = K()

    for j in degree(K):-1:1
        coeff!(res, j-1, Int(data(product[j, 1])))
    end
    return res
end

function embed_pre_mat(x::fq_nmod, K::FqNmodFiniteField, M::nmod_mat)

    d = degree(parent(x))
    S = MatrixSpace(base_ring(M), d, 1)
    col = S()

    for j in 0:(d-1)
        col[j+1, 1] = coeff(x, j)
    end

    product = M*col
    res = K()

    for j in degree(K):-1:1
        coeff!(res, j-1, Int(data(product[j, 1])))
    end

    return res
end

function embedPoly(P::fq_nmod_poly, f::FinFieldMorphism)
    S = PolynomialRing(field(codomain(f)), "T")[1]
    return S([f(coeff(P, j)) for j in 0:degree(P)])
end
