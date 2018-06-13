###############################################################################
#
#   fq_embed.jl : Flint finite fields embeddings
#
###############################################################################

###############################################################################
#
#   Linear factor
#
###############################################################################

function linfactor(x::fq_poly)
    y = parent(x)()
    ccall((:fq_poly_linfactor, :libflint), Void, (Ref{fq_poly},
          Ref{fq_poly}, Ref{FqFiniteField}), y, x, base_ring(x))
    return y
end

###############################################################################
#
#   Naive functions
#
###############################################################################

function embed_gens(k::FqFiniteField, K::FqFiniteField)
    a = k()
    b = K()
    p::fmpz = characteristic(k)
    R = ResidueRing(ZZ, p)
    PR = PolynomialRing(R, "T")[1]
    P = PR()

    ccall((:fq_embed_gens, :libflint), Void, (Ref{fq}, Ref{fq},
    Ref{fmpz_mod_poly}, Ref{FqFiniteField}, Ref{FqFiniteField}), a, b,
    P, k, K)

    return a, b, P
end

function embed_matrices(k::FqFiniteField, K::FqFiniteField)

    m, n = degree(k), degree(K)
    if m == n
        T1, T2 = modulus(k), modulus(K)
        if T1 == T2
            s1 = identity_matrix(base_ring(T1), n)
            s2 = s1
            return s1, s2
        end
    end

    a, b, P = embed_gens(k, K)
    R = base_ring(P)
    S1 = MatrixSpace(R, n, m)
    S2 = MatrixSpace(R, m, n)
    s1 = S1()
    s2 = S2()

    ccall((:fq_embed_matrices, :libflint), Void, (Ref{fmpz_mod_mat},
    Ref{fmpz_mod_mat}, Ref{fq}, Ref{FqFiniteField}, Ref{fq},
    Ref{FqFiniteField}, Ref{fmpz_mod_poly}), s1, s2, a, k, b, K, P)
    return s1, s2
end

function embed_matrices_pre(a::fq, b::fq, P::fmpz_mod_poly)
    k = parent(a)
    K = parent(b)
    m, n = degree(k), degree(K)
    R = base_ring(P)
    S1 = MatrixSpace(R, n, m)
    S2 = MatrixSpace(R, m, n)
    s1 = S1()
    s2 = S2()

    ccall((:fq_embed_matrices, :libflint), Void, (Ref{fmpz_mod_mat},
    Ref{fmpz_mod_mat}, Ref{fq}, Ref{FqFiniteField}, Ref{fq},
    Ref{FqFiniteField}, Ref{fmpz_mod_poly}), s1, s2, a, k, b, K, P)
    return s1, s2
end

function coeff!(x::fq, j::Int, c::Int)
    ccall((:fmpz_mod_poly_set_coeff_ui, :libflint), Void,
          (Ref{fq_nmod}, Int, UInt), x, j, c)
end

function embed_pre_mat(x::fq, K::FqFiniteField, M::fmpz_mod_mat)

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
