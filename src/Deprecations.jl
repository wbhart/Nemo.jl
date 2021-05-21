# Deprecated in 0.16.*

@deprecate ArbField(x::Int, y::Bool) ArbField(x, cached = y)

@deprecate AcbField(x::Int, y::Bool) AcbField(x, cached = y)

@deprecate prec(x::AcbMatSpace) precision(x)

@deprecate prec(x::ArbMatSpace) precision(x)

@deprecate prec(x::AcbPolyRing) precision(x)

@deprecate prec(x::ArbPolyRing) precision(x)

@deprecate prec(x::AcbField) precision(x)

@deprecate prec(x::ArbField) precision(x)

# Deprecated in 0.22.*

@deprecate binom(x::arb, n::UInt) binomial(x, n)

@deprecate binom(n::UInt, k::UInt, r::ArbField) binomial(n, k, r)

# Deprecated in 0.23.*

@deprecate modeta(x::acb) modular_eta(x)

@deprecate modweber_f(x::acb) modular_weber_f(x)

@deprecate modweber_f1(x::acb) modular_weber_f1(x)

@deprecate modweber_f2(x::acb) modular_weber_f2(x)

@deprecate modj(x::acb) modular_j(x)

@deprecate modlambda(x::acb) modular_lambda(x)

@deprecate moddelta(x::acb) modular_delta(x)
