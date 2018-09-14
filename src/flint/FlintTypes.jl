###############################################################################
#
#   FlintTypes.jl : Parent and object types for Flint
#
###############################################################################

###############################################################################
#
#   FlintIntegerRing / fmpz
#
###############################################################################

mutable struct FlintIntegerRing <: Ring
end

const FlintZZ = FlintIntegerRing()

mutable struct fmpz <: RingElem
    d::Int

    function fmpz()
        z = new()
        ccall((:fmpz_init, :libflint), Nothing, (Ref{fmpz},), z)
        finalizer(_fmpz_clear_fn, z)
        return z
    end

    function fmpz(x::Int)
        z = new()
        ccall((:fmpz_init_set_si, :libflint), Nothing, (Ref{fmpz}, Int), z, x)
        finalizer(_fmpz_clear_fn, z)
        return z
    end

    function fmpz(x::UInt)
        z = new()
        ccall((:fmpz_init_set_ui, :libflint), Nothing, (Ref{fmpz}, UInt), z, x)
        finalizer(_fmpz_clear_fn, z)
        return z
    end

    function fmpz(x::BigInt)
        z = new()
        ccall((:fmpz_init, :libflint), Nothing, (Ref{fmpz},), z)
        ccall((:fmpz_set_mpz, :libflint), Nothing, (Ref{fmpz}, Ref{BigInt}), z, x)
        finalizer(_fmpz_clear_fn, z)
        return z
    end

    function fmpz(x::Float64)
        !isinteger(x) && throw(InexactError())
        z = new()
        ccall((:fmpz_init, :libflint), Nothing, (Ref{fmpz},), z)
        ccall((:fmpz_set_d, :libflint), Nothing, (Ref{fmpz}, Cdouble), z, x)
        finalizer(_fmpz_clear_fn, z)
        return z
    end

    fmpz(x::fmpz) = x
end

function _fmpz_clear_fn(a::fmpz)
   ccall((:fmpz_clear, :libflint), Nothing, (Ref{fmpz},), a)
end

mutable struct fmpz_factor
   sign::Cint
   p::Ptr{Nothing} # Array of fmpz_struct's
   exp::Ptr{UInt}
   alloc::Int
   num::Int

   function fmpz_factor()
      z = new()
      ccall((:fmpz_factor_init, :libflint), Nothing, (Ref{fmpz_factor}, ), z)
      finalizer(_fmpz_factor_clear_fn, z)
      return z
   end
end

function _fmpz_factor_clear_fn(a::fmpz_factor)
   ccall((:fmpz_factor_clear, :libflint), Nothing,
         (Ref{fmpz_factor}, ), a)
end

###############################################################################
#
#   FlintRationalField / fmpq
#
###############################################################################

mutable struct FlintRationalField <: FracField{fmpz}
end

const FlintQQ = FlintRationalField()

mutable struct fmpq <: FracElem{fmpz}
   num::Int
   den::Int

   function fmpq()
      z = new()
      ccall((:fmpq_init, :libflint), Nothing, (Ref{fmpq},), z)
      finalizer(_fmpq_clear_fn, z)
      return z
   end

   function fmpq(a::fmpz, b::fmpz)
      iszero(b) && throw(DivideError())
      z = new()
      ccall((:fmpq_init, :libflint), Nothing, (Ref{fmpq},), z)
      ccall((:fmpq_set_fmpz_frac, :libflint), Nothing,
            (Ref{fmpq}, Ref{fmpz}, Ref{fmpz}), z, a, b)
      finalizer(_fmpq_clear_fn, z)
      return z
   end

   function fmpq(a::fmpz)
      z = new()
      ccall((:fmpq_init, :libflint), Nothing, (Ref{fmpq},), z)
      b = fmpz(1)
      ccall((:fmpq_set_fmpz_frac, :libflint), Nothing,
            (Ref{fmpq}, Ref{fmpz}, Ref{fmpz}), z, a, b)
      finalizer(_fmpq_clear_fn, z)
      return z
   end

   function fmpq(a::Int, b::Int)
      b == 0 && throw(DivideError())
      z = new()
      ccall((:fmpq_init, :libflint), Nothing, (Ref{fmpq},), z)
      ccall((:fmpq_set_si, :libflint), Nothing,
            (Ref{fmpq}, Int, Int), z, a, b)
      finalizer(_fmpq_clear_fn, z)
      return z
   end

   function fmpq(a::Int)
      z = new()
      ccall((:fmpq_init, :libflint), Nothing, (Ref{fmpq},), z)
      ccall((:fmpq_set_si, :libflint), Nothing,
            (Ref{fmpq}, Int, Int), z, a, 1)
      finalizer(_fmpq_clear_fn, z)
      return z
   end

   fmpq(a::fmpq) = a
end

_fmpq_clear_fn(a::fmpq) = ccall((:fmpq_clear, :libflint), Nothing, (Ref{fmpq},), a)

###############################################################################
#
#   FmpzPolyRing / fmpz_poly
#
###############################################################################

mutable struct FmpzPolyRing <: PolyRing{fmpz}
   base_ring::FlintIntegerRing
   S::Symbol

   function FmpzPolyRing(s::Symbol, cached::Bool = true)
      if haskey(FmpzPolyID, s)
         return FmpzPolyID[s]
      else
         z = new(FlintZZ, s)
         if cached
            FmpzPolyID[s] = z
         end
         return z
      end
   end
end

const FmpzPolyID = Dict{Symbol, FmpzPolyRing}()

mutable struct fmpz_poly <: PolyElem{fmpz}
   coeffs::Ptr{Nothing}
   alloc::Int
   length::Int
   parent::FmpzPolyRing

   function fmpz_poly()
      z = new()
      ccall((:fmpz_poly_init, :libflint), Nothing, (Ref{fmpz_poly},), z)
      finalizer(_fmpz_poly_clear_fn, z)
      return z
   end

   function fmpz_poly(a::Array{fmpz, 1})
      z = new()
      ccall((:fmpz_poly_init2, :libflint), Nothing,
            (Ref{fmpz_poly}, Int), z, length(a))
      for i = 1:length(a)
         ccall((:fmpz_poly_set_coeff_fmpz, :libflint), Nothing,
                     (Ref{fmpz_poly}, Int, Ref{fmpz}), z, i - 1, a[i])
      end
      finalizer(_fmpz_poly_clear_fn, z)
      return z
   end

   function fmpz_poly(a::Int)
      z = new()
      ccall((:fmpz_poly_init, :libflint), Nothing, (Ref{fmpz_poly},), z)
      ccall((:fmpz_poly_set_si, :libflint), Nothing, (Ref{fmpz_poly}, Int), z, a)
      finalizer(_fmpz_poly_clear_fn, z)
      return z
   end

   function fmpz_poly(a::fmpz)
      z = new()
      ccall((:fmpz_poly_init, :libflint), Nothing, (Ref{fmpz_poly},), z)
      ccall((:fmpz_poly_set_fmpz, :libflint), Nothing,
            (Ref{fmpz_poly}, Ref{fmpz}), z, a)
      finalizer(_fmpz_poly_clear_fn, z)
      return z
   end

   function fmpz_poly(a::fmpz_poly)
      z = new()
      ccall((:fmpz_poly_init, :libflint), Nothing, (Ref{fmpz_poly},), z)
      ccall((:fmpz_poly_set, :libflint), Nothing,
            (Ref{fmpz_poly}, Ref{fmpz_poly}), z, a)
      finalizer(_fmpz_poly_clear_fn, z)
      return z
   end
end

function _fmpz_poly_clear_fn(a::fmpz_poly)
   ccall((:fmpz_poly_clear, :libflint), Nothing, (Ref{fmpz_poly},), a)
end

mutable struct fmpz_poly_factor
  d::Int # fmpz
  p::Ptr{fmpz_poly} # array of flint fmpz_poly_struct's
  exp::Ptr{Int}
  num::Int
  alloc::Int

  function fmpz_poly_factor()
    z = new()
    ccall((:fmpz_poly_factor_init, :libflint), Nothing,
                (Ref{fmpz_poly_factor}, ), z)
    finalizer(_fmpz_poly_factor_clear_fn, z)
    return z
  end
end

function _fmpz_poly_factor_clear_fn(f::fmpz_poly_factor)
  ccall((:fmpz_poly_factor_clear, :libflint), Nothing,
            (Ref{fmpz_poly_factor}, ), f)
  nothing
end

###############################################################################
#
#   FmpqPolyRing / fmpq_poly
#
###############################################################################

mutable struct FmpqPolyRing <: PolyRing{fmpq}
   base_ring::FlintRationalField
   S::Symbol

   function FmpqPolyRing(R::FlintRationalField, s::Symbol, cached::Bool = true)
      if haskey(FmpqPolyID, s)
         return FmpqPolyID[s]
      else
         z = new(R, s)
         if cached
            FmpqPolyID[s] = z
         end
         return z
      end
   end
end

const FmpqPolyID = Dict{Symbol, FmpqPolyRing}()

mutable struct fmpq_poly <: PolyElem{fmpq}
   coeffs::Ptr{Int}
   den::Int
   alloc::Int
   length::Int
   parent::FmpqPolyRing

   function fmpq_poly()
      z = new()
      ccall((:fmpq_poly_init, :libflint), Nothing, (Ref{fmpq_poly},), z)
      finalizer(_fmpq_poly_clear_fn, z)
      return z
   end

   function fmpq_poly(a::Array{fmpq, 1})
      z = new()
      ccall((:fmpq_poly_init2, :libflint), Nothing,
            (Ref{fmpq_poly}, Int), z, length(a))
      for i = 1:length(a)
         ccall((:fmpq_poly_set_coeff_fmpq, :libflint), Nothing,
                     (Ref{fmpq_poly}, Int, Ref{fmpq}), z, i - 1, a[i])
      end
      finalizer(_fmpq_poly_clear_fn, z)
      return z
   end

   function fmpq_poly(a::Int)
      z = new()
      ccall((:fmpq_poly_init, :libflint), Nothing, (Ref{fmpq_poly},), z)
      ccall((:fmpq_poly_set_si, :libflint), Nothing, (Ref{fmpq_poly}, Int), z, a)
      finalizer(_fmpq_poly_clear_fn, z)
      return z
   end

   function fmpq_poly(a::fmpz)
      z = new()
      ccall((:fmpq_poly_init, :libflint), Nothing, (Ref{fmpq_poly},), z)
      ccall((:fmpq_poly_set_fmpz, :libflint), Nothing,
            (Ref{fmpq_poly}, Ref{fmpz}), z, a)
      finalizer(_fmpq_poly_clear_fn, z)
      return z
   end

   function fmpq_poly(a::fmpq)
      z = new()
      ccall((:fmpq_poly_init, :libflint), Nothing, (Ref{fmpq_poly},), z)
      ccall((:fmpq_poly_set_fmpq, :libflint), Nothing,
            (Ref{fmpq_poly}, Ref{fmpq}), z, a)
      finalizer(_fmpq_poly_clear_fn, z)
      return z
   end

   function fmpq_poly(a::fmpz_poly)
      z = new()
      ccall((:fmpq_poly_init, :libflint), Nothing, (Ref{fmpq_poly},), z)
      ccall((:fmpq_poly_set_fmpz_poly, :libflint), Nothing,
            (Ref{fmpq_poly}, Ref{fmpz_poly}), z, a)
      finalizer(_fmpq_poly_clear_fn, z)
      return z
   end

   function fmpq_poly(a::fmpq_poly)
      z = new()
      ccall((:fmpq_poly_init, :libflint), Nothing, (Ref{fmpq_poly},), z)
      ccall((:fmpq_poly_set, :libflint), Nothing,
            (Ref{fmpq_poly}, Ref{fmpq_poly}), z, a)
      finalizer(_fmpq_poly_clear_fn, z)
      return z
   end
end

function _fmpq_poly_clear_fn(a::fmpq_poly)
   ccall((:fmpq_poly_clear, :libflint), Nothing, (Ref{fmpq_poly},), a)
end

###############################################################################
#
#   NmodRing / nmod
#
###############################################################################

mutable struct NmodRing <: Ring
   n::UInt
   ninv::UInt

   function NmodRing(n::UInt, cached::Bool=true)
      if haskey(NmodRingID, n)
         return NmodRingID[n]
      else
         ninv = ccall((:n_preinvert_limb, :libflint), UInt, (UInt,), n)
         z = new(n, ninv)
         if cached
            NmodRingID[n] = z
         end
         return z
      end
   end
end

const NmodRingID = Dict{UInt, NmodRing}()

struct nmod <: ResElem{UInt}
   data::UInt
   parent::NmodRing
end

################################################################################
#
#   GaloisField / gfp
#
###############################################################################

mutable struct GaloisField <: FinField
   n::UInt
   ninv::UInt

   function GaloisField(n::UInt, cached::Bool=true)
      if haskey(GaloisFieldID, n)
         return GaloisFieldID[n]
      else
         ninv = ccall((:n_preinvert_limb, :libflint), UInt, (UInt,), n)
         z = new(n, ninv)
         if cached
            GaloisFieldID[n] = z
         end
         return z
      end
   end
end

const GaloisFieldID = Dict{UInt, GaloisField}()

struct gfp_elem <: FinFieldElem
   data::UInt
   parent::GaloisField
end

###############################################################################
#
#   NmodPolyRing / nmod_poly
#
###############################################################################

mutable struct NmodPolyRing <: PolyRing{nmod}
  base_ring::NmodRing
  S::Symbol
  n::UInt

  function NmodPolyRing(R::NmodRing, s::Symbol, cached::Bool = true)
    m = UInt(modulus(R))
    if haskey(NmodPolyRingID, (m, s))
       return NmodPolyRingID[m, s]
    else
       z = new(R, s, m)
       if cached
          NmodPolyRingID[m, s] = z
       end
       return z
    end
  end
end

const NmodPolyRingID = Dict{Tuple{UInt, Symbol}, NmodPolyRing}()

mutable struct nmod_poly <: PolyElem{nmod}
   coeffs::Ptr{Nothing}
   alloc::Int
   length::Int
   mod_n::UInt
   mod_ninv::UInt
   mod_norm::UInt
   parent::NmodPolyRing

   function nmod_poly(n::UInt)
      z = new()
      ccall((:nmod_poly_init, :libflint), Nothing, (Ref{nmod_poly}, UInt), z, n)
      finalizer(_nmod_poly_clear_fn, z)
      return z
   end

   function nmod_poly(n::UInt, a::UInt)
      z = new()
      ccall((:nmod_poly_init, :libflint), Nothing, (Ref{nmod_poly}, UInt), z, n)
      ccall((:nmod_poly_set_coeff_ui, :libflint), Nothing,
              (Ref{nmod_poly}, Int, UInt), z, 0, a)
      finalizer(_nmod_poly_clear_fn, z)
      return z
   end

   function nmod_poly(n::UInt, a::Int)
      z = new()
      ccall((:nmod_poly_init, :libflint), Nothing, (Ref{nmod_poly}, UInt), z, n)
      ccall((:nmod_poly_set_coeff_ui, :libflint), Nothing,
              (Ref{nmod_poly}, Int, UInt), z, 0, mod(a, n))
      finalizer(_nmod_poly_clear_fn, z)
      return z
   end

   function nmod_poly(n::UInt, arr::Array{fmpz, 1})
      z = new()
      ccall((:nmod_poly_init2, :libflint), Nothing,
            (Ref{nmod_poly}, UInt, Int), z, n, length(arr))
      for i in 1:length(arr)
         tt = ccall((:fmpz_fdiv_ui, :libflint), UInt, (Ref{fmpz}, UInt), arr[i], n)
         ccall((:nmod_poly_set_coeff_ui, :libflint), Nothing,
              (Ref{nmod_poly}, Int, UInt), z, i - 1, tt)
      end
      finalizer(_nmod_poly_clear_fn, z)
      return z
   end

   function nmod_poly(n::UInt, arr::Array{UInt, 1})
      z = new()
      ccall((:nmod_poly_init2, :libflint), Nothing,
            (Ref{nmod_poly}, UInt, Int), z, n, length(arr))
      for i in 1:length(arr)
         ccall((:nmod_poly_set_coeff_ui, :libflint), Nothing,
              (Ref{nmod_poly}, Int, UInt), z, i - 1, arr[i])
      end
      finalizer(_nmod_poly_clear_fn, z)
      return z
   end

   function nmod_poly(n::UInt, arr::Array{nmod, 1})
      z = new()
      ccall((:nmod_poly_init2, :libflint), Nothing,
            (Ref{nmod_poly}, UInt, Int), z, n, length(arr))
      for i in 1:length(arr)
         ccall((:nmod_poly_set_coeff_ui, :libflint), Nothing,
              (Ref{nmod_poly}, Int, UInt), z, i-1, arr[i].data)
      end
      finalizer(_nmod_poly_clear_fn, z)
      return z
   end

   function nmod_poly(n::UInt, f::fmpz_poly)
      z = new()
      ccall((:nmod_poly_init2, :libflint), Nothing,
            (Ref{nmod_poly}, UInt, Int), z, n, length(f))
      ccall((:fmpz_poly_get_nmod_poly, :libflint), Nothing,
            (Ref{nmod_poly}, Ref{fmpz_poly}), z, f)
      finalizer(_nmod_poly_clear_fn, z)
      return z
   end

   function nmod_poly(n::UInt, f::nmod_poly)
      z = new()
      ccall((:nmod_poly_init2, :libflint), Nothing,
            (Ref{nmod_poly}, UInt, Int), z, n, length(f))
      ccall((:nmod_poly_set, :libflint), Nothing,
            (Ref{nmod_poly}, Ref{nmod_poly}), z, f)
      finalizer(_nmod_poly_clear_fn, z)
      return z
   end
end

function _nmod_poly_clear_fn(x::nmod_poly)
  ccall((:nmod_poly_clear, :libflint), Nothing, (Ref{nmod_poly}, ), x)
end

mutable struct nmod_poly_factor
  poly::Ptr{nmod_poly}  # array of flint nmod_poly_struct's
  exp::Ptr{Int}
  num::Int
  alloc::Int
  n::UInt

  function nmod_poly_factor(n::UInt)
    z = new()
    ccall((:nmod_poly_factor_init, :libflint), Nothing,
            (Ref{nmod_poly_factor}, ), z)
    z.n = n
    finalizer(_nmod_poly_factor_clear_fn, z)
    return z
  end
end

function _nmod_poly_factor_clear_fn(a::nmod_poly_factor)
  ccall((:nmod_poly_factor_clear, :libflint), Nothing,
          (Ref{nmod_poly_factor}, ), a)
end

################################################################################
#
#   GFPPolyRing / gfp_poly
#
###############################################################################

mutable struct GFPPolyRing <: PolyRing{gfp_elem}
  base_ring::GaloisField
  S::Symbol
  n::UInt

  function GFPPolyRing(R::GaloisField, s::Symbol, cached::Bool = true)
    m = UInt(modulus(R))
    if haskey(GFPPolyRingID, (m, s))
       return GFPPolyRingID[m, s]
    else
       z = new(R, s, m)
       if cached
          GFPPolyRingID[m, s] = z
       end
       return z
    end
  end
end

const GFPPolyRingID = Dict{Tuple{UInt, Symbol}, GFPPolyRing}()

mutable struct gfp_poly <: PolyElem{gfp_elem}
   coeffs::Ptr{Nothing}
   alloc::Int
   length::Int
   mod_n::UInt
   mod_ninv::UInt
   mod_norm::UInt
   parent::GFPPolyRing

   function gfp_poly(n::UInt)
      z = new()
      ccall((:nmod_poly_init, :libflint), Nothing, (Ref{gfp_poly}, UInt), z, n)
      finalizer(_gfp_poly_clear_fn, z)
      return z
   end

   function gfp_poly(n::UInt, a::UInt)
      z = new()
      ccall((:nmod_poly_init, :libflint), Nothing, (Ref{gfp_poly}, UInt), z, n)
      ccall((:nmod_poly_set_coeff_ui, :libflint), Nothing,
              (Ref{gfp_poly}, Int, UInt), z, 0, a)
      finalizer(_gfp_poly_clear_fn, z)
      return z
   end

   function gfp_poly(n::UInt, a::Int)
      z = new()
      ccall((:nmod_poly_init, :libflint), Nothing, (Ref{gfp_poly}, UInt), z, n)
      ccall((:nmod_poly_set_coeff_ui, :libflint), Nothing,
              (Ref{gfp_poly}, Int, UInt), z, 0, mod(a, n))
      finalizer(_gfp_poly_clear_fn, z)
      return z
   end

   function gfp_poly(n::UInt, arr::Array{fmpz, 1})
      z = new()
      ccall((:nmod_poly_init2, :libflint), Nothing,
            (Ref{gfp_poly}, UInt, Int), z, n, length(arr))
      for i in 1:length(arr)
         tt = ccall((:fmpz_fdiv_ui, :libflint), UInt, (Ref{fmpz}, UInt), arr[i], n)
         ccall((:nmod_poly_set_coeff_ui, :libflint), Nothing,
              (Ref{gfp_poly}, Int, UInt), z, i - 1, tt)
      end
      finalizer(_gfp_poly_clear_fn, z)
      return z
   end

   function gfp_poly(n::UInt, arr::Array{UInt, 1})
      z = new()
      ccall((:nmod_poly_init2, :libflint), Nothing,
            (Ref{gfp_poly}, UInt, Int), z, n, length(arr))
      for i in 1:length(arr)
         ccall((:nmod_poly_set_coeff_ui, :libflint), Nothing,
              (Ref{gfp_poly}, Int, UInt), z, i - 1, arr[i])
      end
      finalizer(_gfp_poly_clear_fn, z)
      return z
   end

   function gfp_poly(n::UInt, arr::Array{gfp_elem, 1})
      z = new()
      ccall((:nmod_poly_init2, :libflint), Nothing,
            (Ref{gfp_poly}, UInt, Int), z, n, length(arr))
      for i in 1:length(arr)
         ccall((:nmod_poly_set_coeff_ui, :libflint), Nothing,
              (Ref{gfp_poly}, Int, UInt), z, i-1, arr[i].data)
      end
      finalizer(_gfp_poly_clear_fn, z)
      return z
   end

   function gfp_poly(n::UInt, f::fmpz_poly)
      z = new()
      ccall((:nmod_poly_init2, :libflint), Nothing,
            (Ref{gfp_poly}, UInt, Int), z, n, length(f))
      ccall((:fmpz_poly_get_nmod_poly, :libflint), Nothing,
            (Ref{gfp_poly}, Ref{fmpz_poly}), z, f)
      finalizer(_gfp_poly_clear_fn, z)
      return z
   end

   function gfp_poly(n::UInt, f::gfp_poly)
      z = new()
      ccall((:nmod_poly_init2, :libflint), Nothing,
            (Ref{gfp_poly}, UInt, Int), z, n, length(f))
      ccall((:nmod_poly_set, :libflint), Nothing,
            (Ref{gfp_poly}, Ref{gfp_poly}), z, f)
      finalizer(_gfp_poly_clear_fn, z)
      return z
   end
end

function _gfp_poly_clear_fn(x::gfp_poly)
  ccall((:nmod_poly_clear, :libflint), Nothing, (Ref{gfp_poly}, ), x)
end

mutable struct gfp_poly_factor
  poly::Ptr{gfp_poly}  # array of flint nmod_poly_struct's
  exp::Ptr{Int}
  num::Int
  alloc::Int
  n::UInt

  function gfp_poly_factor(n::UInt)
    z = new()
    ccall((:nmod_poly_factor_init, :libflint), Nothing,
            (Ref{gfp_poly_factor}, ), z)
    z.n = n
    finalizer(_gfp_poly_factor_clear_fn, z)
    return z
  end
end

function _gfp_poly_factor_clear_fn(a::gfp_poly_factor)
  ccall((:nmod_poly_factor_clear, :libflint), Nothing,
          (Ref{gfp_poly_factor}, ), a)
end

const Zmodn_poly = Union{nmod_poly, gfp_poly}

###############################################################################
#
#   FmpzModPolyRing / fmpz_mod_poly
#
###############################################################################

mutable struct FmpzModPolyRing <: PolyRing{Generic.Res{fmpz}}
  base_ring::Generic.ResRing{fmpz}
  S::Symbol
  n::fmpz

  function FmpzModPolyRing(R::Generic.ResRing{fmpz}, s::Symbol, cached::Bool = true)
    m = modulus(R)
    if haskey(FmpzModPolyRingID, (m, s))
       return FmpzModPolyRingID[m, s]
    else
       z = new(R, s, m)
       if cached
          FmpzModPolyRingID[m ,s] = z
       end
       return z
    end
  end
end

const FmpzModPolyRingID = Dict{Tuple{fmpz, Symbol}, FmpzModPolyRing}()

mutable struct fmpz_mod_poly <: PolyElem{Generic.Res{fmpz}}
   coeffs::Ptr{Nothing}
   alloc::Int
   length::Int
   p::Int
   parent::FmpzModPolyRing

   function fmpz_mod_poly(n::fmpz)
      z = new()
      ccall((:fmpz_mod_poly_init, :libflint), Nothing,
            (Ref{fmpz_mod_poly}, Ref{fmpz}), z, n)
      finalizer(_fmpz_mod_poly_clear_fn, z)
      return z
   end

   function fmpz_mod_poly(n::fmpz, a::UInt)
      z = new()
      ccall((:fmpz_mod_poly_init, :libflint), Nothing,
            (Ref{fmpz_mod_poly}, Ref{fmpz}), z, n)
      ccall((:fmod_poly_set_coeff_ui, :libflint), Nothing,
              (Ref{fmpz_mod_poly}, Int, UInt), z, 0, a)
      finalizer(_fmpz_mod_poly_clear_fn, z)
      return z
   end

   function fmpz_mod_poly(n::fmpz, a::fmpz)
      z = new()
      ccall((:fmpz_mod_poly_init, :libflint), Nothing,
            (Ref{fmpz_mod_poly}, Ref{fmpz}), z, n)
      ccall((:fmpz_mod_poly_set_coeff_fmpz, :libflint), Nothing,
              (Ref{fmpz_mod_poly}, Int, Ref{fmpz}), z, 0, a)
      finalizer(_fmpz_mod_poly_clear_fn, z)
      return z
   end

   function fmpz_mod_poly(n::fmpz, arr::Array{fmpz, 1})
      length(arr) == 0 && error("Array must have length > 0")
      z = new()
      ccall((:fmpz_mod_poly_init2, :libflint), Nothing,
            (Ref{fmpz_mod_poly}, Ref{fmpz}, Int), z, n, length(arr))
      for i in 1:length(arr)
         ccall((:fmpz_mod_poly_set_coeff_fmpz, :libflint), Nothing,
              (Ref{fmpz_mod_poly}, Int, Ref{fmpz}), z, i - 1, arr[i])
      end
      finalizer(_fmpz_mod_poly_clear_fn, z)
      return z
   end

   function fmpz_mod_poly(n::fmpz, arr::Array{Generic.Res{fmpz}, 1})
      z = new()
      ccall((:fmpz_mod_poly_init2, :libflint), Nothing,
            (Ref{fmpz_mod_poly}, Ref{fmpz}, Int), z, n, length(arr))
      for i in 1:length(arr)
         ccall((:fmpz_mod_poly_set_coeff_fmpz, :libflint), Nothing,
              (Ref{fmpz_mod_poly}, Int, Ref{fmpz}), z, i - 1, arr[i].data)
      end
      finalizer(_fmpz_mod_poly_clear_fn, z)
      return z
   end

   function fmpz_mod_poly(n::fmpz, f::fmpz_poly)
      z = new()
      ccall((:fmpz_mod_poly_init2, :libflint), Nothing,
            (Ref{fmpz_mod_poly}, Ref{fmpz}, Int), z, n, length(f))
      ccall((:fmpz_mod_poly_set_fmpz_poly, :libflint), Nothing,
            (Ref{fmpz_mod_poly}, Ref{fmpz_poly}), z, f)
      finalizer(_fmpz_mod_poly_clear_fn, z)
      return z
   end

   function fmpz_mod_poly(n::fmpz, f::fmpz_mod_poly)
      z = new()
      ccall((:fmpz_mod_poly_init2, :libflint), Nothing,
            (Ref{fmpz_mod_poly}, Ref{fmpz}, Int), z, n, length(f))
      ccall((:fmpz_mod_poly_set, :libflint), Nothing,
            (Ref{fmpz_mod_poly}, Ref{fmpz_mod_poly}), z, f)
      finalizer(_fmpz_mod_poly_clear_fn, z)
      return z
   end
end

function _fmpz_mod_poly_clear_fn(x::fmpz_mod_poly)
  ccall((:fmpz_mod_poly_clear, :libflint), Nothing, (Ref{fmpz_mod_poly}, ), x)
end

mutable struct fmpz_mod_poly_factor
  poly::Ptr{fmpz_mod_poly}
  exp::Ptr{Int}
  num::Int
  alloc::Int
  n::fmpz

  function fmpz_mod_poly_factor(n::fmpz)
    z = new()
    ccall((:fmpz_mod_poly_factor_init, :libflint), Nothing,
            (Ref{fmpz_mod_poly_factor}, ), z)
    z.n = n
    finalizer(_fmpz_mod_poly_factor_clear_fn, z)
    return z
  end
end

function _fmpz_mod_poly_factor_clear_fn(a::fmpz_mod_poly_factor)
  ccall((:fmpz_mod_poly_factor_clear, :libflint), Nothing,
          (Ref{fmpz_mod_poly_factor}, ), a)
end

###############################################################################
#
#   FmpzMPolyRing / fmpz_mpoly
#
###############################################################################

# S is a Symbol which can take the values:
# :lex
# :deglex
# :degrevlex
#
# T is an Int which is the number of variables
# (plus one if ordered by total degree)

mutable struct FmpzMPolyRing{S, N} <: PolyRing{fmpz}
   n::Int
   ord::Cint
   base_ring::FlintIntegerRing
   S::Array{Symbol, 1}
   num_vars::Int

   function FmpzMPolyRing{S, N}(s::Array{Symbol, 1}, cached::Bool = true) where {S, N}
      if haskey(FmpzMPolyID, (s, S, N))
         return FmpzMPolyID[s, S, N]
      else
         if S == :lex
            ord = 0
         elseif S == :deglex
            ord = 1
         elseif S == :degrevlex
            ord = 2
         end

         z = new{S, N}()
         ccall((:fmpz_mpoly_ctx_init, :libflint), Nothing,
               (Ref{FmpzMPolyRing}, Int, Int),
               z, length(s), ord)
         z.base_ring = FlintZZ
         z.S = s
         z.num_vars = length(s)
         finalizer(_fmpz_mpoly_ctx_clear_fn, z)
         if cached
            FmpzMPolyID[s, S, N] = z
         end
         return z
      end
   end
end

const FmpzMPolyID = Dict{Tuple{Array{Symbol, 1}, Symbol, Int}, PolyRing{fmpz}}()

mutable struct fmpz_mpoly{S, N} <: PolyElem{fmpz}
   coeffs::Ptr{Nothing}
   exps::Ptr{Nothing}
   alloc::Int
   length::Int
   bits::Int
   parent::FmpzMPolyRing

   function fmpz_mpoly(ctx::FmpzMPolyRing{S, N}) where {S, N}
      z = new{S, N}()
      ccall((:fmpz_mpoly_init, :libflint), Nothing,
            (Ref{fmpz_mpoly}, Ref{FmpzMPolyRing},), z, ctx)
      finalizer(_fmpz_mpoly_clear_fn, z)
      return z
   end

   function fmpz_mpoly(ctx::FmpzMPolyRing{S, N}, a::Array{fmpz, 1}, b::Array{NTuple{N, Int}, 1}) where {S, N}
      z = new{S, N}()
      ccall((:fmpz_mpoly_init, :libflint), Nothing,
            (Ref{fmpz_mpoly}, Ref{FmpzMPolyRing},), z, ctx)
      m = 0
      for i = 1:length(b)
         for j = 1:N
            if b[i][j] > m
               m = b[i][j]
            end
         end
      end
      bits = 8
      while ndigits(m, 2) >= bits
         bits *= 2
      end
      deg = ctx.ord == :deglex || ctx.ord == :degrevlex ? 1 : 0
      ccall((:fmpz_mpoly_fit_length, :libflint), Nothing,
            (Ref{fmpz_mpoly}, Int, Ref{FmpzMPolyRing}), z, length(a), ctx)
      ccall((:fmpz_mpoly_fit_bits, :libflint), Nothing,
            (Ref{fmpz_mpoly}, Int, Ref{FmpzMPolyRing}), z, bits, ctx)
      for i = 1:length(a)
         ccall((:fmpz_mpoly_set_coeff_fmpz, :libflint), Nothing,
            (Ref{fmpz_mpoly}, Int, Ref{fmpz}, Ref{FmpzMPolyRing}),
                                                        z, i - 1, a[i], ctx)
         A = [b[i][j + deg] for j = 1:N - deg]
         ccall((:fmpz_mpoly_set_monomial, :libflint), Nothing,
            (Ref{fmpz_mpoly}, Int, Ptr{Int}, Ref{FmpzMPolyRing}),
                                                            z, i - 1, A, ctx)
      end
      ccall((:_fmpz_mpoly_set_length, :libflint), Nothing,
            (Ref{fmpz_mpoly}, Int, Ref{FmpzMPolyRing}), z, length(a), ctx)
      finalizer(_fmpz_mpoly_clear_fn, z)
      return z
   end

   function fmpz_mpoly(ctx::FmpzMPolyRing{S, N}, a::Int) where {S, N}
      z = new{S, N}()
      ccall((:fmpz_mpoly_init, :libflint), Nothing,
            (Ref{fmpz_mpoly}, Ref{FmpzMPolyRing},), z, ctx)
      ccall((:fmpz_mpoly_set_si, :libflint), Nothing,
            (Ref{fmpz_mpoly}, Int, Ref{FmpzMPolyRing}), z, a, ctx)
      finalizer(_fmpz_mpoly_clear_fn, z)
      return z
   end

   function fmpz_mpoly(ctx::FmpzMPolyRing{S, N}, a::fmpz) where {S, N}
      z = new{S, N}()
      ccall((:fmpz_mpoly_init, :libflint), Nothing,
            (Ref{fmpz_mpoly}, Ref{FmpzMPolyRing},), z, ctx)
      ccall((:fmpz_mpoly_set_fmpz, :libflint), Nothing,
            (Ref{fmpz_mpoly}, Ref{fmpz}, Ref{FmpzMPolyRing}), z, a, ctx)
      finalizer(_fmpz_mpoly_clear_fn, z)
      return z
   end
end

function _fmpz_mpoly_clear_fn(a::fmpz_mpoly)
  ccall((:fmpz_mpoly_clear, :libflint), Nothing,
          (Ref{fmpz_mpoly}, Ref{FmpzMPolyRing}), a, a.parent)
end

function _fmpz_mpoly_ctx_clear_fn(a::FmpzMPolyRing)
  ccall((:fmpz_mpoly_ctx_clear, :libflint), Nothing,
          (Ref{FmpzMPolyRing},), a)
end

###############################################################################
#
#   FmpqMPolyRing / fmpq_mpoly
#
###############################################################################

const flint_orderings = [:lex, :deglex, :degrevlex]

# S is a Symbol which can take the values:
# :lex
# :deglex
# :degrevlex

mutable struct FmpqMPolyRing <: MPolyRing{fmpq}
   nvars::Int
   nfields::Cint
   ord::Int
   deg::Cint
   rev::Cint
   base_ring::FlintRationalField
   S::Array{Symbol, 1}

   function FmpqMPolyRing(s::Array{Symbol, 1}, S::Symbol, cached::Bool = true)
      if cached && haskey(FmpqMPolyID, (s, S))
         return FmpqMPolyID[s, S]
      else
         if S == :lex
            ord = 0
         elseif S == :deglex
            ord = 1
         elseif S == :degrevlex
            ord = 2
         else
            error("$S is not a valid ordering")
         end

         z = new()
         ccall((:fmpq_mpoly_ctx_init, :libflint), Nothing,
               (Ref{FmpqMPolyRing}, Int, Int),
               z, length(s), ord)
         z.base_ring = FlintQQ
         z.S = s
         finalizer(_fmpq_mpoly_ctx_clear_fn, z)
         if cached
            FmpqMPolyID[s, S] = z
         end
         return z
      end
   end
end

function _fmpq_mpoly_ctx_clear_fn(a::FmpqMPolyRing)
  ccall((:fmpq_mpoly_ctx_clear, :libflint), Nothing,
          (Ref{FmpqMPolyRing},), a)
end

const FmpqMPolyID = Dict{Tuple{Array{Symbol, 1}, Symbol}, FmpqMPolyRing}()

mutable struct fmpq_mpoly <: MPolyElem{fmpq}
   content_num::Ptr{Nothing}
   content_den::Ptr{Nothing}
   coeffs::Ptr{Nothing}
   exps::Ptr{Nothing}
   alloc::Int
   length::Int
   bits::Int

   parent::FmpqMPolyRing

   function fmpq_mpoly(ctx::FmpqMPolyRing)
      z = new()
      ccall((:fmpq_mpoly_init, :libflint), Nothing,
            (Ref{fmpq_mpoly}, Ref{FmpqMPolyRing},), z, ctx)
      z.parent = ctx
      finalizer(_fmpq_mpoly_clear_fn, z)
      return z
   end

   function fmpq_mpoly(ctx::FmpqMPolyRing, a::Vector{fmpq}, b::Vector{Vector{UInt}})
      z = new()
      ccall((:fmpq_mpoly_init, :libflint), Nothing,
            (Ref{fmpq_mpoly}, Ref{FmpqMPolyRing},), z, ctx)
      z.parent = ctx
      finalizer(_fmpq_mpoly_clear_fn, z)

      for i in 1:length(a)
        ccall((:fmpq_mpoly_pushterm_fmpq_ui, :libflint), Nothing,
              (Ref{fmpq_mpoly}, Ref{fmpq}, Ptr{UInt}, Ref{FmpqMPolyRing}),
              z, a[i], b[i], ctx)
      end

      ccall((:fmpq_mpoly_sort, :libflint), Nothing,
            (Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}), z, ctx)
      ccall((:fmpq_mpoly_combine_like_terms, :libflint), Nothing,
            (Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}), z, ctx)
      return z
   end

   function fmpq_mpoly(ctx::FmpqMPolyRing, a::Vector{fmpq}, b::Vector{Vector{fmpz}})
      z = new()
      ccall((:fmpq_mpoly_init, :libflint), Nothing,
            (Ref{fmpq_mpoly}, Ref{FmpqMPolyRing},), z, ctx)
      z.parent = ctx
      finalizer(_fmpq_mpoly_clear_fn, z)

      for i in 1:length(a)
        ccall((:fmpq_mpoly_pushterm_fmpq_fmpz, :libflint), Nothing,
              (Ref{fmpq_mpoly}, Ref{fmpq}, Ptr{Ref{fmpz}}, Ref{FmpqMPolyRing}),
              z, a[i], b[i], ctx)
      end

      ccall((:fmpq_mpoly_sort, :libflint), Nothing,
            (Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}), z, ctx)
      ccall((:fmpq_mpoly_combine_like_terms, :libflint), Nothing,
            (Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}), z, ctx)
      return z
   end

   function fmpq_mpoly(ctx::FmpqMPolyRing, a::fmpz)
      z = new()
      ccall((:fmpq_mpoly_init, :libflint), Nothing,
            (Ref{fmpq_mpoly}, Ref{FmpqMPolyRing},), z, ctx)
      ccall((:fmpq_mpoly_set_fmpz, :libflint), Nothing,
            (Ref{fmpq_mpoly}, Ref{fmpz}, Ref{FmpqMPolyRing}), z, a, ctx)
      z.parent = ctx
      finalizer(_fmpq_mpoly_clear_fn, z)
      return z
   end

   function fmpq_mpoly(ctx::FmpqMPolyRing, a::fmpq)
      z = new()
      ccall((:fmpq_mpoly_init, :libflint), Nothing,
            (Ref{fmpq_mpoly}, Ref{FmpqMPolyRing},), z, ctx)
      ccall((:fmpq_mpoly_set_fmpq, :libflint), Nothing,
            (Ref{fmpq_mpoly}, Ref{fmpq}, Ref{FmpqMPolyRing}), z, a, ctx)
      z.parent = ctx
      finalizer(_fmpq_mpoly_clear_fn, z)
      return z
   end
   
   function fmpq_mpoly(ctx::FmpqMPolyRing, a::Int)
      z = new()
      ccall((:fmpq_mpoly_init, :libflint), Nothing,
            (Ref{fmpq_mpoly}, Ref{FmpqMPolyRing},), z, ctx)
      ccall((:fmpq_mpoly_set_si, :libflint), Nothing,
            (Ref{fmpq_mpoly}, Int, Ref{FmpqMPolyRing}), z, a, ctx)
      z.parent = ctx
      finalizer(_fmpq_mpoly_clear_fn, z)
      return z
   end

   function fmpq_mpoly(ctx::FmpqMPolyRing, a::UInt)
      z = new()
      ccall((:fmpq_mpoly_init, :libflint), Nothing,
            (Ref{fmpq_mpoly}, Ref{FmpqMPolyRing},), z, ctx)
      ccall((:fmpq_mpoly_set_ui, :libflint), Nothing,
            (Ref{fmpq_mpoly}, UInt, Ref{FmpqMPolyRing}), z, a, ctx)
      z.parent = ctx
      finalizer(_fmpq_mpoly_clear_fn, z)
      return z
   end
end

function _fmpq_mpoly_clear_fn(a::fmpq_mpoly)
  ccall((:fmpq_mpoly_clear, :libflint), Nothing,
          (Ref{fmpq_mpoly}, Ref{FmpqMPolyRing}), a, a.parent)
end

###############################################################################
#
#   FqNmodFiniteField / fq_nmod
#
###############################################################################

mutable struct FqNmodFiniteField <: FinField
   p :: Int
   n :: Int
   ninv :: Int
   norm :: Int
   sparse_modulus :: Int
   a :: Ptr{Nothing}
   j :: Ptr{Nothing}
   len :: Int
   mod_coeffs :: Ptr{Nothing}
   mod_alloc :: Int
   mod_length :: Int
   mod_n :: Int
   mod_ninv :: Int
   mod_norm :: Int
   inv_coeffs :: Ptr{Nothing}
   inv_alloc :: Int
   inv_length :: Int
   inv_n :: Int
   inv_ninv :: Int
   inv_norm :: Int
   var :: Ptr{Nothing}

   function FqNmodFiniteField(c::fmpz, deg::Int, s::Symbol, cached::Bool = true)
      if haskey(FqNmodFiniteFieldID, (c, deg, s))
         return FqNmodFiniteFieldID[c, deg, s]
      else
         d = new()
         ccall((:fq_nmod_ctx_init, :libflint), Nothing,
               (Ref{FqNmodFiniteField}, Ref{fmpz}, Int, Ptr{UInt8}),
			    d, c, deg, string(s))
         if cached
            FqNmodFiniteFieldID[c, deg, s] = d
         end
         finalizer(_FqNmodFiniteField_clear_fn, d)
         return d
      end
   end

   function FqNmodFiniteField(f::nmod_poly, s::Symbol, cached::Bool = true)
      if haskey(FqNmodFiniteFieldIDPol, (parent(f), f, s))
         return FqNmodFiniteFieldIDPol[parent(f), f, s]
      else
         z = new()
         ccall((:fq_nmod_ctx_init_modulus, :libflint), Nothing,
            (Ref{FqNmodFiniteField}, Ref{nmod_poly}, Ptr{UInt8}),
	      z, f, string(s))
         if cached
            FqNmodFiniteFieldIDPol[parent(f), f, s] = z
         end
         finalizer(_FqNmodFiniteField_clear_fn, z)
         return z
      end
   end
end

const FqNmodFiniteFieldID = Dict{Tuple{fmpz, Int, Symbol}, FqNmodFiniteField}()

const FqNmodFiniteFieldIDPol = Dict{Tuple{NmodPolyRing, nmod_poly, Symbol},
                                    FqNmodFiniteField}()

function _FqNmodFiniteField_clear_fn(a :: FqNmodFiniteField)
   ccall((:fq_nmod_ctx_clear, :libflint), Nothing, (Ref{FqNmodFiniteField},), a)
end

mutable struct fq_nmod <: FinFieldElem
   coeffs :: Ptr{Nothing}
   alloc :: Int
   length :: Int
   n :: Int
   ninv :: Int
   norm :: Int
   parent::FqNmodFiniteField

   function fq_nmod(ctx::FqNmodFiniteField)
      d = new()
      ccall((:fq_nmod_init2, :libflint), Nothing,
            (Ref{fq_nmod}, Ref{FqNmodFiniteField}), d, ctx)
      finalizer(_fq_nmod_clear_fn, d)
      return d
   end

   function fq_nmod(ctx::FqNmodFiniteField, x::Int)
      d = new()
      ccall((:fq_nmod_init2, :libflint), Nothing,
            (Ref{fq_nmod}, Ref{FqNmodFiniteField}), d, ctx)
      finalizer(_fq_nmod_clear_fn, d)
      ccall((:fq_nmod_set_si, :libflint), Nothing,
                (Ref{fq_nmod}, Int, Ref{FqNmodFiniteField}), d, x, ctx)
      return d
   end

   function fq_nmod(ctx::FqNmodFiniteField, x::fmpz)
      d = new()
      ccall((:fq_nmod_init2, :libflint), Nothing,
            (Ref{fq_nmod}, Ref{FqNmodFiniteField}), d, ctx)
      finalizer(_fq_nmod_clear_fn, d)
      ccall((:fq_nmod_set_fmpz, :libflint), Nothing,
            (Ref{fq_nmod}, Ref{fmpz}, Ref{FqNmodFiniteField}), d, x, ctx)
      return d
   end

      function fq_nmod(ctx::FqNmodFiniteField, x::fq_nmod)
      d = new()
      ccall((:fq_nmod_init2, :libflint), Nothing,
            (Ref{fq_nmod}, Ref{FqNmodFiniteField}), d, ctx)
      finalizer(_fq_nmod_clear_fn, d)
      ccall((:fq_nmod_set, :libflint), Nothing,
            (Ref{fq_nmod}, Ref{fq_nmod}, Ref{FqNmodFiniteField}), d, x, ctx)
      return d
   end
end

function _fq_nmod_clear_fn(a::fq_nmod)
   ccall((:fq_nmod_clear, :libflint), Nothing,
         (Ref{fq_nmod}, Ref{FqNmodFiniteField}), a, a.parent)
end

###############################################################################
#
#   FqFiniteField / fq
#
###############################################################################

mutable struct FqFiniteField <: FinField
   p::Int # fmpz
   sparse_modulus::Int
   a::Ptr{Nothing}
   j::Ptr{Nothing}
   len::Int
   mod_coeffs::Ptr{Nothing}
   mod_alloc::Int
   mod_length::Int
   mod_p::Int # fmpz
   inv_coeffs::Ptr{Nothing}
   inv_alloc::Int
   inv_length::Int
   inv_p::Int # fmpz
   var::Ptr{Nothing}

   function FqFiniteField(char::fmpz, deg::Int, s::Symbol, cached::Bool = true)
      if haskey(FqFiniteFieldID, (char, deg, s))
         return FqFiniteFieldID[char, deg, s]
      else
         d = new()
         finalizer(_FqFiniteField_clear_fn, d)
         ccall((:fq_ctx_init, :libflint), Nothing,
               (Ref{FqFiniteField}, Ref{fmpz}, Int, Ptr{UInt8}),
                  d, char, deg, string(s))
         if cached
            FqFiniteFieldID[char, deg, s] = d
         end
         return d
      end
   end

   function FqFiniteField(f::fmpz_mod_poly, s::Symbol, cached::Bool = true)
      if haskey(FqFiniteFieldIDPol, (f, s))
         return FqFiniteFieldIDPol[f, s]
      else
         z = new()
         ccall((:fq_ctx_init_modulus, :libflint), Nothing,
               (Ref{FqFiniteField}, Ref{fmpz_mod_poly}, Ptr{UInt8}),
                  z, f, string(s))
         if cached
            FqFiniteFieldIDPol[f, s] = z
         end
         finalizer(_FqFiniteField_clear_fn, z)
         return z
      end
   end
end

const FqFiniteFieldID = Dict{Tuple{fmpz, Int, Symbol}, FqFiniteField}()

const FqFiniteFieldIDPol = Dict{Tuple{fmpz_mod_poly, Symbol}, FqFiniteField}()

function _FqFiniteField_clear_fn(a :: FqFiniteField)
   ccall((:fq_ctx_clear, :libflint), Nothing, (Ref{FqFiniteField},), a)
end

mutable struct fq <: FinFieldElem
   coeffs :: Ptr{Nothing}
   alloc :: Int
   length :: Int
   parent::FqFiniteField

   function fq(ctx::FqFiniteField)
      d = new()
      ccall((:fq_init2, :libflint), Nothing,
            (Ref{fq}, Ref{FqFiniteField}), d, ctx)
      finalizer(_fq_clear_fn, d)
      d.parent = ctx
      return d
   end

   function fq(ctx::FqFiniteField, x::Int)
      d = new()
      ccall((:fq_init2, :libflint), Nothing,
            (Ref{fq}, Ref{FqFiniteField}), d, ctx)
      finalizer(_fq_clear_fn, d)
      ccall((:fq_set_si, :libflint), Nothing,
                (Ref{fq}, Int, Ref{FqFiniteField}), d, x, ctx)
      d.parent = ctx
      return d
   end

   function fq(ctx::FqFiniteField, x::fmpz)
      d = new()
      ccall((:fq_init2, :libflint), Nothing,
            (Ref{fq}, Ref{FqFiniteField}), d, ctx)
      finalizer(_fq_clear_fn, d)
      ccall((:fq_set_fmpz, :libflint), Nothing,
            (Ref{fq}, Ref{fmpz}, Ref{FqFiniteField}), d, x, ctx)
      d.parent = ctx
      return d
   end

   function fq(ctx::FqFiniteField, x::fq)
      d = new()
      ccall((:fq_init2, :libflint), Nothing,
            (Ref{fq}, Ref{FqFiniteField}), d, ctx)
      finalizer(_fq_clear_fn, d)
      ccall((:fq_set, :libflint), Nothing,
            (Ref{fq}, Ref{fq}, Ref{FqFiniteField}), d, x, ctx)
      d.parent = ctx
      return d
   end
end

function _fq_clear_fn(a::fq)
   ccall((:fq_clear, :libflint), Nothing,
         (Ref{fq}, Ref{FqFiniteField}), a, a.parent)
end

###############################################################################
#
#   FlintPadicField / padic
#
###############################################################################


mutable struct FlintPadicField <: Field
   p::Int
   pinv::Float64
   pow::Ptr{Nothing}
   minpre::Int
   maxpre::Int
   mode::Int
   prec_max::Int

   function FlintPadicField(p::fmpz, prec::Int)
      !isprime(p) && error("Prime base required in FlintPadicField")
      d = new()
      ccall((:padic_ctx_init, :libflint), Nothing,
           (Ref{FlintPadicField}, Ref{fmpz}, Int, Int, Cint),
                                     d, p, 0, 0, 1)
      finalizer(_padic_ctx_clear_fn, d)
      d.prec_max = prec
      return d
   end
end

const PadicBase = Dict{Tuple{fmpz, Int}, FlintPadicField}()

function _padic_ctx_clear_fn(a::FlintPadicField)
   ccall((:padic_ctx_clear, :libflint), Nothing, (Ref{FlintPadicField},), a)
end

mutable struct padic <: FieldElem
   u :: Int
   v :: Int
   N :: Int
   parent::FlintPadicField

   function padic(prec::Int)
      d = new()
      ccall((:padic_init2, :libflint), Nothing, (Ref{padic}, Int), d, prec)
      finalizer(_padic_clear_fn, d)
      return d
   end
end

function _padic_clear_fn(a::padic)
   ccall((:padic_clear, :libflint), Nothing, (Ref{padic},), a)
end

###############################################################################
#
#   FmpzRelSeriesRing / fmpz_rel_series
#
###############################################################################

mutable struct FmpzRelSeriesRing <: SeriesRing{fmpz}
   base_ring::FlintIntegerRing
   prec_max::Int
   S::Symbol

   function FmpzRelSeriesRing(prec::Int, s::Symbol, cached::Bool = true)
      if haskey(FmpzRelSeriesID, (prec, s))
         FmpzRelSeriesID[prec, s]
      else
         z = new(FlintZZ, prec, s)
         if cached
            FmpzRelSeriesID[prec, s] = z
         end
         return z
      end
   end
end

const FmpzRelSeriesID = Dict{Tuple{Int, Symbol}, FmpzRelSeriesRing}()

mutable struct fmpz_rel_series <: RelSeriesElem{fmpz}
   coeffs::Ptr{Nothing}
   alloc::Int
   length::Int
   prec::Int
   val::Int
   parent::FmpzRelSeriesRing

   function fmpz_rel_series()
      z = new()
      ccall((:fmpz_poly_init, :libflint), Nothing,
            (Ref{fmpz_rel_series},), z)
      finalizer(_fmpz_rel_series_clear_fn, z)
      return z
   end

   function fmpz_rel_series(a::Array{fmpz, 1}, len::Int, prec::Int, val::Int)
      z = new()
      ccall((:fmpz_poly_init2, :libflint), Nothing,
            (Ref{fmpz_rel_series}, Int), z, len)
      for i = 1:len
         ccall((:fmpz_poly_set_coeff_fmpz, :libflint), Nothing,
                     (Ref{fmpz_rel_series}, Int, Ref{fmpz}), z, i - 1, a[i])
      end
      z.prec = prec
      z.val = val
      finalizer(_fmpz_rel_series_clear_fn, z)
      return z
   end

   function fmpz_rel_series(a::fmpz_rel_series)
      z = new()
      ccall((:fmpz_poly_init, :libflint), Nothing, (Ref{fmpz_rel_series},), z)
      ccall((:fmpz_poly_set, :libflint), Nothing,
            (Ref{fmpz_rel_series}, Ref{fmpz_rel_series}), z, a)
      finalizer(_fmpz_rel_series_clear_fn, z)
      return z
   end
end

function _fmpz_rel_series_clear_fn(a::fmpz_rel_series)
   ccall((:fmpz_poly_clear, :libflint), Nothing, (Ref{fmpz_rel_series},), a)
end

###############################################################################
#
#   FmpzAbsSeriesRing / fmpz_abs_series
#
###############################################################################

mutable struct FmpzAbsSeriesRing <: SeriesRing{fmpz}
   base_ring::FlintIntegerRing
   prec_max::Int
   S::Symbol

   function FmpzAbsSeriesRing(prec::Int, s::Symbol, cached::Bool = true)
      if haskey(FmpzAbsSeriesID, (prec, s))
         FmpzAbsSeriesID[prec, s]
      else
         z = new(FlintZZ, prec, s)
         if cached
            FmpzAbsSeriesID[prec, s] = z
         end
         return z
      end
   end
end

const FmpzAbsSeriesID = Dict{Tuple{Int, Symbol}, FmpzAbsSeriesRing}()

mutable struct fmpz_abs_series <: AbsSeriesElem{fmpz}
   coeffs::Ptr{Nothing}
   alloc::Int
   length::Int
   prec :: Int
   parent::FmpzAbsSeriesRing

   function fmpz_abs_series()
      z = new()
      ccall((:fmpz_poly_init, :libflint), Nothing,
            (Ref{fmpz_abs_series},), z)
      finalizer(_fmpz_abs_series_clear_fn, z)
      return z
   end

   function fmpz_abs_series(a::Array{fmpz, 1}, len::Int, prec::Int)
      z = new()
      ccall((:fmpz_poly_init2, :libflint), Nothing,
            (Ref{fmpz_abs_series}, Int), z, len)
      for i = 1:len
         ccall((:fmpz_poly_set_coeff_fmpz, :libflint), Nothing,
                     (Ref{fmpz_abs_series}, Int, Ref{fmpz}), z, i - 1, a[i])
      end
      z.prec = prec
      finalizer(_fmpz_abs_series_clear_fn, z)
      return z
   end

   function fmpz_abs_series(a::fmpz_abs_series)
      z = new()
      ccall((:fmpz_poly_init, :libflint), Nothing, (Ref{fmpz_abs_series},), z)
      ccall((:fmpz_poly_set, :libflint), Nothing,
            (Ref{fmpz_abs_series}, Ref{fmpz_abs_series}), z, a)
      finalizer(_fmpz_abs_series_clear_fn, z)
      return z
   end
end

function _fmpz_abs_series_clear_fn(a::fmpz_abs_series)
   ccall((:fmpz_poly_clear, :libflint), Nothing, (Ref{fmpz_abs_series},), a)
end

###############################################################################
#
#   FlintPuiseuxSeriesRing / FlintPuiseuxSeriesRingElem
#
###############################################################################

mutable struct FlintPuiseuxSeriesRing{T <: RingElem} <: Ring where T
   laurent_ring::Ring

   function FlintPuiseuxSeriesRing{T}(R::Ring, cached::Bool = true) where T
      if haskey(FlintPuiseuxSeriesID, R)
         return FlintPuiseuxSeriesID[R]::FlintPuiseuxSeriesRing{T}
      else
         z = new{T}(R)
         if cached
            FlintPuiseuxSeriesID[R] = z
         end
         return z
      end
   end
end

const FlintPuiseuxSeriesID = Dict{Ring, Ring}()

mutable struct FlintPuiseuxSeriesRingElem{T <: RingElem} <: RingElem
   data::T
   scale::Int
   parent::FlintPuiseuxSeriesRing{T}

   function FlintPuiseuxSeriesRingElem{T}(d::T, scale::Int) where T <:
RingElem
      new{T}(d, scale)
   end
end

###############################################################################
#
#   FlintPuiseuxSeriesField / FlintPuiseuxSeriesFieldElem
#
###############################################################################

mutable struct FlintPuiseuxSeriesField{T <: RingElem} <: Field
   laurent_ring::Ring

   function FlintPuiseuxSeriesField{T}(R::Field, cached::Bool = true) where T
      if haskey(FlintPuiseuxSeriesID, R)
         return FlintPuiseuxSeriesID[R]::FlintPuiseuxSeriesField{T}
      else
         z = new{T}(R)
         if cached
            FlintPuiseuxSeriesID[R] = z
         end
         return z
      end
   end
end

mutable struct FlintPuiseuxSeriesFieldElem{T <: RingElem} <: FieldElem
   data::T
   scale::Int
   parent::FlintPuiseuxSeriesField{T}

   function FlintPuiseuxSeriesFieldElem{T}(d::T, scale::Int) where T <:
RingElem
      new{T}(d, scale)
   end
end

const FlintPuiseuxSeriesElem{T} = Union{FlintPuiseuxSeriesRingElem{T}, FlintPuiseuxSeriesFieldElem{T}} where T <: RingElem

###############################################################################
#
#   FmpzLaurentSeriesRing / fmpz_laurent_series
#
###############################################################################

mutable struct FmpzLaurentSeriesRing <: Ring
   base_ring::FlintIntegerRing
   prec_max::Int
   S::Symbol

   function FmpzLaurentSeriesRing(prec::Int, s::Symbol, cached::Bool = true)
      if haskey(FmpzLaurentSeriesID, (prec, s))
         FmpzLaurentSeriesID[prec, s]
      else
         z = new(FlintZZ, prec, s)
         if cached
            FmpzLaurentSeriesID[prec, s] = z
         end
         return z
      end
   end
end

const FmpzLaurentSeriesID = Dict{Tuple{Int, Symbol}, FmpzLaurentSeriesRing}()

mutable struct fmpz_laurent_series <: RingElem
   coeffs::Ptr{Nothing}
   alloc::Int
   length::Int
   prec::Int
   val::Int
   scale::Int
   parent::FmpzLaurentSeriesRing

   function fmpz_laurent_series()
      z = new()
      ccall((:fmpz_poly_init, :libflint), Nothing,
            (Ref{fmpz_laurent_series},), z)
      finalizer(_fmpz_laurent_series_clear_fn, z)
      return z
   end

   function fmpz_laurent_series(a::Array{fmpz, 1}, len::Int, prec::Int, val::Int, scale::Int)
      z = new()
      ccall((:fmpz_poly_init2, :libflint), Nothing,
            (Ref{fmpz_laurent_series}, Int), z, len)
      for i = 1:len
         ccall((:fmpz_poly_set_coeff_fmpz, :libflint), Nothing,
                     (Ref{fmpz_laurent_series}, Int, Ref{fmpz}), z, i - 1, a[i])
      end
      z.prec = prec
      z.val = val
      z.scale = scale
      finalizer(_fmpz_laurent_series_clear_fn, z)
      return z
   end

   function fmpz_laurent_series(a::fmpz_laurent_series)
      z = new()
      ccall((:fmpz_poly_init, :libflint), Nothing, (Ref{fmpz_laurent_series},), z)
      ccall((:fmpz_poly_set, :libflint), Nothing,
            (Ref{fmpz_laurent_series}, Ref{fmpz_laurent_series}), z, a)
      finalizer(_fmpz_laurent_series_clear_fn, z)
      return z
   end
end

function _fmpz_laurent_series_clear_fn(a::fmpz_laurent_series)
   ccall((:fmpz_poly_clear, :libflint), Nothing, (Ref{fmpz_laurent_series},), a)
end

###############################################################################
#
#   FmpqRelSeriesRing / fmpq_rel_series
#
###############################################################################

mutable struct FmpqRelSeriesRing <: SeriesRing{fmpq}
   base_ring::FlintRationalField
   prec_max::Int
   S::Symbol

   function FmpqRelSeriesRing(prec::Int, s::Symbol, cached::Bool = true)
      if haskey(FmpqRelSeriesID, (prec, s))
         return FmpqRelSeriesID[prec, s]
      else
         z = new(FlintQQ, prec, s)
         if cached
            FmpqRelSeriesID[prec, s] = z
         end
         return z
      end
   end
end

const FmpqRelSeriesID = Dict{Tuple{Int, Symbol}, FmpqRelSeriesRing}()

mutable struct fmpq_rel_series <: RelSeriesElem{fmpq}
   coeffs::Ptr{Nothing}
   den::Int
   alloc::Int
   length::Int
   prec::Int
   val::Int
   parent::FmpqRelSeriesRing

   function fmpq_rel_series()
      z = new()
      ccall((:fmpq_poly_init, :libflint), Nothing,
            (Ref{fmpq_rel_series},), z)
      finalizer(_fmpq_rel_series_clear_fn, z)
      return z
   end

   function fmpq_rel_series(a::Array{fmpq, 1}, len::Int, prec::Int, val::Int)
      z = new()
      ccall((:fmpq_poly_init2, :libflint), Nothing,
            (Ref{fmpq_rel_series}, Int), z, len)
      for i = 1:len
         ccall((:fmpq_poly_set_coeff_fmpq, :libflint), Nothing,
                     (Ref{fmpq_rel_series}, Int, Ref{fmpq}), z, i - 1, a[i])
      end
      z.prec = prec
      z.val = val
      finalizer(_fmpq_rel_series_clear_fn, z)
      return z
   end

   function fmpq_rel_series(a::fmpq_rel_series)
      z = new()
      ccall((:fmpq_poly_init, :libflint), Nothing, (Ref{fmpq_rel_series},), z)
      ccall((:fmpq_poly_set, :libflint), Nothing,
            (Ref{fmpq_rel_series}, Ref{fmpq_rel_series}), z, a)
      finalizer(_fmpq_rel_series_clear_fn, z)
      return z
   end
end

function _fmpq_rel_series_clear_fn(a::fmpq_rel_series)
   ccall((:fmpq_poly_clear, :libflint), Nothing, (Ref{fmpq_rel_series},), a)
end

###############################################################################
#
#   FmpqAbsSeriesRing / fmpq_abs_series
#
###############################################################################

mutable struct FmpqAbsSeriesRing <: SeriesRing{fmpq}
   base_ring::FlintRationalField
   prec_max::Int
   S::Symbol

   function FmpqAbsSeriesRing(prec::Int, s::Symbol, cached::Bool = true)
      if haskey(FmpqAbsSeriesID, (prec, s))
         return FmpqAbsSeriesID[prec, s]
      else
         z = new(FlintQQ, prec, s)
         if cached
            FmpqAbsSeriesID[prec, s] = z
         end
         return z
      end
   end
end

const FmpqAbsSeriesID = Dict{Tuple{Int, Symbol}, FmpqAbsSeriesRing}()

mutable struct fmpq_abs_series <: AbsSeriesElem{fmpq}
   coeffs::Ptr{Nothing}
   den::Int
   alloc::Int
   length::Int
   prec :: Int
   parent::FmpqAbsSeriesRing

   function fmpq_abs_series()
      z = new()
      ccall((:fmpq_poly_init, :libflint), Nothing,
            (Ref{fmpq_abs_series},), z)
      finalizer(_fmpq_abs_series_clear_fn, z)
      return z
   end

   function fmpq_abs_series(a::Array{fmpq, 1}, len::Int, prec::Int)
      z = new()
      ccall((:fmpq_poly_init2, :libflint), Nothing,
            (Ref{fmpq_abs_series}, Int), z, len)
      for i = 1:len
         ccall((:fmpq_poly_set_coeff_fmpq, :libflint), Nothing,
                     (Ref{fmpq_abs_series}, Int, Ref{fmpq}), z, i - 1, a[i])
      end
      z.prec = prec
      finalizer(_fmpq_abs_series_clear_fn, z)
      return z
   end

   function fmpq_abs_series(a::fmpq_abs_series)
      z = new()
      ccall((:fmpq_poly_init, :libflint), Nothing, (Ref{fmpq_abs_series},), z)
      ccall((:fmpq_poly_set, :libflint), Nothing,
            (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}), z, a)
      finalizer(_fmpq_abs_series_clear_fn, z)
      return z
   end
end

function _fmpq_abs_series_clear_fn(a::fmpq_abs_series)
   ccall((:fmpq_poly_clear, :libflint), Nothing, (Ref{fmpq_abs_series},), a)
end

###############################################################################
#
#   NmodRelSeriesRing / nmod_rel_series
#
###############################################################################

mutable struct NmodRelSeriesRing <: SeriesRing{nmod}
   base_ring::NmodRing
   prec_max::Int
   S::Symbol

   function NmodRelSeriesRing(R::NmodRing, prec::Int, s::Symbol,
                                 cached::Bool = true)
      if haskey(NmodRelSeriesID, (R, prec, s))
         return NmodRelSeriesID[R, prec, s]
      else
         z = new(R, prec, s)
         if cached
            NmodRelSeriesID[R, prec, s] = z
         end
         return z
      end
   end
end

const NmodRelSeriesID = Dict{Tuple{NmodRing, Int, Symbol},
                                NmodRelSeriesRing}()

mutable struct nmod_rel_series <: RelSeriesElem{nmod}
   coeffs::Ptr{Nothing}
   alloc::Int
   length::Int
   mod_n::UInt
   mod_ninv::UInt
   mod_norm::UInt
   prec::Int
   val::Int
   parent::NmodRelSeriesRing

   function nmod_rel_series(p::UInt)
      z = new()
      ccall((:nmod_poly_init, :libflint), Nothing,
            (Ref{nmod_rel_series}, UInt), z, p)
      finalizer(_nmod_rel_series_clear_fn, z)
      return z
   end

   function nmod_rel_series(p::UInt, a::Array{fmpz, 1}, len::Int, prec::Int, val::Int)
      z = new()
      ccall((:nmod_poly_init2, :libflint), Nothing,
            (Ref{nmod_rel_series}, UInt, Int), z, p, len)
      for i = 1:len
         tt = ccall((:fmpz_fdiv_ui, :libflint), UInt, (Ref{fmpz}, UInt), a[i], p)
         ccall((:nmod_poly_set_coeff_ui, :libflint), Nothing,
                     (Ref{nmod_rel_series}, Int, UInt), z, i - 1, tt)
      end
      z.prec = prec
      z.val = val
      finalizer(_nmod_rel_series_clear_fn, z)
      return z
   end

   function nmod_rel_series(p::UInt, a::Array{UInt, 1}, len::Int, prec::Int, val::Int)
      z = new()
      ccall((:nmod_poly_init2, :libflint), Nothing,
            (Ref{nmod_rel_series}, UInt, Int), z, p, len)
      for i = 1:len
         ccall((:nmod_poly_set_coeff_ui, :libflint), Nothing,
                     (Ref{nmod_rel_series}, Int, UInt), z, i - 1, a[i])
      end
      z.prec = prec
      z.val = val
      finalizer(_nmod_rel_series_clear_fn, z)
      return z
   end

   function nmod_rel_series(p::UInt, a::Array{nmod, 1}, len::Int, prec::Int, val::Int)
      z = new()
      ccall((:nmod_poly_init2, :libflint), Nothing,
            (Ref{nmod_rel_series}, UInt, Int), z, p, len)
      for i = 1:len
         ccall((:nmod_poly_set_coeff_ui, :libflint), Nothing,
                     (Ref{nmod_rel_series}, Int, UInt), z, i - 1, data(a[i]))
      end
      z.prec = prec
      z.val = val
      finalizer(_nmod_rel_series_clear_fn, z)
      return z
   end

   function nmod_rel_series(a::nmod_rel_series)
      z = new()
      p = modulus(base_ring(parent(a)))
      ccall((:nmod_poly_init, :libflint), Nothing,
            (Ref{nmod_rel_series}, UInt), z, p)
      ccall((:nmod_poly_set, :libflint), Nothing,
            (Ref{nmod_rel_series}, Ref{nmod_rel_series}), z, a)
      finalizer(_nmod_rel_series_clear_fn, z)
      return z
   end
end

function _nmod_rel_series_clear_fn(a::nmod_rel_series)
   ccall((:nmod_poly_clear, :libflint), Nothing, (Ref{nmod_rel_series},), a)
end

###############################################################################
#
#   FmpzModRelSeriesRing / fmpz_mod_rel_series
#
###############################################################################

mutable struct FmpzModRelSeriesRing <: SeriesRing{Generic.Res{fmpz}}
   base_ring::Generic.ResRing{fmpz}
   prec_max::Int
   S::Symbol

   function FmpzModRelSeriesRing(R::Ring, prec::Int, s::Symbol,
                                 cached::Bool = true)
      if haskey(FmpzModRelSeriesID, (R, prec, s))
         return FmpzModRelSeriesID[R, prec, s]
      else
         z = new(R, prec, s)
         if cached
            FmpzModRelSeriesID[R, prec, s] = z
         end
         return z
      end
   end
end

const FmpzModRelSeriesID = Dict{Tuple{Generic.ResRing{fmpz}, Int, Symbol},
                                FmpzModRelSeriesRing}()

mutable struct fmpz_mod_rel_series <: RelSeriesElem{Generic.Res{fmpz}}
   coeffs::Ptr{Nothing}
   alloc::Int
   length::Int
   p::Int
   prec::Int
   val::Int
   parent::FmpzModRelSeriesRing

   function fmpz_mod_rel_series(p::fmpz)
      z = new()
      ccall((:fmpz_mod_poly_init, :libflint), Nothing,
            (Ref{fmpz_mod_rel_series}, Ref{fmpz}), z, p)
      finalizer(_fmpz_mod_rel_series_clear_fn, z)
      return z
   end

   function fmpz_mod_rel_series(p::fmpz, a::Array{fmpz, 1}, len::Int, prec::Int, val::Int)
      z = new()
      ccall((:fmpz_mod_poly_init2, :libflint), Nothing,
            (Ref{fmpz_mod_rel_series}, Ref{fmpz}, Int), z, p, len)
      for i = 1:len
         ccall((:fmpz_mod_poly_set_coeff_fmpz, :libflint), Nothing,
                     (Ref{fmpz_mod_rel_series}, Int, Ref{fmpz}), z, i - 1, a[i])
      end
      z.prec = prec
      z.val = val
      finalizer(_fmpz_mod_rel_series_clear_fn, z)
      return z
   end

   function fmpz_mod_rel_series(p::fmpz, a::Array{Generic.Res{fmpz}, 1}, len::Int, prec::Int, val::Int)
      z = new()
      ccall((:fmpz_mod_poly_init2, :libflint), Nothing,
            (Ref{fmpz_mod_rel_series}, Ref{fmpz}, Int), z, p, len)
      for i = 1:len
         ccall((:fmpz_mod_poly_set_coeff_fmpz, :libflint), Nothing,
                     (Ref{fmpz_mod_rel_series}, Int, Ref{fmpz}), z, i - 1, data(a[i]))
      end
      z.prec = prec
      z.val = val
      finalizer(_fmpz_mod_rel_series_clear_fn, z)
      return z
   end

   function fmpz_mod_rel_series(a::fmpz_mod_rel_series)
      z = new()
      p = modulus(base_ring(parent(a)))
      ccall((:fmpz_mod_poly_init, :libflint), Nothing,
            (Ref{fmpz_mod_rel_series}, Ref{fmpz}), z, p)
      ccall((:fmpz_mod_poly_set, :libflint), Nothing,
            (Ref{fmpz_mod_rel_series}, Ref{fmpz_mod_rel_series}), z, a)
      finalizer(_fmpz_mod_rel_series_clear_fn, z)
      return z
   end
end

function _fmpz_mod_rel_series_clear_fn(a::fmpz_mod_rel_series)
   ccall((:fmpz_mod_poly_clear, :libflint), Nothing, (Ref{fmpz_mod_rel_series},), a)
end

###############################################################################
#
#   FmpzModAbsSeriesRing / fmpz_mod_abs_series
#
###############################################################################

mutable struct FmpzModAbsSeriesRing <: SeriesRing{Generic.Res{fmpz}}
   base_ring::Generic.ResRing{fmpz}
   prec_max::Int
   S::Symbol

   function FmpzModAbsSeriesRing(R::Ring, prec::Int, s::Symbol,
                                 cached::Bool = true)
      if haskey(FmpzModAbsSeriesID, (R, prec, s))
         return FmpzModAbsSeriesID[R, prec, s]
      else
         z = new(R, prec, s)
         if cached
            FmpzModAbsSeriesID[R, prec, s]  = z
         end
         return z
      end
   end
end

const FmpzModAbsSeriesID = Dict{Tuple{Generic.ResRing{fmpz}, Int, Symbol},
                                FmpzModAbsSeriesRing}()

mutable struct fmpz_mod_abs_series <: AbsSeriesElem{Generic.Res{fmpz}}
   coeffs::Ptr{Nothing}
   alloc::Int
   length::Int
   p::Int
   prec::Int
   parent::FmpzModAbsSeriesRing

   function fmpz_mod_abs_series(p::fmpz)
      z = new()
      ccall((:fmpz_mod_poly_init, :libflint), Nothing,
            (Ref{fmpz_mod_abs_series}, Ref{fmpz}), z, p)
      finalizer(_fmpz_mod_abs_series_clear_fn, z)
      return z
   end

   function fmpz_mod_abs_series(p::fmpz, a::Array{fmpz, 1}, len::Int, prec::Int)
      z = new()
      ccall((:fmpz_mod_poly_init2, :libflint), Nothing,
            (Ref{fmpz_mod_abs_series}, Ref{fmpz}, Int), z, p, len)
      for i = 1:len
         ccall((:fmpz_mod_poly_set_coeff_fmpz, :libflint), Nothing,
                     (Ref{fmpz_mod_abs_series}, Int, Ref{fmpz}), z, i - 1, a[i])
      end
      z.prec = prec
      finalizer(_fmpz_mod_abs_series_clear_fn, z)
      return z
   end

   function fmpz_mod_abs_series(p::fmpz, a::Array{Generic.Res{fmpz}, 1}, len::Int, prec::Int)
      z = new()
      ccall((:fmpz_mod_poly_init2, :libflint), Nothing,
            (Ref{fmpz_mod_abs_series}, Ref{fmpz}, Int), z, p, len)
      for i = 1:len
         ccall((:fmpz_mod_poly_set_coeff_fmpz, :libflint), Nothing,
                     (Ref{fmpz_mod_abs_series}, Int, Ref{fmpz}), z, i - 1, data(a[i]))
      end
      z.prec = prec
      finalizer(_fmpz_mod_abs_series_clear_fn, z)
      return z
   end

   function fmpz_mod_abs_series(a::fmpz_mod_abs_series)
      z = new()
      p = modulus(base_ring(parent(a)))
      ccall((:fmpz_mod_poly_init, :libflint), Nothing,
            (Ref{fmpz_mod_abs_series}, Ref{fmpz}), z, p)
      ccall((:fmpz_mod_poly_set, :libflint), Nothing,
            (Ref{fmpz_mod_abs_series}, Ref{fmpz_mod_abs_series}), z, a)
      finalizer(_fmpz_mod_abs_series_clear_fn, z)
      return z
   end
end

function _fmpz_mod_abs_series_clear_fn(a::fmpz_mod_abs_series)
   ccall((:fmpz_mod_poly_clear, :libflint), Nothing, (Ref{fmpz_mod_abs_series},), a)
end

###############################################################################
#
#   FqRelSeriesRing / fq_rel_series
#
###############################################################################

mutable struct FqRelSeriesRing <: SeriesRing{fq}
   base_ring::FqFiniteField
   prec_max::Int
   S::Symbol

   function FqRelSeriesRing(R::FqFiniteField, prec::Int, s::Symbol,
                            cached::Bool = true)
      if haskey(FqRelSeriesID, (R, prec, s))
         return FqRelSeriesID[R, prec, s]
      else
         z = new(R, prec, s)
         if cached
            FqRelSeriesID[R, prec, s] = z
         end
         return z
      end
   end
end

const FqRelSeriesID = Dict{Tuple{FqFiniteField, Int, Symbol}, FqRelSeriesRing}()

mutable struct fq_rel_series <: RelSeriesElem{fq}
   coeffs::Ptr{Nothing}
   alloc::Int
   length::Int
   prec::Int
   val::Int
   parent::FqRelSeriesRing

   function fq_rel_series(ctx::FqFiniteField)
      z = new()
      ccall((:fq_poly_init, :libflint), Nothing,
            (Ref{fq_rel_series}, Ref{FqFiniteField}), z, ctx)
      finalizer(_fq_rel_series_clear_fn, z)
      return z
   end

   function fq_rel_series(ctx::FqFiniteField, a::Array{fq, 1}, len::Int, prec::Int, val::Int)
      z = new()
      ccall((:fq_poly_init2, :libflint), Nothing,
            (Ref{fq_rel_series}, Int, Ref{FqFiniteField}), z, len, ctx)
      for i = 1:len
         ccall((:fq_poly_set_coeff, :libflint), Nothing,
               (Ref{fq_rel_series}, Int, Ref{fq}, Ref{FqFiniteField}),
                                               z, i - 1, a[i], ctx)
      end
      z.prec = prec
      z.val = val
      finalizer(_fq_rel_series_clear_fn, z)
      return z
   end

   function fq_rel_series(ctx::FqFiniteField, a::fq_rel_series)
      z = new()
      ccall((:fq_poly_init, :libflint), Nothing,
            (Ref{fq_rel_series}, Ref{FqFiniteField}), z, ctx)
      ccall((:fq_poly_set, :libflint), Nothing,
            (Ref{fq_rel_series}, Ref{fq_rel_series}, Ref{FqFiniteField}), z, a, ctx)
      finalizer(_fq_rel_series_clear_fn, z)
      return z
   end
end

function _fq_rel_series_clear_fn(a::fq_rel_series)
   ctx = base_ring(a)
   ccall((:fq_poly_clear, :libflint), Nothing,
         (Ref{fq_rel_series}, Ref{FqFiniteField}), a, ctx)
end

###############################################################################
#
#   FqAbsSeriesRing / fq_abs_series
#
###############################################################################

mutable struct FqAbsSeriesRing <: SeriesRing{fq}
   base_ring::FqFiniteField
   prec_max::Int
   S::Symbol

   function FqAbsSeriesRing(R::FqFiniteField, prec::Int, s::Symbol,
                            cached::Bool = true)
      if haskey(FqAbsSeriesID, (R, prec, s))
         return FqAbsSeriesID[R, prec, s]
      else
         z = new(R, prec, s)
         if cached
            FqAbsSeriesID[R, prec, s] = z
         end
         return z
      end
   end
end

const FqAbsSeriesID = Dict{Tuple{FqFiniteField, Int, Symbol}, FqAbsSeriesRing}()

mutable struct fq_abs_series <: AbsSeriesElem{fq}
   coeffs::Ptr{Nothing}
   alloc::Int
   length::Int
   prec::Int
   parent::FqAbsSeriesRing

   function fq_abs_series(ctx::FqFiniteField)
      z = new()
      ccall((:fq_poly_init, :libflint), Nothing,
            (Ref{fq_abs_series}, Ref{FqFiniteField}), z, ctx)
      finalizer(_fq_abs_series_clear_fn, z)
      return z
   end

   function fq_abs_series(ctx::FqFiniteField, a::Array{fq, 1}, len::Int, prec::Int)
      z = new()
      ccall((:fq_poly_init2, :libflint), Nothing,
            (Ref{fq_abs_series}, Int, Ref{FqFiniteField}), z, len, ctx)
      for i = 1:len
         ccall((:fq_poly_set_coeff, :libflint), Nothing,
               (Ref{fq_abs_series}, Int, Ref{fq}, Ref{FqFiniteField}),
                                               z, i - 1, a[i], ctx)
      end
      z.prec = prec
      finalizer(_fq_abs_series_clear_fn, z)
      return z
   end

   function fq_abs_series(ctx::FqFiniteField, a::fq_abs_series)
      z = new()
      ccall((:fq_poly_init, :libflint), Nothing,
            (Ref{fq_abs_series}, Ref{FqFiniteField}), z, ctx)
      ccall((:fq_poly_set, :libflint), Nothing,
            (Ref{fq_abs_series}, Ref{fq_abs_series}, Ref{FqFiniteField}), z, a, ctx)
      finalizer(_fq_abs_series_clear_fn, z)
      return z
   end
end

function _fq_abs_series_clear_fn(a::fq_abs_series)
   ctx = base_ring(a)
   ccall((:fq_poly_clear, :libflint), Nothing,
         (Ref{fq_abs_series}, Ref{FqFiniteField}), a, ctx)
end

###############################################################################
#
#   FqNmodRelSeriesRing / fq_nmod_rel_series
#
###############################################################################

mutable struct FqNmodRelSeriesRing <: SeriesRing{fq_nmod}
   base_ring::FqNmodFiniteField
   prec_max::Int
   S::Symbol

   function FqNmodRelSeriesRing(R::FqNmodFiniteField, prec::Int, s::Symbol,
                                cached::Bool = true)
      if haskey(FqNmodRelSeriesID, (R, prec, s))
         return FqNmodRelSeriesID[R, prec, s]
      else
         z = new(R, prec, s)
         if cached
            FqNmodRelSeriesID[R, prec, s] = z
         end
         return z
      end
   end
end

const FqNmodRelSeriesID = Dict{Tuple{FqNmodFiniteField, Int, Symbol},
                               FqNmodRelSeriesRing}()

mutable struct fq_nmod_rel_series <: RelSeriesElem{fq_nmod}
   coeffs::Ptr{Nothing}
   alloc::Int
   length::Int
   prec::Int
   val::Int
   parent::FqNmodRelSeriesRing

   function fq_nmod_rel_series(ctx::FqNmodFiniteField)
      z = new()
      ccall((:fq_nmod_poly_init, :libflint), Nothing,
            (Ref{fq_nmod_rel_series}, Ref{FqNmodFiniteField}), z, ctx)
      finalizer(_fq_nmod_rel_series_clear_fn, z)
      return z
   end

   function fq_nmod_rel_series(ctx::FqNmodFiniteField, a::Array{fq_nmod, 1}, len::Int, prec::Int, val::Int)
      z = new()
      ccall((:fq_nmod_poly_init2, :libflint), Nothing,
            (Ref{fq_nmod_rel_series}, Int, Ref{FqNmodFiniteField}), z, len, ctx)
      for i = 1:len
         ccall((:fq_nmod_poly_set_coeff, :libflint), Nothing,
               (Ref{fq_nmod_rel_series}, Int, Ref{fq_nmod}, Ref{FqNmodFiniteField}),
                                               z, i - 1, a[i], ctx)
      end
      z.prec = prec
      z.val = val
      finalizer(_fq_nmod_rel_series_clear_fn, z)
      return z
   end

   function fq_nmod_rel_series(ctx::FqNmodFiniteField, a::fq_nmod_rel_series)
      z = new()
      ccall((:fq_nmod_poly_init, :libflint), Nothing,
            (Ref{fq_nmod_rel_series}, Ref{FqNmodFiniteField}), z, ctx)
      ccall((:fq_nmod_poly_set, :libflint), Nothing,
            (Ref{fq_nmod_rel_series}, Ref{fq_nmod_rel_series}, Ref{FqNmodFiniteField}), z, a, ctx)
      finalizer(_fq_nmod_rel_series_clear_fn, z)
      return z
   end
end

function _fq_nmod_rel_series_clear_fn(a::fq_nmod_rel_series)
   ctx = base_ring(a)
   ccall((:fq_nmod_poly_clear, :libflint), Nothing,
         (Ref{fq_nmod_rel_series}, Ref{FqNmodFiniteField}), a, ctx)
end

###############################################################################
#
#   FqNmodAbsSeriesRing / fq_nmod_abs_series
#
###############################################################################

mutable struct FqNmodAbsSeriesRing <: SeriesRing{fq_nmod}
   base_ring::FqNmodFiniteField
   prec_max::Int
   S::Symbol

   function FqNmodAbsSeriesRing(R::FqNmodFiniteField, prec::Int, s::Symbol,
                                cached::Bool = true)
      if haskey(FqNmodAbsSeriesID, (R, prec, s))
         return FqNmodAbsSeriesID[R, prec, s]
      else
         z = new(R, prec, s)
         if cached
            FqNmodAbsSeriesID[R, prec, s] = z
         end
         return z
      end
   end
end

const FqNmodAbsSeriesID = Dict{Tuple{FqNmodFiniteField, Int, Symbol},
                               FqNmodAbsSeriesRing}()

mutable struct fq_nmod_abs_series <: AbsSeriesElem{fq_nmod}
   coeffs::Ptr{Nothing}
   alloc::Int
   length::Int
   prec::Int
   parent::FqNmodAbsSeriesRing

   function fq_nmod_abs_series(ctx::FqNmodFiniteField)
      z = new()
      ccall((:fq_nmod_poly_init, :libflint), Nothing,
            (Ref{fq_nmod_abs_series}, Ref{FqNmodFiniteField}), z, ctx)
      finalizer(_fq_nmod_abs_series_clear_fn, z)
      return z
   end

   function fq_nmod_abs_series(ctx::FqNmodFiniteField, a::Array{fq_nmod, 1}, len::Int, prec::Int)
      z = new()
      ccall((:fq_nmod_poly_init2, :libflint), Nothing,
            (Ref{fq_nmod_abs_series}, Int, Ref{FqNmodFiniteField}), z, len, ctx)
      for i = 1:len
         ccall((:fq_nmod_poly_set_coeff, :libflint), Nothing,
               (Ref{fq_nmod_abs_series}, Int, Ref{fq_nmod}, Ref{FqNmodFiniteField}),
                                               z, i - 1, a[i], ctx)
      end
      z.prec = prec
      finalizer(_fq_nmod_abs_series_clear_fn, z)
      return z
   end

   function fq_nmod_abs_series(ctx::FqNmodFiniteField, a::fq_nmod_abs_series)
      z = new()
      ccall((:fq_nmod_poly_init, :libflint), Nothing,
            (Ref{fq_nmod_abs_series}, Ref{FqNmodFiniteField}), z, ctx)
      ccall((:fq_nmod_poly_set, :libflint), Nothing,
            (Ref{fq_nmod_abs_series}, Ref{fq_nmod_abs_series}, Ref{FqNmodFiniteField}), z, a, ctx)
      finalizer(_fq_nmod_abs_series_clear_fn, z)
      return z
   end
end

function _fq_nmod_abs_series_clear_fn(a::fq_nmod_abs_series)
   ctx = base_ring(a)
   ccall((:fq_nmod_poly_clear, :libflint), Nothing,
         (Ref{fq_nmod_abs_series}, Ref{FqNmodFiniteField}), a, ctx)
end

###############################################################################
#
#   FmpqMatSpace / fmpq_mat
#
###############################################################################

# not really a mathematical ring
mutable struct FmpqMatSpace <: MatSpace{fmpq}
   rows::Int
   cols::Int
   base_ring::FlintRationalField

   function FmpqMatSpace(r::Int, c::Int, cached::Bool = true)
      if haskey(FmpqMatID, (r, c))
         return FmpqMatID[r, c]
      else
         z = new(r, c, FlintQQ)
         if cached
            FmpqMatID[r, c] = z
         end
         return z
      end
   end
end

const FmpqMatID = Dict{Tuple{Int, Int}, FmpqMatSpace}()

mutable struct fmpq_mat <: MatElem{fmpq}
   entries::Ptr{Nothing}
   r::Int
   c::Int
   rows::Ptr{Nothing}
   base_ring::FlintRationalField
   view_parent

   # used by windows, not finalised!!
   function fmpq_mat()
      return new()
   end

   function fmpq_mat(r::Int, c::Int)
      z = new()
      ccall((:fmpq_mat_init, :libflint), Nothing,
            (Ref{fmpq_mat}, Int, Int), z, r, c)
      finalizer(_fmpq_mat_clear_fn, z)
      return z
   end

   function fmpq_mat(r::Int, c::Int, arr::Array{fmpq, 2})
      z = new()
      ccall((:fmpq_mat_init, :libflint), Nothing,
            (Ref{fmpq_mat}, Int, Int), z, r, c)
      finalizer(_fmpq_mat_clear_fn, z)
      GC.@preserve z for i = 1:r
         for j = 1:c
            el = ccall((:fmpq_mat_entry, :libflint), Ptr{fmpq},
                       (Ref{fmpq_mat}, Int, Int), z, i - 1, j - 1)
            ccall((:fmpq_set, :libflint), Nothing,
                  (Ptr{fmpq}, Ref{fmpq}), el, arr[i, j])
         end
      end
      return z
   end

   function fmpq_mat(r::Int, c::Int, arr::Array{fmpz, 2})
      z = new()
      ccall((:fmpq_mat_init, :libflint), Nothing,
            (Ref{fmpq_mat}, Int, Int), z, r, c)
      finalizer(_fmpq_mat_clear_fn, z)
      b = fmpz(1)
      GC.@preserve z for i = 1:r
         for j = 1:c
            el = ccall((:fmpq_mat_entry, :libflint), Ptr{fmpq},
                       (Ref{fmpq_mat}, Int, Int), z, i - 1, j - 1)
            ccall((:fmpq_set_fmpz_frac, :libflint), Nothing,
                  (Ptr{fmpq}, Ref{fmpz}, Ref{fmpz}), el, arr[i, j], b)
         end
      end
      return z
   end


   function fmpq_mat(r::Int, c::Int, arr::Array{fmpq, 1})
      z = new()
      ccall((:fmpq_mat_init, :libflint), Nothing,
            (Ref{fmpq_mat}, Int, Int), z, r, c)
      finalizer(_fmpq_mat_clear_fn, z)
      GC.@preserve z for i = 1:r
         for j = 1:c
            el = ccall((:fmpq_mat_entry, :libflint), Ptr{fmpq},
                       (Ref{fmpq_mat}, Int, Int), z, i - 1, j - 1)
            ccall((:fmpq_set, :libflint), Nothing,
                  (Ptr{fmpq}, Ref{fmpq}), el, arr[(i-1)*c+j])
         end
      end
      return z
   end

   function fmpq_mat(r::Int, c::Int, arr::Array{fmpz, 1})
      z = new()
      ccall((:fmpq_mat_init, :libflint), Nothing,
            (Ref{fmpq_mat}, Int, Int), z, r, c)
      finalizer(_fmpq_mat_clear_fn, z)
      b = fmpz(1)
      GC.@preserve z for i = 1:r
         for j = 1:c
            el = ccall((:fmpq_mat_entry, :libflint), Ptr{fmpq},
                       (Ref{fmpq_mat}, Int, Int), z, i - 1, j - 1)
            ccall((:fmpq_set_fmpz_frac, :libflint), Nothing,
                  (Ptr{fmpq}, Ref{fmpz}, Ref{fmpz}), el, arr[(i-1)*c+j], b)
         end
      end
      return z
   end


   function fmpq_mat(r::Int, c::Int, arr::Array{T, 2}) where {T <: Integer}
      z = new()
      ccall((:fmpq_mat_init, :libflint), Nothing,
            (Ref{fmpq_mat}, Int, Int), z, r, c)
      finalizer(_fmpq_mat_clear_fn, z)
      GC.@preserve z for i = 1:r
         for j = 1:c
            el = ccall((:fmpq_mat_entry, :libflint), Ptr{fmpq},
                       (Ref{fmpq_mat}, Int, Int), z, i - 1, j - 1)
            ccall((:fmpq_set, :libflint), Nothing,
                  (Ptr{fmpq}, Ref{fmpq}), el, fmpq(arr[i, j]))
         end
      end
      return z
   end

   function fmpq_mat(r::Int, c::Int, arr::Array{T, 1}) where {T <: Integer}
      z = new()
      ccall((:fmpq_mat_init, :libflint), Nothing,
            (Ref{fmpq_mat}, Int, Int), z, r, c)
      finalizer(_fmpq_mat_clear_fn, z)
      GC.@preserve z for i = 1:r
         for j = 1:c
            el = ccall((:fmpq_mat_entry, :libflint), Ptr{fmpq},
                       (Ref{fmpq_mat}, Int, Int), z, i - 1, j - 1)
            ccall((:fmpq_set, :libflint), Nothing,
                  (Ptr{fmpq}, Ref{fmpq}), el, fmpq(arr[(i-1)*c+j]))
         end
      end
      return z
   end

   function fmpq_mat(r::Int, c::Int, d::fmpq)
      z = new()
      ccall((:fmpq_mat_init, :libflint), Nothing,
            (Ref{fmpq_mat}, Int, Int), z, r, c)
      finalizer(_fmpq_mat_clear_fn, z)
      GC.@preserve z for i = 1:min(r, c)
         el = ccall((:fmpq_mat_entry, :libflint), Ptr{fmpq},
                    (Ref{fmpq_mat}, Int, Int), z, i - 1, i - 1)
         ccall((:fmpq_set, :libflint), Nothing,
               (Ptr{fmpq}, Ref{fmpq}), el, d)
      end
      return z
   end

   function fmpq_mat(m::fmpq_mat)
      z = new()
      ccall((:fmpq_mat_init_set, :libflint), Nothing,
            (Ref{fmpq_mat}, Ref{fmpq_mat}), z, m)
      finalizer(_fmpq_mat_clear_fn, z)
      return z
   end
end

function _fmpq_mat_clear_fn(a::fmpq_mat)
   ccall((:fmpq_mat_clear, :libflint), Nothing, (Ref{fmpq_mat},), a)
end

###############################################################################
#
#   FmpzMatSpace / fmpz_mat
#
###############################################################################

# not really a mathematical ring
mutable struct FmpzMatSpace <: MatSpace{fmpz}
   rows::Int
   cols::Int
   base_ring::FlintIntegerRing

   function FmpzMatSpace(r::Int, c::Int, cached::Bool = true)
      if haskey(FmpzMatID, (r, c))
         return FmpzMatID[r, c]
      else
         z = new(r, c, FlintZZ)
         if cached
            FmpzMatID[r, c] = z
         end
         return z
      end
   end
end

const FmpzMatID = Dict{Tuple{Int, Int}, FmpzMatSpace}()

mutable struct fmpz_mat <: MatElem{fmpz}
   entries::Ptr{Nothing}
   r::Int
   c::Int
   rows::Ptr{Nothing}
   base_ring::FlintIntegerRing
   view_parent

   # Used by view, not finalised!!
   function fmpz_mat()
      return new()
   end

   function fmpz_mat(r::Int, c::Int)
      z = new()
      ccall((:fmpz_mat_init, :libflint), Nothing,
            (Ref{fmpz_mat}, Int, Int), z, r, c)
      finalizer(_fmpz_mat_clear_fn, z)
      return z
   end

   function fmpz_mat(r::Int, c::Int, arr::Array{fmpz, 2})
      z = new()
      ccall((:fmpz_mat_init, :libflint), Nothing,
            (Ref{fmpz_mat}, Int, Int), z, r, c)
      finalizer(_fmpz_mat_clear_fn, z)
      GC.@preserve z for i = 1:r
         for j = 1:c
            el = ccall((:fmpz_mat_entry, :libflint), Ptr{fmpz},
                       (Ref{fmpz_mat}, Int, Int), z, i - 1, j - 1)
            ccall((:fmpz_set, :libflint), Nothing,
                  (Ptr{fmpz}, Ref{fmpz}), el, arr[i, j])
         end
      end
      return z
   end

   function fmpz_mat(r::Int, c::Int, arr::Array{fmpz, 1})
      z = new()
      ccall((:fmpz_mat_init, :libflint), Nothing,
            (Ref{fmpz_mat}, Int, Int), z, r, c)
      finalizer(_fmpz_mat_clear_fn, z)
      GC.@preserve z for i = 1:r
         for j = 1:c
            el = ccall((:fmpz_mat_entry, :libflint), Ptr{fmpz},
                       (Ref{fmpz_mat}, Int, Int), z, i - 1, j - 1)
            ccall((:fmpz_set, :libflint), Nothing,
                  (Ptr{fmpz}, Ref{fmpz}), el, arr[(i-1)*c+j])
         end
      end
      return z
   end

   function fmpz_mat(r::Int, c::Int, arr::Array{T, 2}) where {T <: Integer}
      z = new()
      ccall((:fmpz_mat_init, :libflint), Nothing,
            (Ref{fmpz_mat}, Int, Int), z, r, c)
      finalizer(_fmpz_mat_clear_fn, z)
      GC.@preserve z for i = 1:r
         for j = 1:c
            el = ccall((:fmpz_mat_entry, :libflint), Ptr{fmpz},
                       (Ref{fmpz_mat}, Int, Int), z, i - 1, j - 1)
            ccall((:fmpz_set, :libflint), Nothing,
                  (Ptr{fmpz}, Ref{fmpz}), el, fmpz(arr[i, j]))
         end
      end
      return z
   end

   function fmpz_mat(r::Int, c::Int, arr::Array{T,1}) where {T <: Integer}
      z = new()
      ccall((:fmpz_mat_init, :libflint), Nothing,
            (Ref{fmpz_mat}, Int, Int), z, r, c)
      finalizer(_fmpz_mat_clear_fn, z)
      GC.@preserve z for i = 1:r
         for j = 1:c
            el = ccall((:fmpz_mat_entry, :libflint), Ptr{fmpz},
                       (Ref{fmpz_mat}, Int, Int), z, i - 1, j - 1)
            ccall((:fmpz_set, :libflint), Nothing,
                  (Ptr{fmpz}, Ref{fmpz}), el, fmpz(arr[(i-1)*c+j]))
         end
      end
      return z
   end

   function fmpz_mat(r::Int, c::Int, d::fmpz)
      z = new()
      ccall((:fmpz_mat_init, :libflint), Nothing,
            (Ref{fmpz_mat}, Int, Int), z, r, c)
      finalizer(_fmpz_mat_clear_fn, z)
      GC.@preserve z for i = 1:min(r, c)
         el = ccall((:fmpz_mat_entry, :libflint), Ptr{fmpz},
                    (Ref{fmpz_mat}, Int, Int), z, i - 1, i- 1)
         ccall((:fmpz_set, :libflint), Nothing,
               (Ptr{fmpz}, Ref{fmpz}), el, d)
      end
      return z
   end

   function fmpz_mat(m::fmpz_mat)
      z = new()
      ccall((:fmpz_mat_init_set, :libflint), Nothing,
            (Ref{fmpz_mat}, Ref{fmpz_mat}), z, m)
      finalizer(_fmpz_mat_clear_fn, z)
      return z
   end
end

function _fmpz_mat_clear_fn(a::fmpz_mat)
   ccall((:fmpz_mat_clear, :libflint), Nothing, (Ref{fmpz_mat},), a)
end

###############################################################################
#
#   NmodMatSpace / nmod_mat
#
###############################################################################

mutable struct NmodMatSpace <: MatSpace{nmod}
  base_ring::NmodRing
  n::UInt
  rows::Int
  cols::Int

  function NmodMatSpace(R::NmodRing, r::Int, c::Int,
                        cached::Bool = true)
    (r < 0 || c < 0) && throw(error_dim_negative)
    if haskey(NmodMatID, (R, r, c))
      return NmodMatID[R, r, c]
    else
      z = new(R, R.n, r, c)
      if cached
        NmodMatID[R, r, c] = z
      end
      return z
    end
  end
end

const NmodMatID = Dict{Tuple{NmodRing, Int, Int}, NmodMatSpace}()

mutable struct nmod_mat <: MatElem{nmod}
  entries::Ptr{Nothing}
  r::Int                  # Int
  c::Int                  # Int
  rows::Ptr{Nothing}
  n::UInt                # mp_limb_t / Culong
  ninv::UInt             # mp_limb_t / Culong
  norm::UInt             # mp_limb_t / Culong
  base_ring::NmodRing
  view_parent

  # Used by view, not finalised!!
  function nmod_mat()
    z = new()
    return z
  end

  function nmod_mat(r::Int, c::Int, n::UInt)
    z = new()
    ccall((:nmod_mat_init, :libflint), Nothing,
            (Ref{nmod_mat}, Int, Int, UInt), z, r, c, n)
    finalizer(_nmod_mat_clear_fn, z)
    return z
  end

  function nmod_mat(r::Int, c::Int, n::UInt, arr::Array{UInt, 2}, transpose::Bool = false)
    z = new()
    ccall((:nmod_mat_init, :libflint), Nothing,
            (Ref{nmod_mat}, Int, Int, UInt), z, r, c, n)
    finalizer(_nmod_mat_clear_fn, z)
    if transpose
      se = set_entry_t!
      r,c = c,r
    else
      se = set_entry!
    end
    for i = 1:r
      for j = 1:c
        se(z, i, j, arr[i, j])
      end
    end
    return z
  end

  function nmod_mat(r::Int, c::Int, n::UInt, arr::Array{UInt, 1}, transpose::Bool = false)
    z = new()
    ccall((:nmod_mat_init, :libflint), Nothing,
            (Ref{nmod_mat}, Int, Int, UInt), z, r, c, n)
    finalizer(_nmod_mat_clear_fn, z)
    if transpose
      se = set_entry_t!
      r,c = c,r
    else
      se = set_entry!
    end
    for i = 1:r
      for j = 1:c
        se(z, i, j, arr[(i - 1) * c + j])
      end
    end
    return z
  end

  function nmod_mat(r::Int, c::Int, n::UInt, arr::Array{fmpz, 2}, transpose::Bool = false)
    z = new()
    ccall((:nmod_mat_init, :libflint), Nothing,
            (Ref{nmod_mat}, Int, Int, UInt), z, r, c, n)
    finalizer(_nmod_mat_clear_fn, z)
    if transpose
      se = set_entry_t!
      r,c = c,r
    else
      se = set_entry!
    end
    for i = 1:r
      for j = 1:c
        se(z, i, j, arr[i, j])
      end
    end
    return z
  end

  function nmod_mat(r::Int, c::Int, n::UInt, arr::Array{fmpz, 1}, transpose::Bool = false)
    z = new()
    ccall((:nmod_mat_init, :libflint), Nothing,
            (Ref{nmod_mat}, Int, Int, UInt), z, r, c, n)
    finalizer(_nmod_mat_clear_fn, z)
    if transpose
      se = set_entry_t!
      r,c = c,r
    else
      se = set_entry!
    end
    for i = 1:r
      for j = 1:c
        se(z, i, j, arr[(i - 1) * c + j])
      end
    end
    return z
  end

  function nmod_mat(r::Int, c::Int, n::UInt, arr::Array{T, 2}, transpose::Bool = false) where {T <: Integer}
    arr_fmpz = map(fmpz, arr)
    return nmod_mat(r, c, n, arr_fmpz, transpose)
  end

  function nmod_mat(r::Int, c::Int, n::UInt, arr::Array{T, 1}, transpose::Bool = false) where {T <: Integer}
    arr_fmpz = map(fmpz, arr)
    return nmod_mat(r, c, n, arr_fmpz, transpose)
  end

  function nmod_mat(r::Int, c::Int, n::UInt, arr::Array{nmod, 2}, transpose::Bool = false)
    z = new()
    ccall((:nmod_mat_init, :libflint), Nothing,
            (Ref{nmod_mat}, Int, Int, UInt), z, r, c, n)
    finalizer(_nmod_mat_clear_fn, z)
    if transpose
      se = set_entry_t!
      r,c = c,r
    else
      se = set_entry!
    end
    for i = 1:r
      for j = 1:c
        se(z, i, j, arr[i, j])
      end
    end
    return z
  end

  function nmod_mat(r::Int, c::Int, n::UInt, arr::Array{nmod, 1}, transpose::Bool = false)
    z = new()
    ccall((:nmod_mat_init, :libflint), Nothing,
            (Ref{nmod_mat}, Int, Int, UInt), z, r, c, n)
    finalizer(_nmod_mat_clear_fn, z)
    if transpose
      se = set_entry_t!
      r, c = c, r
    else
      se = set_entry!
    end
    for i = 1:r
      for j = 1:c
        se(z, i, j, arr[(i - 1) * c + j])
      end
    end
    return z
  end

  function nmod_mat(n::UInt, b::fmpz_mat)
    z = new()
    ccall((:nmod_mat_init, :libflint), Nothing,
            (Ref{nmod_mat}, Int, Int, UInt), z, b.r, b.c, n)
    finalizer(_nmod_mat_clear_fn, z)
    ccall((:fmpz_mat_get_nmod_mat, :libflint), Nothing,
            (Ref{nmod_mat}, Ref{fmpz_mat}), z, b)
    return z
  end

  function nmod_mat(n::Int, b::fmpz_mat)
    (n < 0) && error("Modulus must be positive")
    return nmod_mat(UInt(n), b)
  end

  function nmod_mat(n::fmpz, b::fmpz_mat)
    (n < 0) && error("Modulus must be positive")
    (n > typemax(UInt)) &&
          error("Modulus must be smaller than ", fmpz(typemax(UInt)))
    return nmod_mat(UInt(n), b)
  end
end

function _nmod_mat_clear_fn(mat::nmod_mat)
  ccall((:nmod_mat_clear, :libflint), Nothing, (Ref{nmod_mat}, ), mat)
end

################################################################################
#
#   GFPMatSpace / gfp_mat
#
###############################################################################

mutable struct GFPMatSpace <: MatSpace{gfp_elem}
  base_ring::GaloisField
  n::UInt
  rows::Int
  cols::Int

  function GFPMatSpace(R::GaloisField, r::Int, c::Int,
                        cached::Bool = true)
    (r < 0 || c < 0) && throw(error_dim_negative)
    if haskey(GFPMatID, (R, r, c))
      return GFPMatID[R, r, c]
    else
      z = new(R, R.n, r, c)
      if cached
        GFPMatID[R, r, c] = z
      end
      return z
    end
  end
end

const GFPMatID = Dict{Tuple{GaloisField, Int, Int}, GFPMatSpace}()

mutable struct gfp_mat <: MatElem{gfp_elem}
  entries::Ptr{Nothing}
  r::Int                  # Int
  c::Int                  # Int
  rows::Ptr{Nothing}
  n::UInt                # mp_limb_t / Culong
  ninv::UInt             # mp_limb_t / Culong
  norm::UInt             # mp_limb_t / Culong
  base_ring::GaloisField
  view_parent

  # Used by view, not finalised!!
  function gfp_mat()
    z = new()
    return z
  end

  function gfp_mat(r::Int, c::Int, n::UInt)
    z = new()
    ccall((:nmod_mat_init, :libflint), Nothing,
            (Ref{gfp_mat}, Int, Int, UInt), z, r, c, n)
    finalizer(_gfp_mat_clear_fn, z)
    return z
  end

  function gfp_mat(r::Int, c::Int, n::UInt, arr::Array{UInt, 2}, transpose::Bool = false)
    z = new()
    ccall((:nmod_mat_init, :libflint), Nothing,
            (Ref{gfp_mat}, Int, Int, UInt), z, r, c, n)
    finalizer(_gfp_mat_clear_fn, z)
    if transpose
      se = set_entry_t!
      r, c = c, r
    else
      se = set_entry!
    end
    for i = 1:r
      for j = 1:c
        se(z, i, j, arr[i, j])
      end
    end
    return z
  end

  function gfp_mat(r::Int, c::Int, n::UInt, arr::Array{UInt, 1}, transpose::Bool = false)
    z = new()
    ccall((:nmod_mat_init, :libflint), Nothing,
            (Ref{gfp_mat}, Int, Int, UInt), z, r, c, n)
    finalizer(_gfp_mat_clear_fn, z)
    if transpose
      se = set_entry_t!
      r, c = c, r
    else
      se = set_entry!
    end
    for i = 1:r
      for j = 1:c
        se(z, i, j, arr[(i - 1) * c + j])
      end
    end
    return z
  end

  function gfp_mat(r::Int, c::Int, n::UInt, arr::Array{fmpz, 2}, transpose::Bool = false)
    z = new()
    ccall((:nmod_mat_init, :libflint), Nothing,
            (Ref{gfp_mat}, Int, Int, UInt), z, r, c, n)
    finalizer(_gfp_mat_clear_fn, z)
    if transpose
      se = set_entry_t!
      r, c = c, r
    else
      se = set_entry!
    end
    for i = 1:r
      for j = 1:c
        se(z, i, j, arr[i, j])
      end
    end
    return z
  end

  function gfp_mat(r::Int, c::Int, n::UInt, arr::Array{fmpz, 1}, transpose::Bool = false)
    z = new()
    ccall((:nmod_mat_init, :libflint), Nothing,
            (Ref{gfp_mat}, Int, Int, UInt), z, r, c, n)
    finalizer(_gfp_mat_clear_fn, z)
    if transpose
      se = set_entry_t!
      r, c = c, r
    else
      se = set_entry!
    end
    for i = 1:r
      for j = 1:c
        se(z, i, j, arr[(i - 1) * c + j])
      end
    end
    return z
  end

  function gfp_mat(r::Int, c::Int, n::UInt, arr::Array{T, 2}, transpose::Bool = false) where {T <: Integer}
    arr_fmpz = map(fmpz, arr)
    return gfp_mat(r, c, n, arr_fmpz, transpose)
  end

  function gfp_mat(r::Int, c::Int, n::UInt, arr::Array{T, 1}, transpose::Bool = false) where {T <: Integer}
    arr_fmpz = map(fmpz, arr)
    return gfp_mat(r, c, n, arr_fmpz, transpose)
  end

  function gfp_mat(r::Int, c::Int, n::UInt, arr::Array{gfp_elem, 2}, transpose::Bool = false)
    z = new()
    ccall((:nmod_mat_init, :libflint), Nothing,
            (Ref{gfp_mat}, Int, Int, UInt), z, r, c, n)
    finalizer(_gfp_mat_clear_fn, z)
    if transpose
      se = set_entry_t!
      r, c = c, r
    else
      se = set_entry!
    end
    for i = 1:r
      for j = 1:c
        se(z, i, j, arr[i, j])
      end
    end
    return z
  end

  function gfp_mat(r::Int, c::Int, n::UInt, arr::Array{gfp_elem, 1}, transpose::Bool = false)
    z = new()
    ccall((:nmod_mat_init, :libflint), Nothing,
            (Ref{gfp_mat}, Int, Int, UInt), z, r, c, n)
    finalizer(_gfp_mat_clear_fn, z)
    if transpose
      se = set_entry_t!
      r, c = c, r
    else
      se = set_entry!
    end
    for i = 1:r
      for j = 1:c
        se(z, i, j, arr[(i - 1) * c + j])
      end
    end
    return z
  end

  function gfp_mat(n::UInt, b::fmpz_mat)
    z = new()
    ccall((:nmod_mat_init, :libflint), Nothing,
            (Ref{gfp_mat}, Int, Int, UInt), z, b.r, b.c, n)
    finalizer(_gfp_mat_clear_fn, z)
    ccall((:fmpz_mat_get_nmod_mat, :libflint), Nothing,
            (Ref{gfp_mat}, Ref{fmpz_mat}), z, b)
    return z
  end

  function gfp_mat(n::Int, b::fmpz_mat)
    (n < 0) && error("Modulus must be positive")
    return gfp_mat(UInt(n), b)
  end

  function gfp_mat(n::fmpz, b::fmpz_mat)
    (n < 0) && error("Modulus must be positive")
    (n > typemax(UInt)) &&
          error("Modulus must be smaller than ", fmpz(typemax(UInt)))
    return gfp_mat(UInt(n), b)
  end
end

function _gfp_mat_clear_fn(mat::gfp_mat)
  ccall((:nmod_mat_clear, :libflint), Nothing, (Ref{gfp_mat}, ), mat)
end

const Zmodn_mat = Union{nmod_mat, gfp_mat}

###############################################################################
#
#   FqPolyRing / fq_poly
#
###############################################################################

mutable struct FqPolyRing <: PolyRing{fq}
   base_ring::FqFiniteField
   S::Symbol

   function FqPolyRing(R::FqFiniteField, s::Symbol, cached::Bool = true)
      if haskey(FqPolyID, (R, s))
         return FqPolyID[(R, s)]
      else
         z = new(R,s)
         if cached
            FqPolyID[(R, s)] = z
         end
         return z
      end
   end
end

const FqPolyID = Dict{Tuple{FqFiniteField, Symbol}, FqPolyRing}()

mutable struct fq_poly <: PolyElem{fq}
   coeffs::Ptr{Nothing}
   alloc::Int
   length::Int
   parent::FqPolyRing

   function fq_poly()
      z = new()
      ccall((:fq_poly_init, :libflint), Nothing, (Ref{fq_poly},), z)
      finalizer(_fq_poly_clear_fn, z)
      return z
   end

   function fq_poly(a::fq_poly)
      z = new()
      ctx = base_ring(parent(a))
      ccall((:fq_poly_init, :libflint), Nothing,
            (Ref{fq_poly}, Ref{FqFiniteField}), z, ctx)
      ccall((:fq_poly_set, :libflint), Nothing,
            (Ref{fq_poly}, Ref{fq_poly}, Ref{FqFiniteField}),
            z, a, ctx)
      finalizer(_fq_poly_clear_fn, z)
      return z
   end

   function fq_poly(a::fq)
      z = new()
      ctx = parent(a)
      ccall((:fq_poly_init, :libflint), Nothing,
            (Ref{fq_poly}, Ref{FqFiniteField}), z, ctx)
      ccall((:fq_poly_set_fq, :libflint), Nothing,
            (Ref{fq_poly}, Ref{fq}, Ref{FqFiniteField}),
            z, a, ctx)
      finalizer(_fq_poly_clear_fn, z)
      return z
   end

   function fq_poly(a::Array{fq, 1})
      z = new()
      ctx = parent(a[1])
      ccall((:fq_poly_init2, :libflint), Nothing,
            (Ref{fq_poly}, Int, Ref{FqFiniteField}),
            z, length(a), ctx)
      for i = 1:length(a)
         ccall((:fq_poly_set_coeff, :libflint), Nothing,
               (Ref{fq_poly}, Int, Ref{fq}, Ref{FqFiniteField}),
               z, i - 1, a[i], ctx)
      end
      finalizer(_fq_poly_clear_fn, z)
      return z
   end

   function fq_poly(a::Array{fmpz, 1}, ctx::FqFiniteField)
      z = new()
      temp = ctx()
      ccall((:fq_poly_init2, :libflint), Nothing,
            (Ref{fq_poly}, Int, Ref{FqFiniteField}),
            z, length(a), ctx)
      for i = 1:length(a)
         temp = ctx(a[i])
         ccall((:fq_poly_set_coeff, :libflint), Nothing,
               (Ref{fq_poly}, Int, Ref{fq}, Ref{FqFiniteField}),
               z, i - 1, temp, ctx)
      end
      finalizer(_fq_poly_clear_fn, z)
      return z
   end

   function fq_poly(a::fmpz_poly, ctx::FqFiniteField)
      z = new()
      ccall((:fq_poly_init2, :libflint), Nothing,
            (Ref{fq_poly}, Int, Ref{FqFiniteField}),
            z, length(a), ctx)
      for i = 1:length(a)
         temp = ctx(coeff(a, i-1))
         ccall((:fq_poly_set_coeff, :libflint), Nothing,
               (Ref{fq_poly}, Int, Ref{fq}, Ref{FqFiniteField}),
               z, i - 1, temp, ctx)
      end
      finalizer(_fq_poly_clear_fn, z)
      return z
   end
end

function _fq_poly_clear_fn(a::fq_poly)
   ccall((:fq_poly_clear, :libflint), Nothing, (Ref{fq_poly},), a)
end

mutable struct fq_poly_factor
  poly::Ptr{fq_poly}
  exp::Ptr{Int}
  num::Int
  alloc::Int
  base_field::FqFiniteField

  function fq_poly_factor(ctx::FqFiniteField)
    z = new()
    ccall((:fq_poly_factor_init, :libflint), Nothing,
         (Ref{fq_poly_factor}, Ref{FqFiniteField}), z, ctx)
    z.base_field = ctx
    finalizer(_fq_poly_factor_clear_fn, z)
    return z
  end
end

function _fq_poly_factor_clear_fn(a::fq_poly_factor)
   ccall((:fq_poly_factor_clear, :libflint), Nothing,
         (Ref{fq_poly_factor}, Ref{FqFiniteField}),
         a, a.base_field)
end

###############################################################################
#
#   FqNmodPolyRing / fq_nmod_poly
#
###############################################################################

mutable struct FqNmodPolyRing <: PolyRing{fq_nmod}
   base_ring::FqNmodFiniteField
   S::Symbol

   function FqNmodPolyRing(R::FqNmodFiniteField, s::Symbol, cached::Bool = true)
      if haskey(FqNmodPolyID, (R, s))
         return FqNmodPolyID[(R, s)]
      else
         z = new(R,s)
         if cached
            FqNmodPolyID[(R, s)] = z
         end
         return z
      end
   end
end

const FqNmodPolyID = Dict{Tuple{FqNmodFiniteField, Symbol}, FqNmodPolyRing}()

mutable struct fq_nmod_poly <: PolyElem{fq_nmod}
   coeffs::Ptr{Nothing}
   alloc::Int
   length::Int
   parent::FqNmodPolyRing

   function fq_nmod_poly()
      z = new()
      ccall((:fq_nmod_poly_init, :libflint), Nothing, (Ref{fq_nmod_poly},), z)
      finalizer(_fq_nmod_poly_clear_fn, z)
      return z
   end

   function fq_nmod_poly(a::fq_nmod_poly)
      z = new()
      ctx = base_ring(parent(a))
      ccall((:fq_nmod_poly_init, :libflint), Nothing,
            (Ref{fq_nmod_poly}, Ref{FqNmodFiniteField}), z, ctx)
      ccall((:fq_nmod_poly_set, :libflint), Nothing,
            (Ref{fq_nmod_poly}, Ref{fq_nmod_poly}, Ref{FqNmodFiniteField}),
            z, a, ctx)
      finalizer(_fq_nmod_poly_clear_fn, z)
      return z
   end

   function fq_nmod_poly(a::fq_nmod)
      z = new()
      ctx = parent(a)
      ccall((:fq_nmod_poly_init, :libflint), Nothing,
            (Ref{fq_nmod_poly}, Ref{FqNmodFiniteField}), z, ctx)
      ccall((:fq_nmod_poly_set_fq_nmod, :libflint), Nothing,
            (Ref{fq_nmod_poly}, Ref{fq_nmod}, Ref{FqNmodFiniteField}),
            z, a, ctx)
      finalizer(_fq_nmod_poly_clear_fn, z)
      return z
   end

   function fq_nmod_poly(a::Array{fq_nmod, 1})
      z = new()
      ctx = parent(a[1])
      ccall((:fq_nmod_poly_init2, :libflint), Nothing,
            (Ref{fq_nmod_poly}, Int, Ref{FqNmodFiniteField}),
            z, length(a), ctx)
      for i = 1:length(a)
         ccall((:fq_nmod_poly_set_coeff, :libflint), Nothing,
               (Ref{fq_nmod_poly}, Int, Ref{fq_nmod}, Ref{FqNmodFiniteField}),
               z, i - 1, a[i], ctx)
      end
      finalizer(_fq_nmod_poly_clear_fn, z)
      return z
   end

   function fq_nmod_poly(a::Array{fmpz, 1}, ctx::FqNmodFiniteField)
      z = new()
      temp = ctx()
      ccall((:fq_nmod_poly_init2, :libflint), Nothing,
            (Ref{fq_nmod_poly}, Int, Ref{FqNmodFiniteField}),
            z, length(a), ctx)
      for i = 1:length(a)
         temp = ctx(a[i])
         ccall((:fq_nmod_poly_set_coeff, :libflint), Nothing,
               (Ref{fq_nmod_poly}, Int, Ref{fq_nmod}, Ref{FqNmodFiniteField}),
               z, i - 1, temp, ctx)
      end
      finalizer(_fq_nmod_poly_clear_fn, z)
      return z
   end

   function fq_nmod_poly(a::fmpz_poly, ctx::FqNmodFiniteField)
      z = new()
      ccall((:fq_nmod_poly_init2, :libflint), Nothing,
            (Ref{fq_nmod_poly}, Int, Ref{FqNmodFiniteField}),
            z, length(a), ctx)
      for i = 1:length(a)
         temp = ctx(coeff(a,i-1))
         ccall((:fq_nmod_poly_set_coeff, :libflint), Nothing,
               (Ref{fq_nmod_poly}, Int, Ref{fq_nmod}, Ref{FqNmodFiniteField}),
               z, i - 1, temp, ctx)
      end
      finalizer(_fq_nmod_poly_clear_fn, z)
      return z
   end
end

function _fq_nmod_poly_clear_fn(a::fq_nmod_poly)
   ccall((:fq_nmod_poly_clear, :libflint), Nothing, (Ref{fq_nmod_poly},), a)
end

mutable struct fq_nmod_poly_factor
  poly::Ptr{fq_nmod_poly}
  exp::Ptr{Int}
  num::Int
  alloc::Int
  base_field::FqNmodFiniteField

  function fq_nmod_poly_factor(ctx::FqNmodFiniteField)
    z = new()
    ccall((:fq_nmod_poly_factor_init, :libflint), Nothing,
         (Ref{fq_nmod_poly_factor}, Ref{FqNmodFiniteField}), z, ctx)
    z.base_field = ctx
    finalizer(_fq_nmod_poly_factor_clear_fn, z)
    return z
  end
end

function _fq_nmod_poly_factor_clear_fn(a::fq_nmod_poly_factor)
   ccall((:fq_nmod_poly_factor_clear, :libflint), Nothing,
         (Ref{fq_nmod_poly_factor}, Ref{FqNmodFiniteField}),
         a, a.base_field)
end

###############################################################################
#
#   FqMatSpace/fq_mat
#
###############################################################################

mutable struct FqMatSpace <: MatSpace{fq}
  base_ring::FqFiniteField
  rows::Int
  cols::Int

  function FqMatSpace(R::FqFiniteField, r::Int, c::Int, cached::Bool = true)
    (r < 0 || c < 0) && throw(error_dim_negative)
    if haskey(FqMatID, (R, r, c))
      return FqMatID[R, r, c]
    else
      z = new(R, r, c)
      if cached
        FqMatID[R, r, c] = z
      end
      return z
    end
  end
end

const FqMatID = Dict{Tuple{FqFiniteField, Int, Int}, FqMatSpace}()

mutable struct fq_mat <: MatElem{fq}
   entries::Ptr{Nothing}
   r::Int
   c::Int
   rows::Ptr{Nothing}
   base_ring::FqFiniteField
   view_parent
   
   # used by windows, not finalised!!
   function fq_mat()
      return new()
   end

   function fq_mat(r::Int, c::Int, ctx::FqFiniteField)
      z = new()
      ccall((:fq_mat_init, :libflint), Nothing,
            (Ref{fq_mat}, Int, Int, Ref{FqFiniteField}), z, r, c, ctx)
      z.base_ring = ctx
      finalizer(_fq_mat_clear_fn, z)
      return z
   end

   function fq_mat(r::Int, c::Int, arr::Array{fq, 2}, ctx::FqFiniteField)
      z = new()
      ccall((:fq_mat_init, :libflint), Nothing,
            (Ref{fq_mat}, Int, Int, Ref{FqFiniteField}), z, r, c, ctx)
      GC.@preserve z for i = 1:r
         for j = 1:c
            ccall((:fq_mat_entry_set, :libflint), Nothing,
                  (Ref{fq_mat}, Int, Int, Ref{fq}, Ref{FqFiniteField}),
                   z, i - 1, j - 1, arr[i, j], ctx)
         end
      end
      z.base_ring = ctx
      finalizer(_fq_mat_clear_fn, z)
      return z
   end

   function fq_mat(r::Int, c::Int, arr::Array{fq, 1}, ctx::FqFiniteField)
      z = new()
      ccall((:fq_mat_init, :libflint), Nothing,
            (Ref{fq_mat}, Int, Int, Ref{FqFiniteField}), z, r, c, ctx)
      GC.@preserve z for i = 1:r
         for j = 1:c
            ccall((:fq_mat_entry_set, :libflint), Nothing,
                       (Ref{fq_mat}, Int, Int, Ref{fq}, Ref{FqFiniteField}),
                        z, i - 1, j - 1, arr[(i - 1) * c + j], ctx)
         end
      end
      z.base_ring = ctx
      finalizer(_fq_mat_clear_fn, z)
      return z
   end

   function fq_mat(r::Int, c::Int, arr::Array{fmpz, 2}, ctx::FqFiniteField)
      z = new()
      ccall((:fq_mat_init, :libflint), Nothing,
            (Ref{fq_mat}, Int, Int, Ref{FqFiniteField}), z, r, c, ctx)
      GC.@preserve z for i = 1:r
         for j = 1:c
            el = ccall((:fq_mat_entry, :libflint), Ptr{fq},
                       (Ref{fq_mat}, Int, Int), z, i - 1, j - 1)
            ccall((:fq_set_fmpz, :libflint), Nothing,
                  (Ptr{fq}, Ref{fmpz}, Ref{FqFiniteField}), el, arr[i, j], ctx)
         end
      end
      z.base_ring = ctx
      finalizer(_fq_mat_clear_fn, z)
      return z
   end

   function fq_mat(r::Int, c::Int, arr::Array{fmpz, 1}, ctx::FqFiniteField)
      z = new()
      ccall((:fq_mat_init, :libflint), Nothing,
            (Ref{fq_mat}, Int, Int, Ref{FqFiniteField}), z, r, c, ctx)
      GC.@preserve z for i = 1:r
         for j = 1:c
            el = ccall((:fq_mat_entry, :libflint), Ptr{fq},
                       (Ref{fq_mat}, Int, Int), z, i - 1, j - 1)
            ccall((:fq_set_fmpz, :libflint), Nothing,
                  (Ptr{fq}, Ref{fmpz}, Ref{FqFiniteField}), el, arr[(i - 1) * c + j], ctx)
         end
      end
      z.base_ring = ctx
      finalizer(_fq_mat_clear_fn, z)
      return z
   end

   function fq_mat(r::Int, c::Int, arr::Array{T, 2}, ctx::FqFiniteField) where {T <: Integer}
      z = new()
      ccall((:fq_mat_init, :libflint), Nothing,
            (Ref{fq_mat}, Int, Int, Ref{FqFiniteField}), z, r, c, ctx)
      GC.@preserve z for i = 1:r
         for j = 1:c
            ccall((:fq_mat_entry_set, :libflint), Nothing,
                    (Ref{fq_mat}, Int, Int, Ref{fq}, Ref{FqFiniteField}),
                     z, i - 1, j - 1, ctx(arr[i, j]), ctx)
         end
      end
      z.base_ring = ctx
      finalizer(_fq_mat_clear_fn, z)
      return z
   end

   function fq_mat(r::Int, c::Int, arr::Array{T ,1}, ctx::FqFiniteField) where {T <: Integer}
      z = new()
      ccall((:fq_mat_init, :libflint), Nothing,
            (Ref{fq_mat}, Int, Int, Ref{FqFiniteField}), z, r, c, ctx)
      GC.@preserve z for i = 1:r
         for j = 1:c
            ccall((:fq_mat_entry_set, :libflint), Nothing,
                       (Ref{fq_mat}, Int, Int, Ref{fq}, Ref{FqFiniteField}),
                        z, i - 1, j - 1, ctx(arr[(i - 1) * c + j]), ctx)
         end
      end
      z.base_ring = ctx
      finalizer(_fq_mat_clear_fn, z)
      return z
   end

   function fq_mat(r::Int, c::Int, d::fq)
      z = new()
      ctx = parent(d)
      ccall((:fq_mat_init, :libflint), Nothing,
            (Ref{fq_mat}, Int, Int, Ref{FqFiniteField}), z, r, c, ctx)
      for i = 1:min(r, c)
         ccall((:fq_mat_entry_set, :libflint), Nothing,
               (Ref{fq_mat}, Int, Int, Ref{fq}, Ref{FqFiniteField}), z, i - 1, i- 1, d, ctx)
      end
      z.base_ring = ctx
      finalizer(_fq_mat_clear_fn, z)
      return z
   end

   function fq_mat(m::fmpz_mat, ctx::FqFiniteField)
      z = new()
      r = rows(m)
      c = cols(m)
      ccall((:fq_mat_init, :libflint), Nothing,
            (Ref{fq_mat}, Int, Int, Ref{FqFiniteField}), z, r, c, ctx)
      GC.@preserve z for i = 1:r
         for j = 1:c
            el1 = ccall((:fq_mat_entry, :libflint), Ptr{fq},
                        (Ref{fq_mat}, Int, Int), z, i - 1, j - 1)
            el2 = ccall((:fmpz_mat_entry, :libflint), Ptr{fmpz},
                        (Ref{fmpz_mat}, Int, Int), m, i - 1, j - 1)

            ccall((:fq_set_fmpz, :libflint), Nothing,
                  (Ptr{fq}, Ptr{fmpz}, Ref{FqFiniteField}), el1, el2, ctx)
         end
      end
      z.base_ring = ctx
      finalizer(_fq_mat_clear_fn, z)
      return z
   end
end

function _fq_mat_clear_fn(a::fq_mat)
   ccall((:fq_mat_clear, :libflint), Nothing, (Ref{fq_mat}, Ref{FqFiniteField}), a, base_ring(a))
end

###############################################################################
#
#   FqNmodMatSpace/fq_nmod_mat
#
###############################################################################

mutable struct FqNmodMatSpace <: MatSpace{fq_nmod}
  base_ring::FqNmodFiniteField
  rows::Int
  cols::Int

  function FqNmodMatSpace(R::FqNmodFiniteField, r::Int, c::Int, cached::Bool = true)
    (r < 0 || c < 0) && throw(error_dim_negative)
    if haskey(FqNmodMatID, (R, r, c))
      return FqNmodMatID[R, r, c]
    else
      z = new(R, r, c)
      if cached
        FqNmodMatID[R, r, c] = z
      end
      return z
    end
  end
end

const FqNmodMatID = Dict{Tuple{FqNmodFiniteField, Int, Int}, FqNmodMatSpace}()

mutable struct fq_nmod_mat <: MatElem{fq_nmod}
   entries::Ptr{Nothing}
   r::Int
   c::Int
   rows::Ptr{Nothing}
   base_ring::FqNmodFiniteField
   view_parent
   
   # used by windows, not finalised!!
   function fq_nmod_mat()
      return new()
   end

   function fq_nmod_mat(r::Int, c::Int, ctx::FqNmodFiniteField)
      z = new()
      ccall((:fq_nmod_mat_init, :libflint), Nothing,
            (Ref{fq_nmod_mat}, Int, Int, Ref{FqNmodFiniteField}), z, r, c, ctx)
      z.base_ring = ctx
      finalizer(_fq_nmod_mat_clear_fn, z)
      return z
   end

   function fq_nmod_mat(r::Int, c::Int, arr::Array{fq_nmod, 2}, ctx::FqNmodFiniteField)
      z = new()
      ccall((:fq_nmod_mat_init, :libflint), Nothing,
            (Ref{fq_nmod_mat}, Int, Int, Ref{FqNmodFiniteField}), z, r, c, ctx)
      GC.@preserve z for i = 1:r
         for j = 1:c
            ccall((:fq_nmod_mat_entry_set, :libflint), Nothing,
                  (Ref{fq_nmod_mat}, Int, Int, Ref{fq_nmod}, Ref{FqNmodFiniteField}),
                   z, i - 1, j - 1, arr[i, j], ctx)
         end
      end
      z.base_ring = ctx
      finalizer(_fq_nmod_mat_clear_fn, z)
      return z
   end

   function fq_nmod_mat(r::Int, c::Int, arr::Array{fq_nmod, 1}, ctx::FqNmodFiniteField)
      z = new()
      ccall((:fq_nmod_mat_init, :libflint), Nothing,
            (Ref{fq_nmod_mat}, Int, Int, Ref{FqNmodFiniteField}), z, r, c, ctx)
      GC.@preserve z for i = 1:r
         for j = 1:c
            ccall((:fq_nmod_mat_entry_set, :libflint), Nothing,
                       (Ref{fq_nmod_mat}, Int, Int, Ref{fq_nmod}, Ref{FqNmodFiniteField}),
                        z, i - 1, j - 1, arr[(i - 1) * c + j], ctx)
         end
      end
      z.base_ring = ctx
      finalizer(_fq_nmod_mat_clear_fn, z)
      return z
   end

   function fq_nmod_mat(r::Int, c::Int, arr::Array{fmpz, 2}, ctx::FqNmodFiniteField)
      z = new()
      ccall((:fq_nmod_mat_init, :libflint), Nothing,
            (Ref{fq_nmod_mat}, Int, Int, Ref{FqNmodFiniteField}), z, r, c, ctx)
      GC.@preserve z for i = 1:r
         for j = 1:c
            el = ccall((:fq_nmod_mat_entry, :libflint), Ptr{fq_nmod},
                       (Ref{fq_nmod_mat}, Int, Int), z, i - 1, j - 1)
            ccall((:fq_nmod_set_fmpz, :libflint), Nothing,
                  (Ptr{fq_nmod}, Ref{fmpz}, Ref{FqNmodFiniteField}), el, arr[i, j], ctx)
         end
      end
      z.base_ring = ctx
      finalizer(_fq_nmod_mat_clear_fn, z)
      return z
   end

   function fq_nmod_mat(r::Int, c::Int, arr::Array{fmpz, 1}, ctx::FqNmodFiniteField)
      z = new()
      ccall((:fq_nmod_mat_init, :libflint), Nothing,
            (Ref{fq_nmod_mat}, Int, Int, Ref{FqNmodFiniteField}), z, r, c, ctx)
      GC.@preserve z for i = 1:r
         for j = 1:c
            el = ccall((:fq_nmod_mat_entry, :libflint), Ptr{fq_nmod},
                       (Ref{fq_nmod_mat}, Int, Int), z, i - 1, j - 1)
            ccall((:fq_nmod_set_fmpz, :libflint), Nothing,
                  (Ptr{fq_nmod}, Ref{fmpz}, Ref{FqNmodFiniteField}), el, arr[(i - 1) * c + j], ctx)
         end
      end
      z.base_ring = ctx
      finalizer(_fq_nmod_mat_clear_fn, z)
      return z
   end

   function fq_nmod_mat(r::Int, c::Int, arr::Array{T, 2}, ctx::FqNmodFiniteField) where {T <: Integer}
      z = new()
      ccall((:fq_nmod_mat_init, :libflint), Nothing,
            (Ref{fq_nmod_mat}, Int, Int, Ref{FqNmodFiniteField}), z, r, c, ctx)
      GC.@preserve z for i = 1:r
         for j = 1:c
            ccall((:fq_nmod_mat_entry_set, :libflint), Nothing,
                    (Ref{fq_nmod_mat}, Int, Int, Ref{fq_nmod}, Ref{FqNmodFiniteField}),
                     z, i - 1, j - 1, ctx(arr[i, j]), ctx)
         end
      end
      z.base_ring = ctx
      finalizer(_fq_nmod_mat_clear_fn, z)
      return z
   end

   function fq_nmod_mat(r::Int, c::Int, arr::Array{T ,1}, ctx::FqNmodFiniteField) where {T <: Integer}
      z = new()
      ccall((:fq_nmod_mat_init, :libflint), Nothing,
            (Ref{fq_nmod_mat}, Int, Int, Ref{FqNmodFiniteField}), z, r, c, ctx)
      GC.@preserve z for i = 1:r
         for j = 1:c
            ccall((:fq_nmod_mat_entry_set, :libflint), Nothing,
                       (Ref{fq_nmod_mat}, Int, Int, Ref{fq_nmod}, Ref{FqNmodFiniteField}),
                        z, i - 1, j - 1, ctx(arr[(i - 1) * c + j]), ctx)
         end
      end
      z.base_ring = ctx
      finalizer(_fq_nmod_mat_clear_fn, z)
      return z
   end

   function fq_nmod_mat(r::Int, c::Int, d::fq_nmod)
      z = new()
      ctx = parent(d)
      ccall((:fq_nmod_mat_init, :libflint), Nothing,
            (Ref{fq_nmod_mat}, Int, Int, Ref{FqNmodFiniteField}), z, r, c, ctx)
      for i = 1:min(r, c)
         ccall((:fq_nmod_mat_entry_set, :libflint), Nothing,
               (Ref{fq_nmod_mat}, Int, Int, Ref{fq_nmod}, Ref{FqNmodFiniteField}), z, i - 1, i- 1, d, ctx)
      end
      z.base_ring = ctx
      finalizer(_fq_nmod_mat_clear_fn, z)
      return z
   end

   function fq_nmod_mat(m::fmpz_mat, ctx::FqNmodFiniteField)
      z = new()
      r = rows(m)
      c = cols(m)
      ccall((:fq_nmod_mat_init, :libflint), Nothing,
            (Ref{fq_nmod_mat}, Int, Int, Ref{FqNmodFiniteField}), z, r, c, ctx)
      GC.@preserve z for i = 1:r
         for j = 1:c
            el1 = ccall((:fq_nmod_mat_entry, :libflint), Ptr{fq_nmod},
                        (Ref{fq_nmod_mat}, Int, Int), z, i - 1, j - 1)
            el2 = ccall((:fmpz_mat_entry, :libflint), Ptr{fmpz},
                        (Ref{fmpz_mat}, Int, Int), m, i - 1, j - 1)

            ccall((:fq_nmod_set_fmpz, :libflint), Nothing,
                  (Ptr{fq_nmod}, Ptr{fmpz}, Ref{FqNmodFiniteField}), el1, el2, ctx)
         end
      end
      z.base_ring = ctx
      finalizer(_fq_nmod_mat_clear_fn, z)
      return z
   end
end

function _fq_nmod_mat_clear_fn(a::fq_nmod_mat)
   ccall((:fq_nmod_mat_clear, :libflint), Nothing, (Ref{fq_nmod_mat}, Ref{FqNmodFiniteField}), a, base_ring(a))
end

