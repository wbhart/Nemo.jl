###############################################################################
#
#   CalciumTypes.jl : Parent and object types for Calcium
#
###############################################################################

export CalciumQQBarField, qqbar

################################################################################
#
#  Structs for shallow operations
#
################################################################################

mutable struct qqbar_struct
  coeffs::Ptr{Nothing}
  alloc::Int
  length::Int
  real_mid_exp::Int     # fmpz
  real_mid_size::UInt # mp_size_t
  real_mid_d1::UInt    # mantissa_struct
  real_mid_d2::UInt
  real_rad_exp::Int     # fmpz
  real_rad_man::UInt
  imag_mid_exp::Int     # fmpz
  imag_mid_size::UInt # mp_size_t
  imag_mid_d1::UInt    # mantissa_struct
  imag_mid_d2::UInt
  imag_rad_exp::Int     # fmpz
  imag_rad_man::UInt
end

mutable struct fexpr_struct
  data::Ptr{Nothing}
  alloc::Int
end

################################################################################
#
#  Types and memory management for QQBarField
#
################################################################################

mutable struct CalciumQQBarField <: Field
end

const CalciumQQBar = CalciumQQBarField()

mutable struct qqbar <: FieldElem
  coeffs::Ptr{Nothing}
  alloc::Int
  length::Int
  real_mid_exp::Int     # fmpz
  real_mid_size::UInt # mp_size_t
  real_mid_d1::UInt    # mantissa_struct
  real_mid_d2::UInt
  real_rad_exp::Int     # fmpz
  real_rad_man::UInt
  imag_mid_exp::Int     # fmpz
  imag_mid_size::UInt # mp_size_t
  imag_mid_d1::UInt    # mantissa_struct
  imag_mid_d2::UInt
  imag_rad_exp::Int     # fmpz
  imag_rad_man::UInt

  function qqbar()
    z = new()
    ccall((:qqbar_init, libcalcium), Nothing, (Ref{qqbar}, ), z)
    finalizer(_qqbar_clear_fn, z)
    return z
  end

end

function _qqbar_clear_fn(a::qqbar)
   ccall((:qqbar_clear, libcalcium), Nothing, (Ref{qqbar},), a)
end


