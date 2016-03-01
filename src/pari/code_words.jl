###############################################################################
#
#   code_words.jl : functions for low level access to Pari types 
#
###############################################################################

export pari, debug, gclone, gunclone, avma, gen0, gen1, evaltyp, evalsigne,
       evalvarn, setsigne, settyp, varn, typ, setvarn, signe, lg, pari_load,
       pari_print

###############################################################################
#
#   Pari GEN constants
#
###############################################################################

const t_INT      =  1
const t_REAL     =  2
const t_INTMOD   =  3
const t_FRAC     =  4
const t_FFELT    =  5
const t_COMPLEX  =  6
const t_PADIC    =  7
const t_QUAD     =  8
const t_POLMOD   =  9
const t_POL      = 10
const t_SER      = 11
const t_RFRAC    = 13
const t_QFR      = 15
const t_QFI      = 16
const t_VEC      = 17
const t_COL      = 18
const t_MAT      = 19
const t_LIST     = 20
const t_STR      = 21
const t_VECSMALL = 22
const t_CLOSURE  = 23
const t_ERROR    = 24
const t_INFINITY = 25

const BITS_IN_WORD = sizeof(Int)*8

const TYPnumBITS  = 7
const SIGNnumBITS = 2
if BITS_IN_WORD == 64
   const VARNnumBITS = 16
else
   const VARNnumBITS = 14
end
const LGnumBITS = (BITS_IN_WORD - 1 - TYPnumBITS)
const VALPnumBITS = (BITS_IN_WORD - SIGNnumBITS - VARNnumBITS)

const TYPSHIFT = (BITS_IN_WORD - TYPnumBITS)
const SIGNSHIFT = (BITS_IN_WORD - SIGNnumBITS)
const VARNSHIFT = VALPnumBITS

const SIGNBITS = ~((1 << SIGNSHIFT) - 1)
const TYPBITS  = ~((1 << TYPSHIFT) - 1)
const LGBITS = (1 << LGnumBITS) - 1

const CLONEBIT = 1<<LGnumBITS

###############################################################################
#
#   Cloning
#
###############################################################################

gclone(gen::Ptr{Int}) = ccall((:gclone, :libpari), Ptr{Int}, 
                              (Ptr{Int},), gen)

gunclone(gen::Ptr{Int}) = ccall((:gunclone, :libpari), Ptr{Int}, 
                                (Ptr{Int},), gen)

###############################################################################
#
#   Basic manipulation
#
###############################################################################

evaltyp(x::Int) = x << TYPSHIFT

evalsigne(x::Int) = x << SIGNSHIFT

evalvarn(x::Int) = x << VARNSHIFT

function setsigne(x::Ptr{Int}, s::Int)
   unsafe_store!(x, ((unsafe_load(x, 2) & (~SIGNBITS)) | evalsigne(s)), 2)
end

function settyp(x::Ptr{Int}, s::Int)
   unsafe_store!(x, ((unsafe_load(x, 1) & (~TYPBITS)) | evaltyp(s)), 1)
end

function varn(x::Ptr{Int})
   return reinterpret(Int, (reinterpret(UInt, unsafe_load(x, 2))
                                                 & VARNBITS) >> VARNSHIFT)::Int
end

function typ(x::Ptr{Int})
   return reinterpret(Int, reinterpret(UInt, unsafe_load(x, 1))
                                                              >> TYPSHIFT)::Int
end

function setvarn(x::Ptr{Int}, s::Int)
   return unsafe_store(x, reinterpret(Int, reinterpret(UInt, unsafe_load(x, 2))
                                            & ~VARNBITS) | evalvarn(s), 2)::Int
end

function signe(x::Ptr{Int})
   return (reinterpret(Int, unsafe_load(x)) >> SIGNSHIFT)::Int
end

function lg(x::Ptr{Int})
   return reinterpret(Int, (reinterpret(UInt, unsafe_load(x)) & LGBITS))::Int
end

function pari_load(x::Ptr{Int}, n::Int)
   return reinterpret(Ptr{Int}, unsafe_load(x, n))::Ptr{Int}
end

###############################################################################
#
#   Debugging
#
###############################################################################

debug(a::Ptr{Int}) = ccall((:dbgGEN, :libpari), Void, (Ptr{Int}, Int), a, -1)

###############################################################################
#
#   Printing
#
###############################################################################

function pari_print(a::Ptr{Int})
   cstr = ccall((:GENtostr, :libpari), Ptr{UInt8}, (Ptr{Int},), a)
   print(bytestring(cstr))
   ccall((:pari_free, :libpari), Void, (Ptr{UInt8},), cstr)
end
   
function pari_print(io::IO, a::Ptr{Int})
   cstr = ccall((:GENtostr, :libpari), Ptr{UInt8}, (Ptr{Int},), a)
   print(io, bytestring(cstr))
   ccall((:pari_free, :libpari), Void, (Ptr{UInt8},), cstr)
end
