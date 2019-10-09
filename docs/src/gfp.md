```@meta
CurrentModule = Nemo
```

# Galois fields

Nemo allows the creation of Galois fields of the form $\mathbb{Z}/p\mathbb{Z}$ for a
prime $p$. Note that these are not the same as finite fields of degree 1, as Conway
polynomials are not used.

For convenience, the following constructors are provided.

```julia
GF(n::UInt)
GF(n::Int)
```

For example, one can create the Galois field of characteristic $7$ as follows.

```julia
R = GF(7)
```

Elements of the field are then created in the usual way.

```julia
a = R(3)
```

Elements of Galois fields have type `gfp_elem`, and the type of the parent objects is
`GaloisField`.

The modulus $p$ of an element of a Galois field is stored in its parent object.

The `gfp_elem` type belong to the abstract type `FinFieldElem` and the
`GaloisField` parent object type belongs to the abstract type `FinField`.

## Galois field functionality

Galois fields in Nemo implement the residue ring interface of AbstractAlgebra.jl:

[https://nemocas.github.io/AbstractAlgebra.jl/residue_rings.html](https://nemocas.github.io/AbstractAlgebra.jl/residue_rings.html)

In addition, all the functionality for generic residue rings is available:

[https://nemocas.github.io/AbstractAlgebra.jl/residue.html](https://nemocas.github.io/AbstractAlgebra.jl/residue.html)

Below we describe the functionality that is provided in addition to this interface.

## Basic manipulation

```@docs
characteristic(::GaloisField)
```

```@docs
order(::GaloisField)
```

**Examples**

```julia
F = GF(3)

a = characteristic(F)
b = order(F)
```
