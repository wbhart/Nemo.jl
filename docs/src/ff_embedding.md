```@meta
CurrentModule = Nemo
```

# Finite field embeddings

## Introduction

Nemo allows the construction of finite field embeddings making use of the
algorithm of Bosma, Cannon and Steel behind the scenes to ensure compatibility.
Critical routines (e.g. polynomial factorization, matrix computations) are
provided by the C library Flint, whereas high level tasks are written directly in Nemo.

## Embedding functionality

It is possible to explicitly call the embedding `embed` function to create an embedding,
but standard object oriented conversion works too. The resulting embedding is of
type `FinFieldMorphism`. It is also possible to
compute the preimage of an embedding via the `preimage` function, applied to an
embedding or directly to the finite fields (this actually first computes the
embedding), or via conversion. An error is thrown if the element you want to
compute the preimage is not in the image of the embedding.

### Computing an embedding

```@docs
embed(::FqNmodFiniteField, ::FqNmodFiniteField)
```

**Examples**

```julia
k2, x2 = FiniteField(19, 2, "x2")
k4, x4 = FiniteField(19, 4, "x4")

f = embed(k2, k4)

y = f(x2)
z = k4(x2)
```

### Computed the preimage of an embedding

```@docs
preimage(::FqNmodFiniteField, ::FqNmodFiniteField)
preimage(::FinFieldMorphism)
```

**Examples**

```julia
k7, x7 = FiniteField(13, 7, "x7")
k21, x21 = FiniteField(13, 21, "x21")

s = preimage(k7, k21)

y = k21(x7)

z = s(y)
t = k7(y)
```
