```@meta
CurrentModule = Nemo
```

# Future plans

## Ring and CommRing

Currently all commutative ring types belong to `Ring` and their elements to
`RingElem` (and `RingElement`) and we have separate types for noncommutative
rings and elements thereof, i.e. `NCRing` and `NCRingElem` etc.

However, it would be more logical to use `Ring` for not necessarily
commutative rings and CommRing`, `CRing` or `CommutativeRing` (the name has
not been decided on yet) for commutative rings.

This is a big change and should happen with plenty of warning for the
community. It would be convenient if a script could be made available to
automate this.

## Matrices and polynomials over number fields

The Antic C library is used for number field elements, but matrices and
polynomials over these are handled by Julia. There is nothing explicitly
wrong with this, except that the jit compilation is costly for such basic
operations and the code is not as easy for others to use in other projects
if it is not in C.

We plan to write some routines in C in the Antic library to handle at least
univariate polynomials over number fields and matrices as well. These should
use highly optimised routines based on the algorithms with best complexity
for these operations.

## Mono repository

There is currently a proposal to place all Oscar related repositories, or some
subset of them in a single repository called OscarMono.jl. The details are not
finalised and it is not known what impact this will have on Nemo. However,
Nemo developers should be aware that this may happen at some point in the
fairly near future.

Users of Nemo should be unaffected, as Nemo will continue to exist as a
separate package in the OscarMono.jl repository, even if it does become part of
this repository. Julia apparently supports multiple packages in the same
repository nowadays.

The possibility will always exist to separate the repositories again if the
experiment is unsuccessful or serves its purpose and is no longer needed.

## Splitting abstract and generic functionality

Currently we have many generic types in Nemo, e.g. `Poly{T}` and `Mat{T}`. We
also have a lot of generic algorithms implemented. The latter are implemented
for abstract types.

At the present moment, both the specialised code for our generic types and the
generic algorithms for abstract types are all implemented inside the
`src/generic` directory of the AbstractAlgebra package.

We eventually plan to separate the two kinds of implementations to make it
clearer to developers which code is for generic types only and which are the
more general implementations for all types belonging to the AbstractAlgebra
abstract types.

It's not clear at this moment what the new directories will be called when the
separation happens.

## Moving implementations from Hecke

In the Hecke.jl project there are a vast number of implementations that were
intended for AbstractAlgebra.jl and Nemo.jl. They exist in the `src/Misc`
directory of that project.

These implementations will eventually all be moved over to the correct
repositories. Code, documentation and performance improvements will be added.

A number of things must be taken into account when making such moves:

* Substantial chunks of code should be moved at a time. The code can be
  initially placed in a `src/Misc` directory in AbstractAlgebra or Nemo
  until it can finally be integrated fully into the correct place in those
  projects.

* Some of the code calls parent object constructors in generic code. Such calls
  should be removed where possible. If they are essential, they should have
  `AbstractAlgebra` prepended to their calls, as per the developer documentation
  on parent object constructors.

* Some functions such as `exp` and the like requires `Base` to be prepended, as
  we do not import these functions from `Base` into `Generic`.

* Some of the code calls back into convenience functions found only in Hecke.
  These have to be rewritten in terms of AbstractAlgebra/Nemo functions.

* Some of the code relies on `fmpz` being available, but would otherwise be
  suitable for AbstractAlgebra. This code can hopefully be rewritten to be
  agnostic about the integer type.

* Todos, questions and so on should be moved to tickets.

* Sometimes exception types differ between Hecke and Nemo, meaning that tests
  will fail due to the wrong type of exception being raised. Either the tests
  will have to be adjusted, or the Nemo exception types changed.

* `RingElem` is often used where `RingElement` is intended, etc. Also types are
  often unconstrained where `Nemo` would constrain them to `RingElement`.

* Some Hecke functions try to support generic types and specific concrete Nemo
  types in the same implementation. These will unfortunately have to either be
  split between AbstractAlgebra and Nemo or a completely generic implementation
  for abstract types will have to be made.

* Some Hecke implementations assume `sub!` and friends are available in generic
  code. These will have to be rewritten, usually by adding a single unary minus
  outside of a loop and switching to `add!` and friends inside the loops.

* Test code, docstrings and documentation will have to be added where they do
  not already exist.

* Imports and exports of the new functionality will have to be added.

* Some functions should be accompanied by similar functions that don't yet
  exist. For example if there is a `blah_rows` there probably should be a
  `blah_cols` function as well, etc.


