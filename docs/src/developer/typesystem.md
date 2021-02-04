# Types and parents

## The Nemo type system

### Concrete and abstract types

Julia does not provide a traditional class/inheritance approach to programming.
Instead, the basic unit of its object oriented approach is the type definition
(`struct` and `mutable struct`) and inheritance exists only on the function
side of the language rather than data side. Julia provides a rich system of
abstract types and unions on the data side and multimethods on the function
side to effect this.

For example Julia's `Number` type is an abstract type containing all concrete
types that behave like numbers, e.g. `Int64`, `Float64`, and so on.

Abstract types can also belong to other abstract types, forming a tree of
abstract types.

In Nemo the most important abstract types are `Ring` and `Field`, with the
latter belonging to the former so that all fields are rings, and the abstract
types `RingElem` and `FieldElem` for the objects that represent elements of
rings and fields, again with the latter abstract type belonging to the former.

Because this hierarchy of abstract types must form a tree, Julia is strictly
speaking single inheritance, as each concrete and abstract type can belong to
at most one other abstract type. For example, one could not have a diamond of
abstract types with `ExactField` belonging to both `Field` and `ExactRing`.

### Recovering aspects of multiple inheritance in Nemo

Various possibilities exist to get around the limitation that abstract types
must form a 'tree' in Nemo and AbstractAlgebra.

One such possibility is union types. If a function should accept one of a
number of concrete or abstract types that can't all be made to belong to a
single abstract type due to this limitation then one can use a union type.

For example, Nemo defines `RingElement` to be a union of `RingElem` and all
the Julia standard types which behave like ring elements, e.g. all `Integer`
types and types of rationals with `Integer` components.

A second feature we make use of in Nemo is parameterised types. Each concrete
and abstract type can take one or more parameters. These parameter can be any
other type, either concrete or abstract. For example, in Julia `Rational{T}`
is for rationals with numerator and denominator of type `T`.

A great deal of control over parameterised types is possible, e.g. one can
restrict the type parameter `T` using a `where` clause, e.g. to write a
function that accepts all rational types with integer components of the same
type one can use the type `Rational{T} where T <: Integer`.

Nemo makes use of such parameterised types for generic ring constructions
such as generic polynomial rings and matrices over a given base ring. The type
of the elements of the base ring is substituted for the parameter `T` in any
concrete instantiation of the types `Poly{T}` and `Mat{T}`, which are defined
in AbstractAlgebra in `src/generic/GenericTypes.jl`.

The totality of all univariate polynomial types, including those of generic
`Poly{T}` types and those coming from C libraries (such as `fmpz_poly`), is
represented by the abstract type `PolyElem{T}` which in turn belongs to
`RingElem`, both defined in AbstractAlgebra in `src/AbstractTypes.jl`.

Similarly, the totality of all matrix types, including explicit C types
like `fmpz_mat` and the generic `Mat{T}` types is given by the abstract type
`MatElem{T}`, again defined in AbstractAlgebra in `src/AbstractTypes.jl`.

This hierarchy of types allows one to write functions at any level, e.g. for
all univariate polynomial types, just those with a given base type `T`, or
for a specific concrete type corresponding to just one kind of univariate
polynomial.

A third possibility to get around the single inheritance limitation of Julia is
type traits. There is currently no explicit compiler/language support for
traits, however various implementations exist that make use of type parameters
in tricky ways. This allows one to add 'traits' to types, so long as those
traits can be expressed as types. In this way, types can have multiple
'properties' at the same time, instead of belonging to just a single abstract
type.

Nemo does not currently use type traits, though the map types in Nemo do make
use of a custom analogue of this.

Note that unlike class based systems that dispatch on the type of the implicit
`this` parameter, Julia methods dispatch on the type of all arguments. This is
a natural fit for mathematics where all sorts of ad hoc left and right
operations may be required.

### Encapsulation, maps and runtime flags

One limitation of the Julia approach is that the type of an object cannot be
changed at runtime. For example one might like to insist that a given ring is
in fact a field. There are three standard ways to handle this in Julia.

The first approach is to encapsulate the object in another object which does
have the desired type. The second approach is to map the object to a different
one of the required type (e.g. by applying a morphism). The third approach is
to introduce data fields in the original type which can be changed at runtime,
unlike its type. All three approaches come with downsides. 

Encapsulation can be time consuming for the developer as methods which applied
to the original object do not automatically apply to the encapsulated object.
One can write methods which do, but this is not automatic.

Application of a map may come with a performance penalty and may be difficult
for the user to navigate. Moreover, mutation of the resulting object does not
result in mutation of the original object.

The third option of adding runtime data fields essentially takes one back to
writing a (possibly bug ridden) interpreter. It relies on the developer 
implementing outer methods that make use of hand written control statements
to determine which of a range of inner methods should be applied to the object.
This misses the benefits of one of the main defining features of Julia, namely
its multimethod system and can also make introspection more difficult.

Nemo does not apply any of these three approaches widely at present, though
information which can only be known at runtime such as whether a ring is
Euclidean will eventually have to encoded using one of these three methods.

