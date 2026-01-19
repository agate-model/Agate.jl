# Callable dynamics API

Agate expresses tracer tendencies as **callable Julia functions** that are safe to run on
both CPU and GPU architectures.

Each model factory provides a mapping from tracer names to `CompiledEquation` values.
A `CompiledEquation` stores:

- a callable `f` (the tracer tendency), and
- a `Requirements` object describing which model parameters the callable reads.

The constructor uses these requirements to validate that the parameter set is complete before
building the final Oceananigans biogeochemistry object.

## `CompiledEquation` and `Requirements`

The types live in `Agate.Functors`:

- `req(; scalars=..., vectors=..., matrices=...)`
- `CompiledEquation(f, r)`

Requirements are just lists of parameter keys (`Symbol`s). They are intentionally explicit:
there is no expression parsing during construction.

Example:

```julia
using Agate.Functors: CompiledEquation, req

# f must follow the Oceananigans biogeochemistry kernel signature.
#
# Note: `args...` contains tracer values followed by auxiliary fields, in the
# order defined by the constructed model.
@inline function my_tracer(bgc, x, y, z, t, args...)
    k = bgc.parameters.detritus_remineralization
    # ... unpack from args as needed ...
    return -k
end

eq = CompiledEquation(my_tracer, req(scalars=(:detritus_remineralization,)))
```

## Tracer function signature

Callable tracer functions must accept the Oceananigans biogeochemistry signature:

```julia
f(bgc, x, y, z, t, args...)
```

Agate passes the model instance as `bgc`, which stores the resolved parameter NamedTuple at
`bgc.parameters`.

## Building blocks

Agate's reusable building blocks live in `Agate.Library`.
They are small callable structs designed to compose inside tracer functions:

- `Agate.Library.Nutrients`
- `Agate.Library.Photosynthesis`
- `Agate.Library.Predation`
- `Agate.Library.Mortality`
- `Agate.Library.Remineralization`
