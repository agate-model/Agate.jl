# Callable dynamics API

Agate expresses tracer tendencies as **callable Julia functions** that are safe to run on
both CPU and GPU architectures.

Each model factory provides a mapping from tracer names to `CompiledEquation` values.
A `CompiledEquation` stores a callable `f` (the tracer tendency).

Parameter completeness and parameter shape validation are handled by the factory's
**parameter directory** (see `parameter_definitions(factory)` and `ParameterSpec`).
During construction, Agate builds defaults, applies user overrides, computes any derived
interaction matrices, and then errors if any declared parameter is still missing.

## `CompiledEquation`

The type lives in `Agate.Equations`:

  - `CompiledEquation(f)`

Example:

```julia
using Agate.Equations: CompiledEquation

# f must follow the Oceananigans biogeochemistry kernel signature.
#
# Note: `args...` contains tracer values followed by auxiliary fields, in the
# order defined by the constructed model.
@inline function my_tracer(bgc, x, y, z, t, args...)
    k = bgc.parameters.detritus_remineralization
    # ... unpack from args as needed ...
    return -k
end

eq = CompiledEquation(my_tracer)
```

## Tracer function signature

Callable tracer functions must accept the Oceananigans biogeochemistry signature:

```julia
f(bgc, x, y, z, t, args...)
```

Agate passes the model instance as `bgc`, which stores the resolved parameter NamedTuple at
`bgc.parameters`.

**Resolved** means: defaults have been built, user overrides have been applied, derived interaction
matrices (if any) have been computed or recomputed, and role-aware interaction parameters have been
finalized into a canonical representation.

## Convenience unpacking: `tendency_inputs`

Inside tracer tendency functions, Oceananigans passes state values positionally via `args...`.
Agate provides a small helper that converts this positional tuple into explicit, named accessors:

```julia
using Agate.Runtime: tendency_inputs

@inline function my_tracer(bgc, x, y, z, t, args...)
    parameters, tracer_values = tendency_inputs(bgc, args)

    # Scalar parameter (already cast to FT and adapted to the chosen architecture):
    k = parameters.detritus_remineralization

    # Tracer values:
    N = tracer_values.N

    return -k * N
end
```

## Building blocks

Agate's reusable building blocks live in `Agate.Library`.
They are small callable structs designed to compose inside tracer functions:

  - `Agate.Library.Nutrients`
  - `Agate.Library.Photosynthesis`
  - `Agate.Library.Predation`
  - `Agate.Library.Mortality`
  - `Agate.Library.Remineralization`
