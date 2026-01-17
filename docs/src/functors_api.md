# Callable dynamics API

Agate supports *callable tracer functions* as an alternative to construction-time symbolic
[`Agate.Equations.Equation`](@ref) objects.

## Tracer function values

When calling [`Agate.Constructor.define_tracer_functions`](@ref), the `tracers` `NamedTuple` may map tracer names to:

- `Expr` (legacy)
- `Agate.Equations.Equation` (symbolic equation)
- a callable (already compiled)
- `Agate.Functors.CompiledEquation` (callable + explicit requirements)

Callable tracer functions must accept the standard Oceananigans biogeochemistry signature:

```julia
f(bgc, x, y, z, t, tracers..., auxiliary_fields...)
```

## `CompiledEquation`

Wrap a callable with [`Agate.Functors.CompiledEquation`](@ref) to attach a
[`Agate.Equations.Requirements`](@ref) object. This lets the constructor validate that the
required parameters exist without parsing expressions.

## Building blocks

Agate's building blocks are callable functors in `Agate.Library`:

- `Agate.Library.Nutrients` (e.g. `MonodLimitation`)
- `Agate.Library.Photosynthesis` (e.g. `SmithLightLimitation`, `GeiderLightLimitation`)
- `Agate.Library.Predation` (e.g. `HollingTypeII`)
- `Agate.Library.Mortality` (e.g. `LinearLoss`, `QuadraticLoss`)
- `Agate.Library.Remineralization` (e.g. `LinearRemineralization`)
