"""Concrete context for tracer tendency helper functions.

Oceananigans supplies biogeochemistry tracers as a positional `args...` tuple.
Tendency helper functions commonly need three associated values together:

- model parameters (`bgc.parameters`)
- a `Tracers` accessor that maps tracer names to integer indices
- the positional runtime argument tuple (`args`)

`TendencyContext` bundles these values into a concrete, GPU-safe struct so
helper functions can accept a single argument while preserving full type
inference.

Fields
------
- `parameters`: model parameters
- `tracers`: tracer accessor
- `args`: positional tracer and auxiliary values passed at runtime
"""
struct TendencyContext{P,TR,ARGS}
    parameters::P
    tracers::TR
    args::ARGS
end

@inline TendencyContext(bgc, args) = TendencyContext(bgc.parameters, bgc.tracers, args)

"""A value-oriented view of tracer arguments.

`TracerValues` bundles a `Tracers` accessor with the positional runtime state
tuple (`args`) passed by Oceananigans. This lets tracer tendency code read like
mathematical expressions:

```julia
vals = TracerValues(tracers, args)
DIN = vals.DIN
P   = vals.plankton(i)
```

No Symbol lookup is performed at runtime: property access uses the precomputed
indices stored in the `Tracers` accessor.
"""
struct TracerValues{TR,ARGS}
    tracers::TR
    args::ARGS
end

"""Callable view used for `vals.plankton(i)`.

Returned by `getproperty(vals, :plankton)` to avoid closures in kernel-callable
code.
"""
struct BoundPlankton{TR,ARGS}
    tracers::TR
    args::ARGS
end

@inline (bp::BoundPlankton)(i::Int) = bp.tracers.plankton(bp.args, i)

@inline function Base.getproperty(v::TracerValues, name::Symbol)
    if name === :tracers
        return getfield(v, :tracers)
    elseif name === :args
        return getfield(v, :args)
    elseif name === :plankton
        return BoundPlankton(getfield(v, :tracers), getfield(v, :args))
    else
        acc = getproperty(getfield(v, :tracers), name)
        return acc(getfield(v, :args))
    end
end

"""Return model parameters and tracer accessors for a tracer tendency function.

Call this helper at the start of a tracer tendency implementation to unpack
`bgc` and `args` into the values commonly needed by the function:

```julia
parameters, tracer_values = tendency_inputs(bgc, args)
```

The returned `tracer_values` object provides named access to the positional
tracer arguments, including `tracer_values.plankton(i)` for plankton tracers.
"""
@inline function tendency_inputs(bgc, args)
    context = TendencyContext(bgc, args)
    tracer_values = TracerValues(context.tracers, context.args)
    return context.parameters, tracer_values
end