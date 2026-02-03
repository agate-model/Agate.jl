# Tendency-time context for tidy, GPU-friendly code.
#
# Oceananigans calls biogeochemistry tracers with a positional `args...` argument
# list. Within kernel-callable tendency helpers we frequently need to thread
# through the same trio:
#
# - model parameters (`bgc.parameters`)
# - a `Tracers` accessor (name → integer index)
# - the positional argument tuple (`args`)
#
# `TendencyContext` keeps these values together in a concrete struct so helpers
# can accept a single argument without vararg plumbing.

"""Small, concrete context passed through tendency helper functions.

`TendencyContext` is intentionally minimal: it stores only concrete, GPU-safe
objects (typically model parameters, a `Tracers` accessor, and a tuple of
tracer/aux values).

All fields are fully inferred, which helps GPU compilation and reduces
kernel-callsite argument noise.
"""
struct TendencyContext{P,TR,ARGS}
    parameters::P
    tracers::TR
    args::ARGS
end

@inline TendencyContext(bgc, args) = TendencyContext(bgc.parameters, bgc.tracers, args)

# -----------------------------------------------------------------------------
# Convenience views
# -----------------------------------------------------------------------------

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
    if name === :tracers || name === :args
        return getfield(v, name)
    elseif name === :plankton
        return BoundPlankton(v.tracers, v.args)
    else
        accessor = getproperty(v.tracers, name)
        return accessor(v.args)
    end
end

"""Return `(tendency, parameters, vals)` for a tendency function.

This is the recommended prologue for tracer tendency callables:

```julia
tendency, parameters, vals = tendency_views(bgc, args)
```

`tendency` remains available for interaction helpers (predation sums, etc), while
`parameters` and `vals` make the common path terse and intention-revealing.
"""
@inline function tendency_views(bgc, args)
    tendency = TendencyContext(bgc, args)
    return tendency, tendency.parameters, TracerValues(tendency.tracers, tendency.args)
end
