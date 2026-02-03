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
