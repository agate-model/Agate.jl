# Kernel-time bundles for tidy, GPU-friendly code.
#
# Oceananigans calls biogeochemistry tracers with a positional `args...` argument
# list. Within kernel-callable helper functions we frequently need to thread
# through the same trio:
#
# - model parameters (`bgc.parameters`)
# - a `Tracers` accessor (name → integer index)
# - the positional argument tuple (`args`)
#
# `KernelBundle` keeps these values together in a concrete struct so helper
# functions can accept a single argument without vararg plumbing.

"""A small bundle of kernel-time arguments.

`KernelBundle` is intentionally minimal: it stores only concrete, GPU-safe
objects (typically model parameters, a `Tracers` accessor, and a tuple of
tracer/aux values).

All fields are fully inferred, which helps GPU compilation and reduces
kernel-callsite argument noise.
"""
struct KernelBundle{P,TR,ARGS}
    p::P
    tracers::TR
    args::ARGS
end

@inline KernelBundle(bgc, args) = KernelBundle(bgc.parameters, bgc.tracers, args)
