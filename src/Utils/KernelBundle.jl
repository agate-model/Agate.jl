# Kernel-time bundles for tidy, GPU-friendly code.
#
# Oceananigans calls biogeochemistry tracers with a positional `args...` argument
# list. Within kernel-callable helper functions we frequently need to thread
# through the same pair:
#
# - a `Tracers` accessor (name → integer index)
# - the positional argument tuple (`args`)
#
# Agate assumes all tracer/aux values in `args` share the same floating-point
# type `FT` (and uses `Adapt` to enforce this on GPU), so helpers can safely use
# `zero(FT)` internally without additional reduction-seed plumbing.

"""A small bundle of kernel-time arguments.

`KernelBundle` is intentionally minimal: it stores only concrete, GPU-safe
objects (typically a `Tracers` accessor and a tuple of tracer/aux values).

This type is a stepping stone toward the Cycle C "kernel bundle" plan.
"""
struct KernelBundle{TR,ARGS}
    tracers::TR
    args::ARGS
end

