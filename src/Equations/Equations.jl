module Equations

export CompiledEquation

"""Wrap a callable tracer equation in a concrete, type-stable container.

`CompiledEquation` stores the kernel-callable function used for a single tracer
tendency.
"""
struct CompiledEquation{F}
    f::F
end

@inline (eq::CompiledEquation)(args...) = eq.f(args...)

end # module
