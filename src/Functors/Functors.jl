"""Agate.Functors

Small utilities for packaging callables with their construction-time requirements.

The primary use is to attach a `Requirements` object to a callable so the
constructor can validate that all parameters referenced by the callable are
present.
"""

module Functors

using ..Equations: Requirements, req
import ..Equations: requirements

export CompiledEquation

"""
    CompiledEquation(f, req)

A callable wrapper that carries `Agate.Equations.Requirements`.

`f` is invoked directly when the wrapper is called.
"""
struct CompiledEquation{F}
    f::F
    req::Requirements
end

"""Return the stored construction-time requirements."""
@inline requirements(eq::CompiledEquation) = eq.req

@inline (eq::CompiledEquation)(args...; kwargs...) = eq.f(args...; kwargs...)

"""
    CompiledEquation(f)

Wrap `f` with empty requirements.
"""
@inline CompiledEquation(f) = CompiledEquation(f, req())

end # module
