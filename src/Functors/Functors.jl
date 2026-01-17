"""Agate.Functors

Small, construction-time utilities for attaching explicit parameter requirements to
callables.

The runtime biogeochemistry model stores only callables. Requirements exist only
to validate that a parameter bundle provides the keys needed by those callables.
"""

module Functors

export Requirements, req, merge_requirements, requirements
export CompiledEquation

"""A construction-time dependency set for parameter keys."""
struct Requirements
    vectors::Vector{Symbol}
    matrices::Vector{Symbol}
    scalars::Vector{Symbol}
end

"""Construct a `Requirements` object."""
function req(; vectors=(), matrices=(), scalars=())
    return Requirements(
        Symbol[vectors...],
        Symbol[matrices...],
        Symbol[scalars...],
    )
end

@inline function _append_unique!(dest::Vector{Symbol}, src::Vector{Symbol})
    for s in src
        (s in dest) || push!(dest, s)
    end
    return dest
end

"""Merge two `Requirements` with stable ordering (first-seen wins)."""
function merge_requirements(r1::Requirements, r2::Requirements)
    return Requirements(
        _append_unique!(copy(r1.vectors), r2.vectors),
        _append_unique!(copy(r1.matrices), r2.matrices),
        _append_unique!(copy(r1.scalars), r2.scalars),
    )
end

"""Return construction-time requirements for `x`."""
requirements(::Any) = req()

"""A callable wrapper that carries explicit construction-time requirements."""
struct CompiledEquation{F}
    f::F
    req::Requirements
end

"""Return the stored requirements."""
@inline requirements(eq::CompiledEquation) = eq.req

@inline (eq::CompiledEquation)(args...; kwargs...) = eq.f(args...; kwargs...)

"""Wrap `f` with empty requirements."""
@inline CompiledEquation(f) = CompiledEquation(f, req())

end # module Functors
