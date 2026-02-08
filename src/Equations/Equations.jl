module Equations

export CompiledEquation

"""A wrapper around a callable equation.

`CompiledEquation` is used throughout Agate to store GPU-callable tracer tendency
functions in a type-stable way.

Parameter validation is handled by the factory's parameter directory
(`parameter_definitions(factory)`).
"""
struct CompiledEquation{F}
    f::F
end

@inline (eq::CompiledEquation)(args...) = eq.f(args...)

end # module
