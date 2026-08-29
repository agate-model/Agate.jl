"""Component structure, diameter specifications, realized layouts, and interaction defaults."""
module Configuration

using Adapt

export Plankton, Pool
export element, states, reference_state, variable_states, state_element, size_structure
export component_entities, component_state_tracers
export component_tracers, state_tracers, state_tracer
export component_diameters
export AllometricPalatability, ConsumerAssimilation

if !hasmethod(Adapt.adapt_structure, Tuple{Any,NamedTuple})
    @inline function Adapt.adapt_structure(to, nt::NamedTuple{names}) where {names}
        return NamedTuple{names}(map(x -> Adapt.adapt(to, x), values(nt)))
    end
end

include("diameters.jl")
include("components.jl")
include("layout.jl")
include("interaction_derivations.jl")

end # module
