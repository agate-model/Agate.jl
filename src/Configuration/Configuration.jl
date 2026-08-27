module Configuration

using Adapt

export Population, Pool, PopulationStateRef, population_state
export currency, states, state_currency, size_structure
export component_classes, component_state_tracers
export component_tracers, state_tracers, state_tracer
export component_diameters
export PalatabilityAllometric, AssimilationBinary

if !hasmethod(Adapt.adapt_structure, Tuple{Any,NamedTuple})
    @inline function Adapt.adapt_structure(to, nt::NamedTuple{names}) where {names}
        return NamedTuple{names}(map(x -> Adapt.adapt(to, x), values(nt)))
    end
end

include("community.jl")
include("components.jl")
include("interactions_matrices.jl")
include("interactions_derivations.jl")

end # module
