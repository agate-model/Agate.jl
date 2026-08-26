module Configuration

using Adapt
using ..Parameters: parameter_directory

export Population, Pool, PopulationStateRef, population_state
export currency, states, state_currency, size_structure
export realize_components
export component_classes, component_state_tracers
export component_tracers, component_indices, state_tracers, state_tracer, state_indices
export component_diameters, component_class_count
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
