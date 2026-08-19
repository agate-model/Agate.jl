"""Setup-time process contribution realization and static tendency compilation."""
module Compilation

using ..Configuration: ComponentLayout, CommunityContext
using ..Equations: CompiledEquation
using ..Processes:
    AbstractFormulation,
    Mortality,
    NamedProcess,
    NormalizedModelDefinition,
    PartitionRouting,
    formulation,
    parameter_name,
    parameter_requirements,
    process_id,
    process_rate

export AbstractProcessContribution
export MortalityParameterBinding, MortalityTopology
export MortalityLossContribution, MortalityRoutingContribution
export realize_process_topology, process_contributions
export contribution_target, group_contributions
export compile_tendency, compile_tendencies

include("mortality.jl")

end # module Compilation
