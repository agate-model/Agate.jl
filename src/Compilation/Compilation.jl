"""Setup-time process contribution realization and static tendency compilation."""
module Compilation

using ..Configuration: ComponentLayout, CommunityContext
using ..Equations: CompiledEquation
using ..Processes:
    AbstractFormulation,
    Growth,
    Mortality,
    Monod,
    NamedProcess,
    NormalizedModelDefinition,
    NutrientResponse,
    PartitionRouting,
    Smith,
    formulation,
    parameter_name,
    parameter_requirements,
    process_id,
    process_rate

export AbstractProcessContribution
export GrowthParameterBinding, GrowthTopology
export GrowthBiomassContribution, GrowthResourceLossContribution
export MortalityParameterBinding, MortalityTopology
export MortalityLossContribution, MortalityRoutingContribution
export realize_process_topology, process_contributions
export contribution_target, group_contributions
export compile_tendency, compile_tendencies

include("contributions.jl")
include("mortality.jl")
include("growth.jl")

end # module Compilation
