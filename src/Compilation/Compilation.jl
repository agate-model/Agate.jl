"""Setup-time process contribution realization and static tendency compilation."""
module Compilation

using ..Configuration: ComponentLayout, CommunityContext
using ..Equations: CompiledEquation
using ..Processes:
    AbstractFormulation,
    Growth,
    Light,
    Nutrients,
    Grazing,
    Mortality,
    Remineralization,
    Geider,
    Liebig,
    Monod,
    NamedProcess,
    PreferentialGrazing,
    DOMPOMRouting,
    NormalizedModelDefinition,
    NutrientResponse,
    PartitionRouting,
    LinearRemineralization,
    Smith,
    formulation,
    parameter_name,
    parameter_requirements,
    process_id,
    process_rate

export AbstractProcessContribution
export GrowthParameterBinding, GrowthTopology
export GrowthBiomassContribution, GrowthResourceLossContribution
export GrazingParameterBinding, GrazingTopology
export GrazingResourceLossContribution, GrazingConsumerGainContribution
export GrazingUnassimilatedContribution, GrazingRoutedProductContribution
export MortalityParameterBinding, MortalityTopology
export MortalityLossContribution, MortalityRoutingContribution, MortalityRoutedProductContribution
export RemineralizationParameterBinding, RemineralizationTopology
export RemineralizationSourceLossContribution, RemineralizationDestinationGainContribution
export realize_process_topology, process_contributions
export contribution_target, group_contributions
export model_contributions, compile_tendency, compile_tendencies, compile_model_tendencies

include("contributions.jl")
include("routing.jl")
include("mortality.jl")
include("growth.jl")
include("grazing.jl")
include("remineralization.jl")

end # module Compilation
