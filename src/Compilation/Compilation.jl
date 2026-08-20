"""Setup-time process contribution realization and static tendency compilation."""
module Compilation

using ..Configuration: ComponentLayout, CommunityContext
using ..Equations: CompiledEquation
using ..Processes:
    AbstractFormulation,
    Growth,
    Light,
    Nutrients,
    Temperature,
    Grazing,
    Consumption,
    Mortality,
    Remineralization,
    Geider,
    Liebig,
    Monod,
    Q10,
    HeterotrophicConsumption,
    NamedProcess,
    PreferentialGrazing,
    ProductRouting,
    DOMPOMRouting,
    NormalizedModelDefinition,
    NutrientResponse,
    PartitionRouting,
    LinearRemineralization,
    Smith,
    formulation,
    factors,
    factor_value,
    parameter_name,
    parameter_requirements,
    process_id,
    process_rate

export AbstractProcessContribution
export GrowthParameterBinding, GrowthTopology
export GrowthBiomassContribution, GrowthResourceLossContribution
export GrazingParameterBinding, GrazingTopology
export ConsumptionParameterBinding, ConsumptionTopology
export GrazingResourceLossContribution, GrazingConsumerGainContribution
export GrazingUnassimilatedContribution, GrazingRoutedProductContribution
export ConsumptionResourceLossContribution, ConsumptionConsumerGainContribution
export ConsumptionUnassimilatedContribution
export MortalityParameterBinding, MortalityTopology
export MortalityLossContribution, MortalityRoutingContribution, MortalityRoutedProductContribution
export RemineralizationParameterBinding, RemineralizationTopology
export RemineralizationSourceLossContribution, RemineralizationDestinationGainContribution
export realize_process_topology, process_contributions
export contribution_target, group_contributions
export model_contributions, compile_tendency, compile_tendencies, compile_model_tendencies

include("contributions.jl")
include("factors.jl")
include("routing.jl")
include("mortality.jl")
include("growth.jl")
include("grazing.jl")
include("consumption.jl")
include("remineralization.jl")

end # module Compilation
