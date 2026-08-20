"""Setup-time process-flux realization and static tendency compilation."""
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
    ParameterBinding,
    parameter_slot_bindings,
    process_id,
    process_rate

export FluxSpec, RateElement, Weight
export GrowthTopology, GrazingTopology, ConsumptionTopology
export MortalityTopology, RemineralizationTopology
export parameter_operand
export realize_process_topology, process_fluxes
export flux_target, group_fluxes
export model_fluxes, compile_tendency, compile_tendencies, compile_model_tendencies

include("fluxes.jl")
include("factors.jl")
include("routing.jl")
include("mortality.jl")
include("growth.jl")
include("grazing.jl")
include("consumption.jl")
include("remineralization.jl")

end # module Compilation
