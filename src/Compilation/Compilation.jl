"""Setup-time process-flux realization and static tendency compilation."""
module Compilation

using ..Configuration: ComponentLayout, CommunityContext
using ..Equations: CompiledEquation
using ..Processes:
    AbstractFactor,
    Growth,
    Nutrients,
    Grazing,
    Consumption,
    Mortality,
    Remineralization,
    FactorDriver,
    FactorComponent,
    factor_inputs,
    factor_children,
    factor_parameter_context,
    factor_child_path,
    HeterotrophicConsumption,
    NamedProcess,
    IdealizedGrazing,
    PreferentialGrazing,
    ProductRouting,
    DOMPOMRouting,
    NormalizedModelDefinition,
    NutrientResponse,
    DirectRouting,
    PartitionRouting,
    LinearRemineralization,
    formulation,
    factors,
    factor_value,
    ParameterBinding,
    parameter_slots,
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
