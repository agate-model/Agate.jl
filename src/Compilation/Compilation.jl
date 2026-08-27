"""Setup-time process-flux realization and static tendency compilation."""
module Compilation

using ..Configuration:
    ModelLayout, component_classes, state_tracer, state_tracers
using ..Processes:
    AbstractFactor,
    Growth,
    Nutrients,
    Consumption,
    Mortality,
    Remineralization,
    NutrientUptake,
    FactorDriver,
    FactorComponent,
    FactorPopulationState,
    factor_inputs,
    factor_children,
    factor_parameter_context,
    factor_child_path,
    HeterotrophicConsumption,
    NamedProcess,
    SingleResourceRouting,
    QuotaRouting,
    MultiResourceRouting,
    IdealizedGrazing,
    PreferentialGrazing,
    Products,
    NormalizedModelDefinition,
    NutrientResponse,
    formulation,
    factors,
    factor_value,
    ParameterBinding,
    ParameterPlan,
    planned_parameter_slot,
    parameter_slot_bindings,
    product_path,
    process_rate

export compile_model_tendencies

include("runtime_ir.jl")
include("fluxes.jl")
include("factors.jl")
include("products.jl")
include("mortality.jl")
include("growth.jl")
include("nutrient_uptake.jl")
include("consumption.jl")
include("remineralization.jl")

end # module Compilation
