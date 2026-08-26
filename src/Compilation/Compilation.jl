"""Setup-time process-flux realization and static tendency compilation."""
module Compilation

using ..Configuration:
    ModelLayout, PopulationStateRef, component_classes, component_class_indices, state_tracer
using ..Processes:
    AbstractFactor,
    Growth,
    Nutrients,
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

include("fluxes.jl")
include("factors.jl")
include("products.jl")
include("mortality.jl")
include("growth.jl")
include("consumption.jl")
include("remineralization.jl")

end # module Compilation
