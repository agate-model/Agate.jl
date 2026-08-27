"""Setup-time process-flux realization and static tendency compilation."""
module Compilation

using ..Configuration:
    ModelLayout, component_classes, state_tracer, state_tracers
using ..Processes:
    AbstractFactor,
    Growth,
    Light,
    Consumption,
    Mortality,
    Remineralization,
    NutrientUptake,
    FactorDriver,
    FactorComponent,
    FactorPopulationState,
    factor_inputs,
    factor_children,
    HeterotrophicConsumption,
    NamedProcess,
    process_id,
    PreferentialGrazing,
    Products,
    CanonicalModelDefinition,
    formulation,
    factors,
    factor_value,
    ParameterBinding,
    ParameterPlan,
    parameter_storage_index,
    process_products,
    process_rate

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
