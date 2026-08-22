"""Setup-time process-flux realization and static tendency compilation."""
module Compilation

using ..Configuration: ComponentLayout, CommunityContext, axis_indices
using ..Equations: CompiledEquation
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

export compile_model_tendencies

include("fluxes.jl")
include("factors.jl")
include("routing.jl")
include("mortality.jl")
include("growth.jl")
include("consumption.jl")
include("remineralization.jl")

end # module Compilation
