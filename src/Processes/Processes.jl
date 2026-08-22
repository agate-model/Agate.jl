"""Scientific process authoring and setup-time normalization."""
module Processes

using ..Configuration: Population, Pool, ComponentLayout, currency, component_classes
using ..ModelFamilies: AbstractModelFamily, default_components, default_processes
using ..Parameters: ParameterSpec, ParameterDefinition, ParameterProvision, parameter_definitions
using ..Library.Mortality: linear_loss, quadratic_loss
using ..Library.Predation: idealized_predation_loss, preferential_predation_loss
using ..Library.Photosynthesis: geider_light_response, smith_light_limitation
using ..Library.Nutrients: DEFAULT_FRANK_SHARPNESS, FrankTNorm, liebig_minimum, monod_limitation
using ..Library.Temperature: q10_temperature_factor
using ..Library.Remineralization: linear_remineralization

export AbstractProcess, AbstractFormulation, AbstractFactor, AbstractStoichiometry
export Light, Nutrients, Temperature
export FixedStoichiometry
export Growth, NutrientResponse, Consumption, Grazing, Mortality, ProductRouting, Remineralization
export ModelDefinition, NormalizedModelDefinition, NamedProcess
export ParameterSlot, ParameterRequirementIdentity, ParameterRequirement
export ParameterBinding, ParameterApplicability
export process_id, process_kind, formulation, formulation_tag, factors
export participants, drivers, rate_axes, driver_identities, normalize_model
export parameter_slots, parameter_requirements, parameter_bindings, parameter_binding, parameter_name
export parameter_slot_bindings
export resolve_parameter_applicability

include("authoring.jl")
include("rates.jl")
include("normalization.jl")

end # module Processes
