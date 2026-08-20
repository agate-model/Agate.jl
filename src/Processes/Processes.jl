"""Scientific process authoring and setup-time normalization."""
module Processes

using ..Configuration: Population, Pool, ComponentLayout, currency
using ..Factories:
    AbstractBGCFactory,
    ParameterDefinition,
    ParameterProvision,
    default_components,
    default_processes,
    parameter_definitions
using ..Library.Mortality: linear_loss, quadratic_loss
using ..Library.Predation: preferential_predation_loss
using ..Library.Photosynthesis: smith_light_limitation, geider_growth
using ..Library.Nutrients: monod_limitation
using ..Library.Remineralization: linear_remineralization

export AbstractProcess, AbstractFormulation, AbstractFactor, AbstractStoichiometry
export Smith, Geider, Monod, Liebig, Q10, PreferentialGrazing, HeterotrophicConsumption
export Light, Nutrients, Temperature
export LinearMortality, QuadraticMortality, LinearRemineralization
export PartitionRouting, DOMPOMRouting, FixedStoichiometry
export Growth, NutrientResponse, Grazing, Consumption, Mortality, ProductRouting, Remineralization
export ModelDefinition, NormalizedModelDefinition, NamedProcess
export ParameterRequirementIdentity, ParameterRequirement
export ParameterBinding, ParameterApplicability
export process_id, process_kind, formulation, formulation_tag, factors
export participants, drivers, rate_axes, driver_identities, normalize_model
export parameter_requirements, parameter_bindings, parameter_name
export resolve_parameter_applicability
export process_rate

include("authoring.jl")
include("rates.jl")
include("normalization.jl")

end # module Processes
