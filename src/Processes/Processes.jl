"""Scientific process authoring and setup-time normalization."""
module Processes

using ..Configuration: Population, Pool, ComponentLayout
using ..Factories:
    AbstractBGCFactory,
    ParameterDefinition,
    ParameterProvision,
    default_components,
    default_processes,
    parameter_definitions
using ..Library.Mortality: linear_loss, quadratic_loss
using ..Library.Predation: preferential_predation_loss
using ..Library.Photosynthesis: smith_light_limitation
using ..Library.Nutrients: monod_limitation
using ..Library.Remineralization: linear_remineralization

export AbstractProcess, AbstractFormulation, AbstractFactor
export Smith, Monod, PreferentialGrazing
export Light
export LinearMortality, QuadraticMortality, LinearRemineralization, PartitionRouting
export Growth, NutrientResponse, Grazing, Mortality, ProductRouting, Remineralization
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
