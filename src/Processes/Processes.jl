"""Scientific process authoring and setup-time normalization."""
module Processes

using ..Configuration: Population, Pool, ComponentLayout, currency, component_classes
using ..ModelFamilies: AbstractModelFamily, default_components, default_processes
using ..Parameters: Parameter, DerivedDefault, parameter_definitions
using ..Library.Mortality: linear_loss, quadratic_loss
using ..Library.Predation: idealized_predation_loss, preferential_predation_loss
using ..Library.Photosynthesis: geider_light_response, smith_light_limitation
using ..Library.Nutrients: frank_tnorm, liebig_minimum, monod_limitation
using ..Library.Temperature: q10_temperature_factor
using ..Library.Remineralization: linear_remineralization

export AbstractProcess, AbstractFormulation, AbstractFactor, AbstractStoichiometry
export Smith, Geider, Monod, Liebig, FrankTNorm, Q10
export IdealizedGrazing, PreferentialGrazing, HeterotrophicConsumption
export LinearMortality, QuadraticMortality, LinearRemineralization
export Light, Nutrients, Temperature
export FixedStoichiometry
export Growth, NutrientResponse, Consumption, Mortality, Products, Remineralization
export ModelDefinition, NormalizedModelDefinition, NamedProcess
export ParameterSlot
export ParameterBinding, ParameterApplicability
export process_id, formulation, factors
export authored_parameter_bindings
export participants, drivers, driver_identities, normalize_model
export parameter_slots, parameter_bindings
export parameter_slot_bindings
export resolve_parameter_applicability

include("authoring.jl")
include("rates.jl")
include("normalization.jl")

end # module Processes
