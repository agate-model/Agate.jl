"""Scientific process authoring and setup-time normalization."""
module Processes

using ..Configuration: Population, Pool, PopulationStateRef, ModelLayout, currency, state_currency,
    states, size_structure, component_classes
using ..ModelFamilies: AbstractModelFamily, default_components, default_processes
using ..Parameters: Parameter, MetaParameter, DerivedDefault, parameter_definitions
using ..Library.Mortality: linear_loss, quadratic_loss
using ..Library.Predation: preferential_predation_loss
using ..Library.Photosynthesis: geider_light_response, smith_light_limitation
using ..Library.Nutrients:
    frank_tnorm, liebig_minimum, monod_limitation, normalized_droop_limitation,
    quota_uptake_regulation
using ..Library.Temperature: q10_temperature_factor
using ..Library.Remineralization: linear_remineralization

export AbstractProcess, AbstractFormulation, AbstractFactor, AbstractStoichiometry
export Smith, Geider, Monod, NormalizedDroop, QuotaRegulatedMonod, Liebig, FrankTNorm, Q10
export PreferentialGrazing, HeterotrophicConsumption
export LinearMortality, QuadraticMortality, LinearRemineralization
export Light, Nutrients, Temperature
export FixedStoichiometry
export Growth, NutrientResponse, QuotaResponse, NutrientUptake
export Consumption, Mortality, Products, Remineralization
export ModelDefinition, NamedProcess
export ParameterSlot
export process_id, formulation, factors
export authored_parameter_bindings
export participants
export parameter_slots

include("authoring.jl")
include("parameter_schema.jl")
include("rates.jl")
include("normalization.jl")
include("parameter_plan.jl")

end # module Processes
