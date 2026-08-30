"""Scientific process authoring and setup-time validation and canonicalization."""
module Processes

using ..Components: Plankton, Pool, PlanktonStateRef, ModelLayout, element, state_element,
    states, reference_state, component_entities
using ..ModelFamilies: AbstractModelFamily, default_components, default_processes
using ..Parameters: Parameter, ConstructionParameter, DerivedDefault, parameter_definitions
using ..Library.Mortality: linear_loss
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
export Light, NutrientLimitation, Temperature
export FixedStoichiometry
export Growth, NutrientResponse, QuotaResponse, NutrientUptake
export Consumption, Mortality, Products, Remineralization
export ModelDefinition, CanonicalProcess
export ParameterSlot
export process_id, formulation, factors
export authored_parameter_bindings
export participants
export parameter_slots

include("authoring.jl")
include("parameter_schema.jl")
include("rates.jl")
include("validation.jl")
include("canonicalization.jl")
include("parameter_plan.jl")
include("parameter_validation.jl")

end # module Processes
