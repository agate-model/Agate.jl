"""Abstract supertype for scientific process objects."""
abstract type AbstractProcess end

"""Abstract supertype for concrete scientific formulations."""
abstract type AbstractFormulation end

struct Smith <: AbstractFormulation end
struct Monod <: AbstractFormulation end
struct PreferentialGrazing <: AbstractFormulation end
struct LinearMortality <: AbstractFormulation end
struct QuadraticMortality <: AbstractFormulation end
struct LinearRemineralization <: AbstractFormulation end
struct PartitionRouting <: AbstractFormulation end

formulation_tag(::Smith) = :smith
formulation_tag(::Monod) = :monod
formulation_tag(::PreferentialGrazing) = :preferential
formulation_tag(::LinearMortality) = :linear
formulation_tag(::QuadraticMortality) = :quadratic
formulation_tag(::LinearRemineralization) = :linear
formulation_tag(::PartitionRouting) = :partition

_resolve_growth(::Val{:smith}) = Smith()
_resolve_response(::Val{:monod}) = Monod()
_resolve_grazing(::Val{:preferential}) = PreferentialGrazing()
_resolve_mortality(::Val{:linear}) = LinearMortality()
_resolve_mortality(::Val{:quadratic}) = QuadraticMortality()
_resolve_remineralization(::Val{:linear}) = LinearRemineralization()
_resolve_routing(::Val{:partition}) = PartitionRouting()

function _unknown_formulation(kind::Symbol, formulation::Symbol)
    throw(ArgumentError("unknown $(kind) formulation :$(formulation)"))
end

_resolve_growth(::Val{F}) where {F} = _unknown_formulation(:growth, F)
_resolve_response(::Val{F}) where {F} = _unknown_formulation(:nutrient_response, F)
_resolve_grazing(::Val{F}) where {F} = _unknown_formulation(:grazing, F)
_resolve_mortality(::Val{F}) where {F} = _unknown_formulation(:mortality, F)
_resolve_remineralization(::Val{F}) where {F} = _unknown_formulation(:remineralization, F)
_resolve_routing(::Val{F}) where {F} = _unknown_formulation(:product_routing, F)

function _participant_tuple(role::Symbol, singular, plural)
    isnothing(singular) == isnothing(plural) && throw(
        ArgumentError("specify exactly one of `$(role)` or `$(role)s`"),
    )
    values = isnothing(singular) ? plural : (singular,)
    values isa Symbol && (values = (values,))
    values isa Tuple || throw(ArgumentError("participant `$(role)s` must be a Symbol or tuple"))
    isempty(values) && throw(ArgumentError("participant `$(role)s` cannot be empty"))
    all(value -> value isa Symbol, values) || throw(
        ArgumentError("participant `$(role)s` must contain only Symbols"),
    )
    return values
end

function _optional_reference(name::Symbol, value)
    isnothing(value) && return nothing
    value isa Symbol || throw(ArgumentError("`$name` must be a Symbol"))
    return value
end

"""Single-resource response used by a process such as growth."""
struct NutrientResponse{F<:AbstractFormulation}
    formulation::F
    resource::Symbol
end

NutrientResponse(formulation::Symbol; resource::Symbol) =
    NutrientResponse(_resolve_response(Val(formulation)); resource)
NutrientResponse(formulation::Monod; resource::Symbol) = NutrientResponse(formulation, resource)
function NutrientResponse(formulation::AbstractFormulation; resource::Symbol)
    throw(ArgumentError("$(typeof(formulation)) is not a nutrient-response formulation"))
end

"""Product-routing sub-formulation for coupled process effects."""
struct ProductRouting{F<:AbstractFormulation}
    formulation::F
    retained::Symbol
    exported::Symbol
end

ProductRouting(formulation::Symbol; retained::Symbol, exported::Symbol) =
    ProductRouting(_resolve_routing(Val(formulation)); retained, exported)
ProductRouting(formulation::PartitionRouting; retained::Symbol, exported::Symbol) =
    ProductRouting(formulation, retained, exported)
function ProductRouting(
    formulation::AbstractFormulation; retained::Symbol, exported::Symbol
)
    throw(ArgumentError("$(typeof(formulation)) is not a product-routing formulation"))
end

"""Population growth process with optional composable resource response."""
struct Growth{F<:AbstractFormulation,L} <: AbstractProcess
    formulation::F
    populations::Tuple
    light::Symbol
    limitation::L
end

function Growth(
    formulation::Smith;
    population=nothing,
    populations=nothing,
    light=nothing,
    limitation=nothing,
)
    population_refs = _participant_tuple(:population, population, populations)
    light_ref = _optional_reference(:light, light)
    isnothing(light_ref) && throw(ArgumentError("Smith growth requires a `light` driver binding"))
    isnothing(limitation) || limitation isa NutrientResponse || throw(
        ArgumentError("growth `limitation` must be a NutrientResponse"),
    )
    return Growth(formulation, population_refs, light_ref, limitation)
end

Growth(formulation::Symbol; kwargs...) = Growth(_resolve_growth(Val(formulation)); kwargs...)
function Growth(formulation::AbstractFormulation; kwargs...)
    throw(ArgumentError("$(typeof(formulation)) is not a growth formulation"))
end

"""Consumer-resource grazing process."""
struct Grazing{F<:AbstractFormulation} <: AbstractProcess
    formulation::F
    consumers::Tuple
    resources::Tuple
    unassimilated_destination::Union{Nothing,Symbol}
end

function Grazing(
    formulation::PreferentialGrazing;
    consumer=nothing,
    consumers=nothing,
    resource=nothing,
    resources=nothing,
    unassimilated_destination=nothing,
)
    consumer_refs = _participant_tuple(:consumer, consumer, consumers)
    resource_refs = _participant_tuple(:resource, resource, resources)
    destination = _optional_reference(:unassimilated_destination, unassimilated_destination)
    return Grazing(formulation, consumer_refs, resource_refs, destination)
end

Grazing(formulation::Symbol; kwargs...) = Grazing(_resolve_grazing(Val(formulation)); kwargs...)
function Grazing(formulation::AbstractFormulation; kwargs...)
    throw(ArgumentError("$(typeof(formulation)) is not a grazing formulation"))
end

"""Population mortality process with optional product routing."""
struct Mortality{F<:AbstractFormulation,R} <: AbstractProcess
    formulation::F
    populations::Tuple
    routing::R
end

function Mortality(
    formulation::Union{LinearMortality,QuadraticMortality};
    population=nothing,
    populations=nothing,
    routing=nothing,
)
    population_refs = _participant_tuple(:population, population, populations)
    isnothing(routing) || routing isa ProductRouting || throw(
        ArgumentError("mortality `routing` must be a ProductRouting"),
    )
    return Mortality(formulation, population_refs, routing)
end

Mortality(formulation::Symbol; kwargs...) =
    Mortality(_resolve_mortality(Val(formulation)); kwargs...)
function Mortality(formulation::AbstractFormulation; kwargs...)
    throw(ArgumentError("$(typeof(formulation)) is not a mortality formulation"))
end

"""Source-to-destination remineralization process."""
struct Remineralization{F<:AbstractFormulation} <: AbstractProcess
    formulation::F
    sources::Tuple
    destinations::Tuple
end

function Remineralization(
    formulation::LinearRemineralization;
    source=nothing,
    sources=nothing,
    destination=nothing,
    destinations=nothing,
)
    source_refs = _participant_tuple(:source, source, sources)
    destination_refs = _participant_tuple(:destination, destination, destinations)
    return Remineralization(formulation, source_refs, destination_refs)
end

Remineralization(formulation::Symbol; kwargs...) =
    Remineralization(_resolve_remineralization(Val(formulation)); kwargs...)
function Remineralization(formulation::AbstractFormulation; kwargs...)
    throw(ArgumentError("$(typeof(formulation)) is not a remineralization formulation"))
end

formulation(process::AbstractProcess) = process.formulation
formulation(response::NutrientResponse) = response.formulation
formulation(routing::ProductRouting) = routing.formulation

process_kind(::Growth) = :growth
process_kind(::Grazing) = :grazing
process_kind(::Mortality) = :mortality
process_kind(::Remineralization) = :remineralization

participants(process::Growth) = (population=process.populations,)
function participants(process::Grazing)
    base = (consumer=process.consumers, resource=process.resources)
    isnothing(process.unassimilated_destination) && return base
    return merge(base, (unassimilated_destination=(process.unassimilated_destination,),))
end
participants(process::Mortality) = (population=process.populations,)
participants(process::Remineralization) =
    (source=process.sources, destination=process.destinations)

function drivers(process::Growth)
    return (light=process.light,)
end
drivers(::Union{Grazing,Mortality,Remineralization}) = NamedTuple()

rate_axes(::Growth) = (:population,)
rate_axes(::Grazing) = (:consumer, :resource)
rate_axes(::Mortality) = (:population,)
rate_axes(::Remineralization) = (:source,)
