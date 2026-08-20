"""Abstract supertype for scientific process objects."""
abstract type AbstractProcess end

"""Abstract supertype for concrete scientific formulations."""
abstract type AbstractFormulation end

struct Smith <: AbstractFormulation end
struct Monod <: AbstractFormulation end
struct MultiplicativeFactors <: AbstractFormulation end

"""Abstract supertype for named multiplicative process-rate factors."""
abstract type AbstractFactor end
struct PreferentialGrazing <: AbstractFormulation end
struct LinearMortality <: AbstractFormulation end
struct QuadraticMortality <: AbstractFormulation end
struct LinearRemineralization <: AbstractFormulation end
struct PartitionRouting <: AbstractFormulation end

formulation_tag(::Smith) = :smith
formulation_tag(::Monod) = :monod
formulation_tag(::MultiplicativeFactors) = :multiplicative
formulation_tag(::PreferentialGrazing) = :preferential
formulation_tag(::LinearMortality) = :linear
formulation_tag(::QuadraticMortality) = :quadratic
formulation_tag(::LinearRemineralization) = :linear
formulation_tag(::PartitionRouting) = :partition

_resolve_light(::Val{:smith}) = Smith()
_resolve_response(::Val{:monod}) = Monod()
_resolve_grazing(::Val{:preferential}) = PreferentialGrazing()
_resolve_mortality(::Val{:linear}) = LinearMortality()
_resolve_mortality(::Val{:quadratic}) = QuadraticMortality()
_resolve_remineralization(::Val{:linear}) = LinearRemineralization()
_resolve_routing(::Val{:partition}) = PartitionRouting()

function _unknown_formulation(kind::Symbol, formulation::Symbol)
    throw(ArgumentError("unknown $(kind) formulation :$(formulation)"))
end

_resolve_light(::Val{F}) where {F} = _unknown_formulation(:light, F)
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

"""Light-dependent multiplicative factor in a process rate."""
struct Light{F<:AbstractFormulation} <: AbstractFactor
    formulation::F
    driver::Symbol
end

Light(formulation::Symbol; driver::Symbol) = Light(_resolve_light(Val(formulation)); driver)
Light(formulation::Smith; driver::Symbol) = Light(formulation, driver)
function Light(formulation::AbstractFormulation; driver::Symbol)
    throw(ArgumentError("$(typeof(formulation)) is not a light formulation"))
end

"""Single-resource multiplicative nutrient factor used by processes such as growth."""
struct NutrientResponse{F<:AbstractFormulation} <: AbstractFactor
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

function _canonical_factors(factors::NamedTuple)
    isempty(factors) && throw(ArgumentError("growth `factors` cannot be empty"))
    all(factor -> factor isa AbstractFactor, values(factors)) || throw(
        ArgumentError("growth `factors` values must be process factors"),
    )
    names = sort!(collect(keys(factors)); by=String)
    names_tuple = Tuple(names)
    return NamedTuple{names_tuple}(Tuple(getproperty(factors, name) for name in names))
end

"""Population growth process whose named top-level factors multiply."""
struct Growth{F<:NamedTuple} <: AbstractProcess
    formulation::MultiplicativeFactors
    populations::Tuple
    factors::F

    function Growth(populations::Tuple, factors::NamedTuple)
        canonical = _canonical_factors(factors)
        return new{typeof(canonical)}(MultiplicativeFactors(), populations, canonical)
    end
end

function Growth(; population=nothing, populations=nothing, factors::NamedTuple)
    population_refs = _participant_tuple(:population, population, populations)
    return Growth(population_refs, factors)
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
formulation(factor::AbstractFactor) = factor.formulation
formulation(routing::ProductRouting) = routing.formulation

"""Return the named multiplicative factors attached to a process."""
factors(process::Growth) = process.factors
factors(::Union{Grazing,Mortality,Remineralization}) = NamedTuple()

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

factor_drivers(factor::Light) = (driver=factor.driver,)
factor_drivers(::NutrientResponse) = NamedTuple()

function drivers(process::Growth)
    names = Symbol[]
    identities = Symbol[]
    for (name, factor) in pairs(process.factors)
        factor_bindings = factor_drivers(factor)
        isempty(factor_bindings) && continue
        length(factor_bindings) == 1 || throw(
            ArgumentError(
                "growth factor :$name declares multiple drivers; nested driver paths are not implemented"
            ),
        )
        push!(names, name)
        push!(identities, only(values(factor_bindings)))
    end
    return NamedTuple{Tuple(names)}(Tuple(identities))
end
drivers(::Union{Grazing,Mortality,Remineralization}) = NamedTuple()

rate_axes(::Growth) = (:population,)
rate_axes(::Grazing) = (:consumer, :resource)
rate_axes(::Mortality) = (:population,)
rate_axes(::Remineralization) = (:source,)
