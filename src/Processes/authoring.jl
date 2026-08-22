"""Abstract supertype for scientific process objects."""
abstract type AbstractProcess end

"""Abstract supertype for concrete scientific formulations."""
abstract type AbstractFormulation end

struct Smith <: AbstractFormulation end
struct Geider <: AbstractFormulation end
struct Monod <: AbstractFormulation end
struct Liebig <: AbstractFormulation end

"""Differentiable Frank t-norm nutrient combination with configurable sharpness."""
struct Frank{S} <: AbstractFormulation
    sharpness::S
end
Frank() = Frank(DEFAULT_FRANK_SHARPNESS)

struct Q10 <: AbstractFormulation end
struct MultiplicativeFactors <: AbstractFormulation end

"""Abstract supertype for named multiplicative process-rate factors."""
abstract type AbstractFactor end
struct IdealizedGrazing <: AbstractFormulation end
struct PreferentialGrazing <: AbstractFormulation end
struct HeterotrophicConsumption <: AbstractFormulation end
struct LinearMortality <: AbstractFormulation end
struct QuadraticMortality <: AbstractFormulation end
struct LinearRemineralization <: AbstractFormulation end
struct DirectRouting <: AbstractFormulation end
struct PartitionRouting <: AbstractFormulation end
struct DOMPOMRouting <: AbstractFormulation end

abstract type AbstractStoichiometry end

"""Fixed conversion from one reference currency to process target currencies.

Each bound ratio is the amount of its target currency per unit reference currency.
"""
struct FixedStoichiometry <: AbstractStoichiometry
    reference::Symbol
end

FixedStoichiometry(; reference::Symbol) = FixedStoichiometry(reference)

formulation_tag(::Smith) = :smith
formulation_tag(::Geider) = :geider
formulation_tag(::Monod) = :monod
formulation_tag(::Liebig) = :liebig
formulation_tag(::Frank) = :frank
formulation_tag(::Q10) = :q10
formulation_tag(::MultiplicativeFactors) = :multiplicative
formulation_tag(::IdealizedGrazing) = :idealized
formulation_tag(::PreferentialGrazing) = :preferential
formulation_tag(::HeterotrophicConsumption) = :heterotrophic
formulation_tag(::LinearMortality) = :linear
formulation_tag(::QuadraticMortality) = :quadratic
formulation_tag(::LinearRemineralization) = :linear
formulation_tag(::DirectRouting) = :direct
formulation_tag(::PartitionRouting) = :partition
formulation_tag(::DOMPOMRouting) = :dom_pom
formulation_tag(::FixedStoichiometry) = :fixed

function formulation_for(::Type{H}, ::Val{F}) where {H,F}
    throw(ArgumentError("unknown formulation :$F for $H"))
end

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

formulation_for(::Type{Light}, ::Val{:smith}) = Smith()
formulation_for(::Type{Light}, ::Val{:geider}) = Geider()

Light(formulation::Symbol; driver::Symbol) =
    Light(formulation_for(Light, Val(formulation)); driver)
Light(formulation::Union{Smith,Geider}; driver::Symbol) = Light(formulation, driver)

"""Single-resource multiplicative nutrient factor used by processes such as growth."""
struct NutrientResponse{F<:AbstractFormulation} <: AbstractFactor
    formulation::F
    resource::Symbol
end

formulation_for(::Type{NutrientResponse}, ::Val{:monod}) = Monod()

NutrientResponse(formulation::Symbol; resource::Symbol) =
    NutrientResponse(formulation_for(NutrientResponse, Val(formulation)); resource)
NutrientResponse(formulation::Monod; resource::Symbol) = NutrientResponse(formulation, resource)

"""Temperature-dependent multiplicative process-rate factor."""
struct Temperature{F<:AbstractFormulation} <: AbstractFactor
    formulation::F
    driver::Symbol
end

formulation_for(::Type{Temperature}, ::Val{:q10}) = Q10()

Temperature(formulation::Symbol; driver::Symbol=:temperature) =
    Temperature(formulation_for(Temperature, Val(formulation)); driver)
Temperature(formulation::Q10; driver::Symbol=:temperature) = Temperature(formulation, driver)

function _canonical_responses(responses::NamedTuple)
    isempty(responses) && throw(ArgumentError("nutrient `responses` cannot be empty"))
    all(response -> response isa NutrientResponse, values(responses)) || throw(
        ArgumentError("nutrient `responses` values must be NutrientResponse factors"),
    )
    names = sort!(collect(keys(responses)); by=String)
    names_tuple = Tuple(names)
    return NamedTuple{names_tuple}(Tuple(getproperty(responses, name) for name in names))
end

"""Multi-resource nutrient factor with formulation-owned response composition.

When paired with `FixedStoichiometry`, response names identify target currencies.
"""
struct Nutrients{F<:AbstractFormulation,R<:NamedTuple} <: AbstractFactor
    formulation::F
    responses::R
end

formulation_for(::Type{Nutrients}, ::Val{:liebig}) = Liebig()
formulation_for(::Type{Nutrients}, ::Val{:frank}) = Frank()

function Nutrients(formulation::Union{Liebig,Frank}; responses::NamedTuple)
    return Nutrients(formulation, _canonical_responses(responses))
end

function Nutrients(
    formulation::Symbol;
    responses::NamedTuple,
    sharpness=nothing,
)
    if formulation === :frank && !isnothing(sharpness)
        return Nutrients(Frank(sharpness); responses)
    end
    isnothing(sharpness) || throw(
        ArgumentError("`sharpness` is only valid for the :frank nutrient formulation"),
    )
    return Nutrients(formulation_for(Nutrients, Val(formulation)); responses)
end

abstract type AbstractFactorInput end

"""External driver read required by one scientific factor."""
struct FactorDriver <: AbstractFactorInput
    identity::Symbol
end

"""Scalar model-component read required by one scientific factor."""
struct FactorComponent <: AbstractFactorInput
    component::Symbol
end

factor_inputs(::AbstractFactor) = ()
factor_inputs(factor::Light) = (FactorDriver(factor.driver),)
factor_inputs(factor::Temperature) = (FactorDriver(factor.driver),)
factor_inputs(factor::NutrientResponse) = (FactorComponent(factor.resource),)

factor_children(::AbstractFactor) = NamedTuple()
factor_children(factor::Nutrients) = factor.responses

factor_parameter_context(::AbstractFactor) = NamedTuple()
factor_parameter_context(factor::NutrientResponse) = (resource=factor.resource,)

factor_child_path(path::Tuple, ::AbstractFactor, name::Symbol) = (path..., name)
factor_child_path(path::Tuple, ::Nutrients, name::Symbol) = (path..., :responses, name)

"""Product-routing sub-formulation for coupled process effects."""
struct ProductRouting{F<:AbstractFormulation,R,E,P,S}
    formulation::F
    retained::R
    exported::E
    pools::P
    stoichiometry::S
end

formulation_for(::Type{ProductRouting}, ::Val{:direct}) = DirectRouting()
formulation_for(::Type{ProductRouting}, ::Val{:partition}) = PartitionRouting()
formulation_for(::Type{ProductRouting}, ::Val{:dom_pom}) = DOMPOMRouting()

ProductRouting(formulation::Symbol; kwargs...) =
    ProductRouting(formulation_for(ProductRouting, Val(formulation)); kwargs...)
ProductRouting(formulation::DirectRouting; destination::Symbol) =
    ProductRouting(formulation, destination, nothing, nothing, nothing)
ProductRouting(formulation::PartitionRouting; retained::Symbol, exported::Symbol) =
    ProductRouting(formulation, retained, exported, nothing, nothing)
function ProductRouting(
    formulation::DOMPOMRouting; pools::NamedTuple, stoichiometry::FixedStoichiometry
)
    keys(pools) == (:DOM, :POM) || throw(
        ArgumentError("DOM/POM routing `pools` must have exactly the keys (:DOM, :POM)"),
    )
    dom = pools.DOM
    pom = pools.POM
    dom isa NamedTuple && pom isa NamedTuple || throw(
        ArgumentError("DOM/POM routing pools must be named currency-to-component mappings"),
    )
    keys(dom) == keys(pom) || throw(
        ArgumentError("DOM and POM routing pools must declare the same currencies"),
    )
    stoichiometry.reference in keys(dom) || throw(
        ArgumentError(
            "routing pools must include stoichiometric reference currency :$(stoichiometry.reference)"
        ),
    )
    valid_targets = all(value -> value isa Symbol, values(dom)) &&
                    all(value -> value isa Symbol, values(pom))
    valid_targets || throw(
        ArgumentError("DOM/POM routing pool targets must be component Symbols"),
    )
    return ProductRouting(formulation, nothing, nothing, pools, stoichiometry)
end

function _canonical_factors(factors::NamedTuple; allow_empty::Bool=false)
    isempty(factors) && !allow_empty && throw(ArgumentError("process `factors` cannot be empty"))
    all(factor -> factor isa AbstractFactor, values(factors)) || throw(
        ArgumentError("process `factors` values must be process factors"),
    )
    names = sort!(collect(keys(factors)); by=String)
    names_tuple = Tuple(names)
    return NamedTuple{names_tuple}(Tuple(getproperty(factors, name) for name in names))
end

"""Population growth process whose named top-level factors multiply."""
struct Growth{F<:NamedTuple,S,T} <: AbstractProcess
    formulation::MultiplicativeFactors
    populations::Tuple
    factors::F
    source::S
    stoichiometry::T

    function Growth(populations::Tuple, factors::NamedTuple, source, stoichiometry)
        canonical = _canonical_factors(factors)
        isnothing(source) || source isa Symbol ||
            throw(ArgumentError("growth `source` must be a Symbol"))
        isnothing(stoichiometry) || stoichiometry isa AbstractStoichiometry || throw(
            ArgumentError("growth `stoichiometry` must be an AbstractStoichiometry"),
        )
        return new{typeof(canonical),typeof(source),typeof(stoichiometry)}(
            MultiplicativeFactors(), populations, canonical, source, stoichiometry
        )
    end
end

function Growth(;
    population=nothing,
    populations=nothing,
    factors::NamedTuple,
    source=nothing,
    stoichiometry=nothing,
)
    population_refs = _participant_tuple(:population, population, populations)
    return Growth(population_refs, factors, source, stoichiometry)
end

"""Canonical consumer-resource process with optional factors and product routing."""
struct Consumption{F<:AbstractFormulation,A<:NamedTuple,R} <: AbstractProcess
    formulation::F
    consumers::Tuple
    resources::Tuple
    factors::A
    routing::R
end

formulation_for(::Type{Consumption}, ::Val{:idealized}) = IdealizedGrazing()
formulation_for(::Type{Consumption}, ::Val{:preferential}) = PreferentialGrazing()
formulation_for(::Type{Consumption}, ::Val{:heterotrophic}) = HeterotrophicConsumption()

function Consumption(
    formulation::Union{IdealizedGrazing,PreferentialGrazing,HeterotrophicConsumption};
    consumer=nothing,
    consumers=nothing,
    resource=nothing,
    resources=nothing,
    unassimilated_destination=nothing,
    factors::NamedTuple=NamedTuple(),
    routing=nothing,
)
    consumer_refs = _participant_tuple(:consumer, consumer, consumers)
    resource_refs = _participant_tuple(:resource, resource, resources)
    canonical = _canonical_factors(factors; allow_empty=true)
    destination = _optional_reference(:unassimilated_destination, unassimilated_destination)
    isnothing(destination) || isnothing(routing) || throw(
        ArgumentError("consumption cannot specify both `unassimilated_destination` and `routing`"),
    )
    isnothing(routing) || routing isa ProductRouting || throw(
        ArgumentError("consumption `routing` must be a ProductRouting"),
    )
    isnothing(destination) || (routing = ProductRouting(DirectRouting(); destination=destination))
    return Consumption(formulation, consumer_refs, resource_refs, canonical, routing)
end

Consumption(formulation::Symbol; kwargs...) =
    Consumption(formulation_for(Consumption, Val(formulation)); kwargs...)

"""Author-facing grazing convenience that immediately returns canonical `Consumption`."""
function Grazing(formulation::Symbol; kwargs...)
    formulation in (:idealized, :preferential) || throw(
        ArgumentError("unknown grazing formulation :$formulation"),
    )
    return Consumption(formulation; kwargs...)
end
Grazing(formulation::Union{IdealizedGrazing,PreferentialGrazing}; kwargs...) =
    Consumption(formulation; kwargs...)

"""Population mortality process with optional product routing."""
struct Mortality{F<:AbstractFormulation,R} <: AbstractProcess
    formulation::F
    populations::Tuple
    routing::R
end

formulation_for(::Type{Mortality}, ::Val{:linear}) = LinearMortality()
formulation_for(::Type{Mortality}, ::Val{:quadratic}) = QuadraticMortality()

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
    Mortality(formulation_for(Mortality, Val(formulation)); kwargs...)

"""Source-to-destination remineralization process."""
struct Remineralization{F<:AbstractFormulation} <: AbstractProcess
    formulation::F
    sources::Tuple
    destinations::Tuple
end

formulation_for(::Type{Remineralization}, ::Val{:linear}) = LinearRemineralization()

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
    Remineralization(formulation_for(Remineralization, Val(formulation)); kwargs...)

formulation(process::AbstractProcess) = process.formulation
formulation(factor::AbstractFactor) = factor.formulation
formulation(routing::ProductRouting) = routing.formulation

"""Return the named multiplicative factors attached to a process."""
factors(::AbstractProcess) = NamedTuple()
factors(process::Union{Growth,Consumption}) = process.factors

process_routing(::AbstractProcess) = nothing
process_routing(process::Union{Consumption,Mortality}) = process.routing

process_stoichiometry(::AbstractProcess) = nothing
process_stoichiometry(process::Growth) = process.stoichiometry

process_kind(::Growth) = :growth
process_kind(::Consumption) = :consumption
process_kind(::Mortality) = :mortality
process_kind(::Remineralization) = :remineralization

"""Whether a consumer-resource formulation uses living consumer-prey interaction matrices."""
uses_living_interactions(::AbstractFormulation) = false
uses_living_interactions(::Union{IdealizedGrazing,PreferentialGrazing}) = true

function participants(process::Growth)
    base = (population=process.populations,)
    isnothing(process.source) && return base
    return merge(base, (source=(process.source,),))
end
participants(process::Consumption) = (consumer=process.consumers, resource=process.resources)
participants(process::Mortality) = (population=process.populations,)
participants(process::Remineralization) =
    (source=process.sources, destination=process.destinations)

function _factor_driver_identities(factor::AbstractFactor)
    identities = Tuple(
        input.identity for input in factor_inputs(factor) if input isa FactorDriver
    )
    for child in values(factor_children(factor))
        identities = (identities..., _factor_driver_identities(child)...)
    end
    return identities
end

function factor_drivers(factor::AbstractFactor)
    identities = unique(_factor_driver_identities(factor))
    isempty(identities) && return NamedTuple()
    length(identities) == 1 || throw(ArgumentError(
        "one named process factor cannot require multiple external drivers",
    ))
    return (driver=only(identities),)
end

function _factor_drivers(process)
    names = Symbol[]
    identities = Symbol[]
    for (name, factor) in pairs(factors(process))
        factor_bindings = factor_drivers(factor)
        isempty(factor_bindings) && continue
        push!(names, name)
        push!(identities, only(values(factor_bindings)))
    end
    return NamedTuple{Tuple(names)}(Tuple(identities))
end

drivers(process::AbstractProcess) = _factor_drivers(process)

rate_axes(::AbstractProcess) = ()
rate_axes(::Growth) = (:population,)
rate_axes(::Consumption) = (:consumer, :resource)
rate_axes(::Mortality) = (:population,)
rate_axes(::Remineralization) = (:source,)
