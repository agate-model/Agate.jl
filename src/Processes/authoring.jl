"""Abstract supertype for scientific process objects."""
abstract type AbstractProcess end

"""Abstract supertype for concrete scientific formulations."""
abstract type AbstractFormulation end

struct Smith <: AbstractFormulation end
struct Geider <: AbstractFormulation end
struct Monod <: AbstractFormulation end
struct Liebig <: AbstractFormulation end
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

_resolve_light(::Val{:smith}) = Smith()
_resolve_light(::Val{:geider}) = Geider()
_resolve_response(::Val{:monod}) = Monod()
_resolve_nutrients(::Val{:liebig}) = Liebig()
_resolve_temperature(::Val{:q10}) = Q10()
_resolve_grazing(::Val{:idealized}) = IdealizedGrazing()
_resolve_grazing(::Val{:preferential}) = PreferentialGrazing()
_resolve_consumption(::Val{:heterotrophic}) = HeterotrophicConsumption()
_resolve_mortality(::Val{:linear}) = LinearMortality()
_resolve_mortality(::Val{:quadratic}) = QuadraticMortality()
_resolve_remineralization(::Val{:linear}) = LinearRemineralization()
_resolve_routing(::Val{:direct}) = DirectRouting()
_resolve_routing(::Val{:partition}) = PartitionRouting()
_resolve_routing(::Val{:dom_pom}) = DOMPOMRouting()

function _unknown_formulation(kind::Symbol, formulation::Symbol)
    throw(ArgumentError("unknown $(kind) formulation :$(formulation)"))
end

_resolve_light(::Val{F}) where {F} = _unknown_formulation(:light, F)
_resolve_response(::Val{F}) where {F} = _unknown_formulation(:nutrient_response, F)
_resolve_nutrients(::Val{F}) where {F} = _unknown_formulation(:nutrients, F)
_resolve_temperature(::Val{F}) where {F} = _unknown_formulation(:temperature, F)
_resolve_grazing(::Val{F}) where {F} = _unknown_formulation(:grazing, F)
_resolve_consumption(::Val{F}) where {F} = _unknown_formulation(:consumption, F)
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
Light(formulation::Union{Smith,Geider}; driver::Symbol) = Light(formulation, driver)
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

"""Temperature-dependent multiplicative process-rate factor."""
struct Temperature{F<:AbstractFormulation} <: AbstractFactor
    formulation::F
    driver::Symbol
end

Temperature(formulation::Symbol; driver::Symbol=:temperature) =
    Temperature(_resolve_temperature(Val(formulation)); driver)
Temperature(formulation::Q10; driver::Symbol=:temperature) = Temperature(formulation, driver)
function Temperature(formulation::AbstractFormulation; driver::Symbol=:temperature)
    throw(ArgumentError("$(typeof(formulation)) is not a temperature formulation"))
end

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

function Nutrients(formulation::Liebig; responses::NamedTuple)
    return Nutrients(formulation, _canonical_responses(responses))
end

Nutrients(formulation::Symbol; responses::NamedTuple) =
    Nutrients(_resolve_nutrients(Val(formulation)); responses)
function Nutrients(formulation::AbstractFormulation; responses::NamedTuple)
    throw(ArgumentError("$(typeof(formulation)) is not a nutrient-combination formulation"))
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

ProductRouting(formulation::Symbol; kwargs...) =
    ProductRouting(_resolve_routing(Val(formulation)); kwargs...)
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
function ProductRouting(formulation::AbstractFormulation; kwargs...)
    throw(ArgumentError("$(typeof(formulation)) is not a product-routing formulation"))
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

"""Consumer-resource grazing process."""
struct Grazing{F<:AbstractFormulation,R} <: AbstractProcess
    formulation::F
    consumers::Tuple
    resources::Tuple
    unassimilated_destination::Union{Nothing,Symbol}
    routing::R
end

function Grazing(
    formulation::Union{IdealizedGrazing,PreferentialGrazing};
    consumer=nothing,
    consumers=nothing,
    resource=nothing,
    resources=nothing,
    unassimilated_destination=nothing,
    routing=nothing,
)
    consumer_refs = _participant_tuple(:consumer, consumer, consumers)
    resource_refs = _participant_tuple(:resource, resource, resources)
    destination = _optional_reference(:unassimilated_destination, unassimilated_destination)
    isnothing(destination) || isnothing(routing) || throw(
        ArgumentError("grazing cannot specify both `unassimilated_destination` and `routing`"),
    )
    isnothing(routing) || routing isa ProductRouting || throw(
        ArgumentError("grazing `routing` must be a ProductRouting"),
    )
    isnothing(routing) || routing.formulation isa DOMPOMRouting || throw(
        ArgumentError("grazing product routing must use the :dom_pom formulation"),
    )
    return Grazing(formulation, consumer_refs, resource_refs, destination, routing)
end

Grazing(formulation::Symbol; kwargs...) = Grazing(_resolve_grazing(Val(formulation)); kwargs...)
function Grazing(formulation::AbstractFormulation; kwargs...)
    throw(ArgumentError("$(typeof(formulation)) is not a grazing formulation"))
end

"""Consumer-resource heterotrophic consumption with optional multiplicative factors."""
struct Consumption{F<:AbstractFormulation,A<:NamedTuple} <: AbstractProcess
    formulation::F
    consumers::Tuple
    resources::Tuple
    unassimilated_destination::Union{Nothing,Symbol}
    factors::A
end

function Consumption(
    formulation::HeterotrophicConsumption;
    consumer=nothing,
    consumers=nothing,
    resource=nothing,
    resources=nothing,
    unassimilated_destination=nothing,
    factors::NamedTuple=NamedTuple(),
)
    consumer_refs = _participant_tuple(:consumer, consumer, consumers)
    resource_refs = _participant_tuple(:resource, resource, resources)
    destination = _optional_reference(:unassimilated_destination, unassimilated_destination)
    canonical = _canonical_factors(factors; allow_empty=true)
    return Consumption(formulation, consumer_refs, resource_refs, destination, canonical)
end

Consumption(formulation::Symbol; kwargs...) =
    Consumption(_resolve_consumption(Val(formulation)); kwargs...)
function Consumption(formulation::AbstractFormulation; kwargs...)
    throw(ArgumentError("$(typeof(formulation)) is not a consumption formulation"))
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
factors(process::Union{Growth,Consumption}) = process.factors
factors(::Union{Grazing,Mortality,Remineralization}) = NamedTuple()

process_kind(::Growth) = :growth
process_kind(::Grazing) = :grazing
process_kind(::Consumption) = :consumption
process_kind(::Mortality) = :mortality
process_kind(::Remineralization) = :remineralization

function participants(process::Growth)
    base = (population=process.populations,)
    isnothing(process.source) && return base
    return merge(base, (source=(process.source,),))
end
function participants(process::Grazing)
    base = (consumer=process.consumers, resource=process.resources)
    isnothing(process.unassimilated_destination) && return base
    return merge(base, (unassimilated_destination=(process.unassimilated_destination,),))
end
function participants(process::Consumption)
    base = (consumer=process.consumers, resource=process.resources)
    isnothing(process.unassimilated_destination) && return base
    return merge(base, (unassimilated_destination=(process.unassimilated_destination,),))
end
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

drivers(process::Union{Growth,Consumption}) = _factor_drivers(process)
drivers(::Union{Grazing,Mortality,Remineralization}) = NamedTuple()

rate_axes(::Growth) = (:population,)
rate_axes(::Grazing) = (:consumer, :resource)
rate_axes(::Consumption) = (:consumer, :resource)
rate_axes(::Mortality) = (:population,)
rate_axes(::Remineralization) = (:source,)
