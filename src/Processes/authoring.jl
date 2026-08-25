"""Abstract supertype for scientific process objects."""
abstract type AbstractProcess end

"""Abstract supertype for concrete scientific formulations."""
abstract type AbstractFormulation end

struct Smith <: AbstractFormulation end
struct Geider <: AbstractFormulation end
struct Monod <: AbstractFormulation end
struct Liebig <: AbstractFormulation end

"""Differentiable Frank t-norm nutrient-combination formulation."""
struct FrankTNorm <: AbstractFormulation end

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
struct DOMPOMRouting <: AbstractFormulation end

abstract type AbstractStoichiometry end

function _canonical_bindings(bindings::NamedTuple)
    names = sort!(collect(keys(bindings)); by=String)
    names_tuple = Tuple(names)
    values = Tuple(begin
        value = getproperty(bindings, name)
        if value isa Symbol
            value
        elseif value isa NamedTuple
            all(entry -> entry isa Symbol, Base.values(value)) || throw(
                ArgumentError(
                    "binding map for :$name must map qualifier values to parameter Symbols"
                ),
            )
            qualifier_names = sort!(collect(keys(value)); by=String)
            qualifier_names_tuple = Tuple(qualifier_names)
            qualifier_values = Tuple(getproperty(value, key) for key in qualifier_names)
            NamedTuple{qualifier_names_tuple}(qualifier_values)
        else
            throw(ArgumentError(
                "binding :$name must be a parameter Symbol or one-level qualifier NamedTuple",
            ))
        end
    end for name in names)
    return NamedTuple{names_tuple}(values)
end

"""Return setup-only authored model-parameter bindings for one slot-owning node."""
function authored_parameter_bindings end

"""Fixed conversion from one reference currency to process target currencies.

Each bound ratio is the amount of its target currency per unit reference currency.
"""
struct FixedStoichiometry <: AbstractStoichiometry
    reference::Symbol
    bindings::NamedTuple
end

FixedStoichiometry(; reference::Symbol, bindings::NamedTuple=NamedTuple()) =
    FixedStoichiometry(reference, _canonical_bindings(bindings))

authored_parameter_bindings(stoichiometry::FixedStoichiometry) = stoichiometry.bindings

"""Stable semantic formulation identity used by recipes and provenance."""
function formulation_tag(formulation::AbstractFormulation)
    throw(ArgumentError(
        "no semantic formulation tag defined for $(typeof(formulation)); " *
        "extend `formulation_tag` for external formulations that require semantic identity",
    ))
end

formulation_tag(::Smith) = :smith
formulation_tag(::Geider) = :geider
formulation_tag(::Monod) = :monod
formulation_tag(::Liebig) = :liebig
formulation_tag(::FrankTNorm) = :frank_tnorm
formulation_tag(::Q10) = :q10
formulation_tag(::MultiplicativeFactors) = :multiplicative
formulation_tag(::IdealizedGrazing) = :idealized
formulation_tag(::PreferentialGrazing) = :preferential
formulation_tag(::HeterotrophicConsumption) = :heterotrophic
formulation_tag(::LinearMortality) = :linear
formulation_tag(::QuadraticMortality) = :quadratic
formulation_tag(::LinearRemineralization) = :linear
formulation_tag(::DOMPOMRouting) = :dom_pom
formulation_tag(::FixedStoichiometry) = :fixed

"""Return explicit scientific state carried by a formulation for recipe serialization.

Stateless formulations require no method. Formulations with fields must define this hook so
private implementation fields cannot silently become, or disappear from, durable scientific
identity.
"""
function formulation_recipe_fields(formulation::AbstractFormulation)
    fieldcount(typeof(formulation)) == 0 && return NamedTuple()
    throw(ArgumentError(
        "stateful formulation $(typeof(formulation)) must implement `formulation_recipe_fields`",
    ))
end

function _canonical_participants(role::Symbol, values)
    values isa Symbol && (values = (values,))
    values isa Tuple || throw(ArgumentError("participant `$role` must be a Symbol or tuple"))
    isempty(values) && throw(ArgumentError("participant `$role` cannot be empty"))
    all(value -> value isa Symbol, values) || throw(
        ArgumentError("participant `$role` must contain only Symbols"),
    )
    return values
end

"""Light-dependent multiplicative factor in a process rate."""
struct Light{F<:AbstractFormulation} <: AbstractFactor
    formulation::F
    driver::Symbol
    bindings::NamedTuple
end

function Light(
    formulation::Union{Smith,Geider}; driver::Symbol, bindings::NamedTuple=NamedTuple()
)
    return Light(formulation, driver, _canonical_bindings(bindings))
end

authored_parameter_bindings(factor::Light) = factor.bindings

"""Single-resource multiplicative nutrient factor used by processes such as growth."""
struct NutrientResponse{F<:AbstractFormulation} <: AbstractFactor
    formulation::F
    resource::Symbol
    bindings::NamedTuple
end

function NutrientResponse(
    formulation::Monod; resource::Symbol, bindings::NamedTuple=NamedTuple()
)
    return NutrientResponse(formulation, resource, _canonical_bindings(bindings))
end

authored_parameter_bindings(factor::NutrientResponse) = factor.bindings

"""Temperature-dependent multiplicative process-rate factor."""
struct Temperature{F<:AbstractFormulation} <: AbstractFactor
    formulation::F
    driver::Symbol
    bindings::NamedTuple
end

function Temperature(
    formulation::Q10; driver::Symbol=:temperature, bindings::NamedTuple=NamedTuple()
)
    return Temperature(formulation, driver, _canonical_bindings(bindings))
end

authored_parameter_bindings(factor::Temperature) = factor.bindings

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
    bindings::NamedTuple
end

function Nutrients(
    formulation::Union{Liebig,FrankTNorm};
    responses::NamedTuple,
    bindings::NamedTuple=NamedTuple(),
)
    return Nutrients(
        formulation, _canonical_responses(responses), _canonical_bindings(bindings)
    )
end

authored_parameter_bindings(factor::Nutrients) = factor.bindings

"""Return the stable semantic kind of an authored factor."""
factor_kind(::Light) = :light
factor_kind(::Temperature) = :temperature
factor_kind(::NutrientResponse) = :nutrient_response
factor_kind(::Nutrients) = :nutrients

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

"""Legacy DOM/POM product allocation retained until DARWIN migrates to `Products`."""
struct ProductRouting{F<:DOMPOMRouting,P,S}
    formulation::F
    pools::P
    stoichiometry::S
    bindings::NamedTuple
end

function ProductRouting(
    formulation::DOMPOMRouting;
    pools::NamedTuple,
    stoichiometry::FixedStoichiometry,
    bindings::NamedTuple=NamedTuple(),
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
    return ProductRouting(
        formulation, pools, stoichiometry, _canonical_bindings(bindings)
    )
end

authored_parameter_bindings(routing::ProductRouting) = routing.bindings

"""Conservative allocation of one process product flux among named destinations.

For `N` destinations, specify either `N - 1` named fractions or all `N` fractions.
When one fraction is omitted, that product receives the balance
`1 - sum(explicit fractions)`. Fully explicit fractions are validated to sum to one
after parameter values are resolved. A single destination requires no fractions.
"""
struct Products{T,F,I}
    targets::T
    fractions::F
    inferred::I
end

function Products(
    targets::NamedTuple;
    fractions::NamedTuple=NamedTuple(),
)
    isempty(targets) && throw(ArgumentError("products `targets` cannot be empty"))
    all(target -> target isa Symbol, values(targets)) || throw(
        ArgumentError("product targets must be component Symbols"),
    )

    target_names = sort!(collect(keys(targets)); by=String)
    canonical_targets = NamedTuple{Tuple(target_names)}(
        Tuple(getproperty(targets, name) for name in target_names)
    )

    fraction_names = sort!(collect(keys(fractions)); by=String)
    all(name -> name in keys(canonical_targets), fraction_names) || throw(
        ArgumentError("product `fractions` must be keyed by product names declared in `targets`"),
    )
    all(value -> value isa Symbol, values(fractions)) || throw(
        ArgumentError("product `fractions` values must be parameter Symbols"),
    )
    canonical_fractions = NamedTuple{Tuple(fraction_names)}(
        Tuple(getproperty(fractions, name) for name in fraction_names)
    )

    nproducts = length(canonical_targets)
    nfractions = length(canonical_fractions)
    if nproducts == 1
        nfractions == 0 || throw(
            ArgumentError("single-destination products do not take fractions"),
        )
        inferred = only(keys(canonical_targets))
    elseif nfractions == nproducts - 1
        inferred = only(
            name for name in keys(canonical_targets) if !(name in keys(canonical_fractions))
        )
    elseif nfractions == nproducts
        inferred = nothing
    else
        throw(ArgumentError(
            "products has $nproducts destinations but $nfractions fractions; specify either " *
            "$(nproducts - 1) fractions (the omitted product receives the balance) or all " *
            "$nproducts fractions explicitly",
        ))
    end

    return Products(canonical_targets, canonical_fractions, inferred)
end

Products(destination::Symbol) = Products((product=destination,))

authored_parameter_bindings(products::Products) = isempty(products.fractions) ?
    NamedTuple() : (fraction=products.fractions,)

_canonical_products(::Nothing) = nothing
_canonical_products(products::Products) = products
_canonical_products(destination::Symbol) = Products(destination)

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
            populations, canonical, source, stoichiometry
        )
    end
end

function Growth(;
    populations,
    factors::NamedTuple,
    source=nothing,
    stoichiometry=nothing,
)
    population_refs = _canonical_participants(:populations, populations)
    return Growth(population_refs, factors, source, stoichiometry)
end

"""Canonical consumer-resource process with optional factors and unassimilated products."""
struct Consumption{F<:AbstractFormulation,A<:NamedTuple,R} <: AbstractProcess
    formulation::F
    consumers::Tuple
    resources::Tuple
    factors::A
    routing::R
    bindings::NamedTuple
end

function Consumption(
    formulation::Union{IdealizedGrazing,PreferentialGrazing,HeterotrophicConsumption};
    consumers,
    resources,
    factors::NamedTuple=NamedTuple(),
    unassimilated_products=nothing,
    routing=nothing,
    bindings::NamedTuple=NamedTuple(),
)
    consumer_refs = _canonical_participants(:consumers, consumers)
    resource_refs = _canonical_participants(:resources, resources)
    canonical = _canonical_factors(factors; allow_empty=true)
    !isnothing(unassimilated_products) && !isnothing(routing) && throw(
        ArgumentError("specify `unassimilated_products` or legacy `routing`, not both"),
    )
    products = isnothing(unassimilated_products) ? routing : _canonical_products(unassimilated_products)
    isnothing(products) || products isa Union{Products,ProductRouting} || throw(
        ArgumentError("consumption products must be a component Symbol or Products"),
    )
    return Consumption(
        formulation, consumer_refs, resource_refs, canonical, products, _canonical_bindings(bindings)
    )
end

authored_parameter_bindings(process::Consumption) = process.bindings

"""Population mortality process with optional products."""
struct Mortality{F<:AbstractFormulation,R} <: AbstractProcess
    formulation::F
    populations::Tuple
    routing::R
    bindings::NamedTuple
end

function Mortality(
    formulation::Union{LinearMortality,QuadraticMortality};
    populations,
    products=nothing,
    routing=nothing,
    bindings::NamedTuple=NamedTuple(),
)
    population_refs = _canonical_participants(:populations, populations)
    !isnothing(products) && !isnothing(routing) && throw(
        ArgumentError("specify `products` or legacy `routing`, not both"),
    )
    product_spec = isnothing(products) ? routing : _canonical_products(products)
    isnothing(product_spec) || product_spec isa Union{Products,ProductRouting} || throw(
        ArgumentError("mortality products must be a component Symbol or Products"),
    )
    return Mortality(
        formulation, population_refs, product_spec, _canonical_bindings(bindings)
    )
end

authored_parameter_bindings(process::Mortality) = process.bindings

"""Source-to-destination remineralization process."""
struct Remineralization{F<:AbstractFormulation} <: AbstractProcess
    formulation::F
    sources::Tuple
    destinations::Tuple
    bindings::NamedTuple
end

function Remineralization(
    formulation::LinearRemineralization;
    sources,
    destinations,
    bindings::NamedTuple=NamedTuple(),
)
    source_refs = _canonical_participants(:sources, sources)
    destination_refs = _canonical_participants(:destinations, destinations)
    return Remineralization(
        formulation, source_refs, destination_refs, _canonical_bindings(bindings)
    )
end

authored_parameter_bindings(process::Remineralization) = process.bindings

formulation(::Growth) = MultiplicativeFactors()
formulation(process::AbstractProcess) = process.formulation
formulation(factor::AbstractFactor) = factor.formulation
formulation(routing::ProductRouting) = routing.formulation

"""Return the named multiplicative factors attached to a process."""
factors(::AbstractProcess) = NamedTuple()
factors(process::Union{Growth,Consumption}) = process.factors

process_routing(::AbstractProcess) = nothing
process_routing(process::Union{Consumption,Mortality}) =
    process.routing isa ProductRouting ? process.routing : nothing

process_products(::AbstractProcess) = nothing
process_products(process::Union{Consumption,Mortality}) =
    process.routing isa Products ? process.routing : nothing
product_path(::Mortality) = (:products,)
product_path(::Consumption) = (:unassimilated_products,)
product_recipe_key(::Mortality) = :products
product_recipe_key(::Consumption) = :unassimilated_products

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
    identities = Symbol[
        input.identity for input in factor_inputs(factor) if input isa FactorDriver
    ]
    for child in values(factor_children(factor))
        append!(identities, _factor_driver_identities(child))
    end
    return Tuple(identities)
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
