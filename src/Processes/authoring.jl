"""Abstract supertype for scientific process objects."""
abstract type AbstractProcess end

"""Abstract supertype for concrete scientific formulations."""
abstract type AbstractFormulation end

"""Smith light-limitation formulation."""
struct Smith <: AbstractFormulation end

"""Geider light-response formulation."""
struct Geider <: AbstractFormulation end

"""Monod single-resource limitation formulation."""
struct Monod <: AbstractFormulation end

"""Normalized Droop cellular-quota growth-limitation formulation."""
struct NormalizedDroop <: AbstractFormulation end

"""Monod nutrient-uptake formulation regulated by cellular quota capacity."""
struct QuotaRegulatedMonod <: AbstractFormulation end

"""Liebig minimum response-combination formulation."""
struct Liebig <: AbstractFormulation end

"""Differentiable Frank t-norm nutrient-combination formulation."""
struct FrankTNorm <: AbstractFormulation end

"""Q10 temperature-response formulation."""
struct Q10 <: AbstractFormulation end
struct MultiplicativeFactors <: AbstractFormulation end

"""Abstract supertype for named multiplicative process-rate factors."""
abstract type AbstractFactor end

"""Living-prey grazing formulation with consumer-by-prey palatability and assimilation."""
struct PreferentialGrazing <: AbstractFormulation end

"""Heterotrophic resource-consumption formulation without living-prey interaction weights."""
struct HeterotrophicConsumption <: AbstractFormulation end

"""Linear population-mortality formulation."""
struct LinearMortality <: AbstractFormulation end

"""Quadratic population-mortality formulation."""
struct QuadraticMortality <: AbstractFormulation end

"""Linear source-to-destination remineralization formulation."""
struct LinearRemineralization <: AbstractFormulation end

"""Abstract supertype for process stoichiometry mappings."""
abstract type AbstractStoichiometry end

# Shared authoring canonicalizers

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

function _canonical_participants(role::Symbol, values)
    values isa Symbol && (values = (values,))
    values isa Tuple || throw(ArgumentError("participant `$role` must be a Symbol or tuple"))
    isempty(values) && throw(ArgumentError("participant `$role` cannot be empty"))
    all(value -> value isa Symbol, values) || throw(
        ArgumentError("participant `$role` must contain only Symbols"),
    )
    return values
end

"""Light-dependent multiplicative Growth factor using the Growth rate scale."""
struct Light{F<:Union{Smith,Geider}} <: AbstractFactor
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
struct NutrientResponse{F<:Monod} <: AbstractFactor
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

"""Cellular-quota response used by quota-limited growth.

`target` identifies the internal nutrient inventory and `reference` identifies the
biomass inventory used to form the cellular quota.
"""
struct QuotaResponse{F<:NormalizedDroop} <: AbstractFactor
    formulation::F
    target::PopulationStateRef
    reference::PopulationStateRef
    bindings::NamedTuple
end

function QuotaResponse(
    formulation::NormalizedDroop;
    target::PopulationStateRef,
    reference::PopulationStateRef,
    bindings::NamedTuple=NamedTuple(),
)
    return QuotaResponse(formulation, target, reference, _canonical_bindings(bindings))
end

authored_parameter_bindings(factor::QuotaResponse) = factor.bindings

"""Temperature-dependent multiplicative process-rate factor."""
struct Temperature{F<:Q10} <: AbstractFactor
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
    names = sort!(collect(keys(responses)); by=String)
    names_tuple = Tuple(names)
    return NamedTuple{names_tuple}(Tuple(getproperty(responses, name) for name in names))
end

"""Multi-response nutrient factor with formulation-owned response composition.

External `NutrientResponse` children identify resource currencies for fixed-stoichiometry
growth. Internal `QuotaResponse` children identify cellular quota currencies and affect
growth rate without transferring those nutrients. Each `Nutrients` factor uses one response
category consistently.
"""
struct Nutrients{F<:Union{Liebig,FrankTNorm},R<:NamedTuple} <: AbstractFactor
    formulation::F
    responses::R
    bindings::NamedTuple

    function Nutrients(
        formulation::F, responses::R, bindings::NamedTuple
    ) where {F<:Union{Liebig,FrankTNorm},R<:NamedTuple}
        isempty(responses) && throw(ArgumentError("nutrient `responses` cannot be empty"))
        all(
            response -> response isa Union{NutrientResponse,QuotaResponse}, values(responses)
        ) || throw(ArgumentError(
            "nutrient `responses` values must be NutrientResponse or QuotaResponse factors",
        ))
        return new{F,R}(formulation, responses, bindings)
    end
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

abstract type AbstractFactorInput end

"""External driver read required by one scientific factor."""
struct FactorDriver <: AbstractFactorInput
    identity::Symbol
end

"""Scalar model-component read required by one scientific factor."""
struct FactorComponent <: AbstractFactorInput
    component::Symbol
end

"""Setup-only read of one prognostic population state required by a factor."""
struct FactorPopulationState <: AbstractFactorInput
    reference::PopulationStateRef
end

factor_inputs(::AbstractFactor) = ()
factor_inputs(factor::Light) = (FactorDriver(factor.driver),)
factor_inputs(factor::Temperature) = (FactorDriver(factor.driver),)
factor_inputs(factor::NutrientResponse) = (FactorComponent(factor.resource),)
factor_inputs(factor::QuotaResponse) = (
    FactorPopulationState(factor.target),
    FactorPopulationState(factor.reference),
)

factor_children(::AbstractFactor) = NamedTuple()
factor_children(factor::Nutrients) = factor.responses

factor_parameter_context(::AbstractFactor) = NamedTuple()

factor_child_path(path::Tuple, ::AbstractFactor, name::Symbol) = (path..., name)
factor_child_path(path::Tuple, ::Nutrients, name::Symbol) = (path..., :responses, name)

# Product routing

"""Conservative allocation of one process product flux among named destinations.

Each product target may be either one component Symbol or a named currency-to-component
mapping. Multi-currency products require `FixedStoichiometry`; every product then declares
the same currencies and the stoichiometric reference currency must be present.

For `N` products, specify either `N - 1` named fractions or all `N` fractions. When one
fraction is omitted, that product receives the exact conservative remainder
`1 - sum(supplied fractions)`. When all fractions are supplied, every authored value is used
directly and setup requires their sum to equal one within floating-point tolerance. A single
product requires no fractions.
"""
struct Products{T,F,S}
    targets::T
    fractions::F
    stoichiometry::S
end

function _canonical_product_target(target)
    target isa Symbol && return target
    target isa NamedTuple || throw(
        ArgumentError("product targets must be component Symbols or currency-to-component mappings"),
    )
    isempty(target) && throw(ArgumentError("multi-currency product targets cannot be empty"))
    all(component -> component isa Symbol, values(target)) || throw(
        ArgumentError("multi-currency product targets must map currencies to component Symbols"),
    )
    names = sort!(collect(keys(target)); by=String)
    return NamedTuple{Tuple(names)}(Tuple(getproperty(target, name) for name in names))
end

function Products(
    targets::NamedTuple;
    fractions::NamedTuple=NamedTuple(),
    stoichiometry=nothing,
)
    isempty(targets) && throw(ArgumentError("products `targets` cannot be empty"))

    target_names = sort!(collect(keys(targets)); by=String)
    canonical_target_values = Tuple(
        _canonical_product_target(getproperty(targets, name)) for name in target_names
    )
    canonical_targets = NamedTuple{Tuple(target_names)}(canonical_target_values)

    nested = first(canonical_target_values) isa NamedTuple
    all(target -> (target isa NamedTuple) == nested, canonical_target_values) || throw(
        ArgumentError("product targets cannot mix scalar and multi-currency destinations"),
    )
    if nested
        stoichiometry isa FixedStoichiometry || throw(
            ArgumentError("multi-currency products require FixedStoichiometry"),
        )
        currencies = keys(first(canonical_target_values))
        all(target -> keys(target) == currencies, canonical_target_values) || throw(
            ArgumentError("multi-currency products must declare the same currencies"),
        )
        stoichiometry.reference in currencies || throw(
            ArgumentError(
                "multi-currency products must include stoichiometric reference currency :$(stoichiometry.reference)"
            ),
        )
    elseif !isnothing(stoichiometry)
        throw(ArgumentError("scalar product targets do not take stoichiometry"))
    end

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
            ArgumentError("single-product allocations do not take fractions"),
        )
    elseif !(nfractions in (nproducts - 1, nproducts))
        throw(ArgumentError(
            "products has $nproducts destinations but $nfractions fractions; specify either " *
            "$(nproducts - 1) fractions (the omitted product receives the balance) or all " *
            "$nproducts fractions explicitly",
        ))
    end

    return Products(canonical_targets, canonical_fractions, stoichiometry)
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

# Scientific processes

"""Population growth process with one process-owned rate scale and named multiplicative factors.

`bindings.maximum_rate` names the model parameter that sets the growth-rate scale. `source`
identifies the reference resource removed by growth. Fixed-stoichiometry growth draws any
additional external nutrient resources according to its stoichiometric ratios, while quota
nutrients are transferred independently through [`NutrientUptake`](@ref).
"""
struct Growth{F<:NamedTuple,S,T,U} <: AbstractProcess
    populations::Tuple
    factors::F
    source::S
    stoichiometry::T
    state::U
    bindings::NamedTuple

    function Growth(
        populations::Tuple, factors::NamedTuple, source, stoichiometry, state, bindings::NamedTuple
    )
        canonical = _canonical_factors(factors)
        isnothing(source) || source isa Symbol ||
            throw(ArgumentError("growth `source` must be a Symbol"))
        isnothing(stoichiometry) || stoichiometry isa AbstractStoichiometry || throw(
            ArgumentError("growth `stoichiometry` must be an AbstractStoichiometry"),
        )
        isnothing(state) || state isa Symbol ||
            throw(ArgumentError("growth `state` must be a Symbol"))
        return new{typeof(canonical),typeof(source),typeof(stoichiometry),typeof(state)}(
            populations, canonical, source, stoichiometry, state, _canonical_bindings(bindings)
        )
    end
end

function Growth(;
    populations,
    factors::NamedTuple,
    source=nothing,
    stoichiometry=nothing,
    state=nothing,
    bindings::NamedTuple=NamedTuple(),
)
    population_refs = _canonical_participants(:populations, populations)
    return Growth(population_refs, factors, source, stoichiometry, state, bindings)
end

authored_parameter_bindings(process::Growth) = process.bindings

"""Independent external nutrient uptake into one population inventory state.

The `reference` state scales uptake capacity but is not itself transferred. Parameter
bindings are explicit because quota bounds are commonly shared with `QuotaResponse`.
"""
struct NutrientUptake{F<:QuotaRegulatedMonod} <: AbstractProcess
    formulation::F
    population::Symbol
    target_state::Symbol
    reference_state::Symbol
    resource::Symbol
    bindings::NamedTuple
end

function NutrientUptake(
    formulation::QuotaRegulatedMonod;
    population::Symbol,
    target_state::Symbol,
    reference_state::Symbol,
    resource::Symbol,
    bindings::NamedTuple,
)
    return NutrientUptake(
        formulation, population, target_state, reference_state, resource,
        _canonical_bindings(bindings),
    )
end

authored_parameter_bindings(process::NutrientUptake) = process.bindings

"""Consumer-resource process with optional factors and unassimilated products."""
struct Consumption{F<:Union{PreferentialGrazing,HeterotrophicConsumption},A<:NamedTuple,P} <: AbstractProcess
    formulation::F
    consumers::Tuple
    resources::Tuple
    factors::A
    products::P
    bindings::NamedTuple
end

function Consumption(
    formulation::Union{PreferentialGrazing,HeterotrophicConsumption};
    consumers,
    resources,
    factors::NamedTuple=NamedTuple(),
    unassimilated_products=nothing,
    bindings::NamedTuple=NamedTuple(),
)
    consumer_refs = _canonical_participants(:consumers, consumers)
    resource_refs = _canonical_participants(:resources, resources)
    canonical = _canonical_factors(factors; allow_empty=true)
    products = _canonical_products(unassimilated_products)
    return Consumption(
        formulation, consumer_refs, resource_refs, canonical, products, _canonical_bindings(bindings)
    )
end

authored_parameter_bindings(process::Consumption) = process.bindings

"""Population mortality process with optional products."""
struct Mortality{F<:Union{LinearMortality,QuadraticMortality},P} <: AbstractProcess
    formulation::F
    populations::Tuple
    products::P
    bindings::NamedTuple
end

function Mortality(
    formulation::Union{LinearMortality,QuadraticMortality};
    populations,
    products=nothing,
    bindings::NamedTuple=NamedTuple(),
)
    population_refs = _canonical_participants(:populations, populations)
    product_spec = _canonical_products(products)
    return Mortality(
        formulation, population_refs, product_spec, _canonical_bindings(bindings)
    )
end

authored_parameter_bindings(process::Mortality) = process.bindings

"""Source-to-destination remineralization process."""
struct Remineralization{F<:LinearRemineralization} <: AbstractProcess
    formulation::F
    sources::Tuple
    destination::Symbol
    bindings::NamedTuple
end

function Remineralization(
    formulation::LinearRemineralization;
    sources,
    destination,
    bindings::NamedTuple=NamedTuple(),
)
    source_refs = _canonical_participants(:sources, sources)
    destination isa Symbol || throw(
        ArgumentError("remineralization `destination` must be a Symbol"),
    )
    return Remineralization(
        formulation, source_refs, destination, _canonical_bindings(bindings)
    )
end

authored_parameter_bindings(process::Remineralization) = process.bindings

"""Return the scientific formulation carried by a process or factor."""
formulation(::Growth) = MultiplicativeFactors()
formulation(process::AbstractProcess) = process.formulation
formulation(factor::AbstractFactor) = factor.formulation

"""Return the named multiplicative factors attached to a process."""
factors(::AbstractProcess) = NamedTuple()
factors(process::Union{Growth,Consumption}) = process.factors

process_products(::AbstractProcess) = nothing
process_products(process::Union{Consumption,Mortality}) = process.products
product_path(::Mortality) = (:products,)
product_path(::Consumption) = (:unassimilated_products,)

"""Whether a consumer-resource formulation uses living consumer-prey interaction matrices."""
uses_living_interactions(::AbstractFormulation) = false
uses_living_interactions(::PreferentialGrazing) = true

"""Return canonical participant roles for an authored scientific process."""
function participants(process::Growth)
    base = (population=process.populations,)
    isnothing(process.source) && return base
    return merge(base, (source=(process.source,),))
end
participants(process::NutrientUptake) = (
    population=(process.population,), resource=(process.resource,)
)
participants(process::Consumption) = (consumer=process.consumers, resource=process.resources)
participants(process::Mortality) = (population=process.populations,)
participants(process::Remineralization) =
    (source=process.sources, destination=(process.destination,))

# Model definition

"""Author-facing scientific model definition."""
struct ModelDefinition{C,P,A}
    components::C
    processes::P
    parameters::A
end

function ModelDefinition(; components::NamedTuple, processes::NamedTuple, parameters=nothing)
    return ModelDefinition(components, processes, parameters)
end

function ModelDefinition(family::AbstractModelFamily)
    return ModelDefinition(;
        components=default_components(family),
        processes=default_processes(family),
        parameters=parameter_definitions(family),
    )
end
