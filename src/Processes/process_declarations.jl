# Product routing

"""Conservative allocation of one process product flux among named destinations.

Each product destination may be either one component Symbol or a named element-to-component
mapping. A multi-element mapping without `FixedStoichiometry` routes prognostic elemental
states directly. `FixedStoichiometry` derives multi-element products from a one-element source;
every product then declares the same elements and includes the stoichiometric reference element.

For `N` products, specify `N - 1` named fractions. The omitted product receives the exact
conservative remainder `1 - sum(supplied fractions)`. This removes one redundant routing
degree of freedom while preserving exact closure by construction. A single product requires
no fractions.
"""
struct Products{Destinations,Fractions,Stoichiometry}
    destinations::Destinations
    fractions::Fractions
    stoichiometry::Stoichiometry
end

function _canonical_product_destination(destination)
    destination isa Symbol && return destination
    destination isa NamedTuple || throw(
        ArgumentError("product destinations must be component Symbols or element-to-component mappings"),
    )
    isempty(destination) && throw(ArgumentError("multi-element product destinations cannot be empty"))
    all(component -> component isa Symbol, values(destination)) || throw(
        ArgumentError("multi-element product destinations must map elements to component Symbols"),
    )
    names = sort!(collect(keys(destination)); by=String)
    return NamedTuple{Tuple(names)}(Tuple(getproperty(destination, name) for name in names))
end

function Products(
    destinations::NamedTuple;
    fractions::NamedTuple=NamedTuple(),
    stoichiometry=nothing,
)
    isempty(destinations) && throw(ArgumentError("product destinations cannot be empty"))

    destination_names = sort!(collect(keys(destinations)); by=String)
    canonical_destination_values = Tuple(
        _canonical_product_destination(getproperty(destinations, name)) for name in destination_names
    )
    canonical_destinations = NamedTuple{Tuple(destination_names)}(canonical_destination_values)

    nested = first(canonical_destination_values) isa NamedTuple
    all(destination -> (destination isa NamedTuple) == nested, canonical_destination_values) || throw(
        ArgumentError("product destinations cannot mix scalar and multi-element destinations"),
    )
    if nested
        elements = keys(first(canonical_destination_values))
        all(destination -> keys(destination) == elements, canonical_destination_values) || throw(
            ArgumentError("multi-element products must declare the same elements"),
        )
        if !isnothing(stoichiometry)
            stoichiometry isa FixedStoichiometry || throw(
                ArgumentError("product stoichiometry must be FixedStoichiometry"),
            )
            stoichiometry.reference_element in elements || throw(
                ArgumentError(
                    "multi-element products must include stoichiometric reference element :$(stoichiometry.reference_element)"
                ),
            )
        end
    elseif !isnothing(stoichiometry)
        throw(ArgumentError("scalar product destinations do not take stoichiometry"))
    end

    fraction_names = sort!(collect(keys(fractions)); by=String)
    all(name -> name in keys(canonical_destinations), fraction_names) || throw(
        ArgumentError("product `fractions` must be keyed by product names declared in destinations"),
    )
    all(value -> value isa Symbol, values(fractions)) || throw(
        ArgumentError("product `fractions` values must be parameter Symbols"),
    )
    canonical_fractions = NamedTuple{Tuple(fraction_names)}(
        Tuple(getproperty(fractions, name) for name in fraction_names)
    )

    nproducts = length(canonical_destinations)
    nfractions = length(canonical_fractions)
    if nproducts == 1
        nfractions == 0 || throw(
            ArgumentError("single-product allocations do not take fractions"),
        )
    elseif nfractions != nproducts - 1
        throw(ArgumentError(
            "products has $nproducts destinations but $nfractions fractions; specify exactly " *
            "$(nproducts - 1) fractions so the omitted product receives the conservative balance",
        ))
    end

    return Products(canonical_destinations, canonical_fractions, stoichiometry)
end

Products(destination::Symbol) = Products((product=destination,))

authored_parameter_bindings(products::Products) = isempty(products.fractions) ?
    NamedTuple() : (fraction=products.fractions,)

_canonical_products(::Nothing) = nothing
_canonical_products(products::Products) = products
_canonical_products(destination::Symbol) = Products(destination)

function _canonical_factors(factors::NamedTuple)
    all(factor -> factor isa AbstractFactor, values(factors)) || throw(
        ArgumentError("process `factors` values must be process factors"),
    )
    return _canonical_namedtuple(factors)
end

# Scientific processes

"""Plankton growth process with explicit material inputs and optional multiplicative factors.

`bindings.maximum_rate` names the model parameter that sets the growth-rate scale.
`reference_resource` supplies the Element represented by the plankton `reference_state`.
`additional_resources` maps additional Elements to external Pools consumed according to
`FixedStoichiometry`. Factors modify growth rate only; independently prognostic elemental
states are supplied through [`NutrientUptake`](@ref).
"""
struct Growth{Factors<:NamedTuple,AdditionalResources<:NamedTuple,Stoichiometry} <: AbstractProcess
    plankton::Tuple
    factors::Factors
    reference_resource::Symbol
    additional_resources::AdditionalResources
    stoichiometry::Stoichiometry
    bindings::NamedTuple
end

function Growth(;
    plankton,
    reference_resource::Symbol,
    factors::NamedTuple=NamedTuple(),
    additional_resources::NamedTuple=NamedTuple(),
    stoichiometry=nothing,
    bindings::NamedTuple=NamedTuple(),
)
    all(resource -> resource isa Symbol, values(additional_resources)) || throw(
        ArgumentError("growth `additional_resources` values must be Pool Symbols"),
    )
    isnothing(stoichiometry) || stoichiometry isa FixedStoichiometry || throw(
        ArgumentError("growth `stoichiometry` must be FixedStoichiometry"),
    )
    return Growth(
        _canonical_participants(:plankton, plankton),
        _canonical_factors(factors),
        reference_resource,
        _canonical_namedtuple(additional_resources),
        stoichiometry,
        _canonical_bindings(bindings),
    )
end

authored_parameter_bindings(process::Growth) = process.bindings

"""Return whether `factor` is scientifically applicable to `process`."""
factor_applicable(::AbstractProcess, ::AbstractFactor) = true
factor_applicable(::AbstractProcess, ::Union{Light,QuotaResponse}) = false
factor_applicable(::Growth, ::Union{Light,QuotaResponse}) = true

"""Independent external nutrient uptake into one plankton inventory state.

The plankton reference state scales uptake capacity but is not itself transferred. Parameter
bindings are explicit because quota bounds are commonly shared with `QuotaResponse`.
"""
struct NutrientUptake{Formulation<:QuotaRegulatedMonod} <: AbstractProcess
    formulation::Formulation
    plankton::Symbol
    target_state::Symbol
    resource::Symbol
    bindings::NamedTuple
end

function NutrientUptake(
    formulation::QuotaRegulatedMonod;
    plankton::Symbol,
    target_state::Symbol,
    resource::Symbol,
    bindings::NamedTuple,
)
    return NutrientUptake(
        formulation, plankton, target_state, resource, _canonical_bindings(bindings),
    )
end

authored_parameter_bindings(process::NutrientUptake) = process.bindings

"""Consumer-resource process with optional factors and unassimilated products.

The formulation `maximum_rate` is pairwise capacity: one consumer's rate is applied
independently to each declared consumer-resource edge. When one living-prey consumption process
routes multi-element unassimilated products from multiple resources, those resources currently
must expose the same prognostic Element set.
"""
struct Consumption{
    Formulation<:Union{PreferentialGrazing,HeterotrophicConsumption},
    Factors<:NamedTuple,
    ProductRouting,
} <: AbstractProcess
    formulation::Formulation
    consumers::Tuple
    resources::Tuple
    factors::Factors
    products::ProductRouting
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
    canonical = _canonical_factors(factors)
    products = _canonical_products(unassimilated_products)
    return Consumption(
        formulation, consumer_refs, resource_refs, canonical, products, _canonical_bindings(bindings)
    )
end

authored_parameter_bindings(process::Consumption) = process.bindings

"""Plankton mortality process with optional products.

When one mortality process routes multi-element products from multiple plankton, those
plankton currently must expose the same prognostic Element set.
"""
struct Mortality{Formulation<:Union{LinearMortality,QuadraticMortality},ProductRouting} <: AbstractProcess
    formulation::Formulation
    plankton::Tuple
    products::ProductRouting
    bindings::NamedTuple
end

function Mortality(
    formulation::Union{LinearMortality,QuadraticMortality};
    plankton,
    products=nothing,
    bindings::NamedTuple=NamedTuple(),
)
    plankton_refs = _canonical_participants(:plankton, plankton)
    product_spec = _canonical_products(products)
    return Mortality(
        formulation, plankton_refs, product_spec, _canonical_bindings(bindings)
    )
end

authored_parameter_bindings(process::Mortality) = process.bindings

"""Source-to-destination remineralization process."""
struct Remineralization{Formulation<:LinearRemineralization} <: AbstractProcess
    formulation::Formulation
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
formulation(::Growth) = FactorizedGrowth()
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
    resources = (process.reference_resource, Tuple(values(process.additional_resources))...)
    return (plankton=process.plankton, resource=resources)
end
participants(process::NutrientUptake) = (
    plankton=(process.plankton,), resource=(process.resource,)
)
participants(process::Consumption) = (consumer=process.consumers, resource=process.resources)
participants(process::Mortality) = (plankton=process.plankton,)
participants(process::Remineralization) =
    (source=process.sources, destination=(process.destination,))

# Model definition

"""Author-facing scientific model definition."""
struct ModelDefinition{ComponentMap,ProcessMap,ParameterDefinitions}
    components::ComponentMap
    processes::ProcessMap
    parameters::ParameterDefinitions
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
