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

"""Named, validated process instance in canonical normalized model state."""
struct NamedProcess{P<:AbstractProcess}
    id::Symbol
    process::P
end

process_id(process::NamedProcess) = process.id
process_kind(process::NamedProcess) = process_kind(process.process)
formulation(process::NamedProcess) = formulation(process.process)
factors(process::NamedProcess) = factors(process.process)
participants(process::NamedProcess) = participants(process.process)
drivers(process::NamedProcess) = drivers(process.process)
rate_axes(process::NamedProcess) = rate_axes(process.process)

function _canonical_qualifier(qualifier::NamedTuple)
    names = sort!(collect(keys(qualifier)); by=String)
    names_tuple = Tuple(names)
    values = Tuple(getproperty(qualifier, name) for name in names)
    return NamedTuple{names_tuple}(values)
end

"""Formulation-local declaration of one semantic parameter slot.

Dimensionality is structural: zero, one, or two declared axes imply scalar, vector, or
matrix values respectively. `qualify` identifies repeated semantic instances without
changing storage dimensionality. For a scalar slot, a qualifier that is also a process
participant role can still provide ecological applicability without becoming a storage axis.
"""
struct ParameterSlot{A<:Tuple}
    name::Symbol
    axes::A
    qualify::Union{Nothing,Symbol}
end

function ParameterSlot(
    name::Symbol,
    axes::Tuple=();
    qualify::Union{Nothing,Symbol}=nothing,
)
    length(axes) <= 2 || throw(
        ArgumentError("parameter slot axes currently support at most two dimensions"),
    )
    all(axis -> axis isa Symbol, axes) || throw(
        ArgumentError("parameter slot axes must contain only Symbols"),
    )
    return ParameterSlot(name, axes, qualify)
end

parameter_slots(::AbstractFormulation) = ()
parameter_slots(::Smith) = (
    ParameterSlot(:maximum_rate, (:population,)),
    ParameterSlot(:alpha, (:population,)),
)
parameter_slots(::Geider) = (
    ParameterSlot(:maximum_rate, (:population,)),
    ParameterSlot(:alpha, (:population,)),
    ParameterSlot(:chlorophyll_to_carbon_ratio, (:population,)),
)
parameter_slots(::Monod) = (ParameterSlot(:K, (:population,); qualify=:resource),)
parameter_slots(::Liebig) = ()
parameter_slots(::FrankTNorm) = (ParameterSlot(:sharpness),)
parameter_slots(::Q10) = (
    ParameterSlot(:q10),
    ParameterSlot(:reference_temperature),
)
parameter_slots(::IdealizedGrazing) = (
    ParameterSlot(:maximum_rate, (:consumer,)),
    ParameterSlot(:half_saturation, (:consumer,)),
    ParameterSlot(:palatability, (:consumer, :resource)),
    ParameterSlot(:assimilation, (:consumer, :resource)),
)
parameter_slots(::PreferentialGrazing) = (
    ParameterSlot(:maximum_rate, (:consumer,)),
    ParameterSlot(:half_saturation, (:consumer,)),
    ParameterSlot(:palatability, (:consumer, :resource)),
    ParameterSlot(:assimilation, (:consumer, :resource)),
)
parameter_slots(::HeterotrophicConsumption) = (
    ParameterSlot(:maximum_rate, (:consumer,)),
    ParameterSlot(:half_saturation, (:resource,)),
    ParameterSlot(:assimilation, (:consumer, :resource)),
)
parameter_slots(::LinearMortality) = (ParameterSlot(:rate, (:population,); qualify=:population),)
parameter_slots(::QuadraticMortality) = (ParameterSlot(:rate, (:population,); qualify=:population),)
parameter_slots(::LinearRemineralization) = (
    ParameterSlot(:rate; qualify=:source),
)
parameter_slots(products::Products) = length(products.targets) == 1 ? () :
    (ParameterSlot(:fraction; qualify=:product),)
parameter_slots(::FixedStoichiometry) = (ParameterSlot(:ratio; qualify=:currency),)

"""Resolved setup-time mapping from one local parameter slot to model storage.

The scientific-tree location and qualifier identify the local slot during compilation.
`storage_axes=nothing` means slot-local storage; explicit axes identify the global runtime
storage coordinate system used by the bound model parameter.
"""
struct ParameterBinding{P,Q,A,S}
    process::Symbol
    path::P
    slot::Symbol
    qualifier::Q
    axes::A
    parameter::Symbol
    storage_axes::S
end

"""Concrete participant applicability of one resolved parameter binding."""
struct ParameterApplicability{B,C,T}
    binding::B
    axis_components::C
    axis_classes::T
end

_parameter_slot_source(node::Union{AbstractFormulation,AbstractStoichiometry,Products}) = node
_parameter_slot_source(node) = formulation(node)

function _slot_qualifier(slot::ParameterSlot, context::NamedTuple)
    isnothing(slot.qualify) && return NamedTuple()
    hasproperty(context, slot.qualify) || return NamedTuple()
    value = getproperty(context, slot.qualify)
    value isa Symbol || throw(
        ArgumentError("parameter slot qualifier :$(slot.qualify) must identify one Symbol"),
    )
    return NamedTuple{(slot.qualify,)}((value,))
end

function _slot_metadata(
    named::NamedProcess,
    path::Tuple,
    slot::ParameterSlot,
    context::NamedTuple,
)
    all(item -> item isa Symbol, path) || throw(
        ArgumentError("parameter binding path must contain only Symbols"),
    )
    return (;
        process=process_id(named),
        path,
        slot=slot.name,
        qualifier=_canonical_qualifier(_slot_qualifier(slot, context)),
        axes=slot.axes,
    )
end

_binding_key(process::Symbol, path::Tuple, slot::Symbol, qualifier::NamedTuple) =
    (process, path, slot, _canonical_qualifier(qualifier))
_binding_key(binding::ParameterBinding) =
    _binding_key(binding.process, binding.path, binding.slot, binding.qualifier)
_binding_key(metadata::NamedTuple) =
    _binding_key(metadata.process, metadata.path, metadata.slot, metadata.qualifier)

_parameter_node(path::Tuple, node; context=NamedTuple()) =
    (; path, node, context)

_validate_factor_formulation(::AbstractFactor) = nothing
_validate_factor_formulation(factor::Light) =
    formulation(factor) isa Union{Smith,Geider} || throw(ArgumentError("invalid light formulation"))
_validate_factor_formulation(factor::NutrientResponse) =
    formulation(factor) isa Monod || throw(ArgumentError("invalid nutrient-response formulation"))
_validate_factor_formulation(factor::Nutrients) =
    formulation(factor) isa Union{Liebig,FrankTNorm} || throw(ArgumentError("invalid nutrient formulation"))
_validate_factor_formulation(factor::Temperature) =
    formulation(factor) isa Q10 || throw(ArgumentError("invalid temperature formulation"))

function _factor_parameter_nodes(path::Tuple, factor::AbstractFactor)
    _validate_factor_formulation(factor)
    nodes = (_parameter_node(
        path, factor; context=factor_parameter_context(factor),
    ),)
    for (name, child) in pairs(factor_children(factor))
        child isa AbstractFactor || throw(ArgumentError("factor children must be process factors"))
        factor isa Nutrients && !(child isa NutrientResponse) && throw(
            ArgumentError("nutrient responses must be NutrientResponse factors"),
        )
        nodes = (
            nodes...,
            _factor_parameter_nodes(factor_child_path(path, factor, name), child)...,
        )
    end
    return nodes
end

_factor_parameter_nodes(name::Symbol, factor::AbstractFactor) =
    _factor_parameter_nodes((:factors, name), factor)

function _products_parameter_nodes(path::Tuple, products::Products)
    nodes = Any[
        _parameter_node(path, products; context=(product=name,))
        for name in keys(products.fractions)
    ]
    stoichiometry = products.stoichiometry
    if !isnothing(stoichiometry)
        currencies = keys(first(values(products.targets)))
        for currency in currencies
            currency === stoichiometry.reference && continue
            push!(
                nodes,
                _parameter_node(
                    (path..., :stoichiometry),
                    stoichiometry;
                    context=(currency=currency,),
                ),
            )
        end
    end
    return Tuple(nodes)
end

_process_parameter_nodes(process::AbstractProcess) =
    (_parameter_node((), process),)

_process_parameter_nodes(process::Mortality) = Tuple(
    _parameter_node(
        (), process; context=(population=population,)
    )
    for population in process.populations
)

_process_parameter_nodes(process::Remineralization) = Tuple(
    _parameter_node(
        (), process; context=(source=source,)
    )
    for source in process.sources
)

_stoichiometry_parameter_nodes(::AbstractProcess) = ()
function _stoichiometry_parameter_nodes(process::Growth)
    stoichiometry = process_stoichiometry(process)
    isnothing(stoichiometry) && return ()
    nutrients = only(factor for factor in values(process.factors) if factor isa Nutrients)
    return Tuple(
        _parameter_node(
            (:stoichiometry,), stoichiometry; context=(currency=currency,)
        )
        for currency in keys(nutrients.responses)
    )
end

function _parameter_nodes(named::NamedProcess)
    process = named.process
    nodes = Any[_process_parameter_nodes(process)...]
    for (name, factor) in pairs(factors(process))
        append!(nodes, _factor_parameter_nodes(name, factor))
    end
    products = process_products(process)
    isnothing(products) || append!(nodes, _products_parameter_nodes(product_path(process), products))
    append!(nodes, _stoichiometry_parameter_nodes(process))
    return Tuple(nodes)
end

"""Setup-time normalized scientific model definition.

`parameter_bindings` is the canonical ordered contract; `parameter_lookup` is a transient
setup cache used while lowering processes.
"""
struct NormalizedModelDefinition{C,P,A,D,B,L}
    components::C
    processes::P
    parameters::A
    driver_identities::D
    parameter_bindings::B
    parameter_lookup::L
end

"""Return the canonical external-driver identities required by a normalized model."""
driver_identities(definition::NormalizedModelDefinition) = definition.driver_identities

"""Return resolved local-slot-to-model-parameter bindings."""
parameter_bindings(definition::NormalizedModelDefinition) = definition.parameter_bindings

function _parameter_binding(
    definition::NormalizedModelDefinition,
    process::Symbol,
    path::Tuple,
    slot::Symbol,
    qualifier::NamedTuple,
)
    key = _binding_key(process, path, slot, qualifier)
    return get(definition.parameter_lookup, key) do
        throw(ArgumentError(
            "no model parameter is bound to slot :$slot at process :$process path $path qualifier $qualifier",
        ))
    end
end

"""Resolve all parameter slots for one scientific node from its formulation schema."""
function parameter_slot_bindings(
    definition::NormalizedModelDefinition,
    named::NamedProcess,
    path::Tuple,
    node;
    context::NamedTuple=NamedTuple(),
)
    slot_source = _parameter_slot_source(node)
    slots = parameter_slots(slot_source)
    names = Tuple(slot.name for slot in slots)
    bindings = Tuple(
        begin
            qualifier = _slot_qualifier(slot, context)
            binding = _parameter_binding(
                definition, process_id(named), path, slot.name, qualifier
            )
            binding.axes == slot.axes || throw(
                ArgumentError(
                    "resolved parameter binding does not match slot schema for process " *
                    ":$(process_id(named)) path $path slot :$(slot.name) qualifier $qualifier",
                ),
            )
            binding
        end for slot in slots
    )
    return NamedTuple{names}(bindings)
end

function _factor_component_references(factor::AbstractFactor)
    references = Symbol[
        input.component for input in factor_inputs(factor) if input isa FactorComponent
    ]
    for child in values(factor_children(factor))
        append!(references, _factor_component_references(child))
    end
    return Tuple(references)
end

function _product_component_references(products::Products)
    references = Symbol[]
    for target in values(products.targets)
        if target isa Symbol
            push!(references, target)
        else
            append!(references, values(target))
        end
    end
    return Tuple(references)
end

function _process_component_references(process::AbstractProcess)
    references = Symbol[]
    for values_for_role in values(participants(process)), component in values_for_role
        push!(references, component)
    end
    for factor in values(factors(process))
        append!(references, _factor_component_references(factor))
    end
    products = process_products(process)
    isnothing(products) || append!(references, _product_component_references(products))
    return Tuple(references)
end

function _validate_currency_target(
    components::NamedTuple, component::Symbol, expected::Symbol, label::String
)
    actual = currency(getproperty(components, component))
    actual === expected || throw(
        ArgumentError("$label component :$component has currency :$actual, expected :$expected"),
    )
    return nothing
end

_validate_process_science(::AbstractProcess, ::NamedTuple) = nothing

function _validate_process_science(process::Growth, components::NamedTuple)
    single_resource = any(factor -> factor isa NutrientResponse, values(process.factors))
    !single_resource || isnothing(process.source) || throw(
        ArgumentError(
            "single-resource growth derives its source from the nutrient response; omit `source`"
        ),
    )

    stoichiometry = process_stoichiometry(process)
    isnothing(stoichiometry) && return nothing
    stoichiometry isa FixedStoichiometry || throw(
        ArgumentError("unsupported growth stoichiometry $(typeof(stoichiometry))"),
    )
    nutrient_factors = Tuple(
        factor for factor in values(process.factors) if factor isa Nutrients
    )
    length(nutrient_factors) == 1 || throw(
        ArgumentError(
            "fixed-stoichiometry growth requires exactly one multi-resource Nutrients factor"
        ),
    )
    reference = stoichiometry.reference
    isnothing(process.source) && throw(
        ArgumentError("fixed-stoichiometry growth requires a reference-currency source component"),
    )
    _validate_currency_target(components, process.source, reference, "growth source")
    for population in process.populations
        _validate_currency_target(components, population, reference, "growth population")
    end
    for factor in values(process.factors)
        factor isa Nutrients || continue
        for (target_currency, response) in pairs(factor.responses)
            _validate_currency_target(
                components, response.resource, target_currency, "nutrient response"
            )
        end
    end
    return nothing
end

function _validate_products_science(
    products::Products, components::NamedTuple, reference::Symbol
)
    stoichiometry = products.stoichiometry
    if isnothing(stoichiometry)
        for target in values(products.targets)
            target isa Symbol || throw(ArgumentError("scalar products require component targets"))
            getproperty(components, target) isa Pool || throw(
                ArgumentError("product targets must be Pool components"),
            )
            _validate_currency_target(components, target, reference, "product target")
        end
        return nothing
    end

    stoichiometry.reference === reference || throw(ArgumentError(
        "product stoichiometric reference :$(stoichiometry.reference) does not match process currency :$reference",
    ))
    for target_group in values(products.targets), (target_currency, target) in pairs(target_group)
        getproperty(components, target) isa Pool || throw(
            ArgumentError("product targets must be Pool components"),
        )
        _validate_currency_target(components, target, target_currency, "product target")
    end
    return nothing
end

function _validate_process_science(process::Consumption, components::NamedTuple)
    formulation(process) isa Union{
        IdealizedGrazing,PreferentialGrazing,HeterotrophicConsumption
    } || throw(
        ArgumentError("unsupported consumption formulation $(typeof(formulation(process)))"),
    )
    all(consumer -> getproperty(components, consumer) isa Population, process.consumers) || throw(
        ArgumentError("consumption consumers must be Population components"),
    )
    living_resources = uses_living_interactions(process.formulation)
    resource_type = living_resources ? Population : Pool
    valid_resources = all(
        resource -> getproperty(components, resource) isa resource_type, process.resources
    )
    resource_error = living_resources ?
                     "living-interaction resources must be Population components" :
                     "heterotrophic-consumption resources must be Pool components"
    valid_resources || throw(ArgumentError(resource_error))

    reference = currency(getproperty(components, first(process.consumers)))
    for consumer in process.consumers
        _validate_currency_target(components, consumer, reference, "consumption consumer")
    end
    for resource in process.resources
        _validate_currency_target(components, resource, reference, "consumption resource")
    end
    isnothing(process.products) || _validate_products_science(process.products, components, reference)
    return nothing
end

function _validate_process_science(process::Mortality, components::NamedTuple)
    formulation(process) isa Union{LinearMortality,QuadraticMortality} || throw(
        ArgumentError("unsupported mortality formulation $(typeof(formulation(process)))"),
    )
    if !isnothing(process.products)
        references = unique(Tuple(
            currency(getproperty(components, population)) for population in process.populations
        ))
        length(references) == 1 || throw(
            ArgumentError("mortality populations must share one currency when products are declared"),
        )
        _validate_products_science(process.products, components, only(references))
    end
    return nothing
end

function _validate_process_science(process::Remineralization, ::NamedTuple)
    formulation(process) isa LinearRemineralization || throw(
        ArgumentError("unsupported remineralization formulation $(typeof(formulation(process)))"),
    )
    return nothing
end

function _validate_process(id::Symbol, process, components::NamedTuple)
    process isa AbstractProcess || throw(
        ArgumentError("process :$id must be an AbstractProcess; got $(typeof(process))"),
    )
    component_names = keys(components)
    missing = filter(
        reference -> !(reference in component_names), _process_component_references(process)
    )
    isempty(missing) || throw(
        ArgumentError("process :$id references unknown components $missing"),
    )
    _validate_process_science(process, components)
    return NamedProcess(id, process)
end

function _canonical_processes(processes::NamedTuple, components::NamedTuple)
    names = sort!(collect(keys(processes)); by=String)
    names_tuple = Tuple(names)
    values = Tuple(
        _validate_process(name, getproperty(processes, name), components) for name in names
    )
    return NamedTuple{names_tuple}(values)
end

function _canonical_driver_identities(processes::NamedTuple)
    identities = Symbol[]
    for process in values(processes), identity in values(drivers(process))
        identity in identities || push!(identities, identity)
    end
    sort!(identities; by=String)
    return Tuple(identities)
end

function _validate_binding_slot_names(node)
    slots = parameter_slots(_parameter_slot_source(node))
    bindings = if applicable(authored_parameter_bindings, node)
        authored_parameter_bindings(node)
    else
        isempty(slots) || throw(ArgumentError(
            "slot-owning node $(typeof(node)) must implement `authored_parameter_bindings`",
        ))
        NamedTuple()
    end
    slot_names = Tuple(slot.name for slot in slots)
    unknown = Tuple(name for name in keys(bindings) if !(name in slot_names))
    isempty(unknown) || throw(ArgumentError(
        "inline bindings contain unknown slots $unknown; expected $slot_names",
    ))
    return bindings
end

function _binding_value(bindings::NamedTuple, slot::ParameterSlot, qualifier::NamedTuple)
    explicit = hasproperty(bindings, slot.name)
    value = explicit ? getproperty(bindings, slot.name) : slot.name
    value isa Symbol && return value, explicit

    value isa NamedTuple || throw(ArgumentError(
        "binding :$(slot.name) must be a parameter Symbol or one-level qualifier NamedTuple",
    ))
    isnothing(slot.qualify) && throw(ArgumentError(
        "binding :$(slot.name) uses a qualifier map but the slot is unqualified",
    ))
    hasproperty(qualifier, slot.qualify) || throw(ArgumentError(
        "binding :$(slot.name) requires qualifier :$(slot.qualify)",
    ))
    qualifier_value = getproperty(qualifier, slot.qualify)
    hasproperty(value, qualifier_value) || throw(ArgumentError(
        "binding :$(slot.name) has no entry for qualifier :$qualifier_value",
    ))
    parameter = getproperty(value, qualifier_value)
    parameter isa Symbol || throw(ArgumentError(
        "binding :$(slot.name) qualifier :$qualifier_value must map to a parameter Symbol",
    ))
    return parameter, true
end

function _declared_parameter_uses(processes::NamedTuple)
    uses = Any[]
    for named in values(processes), node in _parameter_nodes(named)
        bindings = _validate_binding_slot_names(node.node)
        slot_source = _parameter_slot_source(node.node)
        for slot in parameter_slots(slot_source)
            qualifier = _slot_qualifier(slot, node.context)
            metadata = _slot_metadata(named, node.path, slot, node.context)
            parameter, explicit = _binding_value(bindings, slot, qualifier)
            push!(uses, (; metadata..., parameter, explicit))
        end
    end
    keys = Tuple(_binding_key(use) for use in uses)
    length(unique(keys)) == length(keys) || throw(
        ArgumentError("normalized processes declare duplicate parameter binding keys"),
    )
    return Tuple(uses)
end

function _storage_rank(axes)
    axes === nothing && return nothing
    axes isa Symbol && return 1
    return length(axes)
end

function _normalize_parameter_bindings(processes::NamedTuple, definitions)
    uses = _declared_parameter_uses(processes)
    if isnothing(definitions)
        bindings = Tuple(
            ParameterBinding(
                use.process, use.path, use.slot, use.qualifier,
                use.axes, use.parameter, nothing,
            )
            for use in uses
        )
        return bindings, nothing
    end
    definitions isa NamedTuple || throw(
        ArgumentError("model parameters must be a NamedTuple of Parameter values"),
    )
    all(parameter -> parameter isa Parameter, values(definitions)) || throw(
        ArgumentError("model parameters must contain only Parameter values"),
    )

    definition_names = Set(keys(definitions))
    dependency_names = Set{Symbol}()
    for (name, parameter) in pairs(definitions)
        default = parameter.default
        default isa DerivedDefault || continue
        for dependency in default.deps
            dependency in definition_names || throw(ArgumentError(
                "parameter :$name default depends on undeclared parameter :$dependency",
            ))
            push!(dependency_names, dependency)
        end
    end

    for use in uses
        use.parameter in definition_names || throw(ArgumentError(
            "inline binding for process :$(use.process) path $(use.path) slot :$(use.slot) " *
            "qualifier $(use.qualifier) names undeclared parameter :$(use.parameter)",
        ))
    end

    by_parameter = Dict{Symbol,Vector{Any}}()
    for use in uses
        push!(get!(by_parameter, use.parameter, Any[]), use)
    end
    for (parameter, parameter_uses) in by_parameter
        if length(parameter_uses) > 1 && any(use -> !use.explicit, parameter_uses)
            locations = Tuple(
                (process=use.process, path=use.path, slot=use.slot, qualifier=use.qualifier)
                for use in parameter_uses
            )
            throw(ArgumentError(
                "parameter :$parameter is implicitly bound by multiple parameter slots $locations; " *
                "bind every shared use explicitly",
            ))
        end
        required_ranks = unique(Tuple(length(use.axes) for use in parameter_uses))
        length(required_ranks) == 1 || throw(ArgumentError(
            "parameter :$parameter is bound to slots with incompatible dimensionality $(Tuple(required_ranks))",
        ))
    end

    for (name, parameter) in pairs(definitions)
        parameter_uses = get(by_parameter, name, Any[])
        isempty(parameter_uses) && !(name in dependency_names) && throw(ArgumentError(
            "parameter :$name is neither bound to a process slot nor used by a derived default",
        ))
        isempty(parameter_uses) && continue
        required_rank = only(unique(Tuple(length(use.axes) for use in parameter_uses)))
        storage_rank = _storage_rank(parameter.spec.axes)
        isnothing(storage_rank) || storage_rank == required_rank || throw(ArgumentError(
            "parameter :$name storage axes imply rank $storage_rank but its bound slots require rank $required_rank",
        ))
    end

    bindings = Tuple(
        ParameterBinding(
            use.process,
            use.path,
            use.slot,
            use.qualifier,
            use.axes,
            use.parameter,
            getproperty(definitions, use.parameter).spec.axes,
        )
        for use in uses
    )
    return bindings, definitions
end

"""Normalize process identity and resolve inline parameter bindings.

Process instances are canonicalized by stable process ID, so declaration order does
not change normalized scientific identity. Component ordering is preserved because it
still participates in concrete tracer realization. Local formulation slots bind directly
to stable model parameter names during normalization.
"""
function normalize_model(definition::ModelDefinition)
    all(component -> component isa Union{Population,Pool}, values(definition.components)) ||
        throw(ArgumentError("model components must be Population or Pool values"))
    normalized_processes = _canonical_processes(
        definition.processes, definition.components
    )
    bindings, parameters = _normalize_parameter_bindings(
        normalized_processes, definition.parameters
    )
    lookup = Dict(_binding_key(binding) => binding for binding in bindings)
    return NormalizedModelDefinition(
        definition.components,
        normalized_processes,
        parameters,
        _canonical_driver_identities(normalized_processes),
        bindings,
        lookup,
    )
end

function _axis_components(
    process::NamedProcess, binding::ParameterBinding, axis::Symbol
)
    qualifier = binding.qualifier
    if hasproperty(qualifier, axis)
        value = getproperty(qualifier, axis)
        value isa Symbol || throw(
            ArgumentError("parameter qualifier :$axis must identify one component Symbol"),
        )
        return (value,)
    end
    process_participants = participants(process)
    hasproperty(process_participants, axis) || throw(
        ArgumentError(
            "parameter applicability axis :$axis is not a participant role of process :$(process_id(process))",
        ),
    )
    return getproperty(process_participants, axis)
end

function _axis_classes(layout::ComponentLayout, components::Tuple)
    classes = Symbol[]
    for component in components
        hasproperty(layout.component_classes, component) || throw(
            ArgumentError("parameter applicability references unrealized component :$component"),
        )
        append!(classes, component_classes(layout, component))
    end
    return Tuple(classes)
end

"""Resolve each semantic parameter axis onto ecological component class identities.

The result is setup-time applicability metadata. Population state multiplicity does not
change parameter-axis length: applicability follows ecological classes rather than physical
state tracers. Scalar qualified slots retain participant applicability when their qualifier
names a process participant role; non-participant qualifiers such as currency do not create
ecological axes.
"""
function _applicability_axes(process::NamedProcess, binding::ParameterBinding)
    isempty(binding.axes) || return binding.axes
    process_participants = participants(process)
    return Tuple(
        axis for axis in keys(binding.qualifier) if hasproperty(process_participants, axis)
    )
end

function resolve_parameter_applicability(
    definition::NormalizedModelDefinition, layout::ComponentLayout
)
    return map(definition.parameter_bindings) do binding
        process = getproperty(definition.processes, binding.process)
        applicability_axes = _applicability_axes(process, binding)
        axis_components = map(
            axis -> _axis_components(process, binding, axis), applicability_axes
        )
        axis_classes = map(components -> _axis_classes(layout, components), axis_components)
        ParameterApplicability(binding, axis_components, axis_classes)
    end
end
