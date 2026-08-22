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

"""Stable identity of one semantic formulation parameter requirement.

The identity is scoped by the owning named process, nested factor/formulation path,
formulation tag, slot, and participant qualifier. Applicability axes are resolved from
process participation during setup.
"""
struct ParameterRequirementIdentity{P,Q}
    process::Symbol
    path::P
    formulation::Symbol
    slot::Symbol
    qualifier::Q
end

Base.:(==)(a::ParameterRequirementIdentity, b::ParameterRequirementIdentity) =
    a.process == b.process &&
    a.path == b.path &&
    a.formulation == b.formulation &&
    a.slot == b.slot &&
    a.qualifier == b.qualifier

Base.hash(requirement::ParameterRequirementIdentity, h::UInt) = hash(
    (
        requirement.process,
        requirement.path,
        requirement.formulation,
        requirement.slot,
        requirement.qualifier,
    ),
    h,
)

function ParameterRequirementIdentity(
    process::Symbol,
    path::Tuple,
    formulation_value,
    slot::Symbol;
    qualifier::NamedTuple=NamedTuple(),
)
    all(item -> item isa Symbol, path) || throw(
        ArgumentError("parameter requirement path must contain only Symbols"),
    )
    formulation_name = formulation_value isa Symbol ?
                       formulation_value : formulation_tag(formulation_value)
    return ParameterRequirementIdentity(
        process, path, formulation_name, slot, _canonical_qualifier(qualifier)
    )
end

function _default_requirement_shape(axes::Tuple)
    n_axes = length(axes)
    n_axes == 0 && return :scalar
    n_axes == 1 && return :vector
    n_axes == 2 && return :matrix
    throw(ArgumentError("parameter requirement axes currently support at most two dimensions"))
end

"""Formulation-local declaration of one semantic parameter slot."""
struct ParameterSlot{A<:Tuple}
    name::Symbol
    axes::A
    qualify::Union{Nothing,Symbol}
    shape::Symbol
end

function ParameterSlot(
    name::Symbol,
    axes::Tuple=();
    qualify::Union{Nothing,Symbol}=nothing,
    shape::Symbol=_default_requirement_shape(axes),
)
    length(axes) <= 2 || throw(
        ArgumentError("parameter slot axes currently support at most two dimensions"),
    )
    all(axis -> axis isa Symbol, axes) || throw(
        ArgumentError("parameter slot axes must contain only Symbols"),
    )
    shape in (:scalar, :vector, :matrix) || throw(
        ArgumentError("parameter slot shape must be :scalar, :vector, or :matrix"),
    )
    return ParameterSlot(name, axes, qualify, shape)
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
    ParameterSlot(:rate, (:source,); qualify=:source, shape=:scalar),
)
parameter_slots(::DirectRouting) = ()
parameter_slots(::PartitionRouting) = (ParameterSlot(:export_fraction),)
parameter_slots(::DOMPOMRouting) = (ParameterSlot(:POM_fraction),)
parameter_slots(::FixedStoichiometry) = (ParameterSlot(:ratio; qualify=:currency),)

"""One formulation-declared semantic parameter requirement.

`axes` describe process-local applicability. `shape` describes model-parameter storage,
so a scalar parameter may be shared across an indexed process family.
"""
struct ParameterRequirement{I,A<:Tuple}
    identity::I
    axes::A
    shape::Symbol
end

"""Resolved mapping from one semantic requirement to concrete runtime parameter storage.

`storage_axes=nothing` means requirement-local storage; explicit axes identify the
global runtime storage coordinate system used by the bound parameter.
"""
struct ParameterBinding{R<:ParameterRequirement,A}
    requirement::R
    parameter::Symbol
    storage_axes::A
end

"""Concrete participant applicability of one bound parameter requirement."""
struct ParameterApplicability{B,C,T}
    binding::B
    axis_components::C
    axis_classes::T
end

function _slot_qualifier(slot::ParameterSlot, context::NamedTuple)
    isnothing(slot.qualify) && return NamedTuple()
    hasproperty(context, slot.qualify) || return NamedTuple()
    value = getproperty(context, slot.qualify)
    value isa Symbol || throw(
        ArgumentError("parameter slot qualifier :$(slot.qualify) must identify one Symbol"),
    )
    return NamedTuple{(slot.qualify,)}((value,))
end

function _slot_requirement(
    named::NamedProcess,
    path::Tuple,
    formulation_value,
    slot::ParameterSlot,
    context::NamedTuple,
)
    identity = ParameterRequirementIdentity(
        process_id(named), path, formulation_value, slot.name;
        qualifier=_slot_qualifier(slot, context),
    )
    return ParameterRequirement(identity, slot.axes, slot.shape)
end

function _slot_requirements(
    named::NamedProcess,
    path::Tuple,
    node;
    context::NamedTuple=NamedTuple(),
    formulation_value=node,
)
    return Tuple(
        _slot_requirement(named, path, formulation_value, slot, context)
        for slot in parameter_slots(node)
    )
end

_parameter_node(path::Tuple, node; context=NamedTuple(), formulation_value=node) =
    (; path, node, context, formulation_value)

_validate_factor_formulation(::AbstractFactor) = nothing
_validate_factor_formulation(factor::Light) =
    formulation(factor) isa Union{Smith,Geider} || throw(ArgumentError("invalid light formulation"))
_validate_factor_formulation(factor::NutrientResponse) =
    formulation(factor) isa Monod || throw(ArgumentError("invalid nutrient-response formulation"))
_validate_factor_formulation(factor::Nutrients) =
    formulation(factor) isa Union{Liebig,Frank} || throw(ArgumentError("invalid nutrient formulation"))
_validate_factor_formulation(factor::Temperature) =
    formulation(factor) isa Q10 || throw(ArgumentError("invalid temperature formulation"))

function _factor_parameter_nodes(path::Tuple, factor::AbstractFactor)
    _validate_factor_formulation(factor)
    nodes = (_parameter_node(
        path, formulation(factor); context=factor_parameter_context(factor)
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

function _routing_parameter_nodes(routing::ProductRouting)
    nodes = (_parameter_node((:routing,), routing.formulation),)
    routing.formulation isa DOMPOMRouting || return nodes

    reference = routing.stoichiometry.reference
    for currency in keys(routing.pools.DOM)
        currency === reference && continue
        nodes = (
            nodes...,
            _parameter_node(
                (:routing, :stoichiometry),
                routing.stoichiometry;
                context=(currency=currency,),
            ),
        )
    end
    return nodes
end

_process_parameter_nodes(process::AbstractProcess) =
    (_parameter_node((), formulation(process)),)

function _process_parameter_nodes(process::Mortality)
    context = length(process.populations) == 1 ?
              (population=only(process.populations),) : NamedTuple()
    return (_parameter_node((), formulation(process); context),)
end

_process_parameter_nodes(process::Remineralization) = Tuple(
    _parameter_node((), formulation(process); context=(source=source,))
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
    routing = process_routing(process)
    isnothing(routing) || append!(nodes, _routing_parameter_nodes(routing))
    append!(nodes, _stoichiometry_parameter_nodes(process))
    return Tuple(nodes)
end

"""Return semantic parameter requirements from the named process scientific tree."""
function parameter_requirements(named::NamedProcess)
    requirements = Any[]
    for node in _parameter_nodes(named)
        append!(
            requirements,
            _slot_requirements(
                named,
                node.path,
                node.node;
                context=node.context,
                formulation_value=node.formulation_value,
            ),
        )
    end
    return Tuple(requirements)
end

"""Setup-time normalized scientific model definition.

`parameter_bindings` is the canonical ordered contract; `parameter_lookup` is a transient
setup cache used while lowering processes.
"""
struct NormalizedModelDefinition{C,P,A,D,R,B,L}
    components::C
    processes::P
    parameters::A
    driver_identities::D
    parameter_requirements::R
    parameter_bindings::B
    parameter_lookup::L
end

"""Return the canonical external-driver identities required by a normalized model."""
driver_identities(definition::NormalizedModelDefinition) = definition.driver_identities

"""Return formulation-declared semantic parameter requirements for a normalized model."""
parameter_requirements(definition::NormalizedModelDefinition) = definition.parameter_requirements

"""Return resolved semantic requirement-to-model-parameter bindings."""
parameter_bindings(definition::NormalizedModelDefinition) = definition.parameter_bindings

"""Return the resolved binding for one semantic parameter requirement identity."""
function parameter_binding(
    definition::NormalizedModelDefinition, identity::ParameterRequirementIdentity
)
    return get(definition.parameter_lookup, identity) do
        throw(ArgumentError("no model parameter is bound to requirement $identity"))
    end
end

parameter_binding(definition::NormalizedModelDefinition, requirement::ParameterRequirement) =
    parameter_binding(definition, requirement.identity)

"""Resolve all slots for one scientific node from its formulation slot schema."""
function parameter_slot_bindings(
    definition::NormalizedModelDefinition,
    named::NamedProcess,
    path::Tuple,
    node;
    context::NamedTuple=NamedTuple(),
    formulation_value=node,
)
    requirements = _slot_requirements(named, path, node; context, formulation_value)
    names = Tuple(requirement.identity.slot for requirement in requirements)
    bindings = Tuple(
        begin
            binding = parameter_binding(definition, requirement)
            binding.requirement == requirement || throw(
                ArgumentError(
                    "resolved parameter binding does not match slot schema for $(requirement.identity)",
                ),
            )
            binding
        end for requirement in requirements
    )
    return NamedTuple{names}(bindings)
end

"""Return the model parameter name that supplies `requirement`."""
parameter_name(definition::NormalizedModelDefinition, requirement::ParameterRequirement) =
    parameter_binding(definition, requirement).parameter

parameter_name(
    definition::NormalizedModelDefinition, identity::ParameterRequirementIdentity
) = parameter_binding(definition, identity).parameter

function _factor_component_references(factor::AbstractFactor)
    references = Symbol[
        input.component for input in factor_inputs(factor) if input isa FactorComponent
    ]
    for child in values(factor_children(factor))
        append!(references, _factor_component_references(child))
    end
    return Tuple(references)
end

function _routing_component_references(routing::ProductRouting)
    routing.formulation isa DirectRouting && return (routing.retained,)
    routing.formulation isa PartitionRouting && return (routing.retained, routing.exported)
    routing.formulation isa DOMPOMRouting || throw(
        ArgumentError("unsupported product-routing formulation $(typeof(routing.formulation))"),
    )
    references = Symbol[]
    for pool in values(routing.pools), component in values(pool)
        push!(references, component)
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
    routing = process_routing(process)
    isnothing(routing) || append!(references, _routing_component_references(routing))
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

function _validate_routing_science(routing::ProductRouting, components::NamedTuple)
    routing.formulation isa DOMPOMRouting || return nothing
    for pool in values(routing.pools), (target_currency, component) in pairs(pool)
        _validate_currency_target(components, component, target_currency, "routing target")
    end
    return nothing
end

function _validate_consumption_routing(
    routing::ProductRouting, components::NamedTuple, reference::Symbol
)
    if routing.formulation isa DirectRouting
        getproperty(components, routing.retained) isa Pool || throw(
            ArgumentError("consumption routing destination must be a Pool component"),
        )
        _validate_currency_target(
            components, routing.retained, reference, "consumption destination"
        )
    elseif routing.formulation isa PartitionRouting
        for target in (routing.retained, routing.exported)
            getproperty(components, target) isa Pool || throw(
                ArgumentError("consumption routing targets must be Pool components"),
            )
            _validate_currency_target(
                components, target, reference, "consumption routing target"
            )
        end
    else
        _validate_routing_science(routing, components)
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
    isnothing(process.routing) ||
        _validate_consumption_routing(process.routing, components, reference)
    return nothing
end

function _validate_process_science(process::Mortality, components::NamedTuple)
    formulation(process) isa Union{LinearMortality,QuadraticMortality} || throw(
        ArgumentError("unsupported mortality formulation $(typeof(formulation(process)))"),
    )
    isnothing(process.routing) || _validate_routing_science(process.routing, components)
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

function _declared_parameter_requirements(processes::NamedTuple)
    requirements = Any[]
    for process in values(processes)
        append!(requirements, parameter_requirements(process))
    end
    identities = map(requirement -> requirement.identity, requirements)
    length(unique(identities)) == length(identities) || throw(
        ArgumentError("normalized processes declare duplicate parameter requirement identities"),
    )
    return Tuple(requirements)
end

function _matches_parameter_provision(
    requirement::ParameterRequirement, provision::ParameterProvision
)
    identity = requirement.identity
    identity.process === provision.process || return false
    identity.slot === provision.slot || return false
    isnothing(provision.path) || identity.path == provision.path || return false
    return all(keys(provision.qualifier)) do name
        hasproperty(identity.qualifier, name) &&
            getproperty(identity.qualifier, name) == getproperty(provision.qualifier, name)
    end
end

function _resolve_parameter_requirement(
    provision::ParameterProvision, requirements::Tuple, parameter::Symbol
)
    matches = filter(
        requirement -> _matches_parameter_provision(requirement, provision), requirements
    )
    description = (
        process=provision.process,
        slot=provision.slot,
        qualifier=provision.qualifier,
        path=provision.path,
    )
    isempty(matches) && throw(
        ArgumentError(
            "parameter :$parameter provision $description matches no declared requirement"
        ),
    )
    length(matches) == 1 || throw(
        ArgumentError(
            "parameter :$parameter provision $description is ambiguous; " *
            "matches $(map(r -> r.identity, matches)). Add a qualifier or path.",
        ),
    )
    return only(matches)
end

function _normalize_parameter_bindings(requirements::Tuple, definitions)
    isnothing(definitions) && return (), nothing
    definitions isa Tuple || throw(
        ArgumentError("model parameters must be a tuple of ParameterDefinition values"),
    )
    all(definition -> definition isa ParameterDefinition, definitions) || throw(
        ArgumentError("model parameters must contain only ParameterDefinition values"),
    )

    provided = Dict{ParameterRequirementIdentity,Tuple{Symbol,Union{Nothing,Symbol,NTuple{2,Symbol}}}}()
    resolved_definitions = ParameterDefinition[]

    for definition in definitions
        spec = definition.spec
        matched = Tuple(
            _resolve_parameter_requirement(provision, requirements, spec.name)
            for provision in spec.provides
        )
        required_shapes = unique(Tuple(requirement.shape for requirement in matched))
        shape = spec.shape
        if isnothing(shape)
            isempty(required_shapes) && throw(
                ArgumentError(
                    "parameter :$(spec.name) has no semantic provision or explicit storage axes; declare shape explicitly",
                ),
            )
            length(required_shapes) == 1 || throw(
                ArgumentError(
                    "parameter :$(spec.name) provisions require incompatible shapes $(Tuple(required_shapes))",
                ),
            )
            shape = only(required_shapes)
        end

        for requirement in matched
            shape === requirement.shape || throw(
                ArgumentError(
                    "parameter :$(spec.name) provides $(requirement.identity) with required shape " *
                    "$(requirement.shape), not $shape",
                ),
            )
            identity = requirement.identity
            haskey(provided, identity) && throw(
                ArgumentError(
                    "parameter requirement $identity is provided by both :$(first(provided[identity])) and :$(spec.name)",
                ),
            )
            provided[identity] = (spec.name, spec.axes)
        end

        resolved_spec = ParameterSpec(
            spec.name,
            shape;
            axes=spec.axes,
            provides=spec.provides,
        )
        push!(resolved_definitions, ParameterDefinition(resolved_spec, definition.default))
    end

    missing = filter(requirement -> !haskey(provided, requirement.identity), requirements)
    isempty(missing) || throw(
        ArgumentError(
            "model parameter definitions do not provide requirements $(map(r -> r.identity, missing))",
        ),
    )
    bindings = Tuple(
        ParameterBinding(requirement, provided[requirement.identity]...)
        for requirement in requirements
    )
    return bindings, Tuple(resolved_definitions)
end

"""Normalize process identity, semantic parameter requirements, and model bindings.

Process instances are canonicalized by stable process ID, so declaration order does
not change normalized scientific identity. Component ordering is preserved because it
still participates in concrete tracer realization. Parameter requirements come from the
normalized formulations; model parameter definitions bind their stable names to those
requirements before runtime construction.
"""
function normalize_model(definition::ModelDefinition)
    all(component -> component isa Union{Population,Pool}, values(definition.components)) ||
        throw(ArgumentError("model components must be Population or Pool values"))
    normalized_processes = _canonical_processes(
        definition.processes, definition.components
    )
    requirements = _declared_parameter_requirements(normalized_processes)
    bindings, parameters = _normalize_parameter_bindings(requirements, definition.parameters)
    lookup = Dict{ParameterRequirementIdentity,ParameterBinding}(
        binding.requirement.identity => binding for binding in bindings
    )
    return NormalizedModelDefinition(
        definition.components,
        normalized_processes,
        parameters,
        _canonical_driver_identities(normalized_processes),
        requirements,
        bindings,
        lookup,
    )
end

function _axis_components(
    process::NamedProcess, requirement::ParameterRequirement, axis::Symbol
)
    qualifier = requirement.identity.qualifier
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
state tracers.
"""
function resolve_parameter_applicability(
    definition::NormalizedModelDefinition, layout::ComponentLayout
)
    return map(definition.parameter_bindings) do binding
        requirement = binding.requirement
        process = getproperty(definition.processes, requirement.identity.process)
        axis_components = map(
            axis -> _axis_components(process, requirement, axis), requirement.axes
        )
        axis_classes = map(components -> _axis_classes(layout, components), axis_components)
        ParameterApplicability(binding, axis_components, axis_classes)
    end
end
