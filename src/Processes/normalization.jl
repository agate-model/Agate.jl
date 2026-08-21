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
parameter_slots(::Val{:allometric}) = (
    ParameterSlot(:optimum_predator_prey_ratio, (:consumer,)),
    ParameterSlot(:specificity, (:consumer,)),
    ParameterSlot(:protection, (:resource,)),
)
parameter_slots(::Val{:binary}) = (ParameterSlot(:assimilation_efficiency, (:consumer,)),)

"""One formulation-declared semantic parameter requirement.

`axes` describe process-local applicability. `shape` describes model-parameter storage,
so a scalar parameter may be shared across an indexed process family.
"""
struct ParameterRequirement{I,A<:Tuple}
    identity::I
    axes::A
    shape::Symbol
end

function ParameterRequirement(
    process::Symbol,
    path::Tuple,
    formulation_value,
    slot::Symbol,
    axes::Tuple;
    qualifier::NamedTuple=NamedTuple(),
    shape::Symbol=_default_requirement_shape(axes),
)
    length(axes) <= 2 || throw(
        ArgumentError("parameter requirement axes currently support at most two dimensions"),
    )
    all(axis -> axis isa Symbol, axes) || throw(
        ArgumentError("parameter requirement axes must contain only Symbols"),
    )
    shape in (:scalar, :vector, :matrix) || throw(
        ArgumentError("parameter requirement shape must be :scalar, :vector, or :matrix"),
    )
    identity = ParameterRequirementIdentity(
        process, path, formulation_value, slot; qualifier
    )
    return ParameterRequirement(identity, axes, shape)
end

"""Resolved mapping from one semantic requirement to concrete runtime parameter storage.

`storage_axes=nothing` means requirement-local storage; explicit axes identify the
global runtime storage coordinate system used by the bound parameter.
"""
struct ParameterBinding{R<:ParameterRequirement,P<:Tuple,A}
    requirement::R
    parameter::Symbol
    runtime_path::P
    storage_axes::A
end

"""Concrete participant applicability of one bound parameter requirement."""
struct ParameterApplicability{B,C,T}
    binding::B
    axis_components::C
    axis_tracers::T
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
    return ParameterRequirement(
        process_id(named),
        path,
        formulation_value,
        slot.name,
        slot.axes;
        qualifier=_slot_qualifier(slot, context),
        shape=slot.shape,
    )
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

function _factor_parameter_nodes(path::Tuple, factor::AbstractFactor)
    nodes = (_parameter_node(
        path, formulation(factor); context=factor_parameter_context(factor)
    ),)
    for (name, child) in pairs(factor_children(factor))
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

function _parameter_nodes(named::NamedProcess{P}) where {P<:Growth}
    process = named.process
    nodes = ()
    for (name, factor) in pairs(process.factors)
        nodes = (nodes..., _factor_parameter_nodes(name, factor)...)
    end

    stoichiometry = process.stoichiometry
    isnothing(stoichiometry) && return nodes
    stoichiometry isa FixedStoichiometry || throw(
        ArgumentError("unsupported growth stoichiometry $(typeof(stoichiometry))"),
    )
    nutrient_factors = Tuple(factor for factor in values(process.factors) if factor isa Nutrients)
    length(nutrient_factors) == 1 || throw(
        ArgumentError(
            "fixed-stoichiometry growth requires exactly one multi-resource Nutrients factor"
        ),
    )
    for currency in keys(only(nutrient_factors).responses)
        nodes = (
            nodes...,
            _parameter_node(
                (:stoichiometry,), stoichiometry; context=(currency=currency,)
            ),
        )
    end
    return nodes
end

function _parameter_nodes(named::NamedProcess{P}) where {P<:Grazing}
    process = named.process
    process.formulation isa Union{IdealizedGrazing,PreferentialGrazing} || throw(
        ArgumentError("unsupported grazing formulation $(typeof(process.formulation))"),
    )
    nodes = (
        _parameter_node((), process.formulation),
        _parameter_node(
            (:palatability, :default), Val(:allometric); formulation_value=:allometric
        ),
        _parameter_node(
            (:assimilation, :default), Val(:binary); formulation_value=:binary
        ),
    )
    isnothing(process.routing) && return nodes
    return (nodes..., _routing_parameter_nodes(process.routing)...)
end

function _parameter_nodes(named::NamedProcess{P}) where {P<:Consumption}
    process = named.process
    process.formulation isa HeterotrophicConsumption || throw(
        ArgumentError("unsupported consumption formulation $(typeof(process.formulation))"),
    )
    nodes = (_parameter_node((), process.formulation),)
    for (name, factor) in pairs(process.factors)
        nodes = (nodes..., _factor_parameter_nodes(name, factor)...)
    end
    return nodes
end

function _parameter_nodes(named::NamedProcess{P}) where {P<:Mortality}
    process = named.process
    context = length(process.populations) == 1 ?
              (population=only(process.populations),) : NamedTuple()
    nodes = (_parameter_node((), process.formulation; context),)
    isnothing(process.routing) && return nodes
    return (nodes..., _routing_parameter_nodes(process.routing)...)
end

function _parameter_nodes(named::NamedProcess{P}) where {P<:Remineralization}
    process = named.process
    process.formulation isa LinearRemineralization || throw(
        ArgumentError("unsupported remineralization formulation $(typeof(process.formulation))"),
    )
    return Tuple(
        _parameter_node((), process.formulation; context=(source=source,))
        for source in process.sources
    )
end

"""Return semantic parameter requirements from the named process scientific tree."""
function parameter_requirements(named::NamedProcess)
    requirements = ()
    for node in _parameter_nodes(named)
        requirements = (
            requirements...,
            _slot_requirements(
                named,
                node.path,
                node.node;
                context=node.context,
                formulation_value=node.formulation_value,
            )...,
        )
    end
    return requirements
end

"""Setup-time normalized scientific model definition."""
struct NormalizedModelDefinition{C,P,A,D,R,B}
    components::C
    processes::P
    parameters::A
    driver_identities::D
    parameter_requirements::R
    parameter_bindings::B
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
    for binding in definition.parameter_bindings
        binding.requirement.identity == identity && return binding
    end
    throw(ArgumentError("no model parameter is bound to requirement $identity"))
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
    references = Tuple(
        input.component for input in factor_inputs(factor) if input isa FactorComponent
    )
    for child in values(factor_children(factor))
        references = (references..., _factor_component_references(child)...)
    end
    return references
end

function _routing_component_references(routing::ProductRouting)
    routing.formulation isa DirectRouting && return (routing.retained,)
    routing.formulation isa PartitionRouting && return (routing.retained, routing.exported)
    routing.formulation isa DOMPOMRouting || throw(
        ArgumentError("unsupported product-routing formulation $(typeof(routing.formulation))"),
    )
    references = ()
    for pool in values(routing.pools), component in values(pool)
        references = (references..., component)
    end
    return references
end

function _process_component_references(process::Growth)
    references = process.populations
    for factor in values(process.factors)
        references = (references..., _factor_component_references(factor)...)
    end
    isnothing(process.source) || (references = (references..., process.source))
    return references
end

function _process_component_references(process::Grazing)
    destination = process.unassimilated_destination
    destination_refs = isnothing(destination) ? () : (destination,)
    routing_refs = isnothing(process.routing) ? () : _routing_component_references(process.routing)
    return (process.consumers..., process.resources..., destination_refs..., routing_refs...)
end

function _process_component_references(process::Consumption)
    destination = process.unassimilated_destination
    destination_refs = isnothing(destination) ? () : (destination,)
    references = (process.consumers..., process.resources..., destination_refs...)
    for factor in values(process.factors)
        references = (references..., _factor_component_references(factor)...)
    end
    return references
end

function _process_component_references(process::Mortality)
    routing = process.routing
    routing_refs = isnothing(routing) ? () : _routing_component_references(routing)
    return (process.populations..., routing_refs...)
end

_process_component_references(process::Remineralization) =
    (process.sources..., process.destinations...)

function _validate_currency_target(
    components::NamedTuple, component::Symbol, expected::Symbol, label::String
)
    actual = currency(getproperty(components, component))
    actual === expected || throw(
        ArgumentError("$label component :$component has currency :$actual, expected :$expected"),
    )
    return nothing
end

function _validate_process_science(process::Growth, components::NamedTuple)
    single_resource = any(factor -> factor isa NutrientResponse, values(process.factors))
    !single_resource || isnothing(process.source) || throw(
        ArgumentError(
            "single-resource growth derives its source from the nutrient response; omit `source`"
        ),
    )

    stoichiometry = process.stoichiometry
    stoichiometry isa FixedStoichiometry || return nothing
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

function _validate_process_science(process::Grazing, components::NamedTuple)
    isnothing(process.routing) || _validate_routing_science(process.routing, components)
    return nothing
end

function _validate_process_science(process::Consumption, components::NamedTuple)
    all(consumer -> getproperty(components, consumer) isa Population, process.consumers) || throw(
        ArgumentError("consumption consumers must be Population components"),
    )
    all(resource -> getproperty(components, resource) isa Pool, process.resources) || throw(
        ArgumentError("consumption resources must be Pool components"),
    )
    if !isnothing(process.unassimilated_destination)
        getproperty(components, process.unassimilated_destination) isa Pool || throw(
            ArgumentError("consumption unassimilated destination must be a Pool component"),
        )
    end

    reference = currency(getproperty(components, first(process.consumers)))
    for consumer in process.consumers
        _validate_currency_target(components, consumer, reference, "consumption consumer")
    end
    for resource in process.resources
        _validate_currency_target(components, resource, reference, "consumption resource")
    end
    isnothing(process.unassimilated_destination) || _validate_currency_target(
        components, process.unassimilated_destination, reference, "consumption destination"
    )
    return nothing
end

function _validate_process_science(process::Mortality, components::NamedTuple)
    isnothing(process.routing) || _validate_routing_science(process.routing, components)
    return nothing
end

_validate_process_science(::Remineralization, ::NamedTuple) = nothing

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
    requirements = ()
    for process in values(processes)
        requirements = (requirements..., parameter_requirements(process)...)
    end
    identities = map(requirement -> requirement.identity, requirements)
    length(unique(identities)) == length(identities) || throw(
        ArgumentError("normalized processes declare duplicate parameter requirement identities"),
    )
    return requirements
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
    isnothing(definitions) && return ()
    definitions isa Tuple || throw(
        ArgumentError("model parameters must be a tuple of ParameterDefinition values"),
    )
    all(definition -> definition isa ParameterDefinition, definitions) || throw(
        ArgumentError("model parameters must contain only ParameterDefinition values"),
    )

    provided = Dict{ParameterRequirementIdentity,Tuple{Symbol,Tuple,Union{Nothing,Symbol,NTuple{2,Symbol}}}}()

    for definition in definitions
        spec = definition.spec
        for provision in spec.provides
            requirement = _resolve_parameter_requirement(provision, requirements, spec.name)
            identity = requirement.identity
            spec.shape === requirement.shape || throw(
                ArgumentError(
                    "parameter :$(spec.name) provides $identity with required shape " *
                    "$(requirement.shape), not $(spec.shape)",
                ),
            )
            haskey(provided, identity) && throw(
                ArgumentError(
                    "parameter requirement $identity is provided by both :$(first(provided[identity])) and :$(spec.name)",
                ),
            )
            provided[identity] = (spec.name, spec.runtime_path, spec.axes)
        end
    end

    missing = filter(requirement -> !haskey(provided, requirement.identity), requirements)
    isempty(missing) || throw(
        ArgumentError(
            "model parameter definitions do not provide requirements $(map(r -> r.identity, missing))",
        ),
    )
    return Tuple(
        ParameterBinding(requirement, provided[requirement.identity]...)
        for requirement in requirements
    )
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
    bindings = _normalize_parameter_bindings(requirements, definition.parameters)
    return NormalizedModelDefinition(
        definition.components,
        normalized_processes,
        definition.parameters,
        _canonical_driver_identities(normalized_processes),
        requirements,
        bindings,
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

function _axis_tracers(layout::ComponentLayout, components::Tuple)
    tracers = ()
    for component in components
        hasproperty(layout.component_tracers, component) || throw(
            ArgumentError("parameter applicability references unrealized component :$component"),
        )
        tracers = (tracers..., getproperty(layout.component_tracers, component)...)
    end
    return tracers
end

"""Resolve each semantic parameter axis onto concrete component tracer identities.

The result is setup-time applicability metadata. It derives vector/matrix axes from the
participants of the owning named process rather than from global producer/consumer roles.
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
        axis_tracers = map(components -> _axis_tracers(layout, components), axis_components)
        ParameterApplicability(binding, axis_components, axis_tracers)
    end
end
