"""Author-facing scientific model definition."""
struct ModelDefinition{C,P,A}
    components::C
    processes::P
    parameters::A
end

function ModelDefinition(; components::NamedTuple, processes::NamedTuple, parameters=nothing)
    return ModelDefinition(components, processes, parameters)
end

function ModelDefinition(factory::AbstractBGCFactory)
    return ModelDefinition(;
        components=default_components(factory),
        processes=default_processes(factory),
        parameters=parameter_definitions(factory),
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
participants(process::NamedProcess) = participants(process.process)
drivers(process::NamedProcess) = drivers(process.process)
rate_axes(process::NamedProcess) = rate_axes(process.process)

"""Setup-time normalized scientific model definition."""
struct NormalizedModelDefinition{C,P,A,D}
    components::C
    processes::P
    parameters::A
    driver_identities::D
end

driver_identities(definition::NormalizedModelDefinition) = definition.driver_identities

function _process_component_references(process::Growth)
    limitation = process.limitation
    limitation_refs = isnothing(limitation) ? () : (limitation.resource,)
    return (process.populations..., limitation_refs...)
end

function _process_component_references(process::Grazing)
    destination = process.unassimilated_destination
    destination_refs = isnothing(destination) ? () : (destination,)
    return (process.consumers..., process.resources..., destination_refs...)
end

function _process_component_references(process::Mortality)
    routing = process.routing
    routing_refs = isnothing(routing) ? () : (routing.retained, routing.exported)
    return (process.populations..., routing_refs...)
end

_process_component_references(process::Remineralization) =
    (process.sources..., process.destinations...)

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

"""Normalize process identity and validate model-level references at setup time.

Process instances are canonicalized by stable process ID, so declaration order does
not change normalized scientific identity. Component ordering is preserved because it
still participates in concrete tracer realization.
"""
function normalize_model(definition::ModelDefinition)
    all(component -> component isa Union{Population,Pool}, values(definition.components)) ||
        throw(ArgumentError("model components must be Population or Pool values"))
    normalized_processes = _canonical_processes(
        definition.processes, definition.components
    )
    return NormalizedModelDefinition(
        definition.components,
        normalized_processes,
        definition.parameters,
        _canonical_driver_identities(normalized_processes),
    )
end

function _canonical_qualifier(qualifier::NamedTuple)
    names = sort!(collect(keys(qualifier)); by=String)
    names_tuple = Tuple(names)
    values = Tuple(getproperty(qualifier, name) for name in names)
    return NamedTuple{names_tuple}(values)
end

"""Stable identity of one semantic formulation parameter requirement.

The identity is scoped by the owning named process, nested sub-formulation path,
formulation tag, slot, and participant qualifier. Applicability axes are resolved in
a later construction stage from process participation.
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
