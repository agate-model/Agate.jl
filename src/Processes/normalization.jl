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

"""Resolved mapping from one semantic requirement to a model parameter name."""
struct ParameterBinding{R<:ParameterRequirement}
    requirement::R
    parameter::Symbol
end

"""Concrete participant applicability of one bound parameter requirement."""
struct ParameterApplicability{B,C,T}
    binding::B
    axis_components::C
    axis_tracers::T
end

function _requirement(
    named::NamedProcess,
    path,
    formulation_value,
    slot,
    axes;
    qualifier=NamedTuple(),
    shape=_default_requirement_shape(axes),
)
    return ParameterRequirement(
        process_id(named), path, formulation_value, slot, axes; qualifier, shape
    )
end

function _factor_parameter_requirements(named::NamedProcess, name::Symbol, factor::Light{Smith})
    path = (:factors, name)
    return (
        _requirement(named, path, factor.formulation, :maximum_rate, (:population,)),
        _requirement(named, path, factor.formulation, :alpha, (:population,)),
    )
end

function _factor_parameter_requirements(
    named::NamedProcess, name::Symbol, factor::NutrientResponse{Monod}
)
    return (
        _requirement(
            named,
            (:factors, name),
            factor.formulation,
            :K,
            (:population,);
            qualifier=(resource=factor.resource,),
        ),
    )
end

function _factor_parameter_requirements(named::NamedProcess, name::Symbol, factor::AbstractFactor)
    throw(ArgumentError("unsupported growth factor :$name of type $(typeof(factor))"))
end

"""Return semantic parameter requirements declared by a named process formulation."""
function parameter_requirements(named::NamedProcess{P}) where {P<:Growth}
    requirements = ()
    for (name, factor) in pairs(named.process.factors)
        requirements = (requirements..., _factor_parameter_requirements(named, name, factor)...)
    end
    return requirements
end

function parameter_requirements(named::NamedProcess{P}) where {P<:Grazing}
    process = named.process
    process.formulation isa PreferentialGrazing || throw(
        ArgumentError("unsupported grazing formulation $(typeof(process.formulation))"),
    )
    return (
        _requirement(named, (), process.formulation, :maximum_rate, (:consumer,)),
        _requirement(named, (), process.formulation, :half_saturation, (:consumer,)),
        _requirement(named, (), process.formulation, :palatability, (:consumer, :resource)),
        _requirement(named, (), process.formulation, :assimilation, (:consumer, :resource)),
        _requirement(
            named,
            (:palatability, :default),
            :allometric,
            :optimum_predator_prey_ratio,
            (:consumer,),
        ),
        _requirement(
            named, (:palatability, :default), :allometric, :specificity, (:consumer,)
        ),
        _requirement(
            named, (:palatability, :default), :allometric, :protection, (:resource,)
        ),
        _requirement(
            named,
            (:assimilation, :default),
            :binary,
            :assimilation_efficiency,
            (:consumer,),
        ),
    )
end

function parameter_requirements(named::NamedProcess{P}) where {P<:Mortality}
    process = named.process
    qualifier = length(process.populations) == 1 ?
                (population=only(process.populations),) : NamedTuple()
    requirements = (
        _requirement(
            named, (), process.formulation, :rate, (:population,); qualifier
        ),
    )
    routing = process.routing
    isnothing(routing) && return requirements
    routing.formulation isa PartitionRouting || throw(
        ArgumentError("unsupported product-routing formulation $(typeof(routing.formulation))"),
    )
    return (
        requirements...,
        _requirement(named, (:routing,), routing.formulation, :export_fraction, ()),
    )
end

function parameter_requirements(named::NamedProcess{P}) where {P<:Remineralization}
    process = named.process
    process.formulation isa LinearRemineralization || throw(
        ArgumentError("unsupported remineralization formulation $(typeof(process.formulation))"),
    )
    qualifier = length(process.sources) == 1 ?
                (source=only(process.sources),) : NamedTuple()
    return (
        _requirement(
            named, (), process.formulation, :rate, (:source,); qualifier, shape=:scalar
        ),
    )
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

driver_identities(definition::NormalizedModelDefinition) = definition.driver_identities

"""Return formulation-declared semantic parameter requirements for a normalized model."""
parameter_requirements(definition::NormalizedModelDefinition) = definition.parameter_requirements

"""Return resolved semantic requirement-to-model-parameter bindings."""
parameter_bindings(definition::NormalizedModelDefinition) = definition.parameter_bindings

"""Return the model parameter name that supplies `requirement`."""
function parameter_name(definition::NormalizedModelDefinition, requirement::ParameterRequirement)
    return parameter_name(definition, requirement.identity)
end

function parameter_name(
    definition::NormalizedModelDefinition, identity::ParameterRequirementIdentity
)
    for binding in definition.parameter_bindings
        binding.requirement.identity == identity && return binding.parameter
    end
    throw(ArgumentError("no model parameter is bound to requirement $identity"))
end

_factor_component_references(::Light) = ()
_factor_component_references(factor::NutrientResponse) = (factor.resource,)

function _process_component_references(process::Growth)
    references = process.populations
    for factor in values(process.factors)
        references = (references..., _factor_component_references(factor)...)
    end
    return references
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

function _provision_identity(provision::ParameterProvision)
    return ParameterRequirementIdentity(
        provision.process,
        provision.path,
        provision.formulation,
        provision.slot;
        qualifier=provision.qualifier,
    )
end

function _normalize_parameter_bindings(requirements::Tuple, definitions)
    isnothing(definitions) && return ()
    definitions isa Tuple || throw(
        ArgumentError("model parameters must be a tuple of ParameterDefinition values"),
    )
    all(definition -> definition isa ParameterDefinition, definitions) || throw(
        ArgumentError("model parameters must contain only ParameterDefinition values"),
    )

    requirement_map = Dict(requirement.identity => requirement for requirement in requirements)
    provided = Dict{ParameterRequirementIdentity,Symbol}()

    for definition in definitions
        spec = definition.spec
        for provision in spec.provides
            identity = _provision_identity(provision)
            requirement = get(requirement_map, identity, nothing)
            isnothing(requirement) && throw(
                ArgumentError(
                    "parameter :$(spec.name) provides undeclared requirement $identity",
                ),
            )
            spec.shape === requirement.shape || throw(
                ArgumentError(
                    "parameter :$(spec.name) provides $(requirement.identity) with required shape $(requirement.shape), not $(spec.shape)",
                ),
            )
            haskey(provided, identity) && throw(
                ArgumentError(
                    "parameter requirement $identity is provided by both :$(provided[identity]) and :$(spec.name)",
                ),
            )
            provided[identity] = spec.name
        end
    end

    missing = filter(requirement -> !haskey(provided, requirement.identity), requirements)
    isempty(missing) || throw(
        ArgumentError(
            "model parameter definitions do not provide requirements $(map(r -> r.identity, missing))",
        ),
    )
    return Tuple(ParameterBinding(requirement, provided[requirement.identity]) for requirement in requirements)
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

function _axis_components(process::NamedProcess, axis::Symbol)
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
        axis_components = map(axis -> _axis_components(process, axis), requirement.axes)
        axis_tracers = map(components -> _axis_tracers(layout, components), axis_components)
        ParameterApplicability(binding, axis_components, axis_tracers)
    end
end
