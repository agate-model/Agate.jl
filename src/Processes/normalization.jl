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

function _factor_parameter_requirements(named::NamedProcess, name::Symbol, factor::Light{Geider})
    path = (:factors, name)
    return (
        _requirement(named, path, factor.formulation, :maximum_rate, (:population,)),
        _requirement(named, path, factor.formulation, :alpha, (:population,)),
        _requirement(
            named,
            path,
            factor.formulation,
            :chlorophyll_to_carbon_ratio,
            (:population,),
        ),
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


function _nutrient_response_requirements(
    named::NamedProcess, path::Tuple, response::NutrientResponse{Monod}
)
    return (
        _requirement(
            named,
            path,
            response.formulation,
            :K,
            (:population,);
            qualifier=(resource=response.resource,),
        ),
    )
end


function _factor_parameter_requirements(
    named::NamedProcess, name::Symbol, factor::Nutrients{Liebig}
)
    requirements = ()
    for (response_name, response) in pairs(factor.responses)
        requirements = (
            requirements...,
            _nutrient_response_requirements(
                named, (:factors, name, :responses, response_name), response
            )...,
        )
    end
    return requirements
end

function _factor_parameter_requirements(named::NamedProcess, name::Symbol, factor::Temperature{Q10})
    path = (:factors, name)
    return (
        _requirement(named, path, factor.formulation, :q10, (); shape=:scalar),
        _requirement(
            named, path, factor.formulation, :reference_temperature, (); shape=:scalar
        ),
    )
end

function _factor_parameter_requirements(named::NamedProcess, name::Symbol, factor::AbstractFactor)
    throw(ArgumentError("unsupported process factor :$name of type $(typeof(factor))"))
end

"""Return semantic parameter requirements declared by a named process formulation."""
function parameter_requirements(named::NamedProcess{P}) where {P<:Growth}
    requirements = ()
    for (name, factor) in pairs(named.process.factors)
        requirements = (requirements..., _factor_parameter_requirements(named, name, factor)...)
    end
    stoichiometry = named.process.stoichiometry
    if stoichiometry isa FixedStoichiometry
        nutrient_factors = Tuple(
            factor for factor in values(named.process.factors) if factor isa Nutrients
        )
        length(nutrient_factors) == 1 || throw(
            ArgumentError(
                "fixed-stoichiometry growth requires exactly one multi-resource Nutrients factor"
            ),
        )
        nutrients = only(nutrient_factors)
        for currency in keys(nutrients.responses)
            requirements = (
                requirements...,
                _requirement(
                    named,
                    (:stoichiometry,),
                    stoichiometry,
                    :ratio,
                    ();
                    qualifier=(currency=currency,),
                    shape=:scalar,
                ),
            )
        end
    end
    return requirements
end


function _routing_parameter_requirements(named::NamedProcess, routing::ProductRouting)
    if routing.formulation isa PartitionRouting
        return (
            _requirement(named, (:routing,), routing.formulation, :export_fraction, ()),
        )
    elseif routing.formulation isa DOMPOMRouting
        requirements = (
            _requirement(named, (:routing,), routing.formulation, :POM_fraction, ()),
        )
        reference = routing.stoichiometry.reference
        for currency in keys(routing.pools.DOM)
            currency === reference && continue
            requirements = (
                requirements...,
                _requirement(
                    named,
                    (:routing, :stoichiometry),
                    routing.stoichiometry,
                    :ratio,
                    ();
                    qualifier=(currency=currency,),
                    shape=:scalar,
                ),
            )
        end
        return requirements
    end
    throw(ArgumentError("unsupported product-routing formulation $(typeof(routing.formulation))"))
end

function parameter_requirements(named::NamedProcess{P}) where {P<:Grazing}
    process = named.process
    process.formulation isa PreferentialGrazing || throw(
        ArgumentError("unsupported grazing formulation $(typeof(process.formulation))"),
    )
    requirements = (
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
    isnothing(process.routing) && return requirements
    return (requirements..., _routing_parameter_requirements(named, process.routing)...)
end

function parameter_requirements(named::NamedProcess{P}) where {P<:Consumption}
    process = named.process
    process.formulation isa HeterotrophicConsumption || throw(
        ArgumentError("unsupported consumption formulation $(typeof(process.formulation))"),
    )
    requirements = (
        _requirement(named, (), process.formulation, :maximum_rate, (:consumer,)),
        _requirement(named, (), process.formulation, :half_saturation, (:resource,)),
        _requirement(named, (), process.formulation, :assimilation, (:consumer, :resource)),
    )
    for (name, factor) in pairs(process.factors)
        requirements = (requirements..., _factor_parameter_requirements(named, name, factor)...)
    end
    return requirements
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
    return (requirements..., _routing_parameter_requirements(named, routing)...)
end

function parameter_requirements(named::NamedProcess{P}) where {P<:Remineralization}
    process = named.process
    process.formulation isa LinearRemineralization || throw(
        ArgumentError("unsupported remineralization formulation $(typeof(process.formulation))"),
    )
    return Tuple(
        _requirement(
            named,
            (),
            process.formulation,
            :rate,
            (:source,);
            qualifier=(source=source,),
            shape=:scalar,
        ) for source in process.sources
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
_factor_component_references(::Temperature) = ()
_factor_component_references(factor::NutrientResponse) = (factor.resource,)
function _factor_component_references(factor::Nutrients)
    references = ()
    for response in values(factor.responses)
        references = (references..., _factor_component_references(response)...)
    end
    return references
end

function _routing_component_references(routing::ProductRouting)
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
    return (process.consumers..., process.resources..., destination_refs...)
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
