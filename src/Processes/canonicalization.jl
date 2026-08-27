"""Named, validated process instance in canonical model state.

`facts` contains setup-only decisions that compilation may trust, such as resolved
population-state references and Growth routing ownership.
"""
struct NamedProcess{P<:AbstractProcess,F}
    id::Symbol
    process::P
    facts::F
end

NamedProcess(id::Symbol, process::P) where {P<:AbstractProcess} =
    NamedProcess(id, process, NamedTuple())

abstract type AbstractGrowthRouting end

struct SingleResourceRouting <: AbstractGrowthRouting
    factor::Symbol
end

struct QuotaRouting <: AbstractGrowthRouting
    factor::Symbol
end

struct MultiResourceRouting <: AbstractGrowthRouting
    factor::Symbol
end

process_id(process::NamedProcess) = process.id
formulation(process::NamedProcess) = formulation(process.process)
factors(process::NamedProcess) = factors(process.process)
participants(process::NamedProcess) = participants(process.process)

function _resolve_slot_qualifier(slot::ParameterSlot, context::NamedTuple)
    isnothing(slot.qualify) && return nothing
    hasproperty(context, slot.qualify) || return nothing
    value = getproperty(context, slot.qualify)
    value isa Symbol || throw(
        ArgumentError("parameter slot qualifier :$(slot.qualify) must identify one Symbol"),
    )
    return Qualifier(slot.qualify, value)
end

function _parameter_slot_metadata(
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
        qualifier=_resolve_slot_qualifier(slot, context),
        axes=slot.axes,
    )
end

_qualifier_key(::Nothing) = nothing
_qualifier_key(qualifier::Qualifier) = (qualifier.axis, qualifier.value)
_binding_key(process::Symbol, path::Tuple, slot::Symbol, qualifier) =
    (process, path, slot, _qualifier_key(qualifier))
_binding_key(binding::ParameterBinding) =
    _binding_key(binding.process, binding.path, binding.slot, binding.qualifier)
_binding_key(metadata::NamedTuple) =
    _binding_key(metadata.process, metadata.path, metadata.slot, metadata.qualifier)

function _emit_parameter_slots!(
    uses::Vector{Any},
    named::NamedProcess,
    path::Tuple,
    node;
    context::NamedTuple=NamedTuple(),
)
    bindings = _resolve_authored_bindings(node)
    for slot in parameter_slots(_parameter_slot_source(node))
        qualifier = _resolve_slot_qualifier(slot, context)
        metadata = _parameter_slot_metadata(named, path, slot, context)
        parameter, explicit = _resolve_binding_value(bindings, slot, qualifier)
        push!(uses, (; metadata..., parameter, explicit))
    end
    return nothing
end

function _visit_factor_slots!(
    uses::Vector{Any}, named::NamedProcess, path::Tuple, factor::AbstractFactor
)
    _emit_parameter_slots!(
        uses, named, path, factor; context=factor_parameter_context(factor)
    )
    for (name, child) in pairs(_resolve_factor_children(factor))
        _visit_factor_slots!(
            uses, named, factor_child_path(path, factor, name), child
        )
    end
    return nothing
end

function _visit_product_slots!(
    uses::Vector{Any}, named::NamedProcess, path::Tuple, products::Products
)
    for product in keys(products.fractions)
        _emit_parameter_slots!(uses, named, path, products; context=(product=product,))
    end
    stoichiometry = products.stoichiometry
    isnothing(stoichiometry) && return nothing
    currencies = keys(first(values(products.targets)))
    for currency in currencies
        currency === stoichiometry.reference && continue
        _emit_parameter_slots!(
            uses,
            named,
            (path..., :stoichiometry),
            stoichiometry;
            context=(currency=currency,),
        )
    end
    return nothing
end

function _visit_process_slots!(uses::Vector{Any}, named::NamedProcess)
    process = named.process
    if process isa Mortality
        for population in process.populations
            _emit_parameter_slots!(uses, named, (), process; context=(population=population,))
        end
    elseif process isa Remineralization
        for source in process.sources
            _emit_parameter_slots!(uses, named, (), process; context=(source=source,))
        end
    else
        _emit_parameter_slots!(uses, named, (), process)
    end

    for (name, factor) in pairs(factors(process))
        _visit_factor_slots!(uses, named, (:factors, name), factor)
    end

    products = process_products(process)
    isnothing(products) || _visit_product_slots!(uses, named, product_path(process), products)

    if process isa Growth && !isnothing(process.stoichiometry)
        nutrients = getproperty(process.factors, named.facts.routing.factor)
        for currency in keys(nutrients.responses)
            _emit_parameter_slots!(
                uses,
                named,
                (:stoichiometry,),
                process.stoichiometry;
                context=(currency=currency,),
            )
        end
    end
    return nothing
end

"""Setup-time canonical scientific model definition.

`parameter_bindings` is the canonical ordered contract; `parameter_lookup` is a transient
setup cache used while lowering processes.
"""
struct CanonicalModelDefinition{C,P,A,D,B,L}
    components::C
    processes::P
    parameters::A
    driver_identities::D
    parameter_bindings::B
    parameter_lookup::L
end

"""Return the canonical external-driver identities required by a canonical model."""
driver_identities(definition::CanonicalModelDefinition) = definition.driver_identities

function _parameter_binding(
    definition::CanonicalModelDefinition,
    process::Symbol,
    path::Tuple,
    slot::Symbol,
    qualifier,
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
    definition::CanonicalModelDefinition,
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
            qualifier = _resolve_slot_qualifier(slot, context)
            _parameter_binding(
                definition, process_id(named), path, slot.name, qualifier
            )
        end for slot in slots
    )
    return NamedTuple{names}(bindings)
end


"""Attach setup-validated facts to a process before compilation.

Custom process implementations may extend this hook when lowering needs setup-resolved facts
beyond the authored process object.
"""
process_facts(::AbstractProcess, ::Symbol, ::NamedTuple) = NamedTuple()

function process_facts(process::Growth, id::Symbol, components::NamedTuple)
    routing = Tuple(
        (name, factor) for (name, factor) in pairs(process.factors)
        if factor isa Union{NutrientResponse,Nutrients}
    )
    length(routing) == 1 || throw(ArgumentError(
        "process :$id growth must declare exactly one NutrientResponse or Nutrients routing factor",
    ))
    factor_name, nutrient_factor = only(routing)

    rate_owners = Tuple(
        name for (name, factor) in pairs(process.factors)
        if any(slot -> slot.name === :maximum_rate, parameter_slots(formulation(factor)))
    )
    length(rate_owners) == 1 || throw(ArgumentError(
        "process :$id growth must declare exactly one factor that owns the maximum_rate slot",
    ))
    population_states = Tuple(
        _resolve_population_state(components, population, process.state, id, "population")
        for population in process.populations
    )

    if nutrient_factor isa NutrientResponse
        isnothing(process.source) || throw(ArgumentError(
            "process :$id single-resource growth derives its source from the nutrient response; omit `source`",
        ))
        isnothing(process.stoichiometry) || throw(ArgumentError(
            "process :$id single-resource growth does not take fixed stoichiometry",
        ))
        resource = _resolve_scalar_pool(components, nutrient_factor.resource, id, "nutrient factor resource")
        reference = currency(resource)
        _validate_state_currencies(components, population_states, reference, id, "population state")
        route = SingleResourceRouting(factor_name)
    else
        responses = values(_resolve_factor_children(nutrient_factor))
        quota_routing = all(response -> response isa QuotaResponse, responses)
        if quota_routing
            length(process.populations) == 1 || throw(ArgumentError(
                "process :$id quota growth requires exactly one logical population",
            ))
            process.source isa Symbol || throw(
                ArgumentError("process :$id quota growth requires a source component"),
            )
            isnothing(process.stoichiometry) || throw(ArgumentError(
                "process :$id quota growth uses independent NutrientUptake processes; omit fixed stoichiometry",
            ))
            reference_state = only(population_states)
            reference_currency = _state_currency(components, reference_state)
            source = _resolve_scalar_pool(components, process.source, id, "growth source")
            _validate_currency(currency(source), reference_currency, id, "growth source :$(process.source)")
            population = only(process.populations)
            for (target_currency, response) in pairs(nutrient_factor.responses)
                target = _resolve_population_state(components, response.target, id, "quota response target")
                reference = _resolve_population_state(components, response.reference, id, "quota response reference")
                target.population === population && reference.population === population || throw(
                    ArgumentError(
                        "process :$id quota response :$target_currency must reference growth population :$population",
                    ),
                )
                actual_target_currency = _state_currency(components, target)
                _validate_currency(
                    actual_target_currency, target_currency, id,
                    "quota response :$target_currency target state",
                )
                _validate_currency(
                    _state_currency(components, reference), reference_currency, id,
                    "quota response :$target_currency reference state",
                )
            end
            route = QuotaRouting(factor_name)
        else
            process.source isa Symbol || throw(
                ArgumentError("process :$id multi-resource growth requires a source component"),
            )
            process.stoichiometry isa FixedStoichiometry || throw(
                ArgumentError("process :$id multi-resource growth requires FixedStoichiometry"),
            )
            reference = process.stoichiometry.reference
            source = _resolve_scalar_pool(components, process.source, id, "growth source")
            _validate_currency(currency(source), reference, id, "growth source :$(process.source)")
            _validate_state_currencies(components, population_states, reference, id, "population state")
            for (target_currency, response) in pairs(nutrient_factor.responses)
                resource = _resolve_scalar_pool(
                    components, response.resource, id, "nutrient response :$target_currency resource",
                )
                _validate_currency(
                    currency(resource), target_currency, id,
                    "nutrient response :$target_currency resource :$(response.resource)",
                )
            end
            route = MultiResourceRouting(factor_name)
        end
    end
    return (;
        population_states, routing=route, maximum_rate_factor=only(rate_owners),
    )
end

function process_facts(
    process::NutrientUptake, id::Symbol, components::NamedTuple
)
    target = _resolve_population_state(
        components, process.population, process.target_state, id, "uptake target"
    )
    reference = _resolve_population_state(
        components, process.population, process.reference_state, id, "uptake reference"
    )
    target.state === reference.state && throw(ArgumentError(
        "process :$id nutrient uptake target and reference states must be distinct",
    ))
    resource = _resolve_scalar_pool(components, process.resource, id, "nutrient uptake resource")
    _validate_currency(
        currency(resource), _state_currency(components, target), id,
        "nutrient uptake resource :$(process.resource)",
    )
    required = Tuple(slot.name for slot in parameter_slots(process.formulation))
    missing = Tuple(name for name in required if !hasproperty(process.bindings, name))
    isempty(missing) || throw(ArgumentError(
        "process :$id nutrient uptake requires explicit bindings for slots $missing",
    ))
    return (; target, reference, resource=process.resource)
end

function process_facts(process::Consumption, id::Symbol, components::NamedTuple)
    consumer_states = Tuple(
        _resolve_population_state(components, consumer, id, "consumer") for consumer in process.consumers
    )
    reference = _state_currency(components, first(consumer_states))
    _validate_state_currencies(components, consumer_states, reference, id, "consumer state")

    if uses_living_interactions(process.formulation)
        resources = Tuple(
            _resolve_population_state(components, resource, id, "living resource")
            for resource in process.resources
        )
        _validate_state_currencies(components, resources, reference, id, "resource state")
    else
        resources = process.resources
        for resource in resources
            pool = getproperty(components, resource)
            pool isa Pool || throw(
                ArgumentError("process :$id heterotrophic resource :$resource must be a Pool"),
            )
            _validate_currency(currency(pool), reference, id, "heterotrophic resource :$resource")
        end
    end
    isnothing(process.products) || _validate_products(
        id, process.products, components, reference, "unassimilated products"
    )
    return (; consumer_states, resources)
end

function process_facts(process::Mortality, id::Symbol, components::NamedTuple)
    population_states = Tuple(
        _resolve_population_state(components, population, id, "mortality population")
        for population in process.populations
    )
    if !isnothing(process.products)
        reference = _state_currency(components, first(population_states))
        _validate_state_currencies(
            components, population_states, reference, id, "mortality population state"
        )
        _validate_products(id, process.products, components, reference, "mortality products")
    end
    return (; population_states)
end

function process_facts(process::Remineralization, id::Symbol, components::NamedTuple)
    destination = _resolve_scalar_pool(components, process.destination, id, "remineralization destination")
    reference = currency(destination)
    for source in process.sources
        pool = _resolve_scalar_pool(components, source, id, "remineralization source")
        _validate_currency(currency(pool), reference, id, "remineralization source :$source")
    end
    return NamedTuple()
end

function _canonicalize_process(id::Symbol, process, components::NamedTuple)
    process isa AbstractProcess || throw(
        ArgumentError("process :$id must be an AbstractProcess; got $(typeof(process))"),
    )
    component_names = keys(components)
    missing = unique(Tuple(
        reference for reference in _collect_process_component_references(process)
        if !(reference in component_names)
    ))
    isempty(missing) || throw(
        ArgumentError("process :$id references unknown components $missing"),
    )
    facts = process_facts(process, id, components)
    for (name, factor) in pairs(factors(process))
        _validate_factor_component_inputs(id, factor, components, (:factors, name))
    end
    return NamedProcess(id, process, facts)
end

function _canonical_processes(processes::NamedTuple, components::NamedTuple)
    names = sort!(collect(keys(processes)); by=String)
    names_tuple = Tuple(names)
    process_values = Tuple(
        _canonicalize_process(name, getproperty(processes, name), components) for name in names
    )
    return NamedTuple{names_tuple}(process_values)
end

function _collect_driver_identities!(identities::Vector{Symbol}, factor::AbstractFactor)
    for input in factor_inputs(factor)
        input isa FactorDriver || continue
        input.identity in identities || push!(identities, input.identity)
    end
    for child in values(_resolve_factor_children(factor))
        _collect_driver_identities!(identities, child)
    end
    return nothing
end

function _canonical_driver_identities(processes::NamedTuple)
    identities = Symbol[]
    for process in values(processes), factor in values(factors(process))
        _collect_driver_identities!(identities, factor)
    end
    sort!(identities; by=String)
    return Tuple(identities)
end


function _resolve_binding_value(bindings::NamedTuple, slot::ParameterSlot, qualifier)
    explicit = hasproperty(bindings, slot.name)
    value = explicit ? getproperty(bindings, slot.name) : slot.name
    value isa Symbol && return value, explicit

    value isa NamedTuple || throw(ArgumentError(
        "binding :$(slot.name) must be a parameter Symbol or one-level qualifier NamedTuple",
    ))
    isnothing(slot.qualify) && throw(ArgumentError(
        "binding :$(slot.name) uses a qualifier map but the slot is unqualified",
    ))
    qualifier isa Qualifier && qualifier.axis === slot.qualify || throw(ArgumentError(
        "binding :$(slot.name) requires qualifier :$(slot.qualify)",
    ))
    qualifier_value = qualifier.value
    hasproperty(value, qualifier_value) || throw(ArgumentError(
        "binding :$(slot.name) has no entry for qualifier :$qualifier_value",
    ))
    parameter = getproperty(value, qualifier_value)
    parameter isa Symbol || throw(ArgumentError(
        "binding :$(slot.name) qualifier :$qualifier_value must map to a parameter Symbol",
    ))
    return parameter, true
end

function _collect_parameter_uses(processes::NamedTuple)
    uses = Any[]
    for named in values(processes)
        _visit_process_slots!(uses, named)
    end
    keys = Tuple(_binding_key(use) for use in uses)
    length(unique(keys)) == length(keys) || throw(
        ArgumentError("canonical processes declare duplicate parameter binding keys"),
    )
    return Tuple(uses)
end

const RESERVED_PARAMETER_KEYS = (:x, :y, :z, :t)

function _resolve_parameter_definitions(definitions)
    isnothing(definitions) && return nothing, Set{Symbol}()
    definitions isa NamedTuple || throw(
        ArgumentError("model parameters must be a NamedTuple of Parameter or MetaParameter values"),
    )
    all(parameter -> parameter isa Union{Parameter,MetaParameter}, values(definitions)) || throw(
        ArgumentError("model parameters must contain only Parameter or MetaParameter values"),
    )

    definition_names = Tuple(keys(definitions))
    definition_set = Set(definition_names)
    for name in definition_names
        name in RESERVED_PARAMETER_KEYS && throw(
            ArgumentError("model parameters declare reserved parameter name :$name"),
        )
    end

    dependency_names = Set{Symbol}()
    for (name, parameter) in pairs(definitions)
        default = parameter.default
        default isa DerivedDefault || continue
        for dependency in default.deps
            dependency in definition_set || throw(
                ArgumentError(
                    "parameter :$name default depends on undeclared parameter :$dependency",
                ),
            )
            dependency_default = getproperty(definitions, dependency).default
            dependency_default isa DerivedDefault && throw(
                ArgumentError(
                    "parameter :$name default depends on derived parameter :$dependency; " *
                    "DerivedDefault dependencies must be non-derived parameters",
                ),
            )
            push!(dependency_names, dependency)
        end
    end
    return definitions, dependency_names
end

function _resolve_parameter_bindings(
    processes::NamedTuple, definitions, dependency_names::Set{Symbol}
)
    uses = _collect_parameter_uses(processes)
    if isnothing(definitions)
        bindings = Tuple(
            ParameterBinding(
                use.process, use.path, use.slot, use.qualifier, use.axes, use.parameter,
            )
            for use in uses
        )
        return bindings
    end

    definition_names = Set(keys(definitions))
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
        required_axes = unique(Tuple(use.axes for use in parameter_uses))
        length(required_axes) == 1 || throw(ArgumentError(
            "parameter :$parameter is bound to slots with incompatible semantic axes $(Tuple(required_axes))",
        ))
    end

    for (name, parameter) in pairs(definitions)
        parameter_uses = get(by_parameter, name, Any[])
        if parameter isa MetaParameter
            isempty(parameter_uses) || throw(ArgumentError(
                "MetaParameter :$name is construction-only and cannot bind to a process slot; use Parameter for runtime process parameters",
            ))
            name in dependency_names || throw(ArgumentError(
                "MetaParameter :$name must be used by at least one DerivedDefault",
            ))
        else
            isempty(parameter_uses) && throw(ArgumentError(
                "Parameter :$name is not bound to a process slot; use MetaParameter for construction-only DerivedDefault inputs",
            ))
        end
    end

    return Tuple(
        ParameterBinding(
            use.process,
            use.path,
            use.slot,
            use.qualifier,
            use.axes,
            use.parameter,
        )
        for use in uses
    )
end

"""Canonicalize process identity and resolve inline parameter bindings.

Process instances are canonicalized by stable process ID, so declaration order does
not change canonical scientific identity. Component ordering is preserved because it
still participates in concrete tracer realization. Local formulation slots bind directly
to stable model parameter names during canonicalization.
"""
function canonicalize_model(definition::ModelDefinition)
    all(component -> component isa Union{Population,Pool}, values(definition.components)) ||
        throw(ArgumentError("model components must be Population or Pool values"))
    canonical_processes = _canonical_processes(
        definition.processes, definition.components
    )
    parameters, dependency_names = _resolve_parameter_definitions(definition.parameters)
    bindings = _resolve_parameter_bindings(
        canonical_processes, parameters, dependency_names
    )
    lookup = Dict(_binding_key(binding) => binding for binding in bindings)
    return CanonicalModelDefinition(
        definition.components,
        canonical_processes,
        parameters,
        _canonical_driver_identities(canonical_processes),
        bindings,
        lookup,
    )
end
