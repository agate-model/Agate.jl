"""Named, validated process instance in canonical model state.

`facts` contains setup-only scientific decisions that compilation may trust. `binding_refs`
contains dense references into the canonical model's ordered parameter-binding tuple, arranged
alongside the process/factor/product structure so lowering never reconstructs scientific paths.
"""
struct NamedProcess{P<:AbstractProcess,F,R}
    id::Symbol
    process::P
    facts::F
    binding_refs::R
end

NamedProcess(id::Symbol, process::P) where {P<:AbstractProcess} =
    NamedProcess(id, process, NamedTuple(), NamedTuple())
NamedProcess(id::Symbol, process::P, facts::F) where {P<:AbstractProcess,F} =
    NamedProcess(id, process, facts, NamedTuple())

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

function _binding_axis_components(
    named::NamedProcess, axes::Tuple, qualifier
)
    process_participants = participants(named)
    return map(axes) do axis
        qualifier isa Qualifier && qualifier.axis === axis && return (qualifier.value,)
        hasproperty(process_participants, axis) || throw(ArgumentError(
            "parameter applicability axis :$axis is not a participant role of process :$(process_id(named))",
        ))
        getproperty(process_participants, axis)
    end
end

function _parameter_slot_metadata(
    named::NamedProcess,
    path::Tuple,
    slot::ParameterSlot,
    qualifier,
)
    all(item -> item isa Symbol, path) || throw(
        ArgumentError("parameter binding path must contain only Symbols"),
    )
    return (;
        process=process_id(named),
        path,
        slot=slot.name,
        qualifier,
        axes=slot.axes,
        axis_components=_binding_axis_components(named, slot.axes, qualifier),
    )
end

function _emit_parameter_slots!(
    uses::Vector{Any},
    seen::Set{Any},
    named::NamedProcess,
    path::Tuple,
    node;
    context::NamedTuple=NamedTuple(),
    runtime_bound::Bool=true,
)
    slots = parameter_slots(_parameter_slot_source(node))
    names = Tuple(slot.name for slot in slots)
    bindings = _resolve_authored_bindings(node)
    refs = ntuple(length(slots)) do i
        slot = slots[i]
        qualifier = _resolve_slot_qualifier(slot, context)
        metadata = _parameter_slot_metadata(named, path, slot, qualifier)
        qualifier_key = isnothing(qualifier) ? nothing : (qualifier.axis, qualifier.value)
        identity = (metadata.process, metadata.path, metadata.slot, qualifier_key)
        identity in seen && throw(
            ArgumentError("canonical processes declare duplicate parameter binding key $identity"),
        )
        push!(seen, identity)
        parameter, explicit = _resolve_binding_value(bindings, slot, qualifier)
        push!(uses, (; metadata..., parameter, explicit, runtime_bound))
        length(uses)
    end
    return NamedTuple{names}(refs)
end

function _visit_factor_slots!(
    uses::Vector{Any}, seen::Set{Any}, named::NamedProcess, path::Tuple, factor::AbstractFactor
)
    slots = _emit_parameter_slots!(
        uses, seen, named, path, factor; context=factor_parameter_context(factor)
    )
    children = _resolve_factor_children(factor)
    names = keys(children)
    child_refs = NamedTuple{names}(Tuple(
        _visit_factor_slots!(
            uses, seen, named, factor_child_path(path, factor, name), child
        )
        for (name, child) in pairs(children)
    ))
    return (; slots, children=child_refs)
end

function _visit_product_slots!(
    uses::Vector{Any}, seen::Set{Any}, named::NamedProcess, path::Tuple, products::Products
)
    fraction_names = keys(products.fractions)
    fractions = NamedTuple{fraction_names}(Tuple(
        _emit_parameter_slots!(
            uses,
            seen,
            named,
            path,
            products;
            context=(product=product,),
            runtime_bound=product !== products.balanced,
        )
        for product in fraction_names
    ))

    stoichiometry = products.stoichiometry
    isnothing(stoichiometry) && return (; fractions, stoichiometry=NamedTuple())
    currencies = Tuple(
        currency for currency in keys(first(values(products.targets)))
        if currency !== stoichiometry.reference
    )
    ratios = NamedTuple{currencies}(Tuple(
        _emit_parameter_slots!(
            uses,
            seen,
            named,
            (path..., :stoichiometry),
            stoichiometry;
            context=(currency=currency,),
        )
        for currency in currencies
    ))
    return (; fractions, stoichiometry=ratios)
end

function _visit_process_slots!(uses::Vector{Any}, seen::Set{Any}, named::NamedProcess)
    process = named.process
    process_refs = if process isa Mortality
        Tuple(
            _emit_parameter_slots!(
                uses, seen, named, (), process; context=(population=population,)
            )
            for population in process.populations
        )
    elseif process isa Remineralization
        Tuple(
            _emit_parameter_slots!(
                uses, seen, named, (), process; context=(source=source,)
            )
            for source in process.sources
        )
    else
        _emit_parameter_slots!(uses, seen, named, (), process)
    end

    process_factors = factors(process)
    factor_names = keys(process_factors)
    factor_refs = NamedTuple{factor_names}(Tuple(
        _visit_factor_slots!(uses, seen, named, (:factors, name), factor)
        for (name, factor) in pairs(process_factors)
    ))

    products = process_products(process)
    product_refs = isnothing(products) ? nothing :
        _visit_product_slots!(uses, seen, named, product_path(process), products)

    stoichiometry_refs = NamedTuple()
    if process isa Growth && !isnothing(process.stoichiometry)
        currencies = keys(named.facts.additional_resources)
        stoichiometry_refs = NamedTuple{currencies}(Tuple(
            _emit_parameter_slots!(
                uses,
                seen,
                named,
                (:stoichiometry,),
                process.stoichiometry;
                context=(currency=currency,),
            )
            for currency in currencies
        ))
    end
    return (;
        process=process_refs,
        factors=factor_refs,
        products=product_refs,
        stoichiometry=stoichiometry_refs,
    )
end

"""Setup-time canonical scientific model definition.

`parameter_bindings` is the single ordered parameter-binding representation. Canonical
processes carry dense references into this tuple so setup and lowering do not maintain a
second path-key lookup representation.
"""
struct CanonicalModelDefinition{C,P,A,D,B}
    components::C
    processes::P
    parameters::A
    driver_identities::D
    parameter_bindings::B
end

"""Return the canonical external-driver identities required by a canonical model."""
driver_identities(definition::CanonicalModelDefinition) = definition.driver_identities

function _resolved_slot_bindings(definition::CanonicalModelDefinition, refs::NamedTuple)
    names = keys(refs)
    return NamedTuple{names}(Tuple(
        definition.parameter_bindings[ref] for ref in values(refs)
    ))
end

_canonical_target_currencies(target::Symbol, reference::Symbol) =
    NamedTuple{(reference,)}((target,))
_canonical_target_currencies(target::NamedTuple, ::Symbol) = target

function _canonical_product_targets(id, products, components, reference, label)
    names = keys(products.targets)
    targets = NamedTuple{names}(Tuple(
        _canonical_target_currencies(getproperty(products.targets, name), reference)
        for name in names
    ))

    stoichiometry = products.stoichiometry
    isnothing(stoichiometry) || _validate_currency(
        stoichiometry.reference, reference, id, "$label stoichiometric reference"
    )
    for (name, currencies) in pairs(targets), (target_currency, target) in pairs(currencies)
        pool = _resolve_scalar_pool(components, target, id, "$label product :$name target")
        _validate_currency(
            currency(pool), target_currency, id, "$label product :$name target :$target",
        )
    end
    return targets
end

"""Attach setup-validated facts to a process before compilation.

Custom process implementations may extend this hook when lowering needs setup-resolved facts
beyond the authored process object.
"""
process_facts(::AbstractProcess, ::Symbol, ::NamedTuple) = NamedTuple()

function process_facts(process::Growth, id::Symbol, components::NamedTuple)
    nutrient_factors = Tuple(
        factor for factor in values(process.factors)
        if factor isa Union{NutrientResponse,Nutrients}
    )
    length(nutrient_factors) == 1 || throw(ArgumentError(
        "process :$id growth must declare exactly one NutrientResponse or Nutrients factor",
    ))
    nutrient_factor = only(nutrient_factors)
    population_states = Tuple(
        _resolve_population_state(components, population, process.state, id, "population")
        for population in process.populations
    )

    reference_source, additional_resources = if nutrient_factor isa NutrientResponse
        isnothing(process.source) || throw(ArgumentError(
            "process :$id single-resource growth derives its source from the nutrient response; omit `source`",
        ))
        isnothing(process.stoichiometry) || throw(ArgumentError(
            "process :$id single-resource growth does not take fixed stoichiometry",
        ))
        resource = _resolve_scalar_pool(
            components, nutrient_factor.resource, id, "nutrient factor resource"
        )
        reference = currency(resource)
        _validate_state_currencies(
            components, population_states, reference, id, "population state"
        )
        nutrient_factor.resource, NamedTuple()
    else
        responses = values(_resolve_factor_children(nutrient_factor))
        quota_growth = all(response -> response isa QuotaResponse, responses)
        if quota_growth
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
            _validate_currency(
                currency(source), reference_currency, id, "growth source :$(process.source)"
            )
            population = only(process.populations)
            for (target_currency, response) in pairs(nutrient_factor.responses)
                target = _resolve_population_state(
                    components, response.target, id, "quota response target"
                )
                reference = _resolve_population_state(
                    components, response.reference, id, "quota response reference"
                )
                target.population === population && reference.population === population || throw(
                    ArgumentError(
                        "process :$id quota response :$target_currency must reference growth population :$population",
                    ),
                )
                _validate_currency(
                    _state_currency(components, target), target_currency, id,
                    "quota response :$target_currency target state",
                )
                _validate_currency(
                    _state_currency(components, reference), reference_currency, id,
                    "quota response :$target_currency reference state",
                )
            end
            process.source, NamedTuple()
        else
            process.source isa Symbol || throw(
                ArgumentError("process :$id multi-resource growth requires a source component"),
            )
            process.stoichiometry isa FixedStoichiometry || throw(
                ArgumentError("process :$id multi-resource growth requires FixedStoichiometry"),
            )
            reference = process.stoichiometry.reference
            source = _resolve_scalar_pool(components, process.source, id, "growth source")
            _validate_currency(
                currency(source), reference, id, "growth source :$(process.source)"
            )
            _validate_state_currencies(
                components, population_states, reference, id, "population state"
            )
            for (target_currency, response) in pairs(nutrient_factor.responses)
                resource = _resolve_scalar_pool(
                    components, response.resource, id,
                    "nutrient response :$target_currency resource",
                )
                _validate_currency(
                    currency(resource), target_currency, id,
                    "nutrient response :$target_currency resource :$(response.resource)",
                )
            end
            names = Tuple(
                currency for currency in keys(nutrient_factor.responses)
                if currency !== reference
            )
            additional = NamedTuple{names}(Tuple(
                getproperty(nutrient_factor.responses, currency).resource for currency in names
            ))
            process.source, additional
        end
    end

    return (; population_states, reference_source, additional_resources)
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
    product_targets = isnothing(process.products) ? nothing : _canonical_product_targets(
        id, process.products, components, reference, "unassimilated products"
    )
    return (; consumer_states, resources, product_targets)
end

function process_facts(process::Mortality, id::Symbol, components::NamedTuple)
    population_states = Tuple(
        _resolve_population_state(components, population, id, "mortality population")
        for population in process.populations
    )
    product_targets = if isnothing(process.products)
        nothing
    else
        reference = _state_currency(components, first(population_states))
        _validate_state_currencies(
            components, population_states, reference, id, "mortality population state"
        )
        _canonical_product_targets(
            id, process.products, components, reference, "mortality products"
        )
    end
    return (; population_states, product_targets)
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
    seen = Set{Any}()
    names = keys(processes)
    refs = NamedTuple{names}(Tuple(
        _visit_process_slots!(uses, seen, named) for named in values(processes)
    ))
    return Tuple(uses), refs
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
    uses::Tuple, definitions, dependency_names::Set{Symbol}
)
    if isnothing(definitions)
        bindings = Tuple(
            ParameterBinding(
                use.axes, use.axis_components, use.parameter, use.runtime_bound,
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
            use.axes,
            use.axis_components,
            use.parameter,
            use.runtime_bound,
        )
        for use in uses
    )
end

function _attach_binding_refs(processes::NamedTuple, refs::NamedTuple)
    names = keys(processes)
    return NamedTuple{names}(Tuple(
        begin
            named = getproperty(processes, name)
            NamedProcess(named.id, named.process, named.facts, getproperty(refs, name))
        end
        for name in names
    ))
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
    uses, binding_refs = _collect_parameter_uses(canonical_processes)
    bindings = _resolve_parameter_bindings(uses, parameters, dependency_names)
    bound_processes = _attach_binding_refs(canonical_processes, binding_refs)
    return CanonicalModelDefinition(
        definition.components,
        bound_processes,
        parameters,
        _canonical_driver_identities(canonical_processes),
        bindings,
    )
end
