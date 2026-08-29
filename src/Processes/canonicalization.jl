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

"""Return the stable identity of a canonical named process."""
process_id(process::NamedProcess) = process.id
formulation(process::NamedProcess) = formulation(process.process)
factors(process::NamedProcess) = factors(process.process)
participants(process::NamedProcess) = participants(process.process)

# Dense parameter-binding reference collection

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
        push!(uses, (; metadata..., parameter, explicit))
        length(uses)
    end
    return NamedTuple{names}(refs)
end

function _visit_factor_slots!(
    uses::Vector{Any}, seen::Set{Any}, named::NamedProcess, path::Tuple, factor::AbstractFactor
)
    slots = _emit_parameter_slots!(uses, seen, named, path, factor)
    subfactors = _resolve_factor_subfactors(factor)
    names = keys(subfactors)
    subfactor_refs = NamedTuple{names}(Tuple(
        _visit_factor_slots!(
            uses, seen, named, factor_subfactor_path(path, factor, name), subfactor
        )
        for (name, subfactor) in pairs(subfactors)
    ))
    return (; slots, subfactors=subfactor_refs)
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
        )
        for product in fraction_names
    ))

    stoichiometry = products.stoichiometry
    isnothing(stoichiometry) && return (; fractions, stoichiometry=NamedTuple())
    elements = Tuple(
        element for element in keys(first(values(products.targets)))
        if element !== stoichiometry.reference_element
    )
    ratios = NamedTuple{elements}(Tuple(
        _emit_parameter_slots!(
            uses,
            seen,
            named,
            (path..., :stoichiometry),
            stoichiometry;
            context=(element=element,),
        )
        for element in elements
    ))
    return (; fractions, stoichiometry=ratios)
end

function _visit_process_slots!(uses::Vector{Any}, seen::Set{Any}, named::NamedProcess)
    process = named.process
    process_refs = if process isa Mortality
        names = process.plankton
        NamedTuple{names}(Tuple(
            _emit_parameter_slots!(
                uses, seen, named, (), process; context=(plankton=plankton,)
            )
            for plankton in names
        ))
    elseif process isa Remineralization
        names = process.sources
        NamedTuple{names}(Tuple(
            _emit_parameter_slots!(
                uses, seen, named, (), process; context=(source=source,)
            )
            for source in names
        ))
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
        elements = keys(process.additional_resources)
        stoichiometry_refs = NamedTuple{elements}(Tuple(
            _emit_parameter_slots!(
                uses,
                seen,
                named,
                (:stoichiometry,),
                process.stoichiometry;
                context=(element=element,),
            )
            for element in elements
        ))
    end
    return (;
        process=process_refs,
        factors=factor_refs,
        products=product_refs,
        stoichiometry=stoichiometry_refs,
    )
end

# Canonical model state and process facts

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

_canonical_target_elements(target::Symbol, reference_element::Symbol) =
    NamedTuple{(reference_element,)}((target,))
_canonical_target_elements(target::NamedTuple, _) = target

function _canonical_product_targets(id, products, components, reference_element, label)
    isnothing(reference_element) && any(target -> target isa Symbol, values(products.targets)) &&
        throw(ArgumentError(
            "process :$id $label uses scalar product targets but its reference state has no element",
        ))
    names = keys(products.targets)
    targets = NamedTuple{names}(Tuple(
        _canonical_target_elements(getproperty(products.targets, name), reference_element)
        for name in names
    ))

    for (name, elements) in pairs(targets), (target_element, target) in pairs(elements)
        pool = _resolve_pool(components, target, id, "$label product :$name target")
        _validate_element(
            element(pool), target_element, id, "$label product :$name target :$target",
        )
    end
    return targets
end

function _product_transfer_mode(
    id::Symbol,
    products::Products,
    product_targets::NamedTuple,
    source_elements::Tuple,
    reference_element,
    label::AbstractString,
)
    target_elements = Tuple(keys(first(values(product_targets))))
    stoichiometry = products.stoichiometry
    if isnothing(stoichiometry)
        sort(collect(target_elements); by=String) == sort(collect(source_elements); by=String) ||
            throw(ArgumentError(
                "process :$id $label product elements $target_elements must match the source " *
                "elemental states $source_elements when FixedStoichiometry is omitted",
            ))
        return :state
    end

    isnothing(reference_element) && throw(ArgumentError(
        "process :$id $label cannot use FixedStoichiometry because its reference state has no element",
    ))
    _validate_element(
        stoichiometry.reference_element,
        reference_element,
        id,
        "$label stoichiometric reference",
    )
    source_elements == (reference_element,) || throw(ArgumentError(
        "process :$id $label with multiple elemental states uses their prognostic inventories " *
        "directly; omit FixedStoichiometry and route each source element explicitly",
    ))
    return :stoichiometric
end

function _matching_element_sets(
    element_states::NamedTuple, id::Symbol, label::AbstractString
)
    isempty(element_states) && return ()
    first_elements = Tuple(keys(first(values(element_states))))
    for (name, states_for_element) in pairs(element_states)
        elements = Tuple(keys(states_for_element))
        sort(collect(elements); by=String) == sort(collect(first_elements); by=String) || throw(ArgumentError(
            "process :$id $label :$name has elemental states $elements; expected $first_elements " *
            "so one shared Products mapping can route every participant",
        ))
    end
    return first_elements
end

"""Attach setup-validated facts to a process before compilation.

Custom process implementations may extend this hook when lowering needs setup-resolved facts
beyond the authored process object.
"""
process_facts(::AbstractProcess, ::Symbol, ::NamedTuple) = NamedTuple()

function process_facts(process::Growth, id::Symbol, components::NamedTuple)
    plankton_states = Tuple(
        _resolve_reference_state(components, plankton, id, "plankton")
        for plankton in process.plankton
    )
    reference_element = _state_element(components, first(plankton_states))
    isnothing(reference_element) && throw(ArgumentError(
        "process :$id growth reference state must represent an Element",
    ))
    _validate_state_elements(
        components, plankton_states, reference_element, id, "plankton reference state"
    )

    reference_resource = _resolve_pool(
        components, process.reference_resource, id, "growth reference resource"
    )
    _validate_element(
        element(reference_resource),
        reference_element,
        id,
        "growth reference resource :$(process.reference_resource)",
    )

    has_stoichiometry = !isnothing(process.stoichiometry)
    has_stoichiometry == !isempty(process.additional_resources) || throw(ArgumentError(
        "process :$id `additional_resources` and FixedStoichiometry must be provided together",
    ))

    if has_stoichiometry
        _validate_element(
            process.stoichiometry.reference_element, reference_element, id,
            "growth stoichiometric reference",
        )
        for (target_element, resource_name) in pairs(process.additional_resources)
            target_element === reference_element && throw(ArgumentError(
                "process :$id `additional_resources` must not repeat reference element :$reference_element",
            ))
            resource = _resolve_pool(
                components, resource_name, id, "growth additional resource :$target_element"
            )
            _validate_element(
                element(resource), target_element, id,
                "growth additional resource :$target_element component :$resource_name",
            )
            for plankton_name in process.plankton
                explicit_elements = _plankton_element_states(
                    components, plankton_name, id, "growth plankton"
                )
                target_element in keys(explicit_elements) && throw(ArgumentError(
                    "process :$id growth additional resource :$target_element conflicts with " *
                    "explicit prognostic state :$(getproperty(explicit_elements, target_element).state) " *
                    "for plankton :$plankton_name; explicit elemental states must be updated by " *
                    "explicit state-changing processes",
                ))
            end
        end
    end

    return (; plankton_states)
end

function process_facts(
    process::NutrientUptake, id::Symbol, components::NamedTuple
)
    target = _resolve_plankton_state(
        components, process.plankton, process.target_state, id, "uptake target"
    )
    reference = _resolve_reference_state(
        components, process.plankton, id, "uptake reference"
    )
    target.state === reference.state && throw(ArgumentError(
        "process :$id nutrient uptake target and reference states must be distinct",
    ))
    resource = _resolve_pool(components, process.resource, id, "nutrient uptake resource")
    _validate_element(
        element(resource), _state_element(components, target), id,
        "nutrient uptake resource :$(process.resource)",
    )
    required = Tuple(slot.name for slot in parameter_slots(process.formulation))
    missing_names = Tuple(name for name in required if !hasproperty(process.bindings, name))
    isempty(missing_names) || throw(ArgumentError(
        "process :$id nutrient uptake requires explicit bindings for slots $missing_names",
    ))
    return (; target, reference, resource=process.resource)
end

function process_facts(process::Consumption, id::Symbol, components::NamedTuple)
    consumer_states = Tuple(
        _resolve_reference_state(components, consumer, id, "consumer")
        for consumer in process.consumers
    )
    reference_element = _state_element(components, first(consumer_states))
    isnothing(reference_element) && throw(ArgumentError(
        "process :$id consumption requires an elemental consumer reference state",
    ))
    _validate_state_elements(
        components, consumer_states, reference_element, id, "consumer reference state"
    )
    consumer_element_states = NamedTuple{process.consumers}(Tuple(
        _plankton_element_states(components, consumer, id, "consumer")
        for consumer in process.consumers
    ))

    if uses_living_interactions(process.formulation)
        resources = Tuple(
            _resolve_reference_state(components, resource, id, "living resource")
            for resource in process.resources
        )
        _validate_state_elements(
            components, resources, reference_element, id, "resource reference state"
        )
        resource_state_sets = NamedTuple{process.resources}(Tuple(
            _plankton_state_refs(components, resource, id, "living resource")
            for resource in process.resources
        ))
        resource_state_elements = NamedTuple{process.resources}(Tuple(
            _plankton_state_elements(components, resource, id, "living resource")
            for resource in process.resources
        ))
        resource_element_states = NamedTuple{process.resources}(Tuple(
            _plankton_element_states(components, resource, id, "living resource")
            for resource in process.resources
        ))

        for (resource, element_states) in pairs(resource_element_states)
            for source_element in keys(element_states), (consumer, consumer_states) in pairs(consumer_element_states)
                hasproperty(consumer_states, source_element) || throw(ArgumentError(
                    "process :$id consumer :$consumer has no state for resource :$resource " *
                    "element :$source_element required by elemental assimilation",
                ))
            end
        end

        product_targets, product_mode = if isnothing(process.products)
            nothing, nothing
        else
            source_elements = _matching_element_sets(
                resource_element_states, id, "living resource"
            )
            targets = _canonical_product_targets(
                id, process.products, components, reference_element, "unassimilated products"
            )
            mode = _product_transfer_mode(
                id,
                process.products,
                targets,
                source_elements,
                reference_element,
                "unassimilated products",
            )
            targets, mode
        end
        return (;
            consumer_states,
            consumer_element_states,
            resources,
            resource_state_sets,
            resource_state_elements,
            resource_element_states,
            reference_element,
            product_targets,
            product_mode,
        )
    end

    resources = process.resources
    for resource in resources
        pool = getproperty(components, resource)
        pool isa Pool || throw(
            ArgumentError("process :$id heterotrophic resource :$resource must be a Pool"),
        )
        _validate_element(
            element(pool), reference_element, id, "heterotrophic resource :$resource"
        )
    end
    product_targets, product_mode = if isnothing(process.products)
        nothing, nothing
    else
        targets = _canonical_product_targets(
            id, process.products, components, reference_element, "unassimilated products"
        )
        mode = _product_transfer_mode(
            id,
            process.products,
            targets,
            (reference_element,),
            reference_element,
            "unassimilated products",
        )
        targets, mode
    end
    return (;
        consumer_states,
        consumer_element_states,
        resources,
        resource_state_sets=NamedTuple(),
        resource_state_elements=NamedTuple(),
        resource_element_states=NamedTuple(),
        reference_element,
        product_targets,
        product_mode,
    )
end

function process_facts(process::Mortality, id::Symbol, components::NamedTuple)
    plankton_states = Tuple(
        _resolve_reference_state(components, plankton, id, "mortality plankton")
        for plankton in process.plankton
    )
    state_sets = NamedTuple{process.plankton}(Tuple(
        _plankton_state_refs(components, plankton, id, "mortality plankton")
        for plankton in process.plankton
    ))
    state_elements = NamedTuple{process.plankton}(Tuple(
        _plankton_state_elements(components, plankton, id, "mortality plankton")
        for plankton in process.plankton
    ))
    element_states = NamedTuple{process.plankton}(Tuple(
        _plankton_element_states(components, plankton, id, "mortality plankton")
        for plankton in process.plankton
    ))

    product_targets, product_mode = if isnothing(process.products)
        nothing, nothing
    else
        reference_element = _state_element(components, first(plankton_states))
        isnothing(reference_element) && throw(ArgumentError(
            "process :$id mortality products require an elemental reference state",
        ))
        _validate_state_elements(
            components, plankton_states, reference_element, id, "mortality reference state"
        )
        source_elements = _matching_element_sets(element_states, id, "mortality plankton")
        targets = _canonical_product_targets(
            id, process.products, components, reference_element, "mortality products"
        )
        mode = _product_transfer_mode(
            id,
            process.products,
            targets,
            source_elements,
            reference_element,
            "mortality products",
        )
        targets, mode
    end
    return (;
        plankton_states,
        state_sets,
        state_elements,
        element_states,
        product_targets,
        product_mode,
    )
end

function process_facts(process::Remineralization, id::Symbol, components::NamedTuple)
    destination = _resolve_pool(components, process.destination, id, "remineralization destination")
    reference = element(destination)
    for source in process.sources
        pool = _resolve_pool(components, source, id, "remineralization source")
        _validate_element(element(pool), reference, id, "remineralization source :$source")
    end
    return NamedTuple()
end

# Canonical process identity and drivers

function _canonicalize_process(id::Symbol, process, components::NamedTuple)
    process isa AbstractProcess || throw(
        ArgumentError("process :$id must be an AbstractProcess; got $(typeof(process))"),
    )
    component_names = keys(components)
    missing_names = unique(Tuple(
        reference for reference in _collect_process_component_references(process)
        if !(reference in component_names)
    ))
    isempty(missing_names) || throw(
        ArgumentError("process :$id references unknown components $missing_names"),
    )
    facts = process_facts(process, id, components)
    for (name, factor) in pairs(factors(process))
        _validate_factor_for_process(id, process, factor, components, (:factors, name))
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
    for subfactor in values(_resolve_factor_subfactors(factor))
        _collect_driver_identities!(identities, subfactor)
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


# Parameter definition and binding resolution

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

function _resolve_parameter_definitions(definitions)
    isnothing(definitions) && return nothing, Set{Symbol}()
    definitions isa NamedTuple || throw(
        ArgumentError("model parameters must be a NamedTuple of Parameter or MetaParameter values"),
    )
    all(parameter -> parameter isa Union{Parameter,MetaParameter}, values(definitions)) || throw(
        ArgumentError("model parameters must contain only Parameter or MetaParameter values"),
    )

    definition_set = Set(keys(definitions))

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
            ParameterBinding(use.axes, use.axis_components, use.parameter)
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
        ParameterBinding(use.axes, use.axis_components, use.parameter)
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

# Complete model canonicalization

"""Canonicalize process identity and resolve inline parameter bindings.

Process instances are canonicalized by stable process ID, so declaration order does
not change canonical scientific identity. Component order is preserved within component kinds;
`ModelLayout` realizes scalar Pools first and then Plankton in authored component order. Local
formulation slots bind directly to stable model parameter names during canonicalization.
"""
function canonicalize_model(definition::ModelDefinition)
    all(component -> component isa Union{Plankton,Pool}, values(definition.components)) ||
        throw(ArgumentError("model components must be Plankton or Pool values"))
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
