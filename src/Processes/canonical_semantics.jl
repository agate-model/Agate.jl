# Canonical model state and process semantic facts

"""Setup-time canonical scientific model definition.

`parameter_bindings` is the single ordered parameter-binding representation. Canonical
processes carry dense references into this tuple so setup and lowering do not maintain a
second path-key lookup representation.
"""
struct CanonicalModelDefinition{
    ComponentMap,ProcessMap,ParameterDefinitions,DriverIdentities,ParameterBindings
}
    components::ComponentMap
    processes::ProcessMap
    parameters::ParameterDefinitions
    driver_identities::DriverIdentities
    parameter_bindings::ParameterBindings
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
    isnothing(reference_element) && any(destination -> destination isa Symbol, values(products.destinations)) &&
        throw(ArgumentError(
            "process :$id $label uses scalar product destinations but its reference state has no element",
        ))
    names = keys(products.destinations)
    targets = NamedTuple{names}(Tuple(
        _canonical_target_elements(getproperty(products.destinations, name), reference_element)
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

"""Attach setup-validated semantic facts to a process before compilation.

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
