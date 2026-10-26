# Canonical model state and process semantic facts

"""Setup-time canonical scientific model definition.

`parameter_bindings` is the single ordered parameter-binding representation. Canonical
processes carry dense references into this tuple so setup and lowering do not maintain a
second path-key lookup representation.
"""
struct CanonicalModelDefinition{
    ComponentMap,ProcessMap,ParameterDefinitions,DriverIdentities,ParameterBindings,CompositionRefs
}
    components::ComponentMap
    processes::ProcessMap
    parameters::ParameterDefinitions
    driver_identities::DriverIdentities
    parameter_bindings::ParameterBindings
    composition_refs::CompositionRefs
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

"""Resolve represented elemental inventories to states and optional composition parameters."""
function _material_sources(components, name::Symbol, id::Symbol, label::AbstractString)
    plankton = _plankton_component(components, name, id, label)
    explicit = _plankton_element_states(components, name, id, label)
    implicit = composition(plankton)
    names = Tuple(sort!(collect((keys(explicit)..., keys(implicit)...)); by=String))
    reference = PlanktonStateRef(name, reference_state(plankton))
    return NamedTuple{names}(Tuple(
        hasproperty(explicit, element) ? getproperty(explicit, element) : reference
        for element in names
    ))
end

function _check_product_elements(id, targets, represented, label)
    actual = Set(keys(first(values(targets))))
    expected = Set(represented)
    actual == expected || throw(ArgumentError(
        "process :$id $label product elements $actual must match represented elements $expected",
    ))
    return nothing
end

"""Elements represented by any participant; products route only each participant's inventory."""
_represented_elements(sources::NamedTuple) = Tuple(sort!(unique([
    element for inventory in values(sources) for element in keys(inventory)
]); by=String))

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

    required = Set(keys(process.additional_resources))
    implicit = Set{Symbol}()
    for plankton_name in process.plankton
        union!(implicit, keys(composition(getproperty(components, plankton_name))))
    end
    required == implicit || throw(ArgumentError(
        "process :$id growth additional_resources $required must match " *
        "implicit composition elements $implicit across participating plankton",
    ))
    for (target_element, resource_name) in pairs(process.additional_resources)
        resource = _resolve_pool(
            components, resource_name, id, "growth additional resource :$target_element"
        )
        _validate_element(
            element(resource), target_element, id,
            "growth additional resource :$target_element component :$resource_name",
        )
    end

    product_targets = if isnothing(process.products)
        nothing
    else
        targets = _canonical_product_targets(
            id, process.products, components, reference_element, "growth products"
        )
        _check_product_elements(
            id, targets, (reference_element, keys(process.additional_resources)...),
            "growth products",
        )
        targets
    end

    return (; plankton_states, reference_element, product_targets)
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
        _material_sources(components, consumer, id, "consumer")
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
        resource_element_states = NamedTuple{process.resources}(Tuple(
            _material_sources(components, resource, id, "living resource")
            for resource in process.resources
        ))

        for (resource, source) in pairs(resource_element_states)
            for (consumer, destination) in pairs(consumer_element_states)
                all(element -> hasproperty(destination, element), keys(source)) || throw(ArgumentError(
                    "process :$id consumer :$consumer lacks an elemental state " *
                    "required to assimilate prey :$resource",
                ))
            end
        end

        product_targets = if isnothing(process.products)
            nothing
        else
            source_elements = _represented_elements(resource_element_states)
            targets = _canonical_product_targets(
                id, process.products, components, reference_element, "unassimilated products"
            )
            _check_product_elements(id, targets, source_elements, "unassimilated products")
            targets
        end
        return (;
            consumer_states,
            consumer_element_states,
            resources,
            resource_state_sets,
            resource_element_states,
            reference_element,
            product_targets,
        )
    end

    any(!isempty(composition(getproperty(components, consumer))) for consumer in process.consumers) &&
        throw(ArgumentError(
            "process :$id heterotrophic uptake requires elemental balancing transfers for " *
            "consumers with implicit composition",
        ))
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
    product_targets = if isnothing(process.products)
        nothing
    else
        targets = _canonical_product_targets(
            id, process.products, components, reference_element, "unassimilated products"
        )
        _check_product_elements(id, targets, (reference_element,), "unassimilated products")
        targets
    end
    return (;
        consumer_states,
        consumer_element_states,
        resources,
        resource_state_sets=NamedTuple(),
        resource_element_states=NamedTuple(),
        reference_element,
        product_targets,
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
    material_sources = NamedTuple{process.plankton}(Tuple(
        _material_sources(components, plankton, id, "mortality plankton")
        for plankton in process.plankton
    ))
    product_targets = if isnothing(process.products)
        nothing
    else
        reference_element = _state_element(components, first(plankton_states))
        isnothing(reference_element) && throw(ArgumentError(
            "process :$id mortality products require an elemental reference state",
        ))
        _validate_state_elements(
            components, plankton_states, reference_element, id, "mortality reference state"
        )
        source_elements = _represented_elements(material_sources)
        targets = _canonical_product_targets(
            id, process.products, components, reference_element, "mortality products"
        )
        _check_product_elements(id, targets, source_elements, "mortality products")
        targets
    end
    return (;
        plankton_states,
        state_sets,
        material_sources,
        product_targets,
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
