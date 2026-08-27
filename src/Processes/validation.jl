function _resolve_factor_children(factor::AbstractFactor)
    children = factor_children(factor)
    children isa NamedTuple || throw(
        ArgumentError("factor children for $(typeof(factor)) must be a NamedTuple"),
    )
    all(child -> child isa AbstractFactor, values(children)) || throw(
        ArgumentError("factor children for $(typeof(factor)) must be process factors"),
    )
    if factor isa Nutrients
        responses = values(children)
        all_external = all(response -> response isa NutrientResponse, responses)
        all_internal = all(response -> response isa QuotaResponse, responses)
        all_external || all_internal || throw(ArgumentError(
            "nutrient `responses` must be either all NutrientResponse or all QuotaResponse factors",
        ))
    end
    return children
end

function _collect_factor_component_references(factor::AbstractFactor)
    references = Symbol[]
    for input in factor_inputs(factor)
        input isa FactorComponent && push!(references, input.component)
        input isa FactorPopulationState && push!(references, input.reference.population)
    end
    for child in values(_resolve_factor_children(factor))
        append!(references, _collect_factor_component_references(child))
    end
    return Tuple(references)
end

function _collect_product_component_references(products::Products)
    return Tuple(Iterators.flatten(
        target isa Symbol ? (target,) : Tuple(values(target)) for target in values(products.targets)
    ))
end

function _collect_process_component_references(process::AbstractProcess)
    references = Symbol[]
    for values_for_role in values(participants(process))
        append!(references, values_for_role)
    end
    for factor in values(factors(process))
        append!(references, _collect_factor_component_references(factor))
    end
    products = process_products(process)
    isnothing(products) || append!(references, _collect_product_component_references(products))
    return Tuple(references)
end

_is_scalar_component(component::Pool) = isnothing(size_structure(component))
_is_scalar_component(component::Population) =
    isnothing(size_structure(component)) && length(states(component)) == 1

function _resolve_scalar_pool(components, name::Symbol, id::Symbol, label::AbstractString)
    pool = getproperty(components, name)
    pool isa Pool || throw(ArgumentError("process :$id $label :$name must be a Pool"))
    _is_scalar_component(pool) || throw(
        ArgumentError("process :$id $label :$name must be a scalar Pool"),
    )
    return pool
end

function _validate_factor_component_inputs(
    id::Symbol,
    factor::AbstractFactor,
    components::NamedTuple,
    path::Tuple,
)
    for input in factor_inputs(factor)
        input isa FactorComponent || continue
        component = getproperty(components, input.component)
        _is_scalar_component(component) || throw(ArgumentError(
            "process :$id factor path $path component :$(input.component) must be a scalar component",
        ))
    end
    for (name, child) in pairs(_resolve_factor_children(factor))
        _validate_factor_component_inputs(
            id, child, components, factor_child_path(path, factor, name)
        )
    end
    return nothing
end

function _resolve_population_state(
    components, name::Symbol, state, id::Symbol, label::AbstractString
)
    hasproperty(components, name) || throw(
        ArgumentError("process :$id $label references unknown population :$name"),
    )
    population = getproperty(components, name)
    population isa Population || throw(
        ArgumentError("process :$id $label :$name must be a Population"),
    )
    state_names = keys(states(population))
    if isnothing(state)
        length(state_names) == 1 || throw(
            ArgumentError("process :$id $label :$name requires explicit state selection"),
        )
        state = only(state_names)
    end
    state isa Symbol || throw(ArgumentError("process :$id $label state must be a Symbol"))
    state in state_names || throw(
        ArgumentError("process :$id $label :$name has no state :$state"),
    )
    return PopulationStateRef(name, state)
end

_resolve_population_state(components, name::Symbol, id::Symbol, label::AbstractString) =
    _resolve_population_state(components, name, nothing, id, label)

function _resolve_population_state(
    components, reference::PopulationStateRef, id::Symbol, label::AbstractString
)
    return _resolve_population_state(
        components, reference.population, reference.state, id, label
    )
end

_state_currency(components, ref::PopulationStateRef) =
    state_currency(getproperty(components, ref.population), ref.state)

function _validate_currency(actual, expected, id::Symbol, label::AbstractString)
    actual === expected || throw(
        ArgumentError("process :$id $label has currency :$actual, expected :$expected"),
    )
end

function _validate_state_currencies(components, refs, expected, id, label)
    for ref in refs
        _validate_currency(
            _state_currency(components, ref), expected, id,
            "$label :$(ref.population).$(ref.state)",
        )
    end
    return nothing
end

function _resolve_authored_bindings(node)
    slots = parameter_slots(_parameter_slot_source(node))
    bindings = if applicable(authored_parameter_bindings, node)
        authored_parameter_bindings(node)
    else
        isempty(slots) || throw(ArgumentError(
            "slot-owning node $(typeof(node)) must implement `authored_parameter_bindings`",
        ))
        NamedTuple()
    end
    slot_names = Tuple(slot.name for slot in slots)
    unknown = Tuple(name for name in keys(bindings) if !(name in slot_names))
    isempty(unknown) || throw(ArgumentError(
        "inline bindings contain unknown slots $unknown; expected $slot_names",
    ))
    return bindings
end
