function _resolve_factor_subfactors(factor::AbstractFactor)
    subfactors = factor_subfactors(factor)
    subfactors isa NamedTuple || throw(
        ArgumentError("factor subfactors for $(typeof(factor)) must be a NamedTuple"),
    )
    all(subfactor -> subfactor isa AbstractFactor, values(subfactors)) || throw(
        ArgumentError("factor subfactors for $(typeof(factor)) must be process factors"),
    )
    if factor isa Nutrients
        responses = values(subfactors)
        all_external = all(response -> response isa NutrientResponse, responses)
        all_internal = all(response -> response isa QuotaResponse, responses)
        all_external || all_internal || throw(ArgumentError(
            "nutrient `responses` must be either all NutrientResponse or all QuotaResponse factors",
        ))
    end
    return subfactors
end

function _factor_contains(factor::AbstractFactor, ::Type{T}) where {T<:AbstractFactor}
    factor isa T && return true
    return any(
        subfactor -> _factor_contains(subfactor, T), values(_resolve_factor_subfactors(factor))
    )
end

function _collect_factor_component_references(factor::AbstractFactor)
    references = Symbol[]
    for input in factor_inputs(factor)
        input isa FactorComponent && push!(references, input.component)
        input isa FactorPlanktonState && push!(references, input.reference.plankton)
    end
    for subfactor in values(_resolve_factor_subfactors(factor))
        append!(references, _collect_factor_component_references(subfactor))
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

_is_scalar_component(::Pool) = true
_is_scalar_component(component::Plankton) =
    isnothing(size_structure(component)) && length(states(component)) == 1

function _resolve_pool(components, name::Symbol, id::Symbol, label::AbstractString)
    pool = getproperty(components, name)
    pool isa Pool || throw(ArgumentError("process :$id $label :$name must be a Pool"))
    return pool
end

function _factor_element(
    id::Symbol, process::AbstractProcess, factor::NutrientResponse, components, path::Tuple
)
    resource = _resolve_pool(
        components, factor.resource, id, "factor path $path nutrient resource"
    )
    return element(resource)
end

function _factor_element(
    id::Symbol, process::AbstractProcess, factor::QuotaResponse, components, path::Tuple
)
    process isa Growth || throw(ArgumentError(
        "process :$id factor path $path uses QuotaResponse, which is only valid for Growth processes",
    ))
    length(process.plankton) == 1 || throw(ArgumentError(
        "process :$id quota growth requires exactly one logical plankton",
    ))
    target = _resolve_plankton_state(
        components, only(process.plankton), factor.variable_state, id,
        "quota response variable state",
    )
    target_element = _state_element(components, target)
    isnothing(target_element) && throw(ArgumentError(
        "process :$id quota response variable state :$(factor.variable_state) must represent an Element",
    ))
    return target_element
end

function _validate_factor_for_process(
    id::Symbol,
    process::AbstractProcess,
    factor::AbstractFactor,
    components::NamedTuple,
    path::Tuple;
    expected_element=nothing,
)
    factor isa Light && !(process isa Growth) && throw(ArgumentError(
        "process :$id factor path $path uses Light, which is only valid for Growth processes",
    ))
    subfactors = _resolve_factor_subfactors(factor)
    if factor isa Union{NutrientResponse,QuotaResponse}
        factor_element = _factor_element(id, process, factor, components, path)
        isnothing(expected_element) || _validate_element(
            factor_element, expected_element, id, "factor path $path response element"
        )
    end
    for input in factor_inputs(factor)
        input isa FactorComponent || continue
        component = getproperty(components, input.component)
        _is_scalar_component(component) || throw(ArgumentError(
            "process :$id factor path $path component :$(input.component) must be a scalar component",
        ))
    end
    for (name, subfactor) in pairs(subfactors)
        _validate_factor_for_process(
            id,
            process,
            subfactor,
            components,
            factor_subfactor_path(path, factor, name);
            expected_element=factor isa Nutrients ? name : nothing,
        )
    end
    return nothing
end

function _plankton_component(components, name::Symbol, id::Symbol, label::AbstractString)
    hasproperty(components, name) || throw(
        ArgumentError("process :$id $label references unknown plankton :$name"),
    )
    plankton = getproperty(components, name)
    plankton isa Plankton || throw(
        ArgumentError("process :$id $label :$name must be a Plankton"),
    )
    return plankton
end

function _resolve_plankton_state(
    components, name::Symbol, state::Symbol, id::Symbol, label::AbstractString
)
    plankton = _plankton_component(components, name, id, label)
    state in states(plankton) || throw(
        ArgumentError("process :$id $label :$name has no state :$state"),
    )
    return PlanktonStateRef(name, state)
end

function _resolve_reference_state(components, name::Symbol, id::Symbol, label::AbstractString)
    plankton = _plankton_component(components, name, id, label)
    return PlanktonStateRef(name, reference_state(plankton))
end

function _plankton_state_refs(components, name::Symbol, id::Symbol, label::AbstractString)
    plankton = _plankton_component(components, name, id, label)
    return Tuple(PlanktonStateRef(name, state) for state in states(plankton))
end

function _plankton_state_elements(components, name::Symbol, id::Symbol, label::AbstractString)
    plankton = _plankton_component(components, name, id, label)
    names = states(plankton)
    return NamedTuple{names}(Tuple(state_element(plankton, state) for state in names))
end

function _plankton_element_states(components, name::Symbol, id::Symbol, label::AbstractString)
    plankton = _plankton_component(components, name, id, label)
    element_names = Symbol[]
    refs = PlanktonStateRef[]
    for state in states(plankton)
        state_element_value = state_element(plankton, state)
        isnothing(state_element_value) && continue
        state_element_value in element_names && throw(ArgumentError(
            "process :$id $label :$name has multiple states for element :$state_element_value; " *
            "explicit same-element state mappings are not yet enabled",
        ))
        push!(element_names, state_element_value)
        push!(refs, PlanktonStateRef(name, state))
    end
    return NamedTuple{Tuple(element_names)}(Tuple(refs))
end

_state_element(components, ref::PlanktonStateRef) =
    state_element(getproperty(components, ref.plankton), ref.state)

function _validate_element(actual, expected, id::Symbol, label::AbstractString)
    actual === expected || throw(
        ArgumentError("process :$id $label has element $(repr(actual)), expected $(repr(expected))"),
    )
end

function _validate_state_elements(components, refs, expected, id, label)
    for ref in refs
        _validate_element(
            _state_element(components, ref), expected, id,
            "$label :$(ref.plankton).$(ref.state)",
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
