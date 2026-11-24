using ..Parameters:
    ConstantDefault,
    DerivedDefault,
    NoDefault,
    DiameterIndexedVectorDefault,
    _prepare_parameter_default

using ..Processes:
    planned_parameter, runtime_parameter_values, parameter_constraints,
    validate_realized_parameters, ParameterValueRef, _parameter_value

using ..Library.Allometry: AbstractParamDef, AllometricParam, resolve_diameter_indexed_vector

function validate_parameter_value(parameter, value, ::Type{T}; derived::Bool=false) where {T<:Real}
    rank = parameter.rank
    name = parameter.name
    shape = parameter.storage_shape

    if rank == 0
        value isa Bool && return nothing
        value isa Number || throw(ArgumentError(
            "parameter :$name must be scalar; got $(typeof(value)).",
        ))
        if derived && typeof(value) !== T
            throw(ArgumentError(
                "derived default :$name must have type $(T); got $(typeof(value)). No implicit casting is performed.",
            ))
        end
        return nothing
    elseif rank == 1
        value isa AbstractVector || throw(ArgumentError(
            "parameter :$name must be a vector; got $(typeof(value)).",
        ))
        length(value) == only(shape) || throw(ArgumentError(
            "parameter :$name must have length $(only(shape)) (got $(length(value))).",
        ))
        if derived && eltype(value) !== T
            throw(ArgumentError(
                "derived default :$name must have eltype $(T); got eltype $(eltype(value)). No implicit casting is performed.",
            ))
        end
        return nothing
    elseif rank == 2
        value isa AbstractMatrix || throw(ArgumentError(
            "parameter :$name must be a matrix; got $(typeof(value)).",
        ))
        size(value) == shape || throw(ArgumentError(
            "parameter :$name must have size $shape (got $(size(value))).",
        ))
        if derived && eltype(value) !== T
            throw(ArgumentError(
                "derived default :$name must have eltype $(T); got eltype $(eltype(value)). No implicit casting is performed.",
            ))
        end
        return nothing
    end

    throw(ArgumentError("parameter :$name has unsupported rank $rank"))
end

function validate_parameter_storage(plan, values::NamedTuple, ::Type{T}) where {T<:Real}
    for name in keys(plan.parameters)
        validate_parameter_value(getproperty(plan.parameters, name), getproperty(values, name), T)
    end
    return nothing
end

function parameter_axis_names(parameter)
    parameter.rank == 1 || throw(ArgumentError(
        "parameter :$(parameter.name) does not have one vector storage axis",
    ))
    return only(parameter.storage_labels)
end

function expand_named_vector_override(
    parameter, default_value, user_value::NamedTuple, ::Type{T}
) where {T<:Real}
    name = parameter.name
    names = parameter_axis_names(parameter)
    length(default_value) == length(names) || throw(ArgumentError(
        "parameter :$name default length $(length(default_value)) does not match its planned axis length $(length(names)).",
    ))

    expanded = copy(default_value)
    for (key, value) in pairs(user_value)
        idx = findfirst(==(key), names)
        if idx === nothing
            expected = join(string.(names), ", ")
            throw(ArgumentError(
                "Unknown key `$(key)` for parameter `$name`. Expected one of: $(expected).",
            ))
        end
        expanded[idx] = value isa Bool ? value : T(value)
    end
    return expanded
end

function materialize_parameter_value(parameter, value, ::Type{T}) where {T<:Real}
    rank = parameter.rank
    if rank == 0
        return value isa Bool ? value : T(value)
    elseif rank == 1
        value isa AbstractVector || return value
        eltype(value) === T && return copy(value)
        out = similar(value, T, axes(value))
        copyto!(out, value)
        return out
    elseif rank == 2
        value isa AbstractMatrix || return value
        eltype(value) === T && return copy(value)
        out = similar(value, T, axes(value))
        copyto!(out, value)
        return out
    end
    throw(ArgumentError("parameter :$(parameter.name) has unsupported rank $rank"))
end

function validate_override_keys(plan, overrides::NamedTuple)
    for (key, value) in pairs(overrides)
        hasproperty(plan.parameters, key) || throw(
            ArgumentError("parameters: unknown parameter key :$key."),
        )
        value === nothing && throw(
            ArgumentError("parameter :$key cannot be overridden with `nothing`."),
        )
    end
    return nothing
end

@inline contains_missing(x) = x === missing

function contains_missing(x::NamedTuple)
    for v in values(x)
        contains_missing(v) && return true
    end
    return false
end

function contains_missing(x::AbstractArray)
    return any(ismissing, x)
end

function reject_missing_parameter_values(params::NamedTuple)
    for (k, v) in pairs(params)
        contains_missing(v) && throw(
            ArgumentError(
                "parameter :$k contains `missing`; all required parameters must be explicitly defined.",
            ),
        )
    end
    return nothing
end

function materialize_parameter_default(
    provider::ConstantDefault, parameter, ::Type{T}
) where {T<:Real}
    value = provider.value
    rank = parameter.rank
    rank == 0 && return value isa Bool ? value : T(value)

    expected = parameter.storage_shape
    if rank == 1
        value isa AbstractVector && return materialize_parameter_value(parameter, value, T)
        return fill(value isa Bool ? value : T(value), only(expected))
    elseif rank == 2
        value isa AbstractMatrix && return materialize_parameter_value(parameter, value, T)
        return fill(value isa Bool ? value : T(value), expected...)
    end
    throw(ArgumentError("parameter :$(parameter.name) has unsupported rank $rank"))
end

function _require_allometric_diameters(parameter, diameters, value, indices)
    value isa AllometricParam || return nothing
    labels = only(parameter.storage_labels)
    missing_labels = Tuple(
        labels[i] for i in indices
        if !(
            diameters[i] isa Real && isfinite(diameters[i]) &&
            diameters[i] > zero(diameters[i])
        )
    )
    isempty(missing_labels) || throw(ArgumentError(
        "parameter :$(parameter.name) uses an allometric definition but has no diameter metadata for SizeClasses $missing_labels; provide an explicit size structure or explicit values for those SizeClasses",
    ))
    return nothing
end

function materialize_parameter_default(
    provider::DiameterIndexedVectorDefault,
    parameter,
    ::Type{T};
    overridden_labels::Tuple=(),
) where {T<:Real}
    parameter.rank == 1 || throw(
        ArgumentError("DiameterIndexedVectorDefault requires vector parameter storage."),
    )
    diameters = parameter.storage_diameters
    diameters === nothing && throw(ArgumentError(
        "parameter :$(parameter.name) has no realized diameter axis for DiameterIndexedVectorDefault",
    ))
    labels = only(parameter.storage_labels)
    indices = Tuple(i for i in eachindex(labels) if labels[i] ∉ overridden_labels)
    _require_allometric_diameters(parameter, diameters, provider.value, indices)
    default = T(provider.default)
    return resolve_diameter_indexed_vector(T, diameters, indices, provider.value; default)
end

function materialize_parameter_defaults(
    plan, ::Type{T}, overrides::NamedTuple=(;)
) where {T<:Real}
    entries = Pair{Symbol,Any}[]
    for (name, parameter) in pairs(plan.parameters)
        provider = parameter.definition.default
        (provider isa NoDefault || provider isa DerivedDefault) && continue
        override = hasproperty(overrides, name) ? getproperty(overrides, name) : nothing
        override !== nothing && !(override isa NamedTuple) && continue
        if provider isa DiameterIndexedVectorDefault && override isa NamedTuple
            labels = parameter_axis_names(parameter)
            for key in keys(override)
                key in labels || throw(ArgumentError(
                    "Unknown key `$(key)` for parameter `$name`. Expected one of: $(join(string.(labels), ", ")).",
                ))
            end
            value = materialize_parameter_default(
                provider, parameter, T; overridden_labels=Tuple(keys(override))
            )
            push!(entries, name => value)
        else
            push!(entries, name => materialize_parameter_default(provider, parameter, T))
        end
    end
    return (; entries...)
end

function materialize_parameter_law_override(
    parameter, value::AbstractParamDef, ::Type{T}
) where {T<:Real}
    provider = parameter.definition.default
    provider isa DiameterIndexedVectorDefault || throw(ArgumentError(
        "parameter :$(parameter.name) only supports parameter-law overrides with a diameter-indexed vector default provider (DiameterIndexedVectorDefault).",
    ))
    parameter.rank == 1 || throw(ArgumentError(
        "parameter :$(parameter.name) diameter-indexed override requires vector storage",
    ))
    diameters = parameter.storage_diameters
    diameters === nothing && throw(ArgumentError(
        "parameter :$(parameter.name) has no realized diameter axis for a diameter-indexed override",
    ))
    indices = Tuple(eachindex(diameters))
    _require_allometric_diameters(parameter, diameters, value, indices)
    return resolve_diameter_indexed_vector(
        T,
        diameters,
        indices,
        value;
        default=T(provider.default),
    )
end

function materialize_parameter_overrides(
    plan,
    defaults::NamedTuple,
    overrides::NamedTuple,
    ::Type{T},
) where {T<:Real}
    isempty(overrides) && return overrides
    entries = Pair{Symbol,Any}[]
    for (key, value) in pairs(overrides)
        parameter = planned_parameter(plan, key)
        if value isa AbstractParamDef
            push!(entries, key => materialize_parameter_law_override(
                parameter, value, T
            ))
        elseif value isa NamedTuple
            parameter.rank == 1 || throw(ArgumentError(
                "parameter :$key does not support NamedTuple overrides because it is not vector-valued.",
            ))
            hasproperty(defaults, key) || throw(ArgumentError(
                "parameter :$key has no direct default for partial overrides.",
            ))
            push!(entries, key => expand_named_vector_override(
                parameter, getproperty(defaults, key), value, T
            ))
        else
            push!(entries, key => materialize_parameter_value(parameter, value, T))
        end
    end
    return (; entries...)
end

# A dependency tuple is fixed during preparation, rather than discovered inside AD.
function _derived_evaluator(kernel, ::Val{Deps}) where {Deps}
    return inputs -> kernel(NamedTuple{Deps}(map(dep -> getproperty(inputs, dep), Deps)))
end

function _prepare_derived_defaults(plan, layout, inputs; derivation_owner)
    names = Tuple(
        name for (name, parameter) in pairs(plan.parameters)
        if parameter.definition.default isa DerivedDefault && !hasproperty(inputs, name)
    )
    return NamedTuple{names}(map(names) do name
        parameter = planned_parameter(plan, name)
        provider = parameter.definition.default
        missing_deps = Tuple(dep for dep in provider.deps if !hasproperty(inputs, dep))
        isempty(missing_deps) || throw(ArgumentError(
            "derived default :$name is missing dependencies: " * join(string.(missing_deps), ", "),
        ))
        kernel = _prepare_parameter_default(provider.deriver, derivation_owner, layout, parameter)
        _derived_evaluator(kernel, Val(provider.deps))
    end)
end

@inline function _evaluate_derived_defaults(inputs, derivations)
    return merge(inputs, map(evaluate -> evaluate(inputs), derivations))
end

function _checked_parameter_values(canonical, layout, plan, resolved, derivations, constraints)
    T = layout.scalar_type
    for name in keys(derivations)
        validate_parameter_value(planned_parameter(plan, name), getproperty(resolved, name), T;
                                 derived=true)
    end
    required = keys(plan.parameters)
    missing_names = Tuple(name for name in required if !hasproperty(resolved, name))
    isempty(missing_names) || throw(
        ArgumentError("missing required parameters: $(join(string.(missing_names), ", "))")
    )
    ordered = NamedTuple{required}(map(name -> getproperty(resolved, name), required))
    reject_missing_parameter_values(ordered)
    validate_parameter_storage(plan, ordered, T)
    validate_realized_parameters(canonical, ordered, constraints)
    return ordered
end

"""Prepare ordinary parameter realization, retaining the numerical derivations and override provenance."""
function _prepare_parameter_realization(canonical, layout, plan, overrides; derivation_owner)
    T = layout.scalar_type
    validate_override_keys(plan, overrides)
    defaults = materialize_parameter_defaults(plan, T, overrides)
    inputs = merge(defaults, materialize_parameter_overrides(plan, defaults, overrides, T))
    derivations = _prepare_derived_defaults(plan, layout, inputs; derivation_owner)
    constraints = parameter_constraints(canonical, layout, plan)
    resolved = _evaluate_derived_defaults(inputs, derivations)
    values = _checked_parameter_values(canonical, layout, plan, resolved, derivations, constraints)
    return (; canonical, layout, plan, constraints, values, derivations, derivation_owner,
            explicit_override_keys=keys(overrides))
end

_copy_parameter_value(x) = x isa AbstractArray ? copy(x) : x

# Select whole materialized trait parameters with a fixed runtime active set.
# Authoring-level entry selection is separate from this numerical preparation.
function _prepare_trait_realization(prepared, traits::Tuple, active)
    T = prepared.layout.scalar_type
    T <: AbstractFloat || throw(ArgumentError("trait realization requires floating-point parameters"))
    allunique(traits) || throw(ArgumentError("duplicate construction traits"))
    for name in traits
        parameter = planned_parameter(prepared.plan, name)
        parameter.runtime_bound && throw(ArgumentError("parameter :$name is not a construction trait"))
        _validate_trait_input(getproperty(prepared.values, name), name)
    end
    targets = Tuple(
        name for name in keys(prepared.derivations)
        if any(dep -> dep in traits, planned_parameter(prepared.plan, name).definition.default.deps)
    )
    derivations = NamedTuple{targets}(map(name -> getproperty(prepared.derivations, name), targets))
    # Explicit output overrides were excluded at preparation. Other derived
    # outputs retain their baseline values, even when traits are updated repeatedly.
    numerics = (; baseline=map(_copy_parameter_value, prepared.values), derivations)
    return (; prepared, traits, numerics, active, gather=_prepare_runtime_gather(active))
end

function _validate_trait_input(value, name)
    entries = value isa Real ? (value,) : value
    entries isa Union{Tuple,NamedTuple,AbstractVector} &&
        all(v -> v isa Real && !(v isa Bool) && isfinite(v), entries) && return nothing
    throw(ArgumentError("construction trait :$name must contain finite real values"))
end

# Reuse the layout-resolved references used by coupled parameter constraints.
# The runtime active indices, not NamedTuple or matrix storage order, define p.
function _prepare_runtime_gather(active)
    refs = Vector{Any}(undef, length(active))
    for (name, selection) in pairs(active.map)
        if selection isa Integer
            refs[selection] = ParameterValueRef(name, ())
        else
            for slot in selection
                refs[slot.active_index] = ParameterValueRef(name, slot.indices)
            end
        end
    end
    return Tuple(refs)
end

# Numerical functions accept already-materialized updates, and only evaluate
# prepared derivations. The inputs and baseline are read-only during evaluation.
@inline function _rematerialize_traits(numerics, updates::NamedTuple)
    inputs = merge(numerics.baseline, updates)
    return _evaluate_derived_defaults(inputs, numerics.derivations)
end

@inline function _gather_runtime_vector(::Type{T}, resolved, gather) where {T<:Real}
    return collect(T, map(ref -> _parameter_value(resolved, ref), gather))
end

@inline function _realize_trait_vector(::Type{T}, numerics, updates, gather) where {T<:Real}
    return _gather_runtime_vector(T, _rematerialize_traits(numerics, updates), gather)
end

# Host-side checked realization accepts the ordinary named-vector override
# vocabulary. Keep validation and authoring-input expansion outside Enzyme.
function _realize_trait_parameters(selection, updates::NamedTuple)
    prepared = selection.prepared
    T = prepared.layout.scalar_type
    Set(keys(updates)) == Set(selection.traits) || throw(
        ArgumentError("updates must match the prepared construction traits $(selection.traits)")
    )
    for (name, value) in pairs(updates)
        _validate_trait_input(value, name)
    end
    materialized = materialize_parameter_overrides(
        prepared.plan, selection.numerics.baseline, updates, T
    )
    for (name, value) in pairs(materialized)
        validate_parameter_value(planned_parameter(prepared.plan, name), value, T)
        _validate_trait_input(value, name)
    end
    resolved = _rematerialize_traits(selection.numerics, materialized)
    checked = _checked_parameter_values(
        prepared.canonical, prepared.layout, prepared.plan, resolved,
        selection.numerics.derivations, prepared.constraints
    )
    return map(_copy_parameter_value, runtime_parameter_values(prepared.plan, checked))
end
