using ..Parameters:
    ConstantDefault,
    DerivedDefault,
    NoDefault,
    DiameterIndexedVectorDefault,
    _derive_parameter_default

using ..Processes: planned_parameter

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
    for key in keys(overrides)
        hasproperty(plan.parameters, key) || throw(
            ArgumentError("parameters: unknown parameter key :$key."),
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
    value = value isa Bool ? value : T(value)
    rank = parameter.rank
    rank == 0 && return value

    expected = parameter.storage_shape
    if rank == 1
        return fill(value, only(expected))
    elseif rank == 2
        return fill(value, expected...)
    end
    throw(ArgumentError("parameter :$(parameter.name) has unsupported rank $rank"))
end

function _require_allometric_diameters(parameter, diameters, value, indices)
    value isa AllometricParam || return nothing
    labels = only(parameter.storage_labels)
    missing = Tuple(
        labels[i] for i in indices
        if !(diameters[i] isa Real && isfinite(diameters[i]) && diameters[i] > zero(diameters[i]))
    )
    isempty(missing) || throw(ArgumentError(
        "parameter :$(parameter.name) uses an allometric definition but has no diameter metadata for unsized PFT entities $missing; provide a size structure or an explicit value for those PFTs",
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

"""Resolve one-level `DerivedDefault` values from already-materialized parameters."""
function resolve_derived_parameter_defaults(
    plan,
    layout,
    params::NamedTuple;
    derivation_owner,
)
    resolved = params
    T = layout.scalar_type

    for (key, parameter) in pairs(plan.parameters)
        provider = parameter.definition.default
        provider isa DerivedDefault || continue
        hasproperty(resolved, key) && continue

        missing_deps = Tuple(dep for dep in provider.deps if !hasproperty(params, dep))
        isempty(missing_deps) || throw(ArgumentError(
            "derived default :$key is missing dependencies: " * join(string.(missing_deps), ", "),
        ))
        dependencies = NamedTuple{provider.deps}(
            Tuple(getproperty(params, dep) for dep in provider.deps)
        )

        value = _derive_parameter_default(
            provider.deriver, derivation_owner, layout, parameter, dependencies
        )
        validate_parameter_value(parameter, value, T; derived=true)
        resolved = merge(resolved, NamedTuple{(key,)}((value,)))
    end
    return resolved
end
