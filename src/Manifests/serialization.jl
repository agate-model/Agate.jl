module Serialization

using ...Factories: AbstractBGCFactory, parameter_spec

function manifest_axis(axis)
    axis === nothing && return nothing
    axis isa Symbol && return string(axis)
    axis isa Tuple && return Any[string(a) for a in axis]
    return string(axis)
end

function manifest_value(x, name=nothing)
    if x === nothing || x isa Bool || x isa Integer || x isa AbstractString
        return x
    elseif x isa AbstractFloat
        return isfinite(x) ? x : string(x)
    elseif x isa Symbol
        return string(x)
    elseif x isa NamedTuple
        return Dict{String,Any}(string(k) => manifest_value(v, k) for (k, v) in pairs(x))
    elseif x isa AbstractDict
        return Dict{String,Any}(string(k) => manifest_value(v, k) for (k, v) in pairs(x))
    elseif x isa AbstractMatrix
        return Any[Any[manifest_value(x[i, j], name) for j in axes(x, 2)] for i in axes(x, 1)]
    elseif x isa AbstractVector || x isa Tuple
        return Any[manifest_value(v, name) for v in x]
    else
        label = isnothing(name) ? "" : " $(repr(name))"
        throw(ArgumentError("Cannot serialize manifest value$(label) of type $(typeof(x))."))
    end
end

function manifest_ordered_pairs(x)
    x === nothing && return nothing
    if x isa NamedTuple || x isa AbstractDict
        return Any[
            Dict{String,Any}("name" => string(k), "value" => manifest_value(v, k))
            for (k, v) in pairs(x)
        ]
    else
        return manifest_value(x)
    end
end

function parameter_record(spec, value)
    return Dict{String,Any}(
        "shape" => string(spec.shape),
        "axes" => manifest_axis(spec.axes),
        "doc" => spec.doc,
        "value" => manifest_value(value),
    )
end

function parameter_manifest(factory::AbstractBGCFactory, params::NamedTuple, required::Tuple)
    records = Dict{String,Any}()
    for key in required
        spec = parameter_spec(factory, key)
        spec === nothing && throw(
            ArgumentError("Factory $(typeof(factory)) is missing a ParameterSpec for parameter :$key.")
        )
        records[string(key)] = parameter_record(spec, getproperty(params, key))
    end
    return records
end

function parameter_values(params::NamedTuple, required::Tuple)
    return Dict{String,Any}(
        string(key) => manifest_value(getproperty(params, key)) for key in required
    )
end

function plankton_diameter_groups(context)
    return Dict{String,Any}(
        string(group) => Any[manifest_value(context.diameters[i]) for i in indices]
        for (group, indices) in context.group_indices
    )
end

end
