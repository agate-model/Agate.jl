module Serialization

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
        throw(ArgumentError("Cannot serialize model setup value$(label) of type $(typeof(x))."))
    end
end

function manifest_ordered_pairs(x)
    x === nothing && return nothing
    return Any[
        Dict{String,Any}("name" => string(k), "value" => manifest_value(v, k))
        for (k, v) in pairs(x)
    ]
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
