import Adapt

"""Internal array view for one active parameter leaf.

`base` is the stored parameter value and `p` is the external active parameter
vector. `slots` maps indices in `base` to entries in `p`.
"""
struct ActiveParameterArray{B,P,S}
    base::B
    p::P
    slots::S
end

@inline Base.length(a::ActiveParameterArray) = length(a.base)
@inline Base.size(a::ActiveParameterArray) = size(a.base)
@inline Base.axes(a::ActiveParameterArray) = axes(a.base)
@inline Base.eachindex(a::ActiveParameterArray) = eachindex(a.base)
@inline Base.IndexStyle(::Type{<:ActiveParameterArray}) = IndexLinear()

@inline function Base.eltype(::Type{<:ActiveParameterArray{B,P}}) where {B,P}
    return promote_type(eltype(B), eltype(P))
end

@inline function Base.getindex(a::ActiveParameterArray, indices::Vararg{Int,N}) where {N}
    @inbounds for slot in a.slots
        slot.indices == indices && return a.p[slot.active_index]
    end
    return a.base[indices...]
end

"""Internal parameter container that overrides selected parameter fields from `p`."""
struct ActiveParameters{B,P,M}
    base::B
    p::P
    map::M
end

@inline Base.propertynames(ap::ActiveParameters) = propertynames(ap.base)

@inline function Base.getproperty(ap::ActiveParameters, name::Symbol)
    name === :base && return getfield(ap, :base)
    name === :p && return getfield(ap, :p)
    name === :map && return getfield(ap, :map)

    base = getfield(ap, :base)
    p = getfield(ap, :p)
    map = getfield(ap, :map)

    if hasproperty(map, name)
        selector = getproperty(map, name)
        value = getproperty(base, name)
        selector isa Integer && return p[selector]
        selector isa NamedTuple && return ActiveParameters(value, p, selector)
        return ActiveParameterArray(value, p, selector)
    end

    return getproperty(base, name)
end

"""Oceananigans/OceanBioME-compatible BGC wrapper with external active parameters."""
struct ParameterizedBGC{B,P} <: AbstractContinuousFormBiogeochemistry
    bgc::B
    parameters::P
end

@inline Base.getproperty(bgc_p::ParameterizedBGC, name::Symbol) = begin
    name === :bgc && return getfield(bgc_p, :bgc)
    name === :parameters && return getfield(bgc_p, :parameters)
    return getproperty(getfield(bgc_p, :bgc), name)
end

@inline Adapt.adapt_structure(to, a::ActiveParameterArray) =
    ActiveParameterArray(Adapt.adapt(to, a.base), Adapt.adapt(to, a.p), Adapt.adapt(to, a.slots))

@inline Adapt.adapt_structure(to, ap::ActiveParameters) =
    ActiveParameters(Adapt.adapt(to, ap.base), Adapt.adapt(to, ap.p), Adapt.adapt(to, ap.map))

@inline Adapt.adapt_structure(to, bgc_p::ParameterizedBGC) =
    ParameterizedBGC(Adapt.adapt(to, bgc_p.bgc), Adapt.adapt(to, bgc_p.parameters))


@inline function evaluate_tendency(bgc_p::ParameterizedBGC, ::Val{tracer}, args...) where {tracer}
    equation = getfield(bgc_p.bgc.equations, tracer)
    return equation(bgc_p, args...)
end

@inline function evaluate_tendency(bgc, ::Val{tracer}, args...) where {tracer}
    equation = getfield(bgc.equations, tracer)
    return equation(bgc, args...)
end

@inline function (bgc_p::ParameterizedBGC)(val_name::Val, args...)
    return evaluate_tendency(bgc_p, val_name, args...)
end

"""Selected active parameters for a BGC.

Returned by [`active_parameters`](@ref). `labels` names the flat-vector entries,
and `values` stores the corresponding values from the BGC used to create the set.
"""
struct ActiveParameterSet{M,V}
    map::M
    labels::Tuple{Vararg{String}}
    values::V
end

Base.length(active::ActiveParameterSet) = length(active.values)
Base.isempty(active::ActiveParameterSet) = isempty(active.values)


"""Return selected BGC parameters as a labelled flat-vector parameter set.

The keyword arguments follow the structure of `bgc.parameters`. Scalar
parameters are selected with `true`, vector parameters with tracer symbols, and
matrix parameters with row-column tracer pairs.

```julia
active = active_parameters(
    bgc;
    maximum_growth_rate = (:P_1, :P_2),
    detritus_remineralization = true,
    palatability_matrix = ((:Z_1, :P_1), (:Z_1, :P_2)),
)

θ = copy(active.values)
parameterized(bgc, θ; active_parameters = active)
```
"""
function active_parameters(bgc; kwargs...)
    active_index = Ref(1)
    labels = String[]
    values = Any[]
    map = active_parameter_map!(labels, values, bgc, (), bgc.parameters, (; kwargs...), active_index)
    active_values = isempty(values) ? Float64[] : collect(promote(values...))
    return ActiveParameterSet(map, Tuple(labels), active_values)
end

"""
    parameterized(bgc, p; active_parameters=nothing)

Return an Oceananigans/OceanBioME-compatible BGC wrapper whose selected
runtime parameters are read from the flat vector `p`.

This keeps the BGC structure fixed while allowing AD backends and external
solvers to differentiate with respect to `p`. When `active_parameters` is
provided, it should be the [`ActiveParameterSet`](@ref) returned by
[`active_parameters`](@ref).
"""
function parameterized(bgc, p; active_parameters=nothing)
    map = active_parameters === nothing ? (;) : active_parameters.map
    parameters = ActiveParameters(bgc.parameters, p, map)
    return ParameterizedBGC(bgc, parameters)
end

function active_parameter_map!(labels, values, bgc, path::Tuple, container, selections::NamedTuple, active_index)
    entries = Pair{Symbol, Any}[]

    for (name, selection) in pairs(selections)
        entry_path = (path..., name)
        validate_runtime_active_parameter(bgc, entry_path)
        hasproperty(container, name) || throw(ArgumentError(
            "Runtime parameter storage is missing $(path_label(entry_path))."
        ))
        value = getproperty(container, name)
        slots = active_parameter_entry!(labels, values, bgc, entry_path, value, selection, active_index)
        push!(entries, name => slots)
    end

    return (; entries...)
end

function active_parameter_entry!(labels, values, bgc, path::Tuple, value, selected::Bool, active_index)
    selected || throw(ArgumentError("Boolean active parameter selections must be true."))
    value isa Number || throw(ArgumentError(
        "Scalar active selector for $(path_label(path)) is not supported because the stored parameter is not scalar."
    ))
    return push_active_value!(labels, values, active_index, path_label(path), value)
end

function active_parameter_entry!(labels, values, bgc, path::Tuple, value, selection::NamedTuple, active_index)
    return active_parameter_map!(labels, values, bgc, path, value, selection, active_index)
end

function active_parameter_entry!(labels, values, bgc, path::Tuple, value, selection::Tuple, active_index)
    isempty(selection) && return ()

    if all(item -> item isa Symbol, selection)
        return vector_active_parameter_entry!(labels, values, bgc, path, value, selection, active_index)
    elseif all(is_pair_selection, selection)
        return matrix_active_parameter_entry!(labels, values, bgc, path, value, selection, active_index)
    end

    throw(ArgumentError(
        "Active parameter tuple selections must contain tracer symbols, such as (:P_1, :P_2), " *
        "or row-column tracer pairs, such as ((:Z_1, :P_1), (:Z_1, :P_2))."
    ))
end

function push_active_value!(labels, values, active_index, label, value)
    index = active_index[]
    active_index[] += 1
    push!(labels, label)
    push!(values, value)
    return index
end

function push_active_slot!(entries, labels, values, active_index, label, value, indices)
    index = push_active_value!(labels, values, active_index, label, value)
    push!(entries, (; indices, active_index=index))
    return nothing
end

is_pair_selection(item) = item isa Tuple && length(item) == 2 && item[1] isa Symbol && item[2] isa Symbol
path_label(path::Tuple) = join(string.(path), ".")

function vector_active_parameter_entry!(labels, values, bgc, path::Tuple, value, selection::Tuple, active_index)
    value isa AbstractVector || throw(ArgumentError(
        "Vector active selector for $(path_label(path)) is not supported because the stored parameter is not a vector."
    ))

    entries = []
    for tracer in selection
        index = parameter_label_index(bgc, only(path), 1, tracer)
        push_active_slot!(entries, labels, values, active_index, "$(path_label(path)).$(tracer)", value[index], (index,))
    end

    return Tuple(entries)
end

function matrix_active_parameter_entry!(labels, values, bgc, path::Tuple, value, selection::Tuple, active_index)
    value isa AbstractMatrix || throw(ArgumentError(
        "Matrix active selector for $(path_label(path)) is not supported because the stored parameter is not a matrix."
    ))

    entries = []
    for (row, column) in selection
        indices = (parameter_label_index(bgc, only(path), 1, row), parameter_label_index(bgc, only(path), 2, column))
        push_active_slot!(entries, labels, values, active_index, "$(path_label(path))[$row, $column]", value[indices...], indices)
    end

    return Tuple(entries)
end

function _parameter_metadata(bgc, name::Symbol)
    metadata = getproperty(bgc, :metadata)
    metadata === nothing && throw(ArgumentError("Model has no parameter metadata."))
    axes = metadata.parameter_axes
    hasproperty(axes, name) || throw(
        ArgumentError("Unknown active parameter path: $name."),
    )
    return getproperty(axes, name)
end

function validate_runtime_active_parameter(bgc, path::Tuple)
    length(path) == 1 || return nothing
    name = only(path)
    metadata = _parameter_metadata(bgc, name)
    metadata.runtime_bound && return nothing

    targets = metadata.derived_runtime_parameters
    detail = isempty(targets) ? "" :
        " It is used to derive runtime parameter" * (length(targets) == 1 ? " " : "s ") *
        join(":" .* string.(targets), ", ") * "."
    throw(ArgumentError(
        "Active parameter :$name is construction-only and is not stored in the runtime BGC." * detail,
    ))
end

function parameter_label_index(bgc, parameter::Symbol, dimension::Int, label::Symbol)
    metadata = _parameter_metadata(bgc, parameter)
    1 <= dimension <= length(metadata.labels) || throw(ArgumentError(
        "parameter :$parameter has no storage dimension $dimension",
    ))
    labels = metadata.labels[dimension]
    index = findfirst(==(label), labels)
    index === nothing && throw(ArgumentError(
        "Unknown class :$label for parameter :$parameter; expected one of $(labels).",
    ))
    return index
end
