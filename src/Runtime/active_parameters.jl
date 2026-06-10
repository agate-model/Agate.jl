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

@inline function tendency_inputs(bgc_p::ParameterizedBGC, args)
    tracer_values = TracerValues(bgc_p.bgc.tracers, args)
    return bgc_p.parameters, tracer_values
end

@inline Adapt.adapt_structure(to, a::ActiveParameterArray) =
    ActiveParameterArray(Adapt.adapt(to, a.base), Adapt.adapt(to, a.p), Adapt.adapt(to, a.slots))

@inline Adapt.adapt_structure(to, ap::ActiveParameters) =
    ActiveParameters(Adapt.adapt(to, ap.base), Adapt.adapt(to, ap.p), Adapt.adapt(to, ap.map))

@inline Adapt.adapt_structure(to, bgc_p::ParameterizedBGC) =
    ParameterizedBGC(Adapt.adapt(to, bgc_p.bgc), Adapt.adapt(to, bgc_p.parameters))


@inline function evaluate_tendency(bgc, parameters, ::Val{tracer}, args...) where {tracer}
    f = getfield(bgc.tracer_functions, tracer)
    return f(ParameterizedBGC(bgc, parameters), args...)
end

@inline function (bgc_p::ParameterizedBGC)(val_name::Val, args...)
    return evaluate_tendency(bgc_p.bgc, bgc_p.parameters, val_name, args...)
end

"""Selected active parameters for a BGC.

Returned by [`active_parameters`](@ref). `map` is the internal slot map used by
[`parameterized`](@ref) and [`ode_problem`](@ref), `labels` names the flat-vector
entries, and `values` stores the corresponding values from the BGC used to
create the set.
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
    maximum_growth_rate = (:P1, :P2),
    detritus_remineralization = true,
    interactions = (;
        palatability = ((:Z1, :P1), (:Z1, :P2)),
    ),
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
        hasproperty(container, name) || throw(ArgumentError("Unknown active parameter path: $(path_label((path..., name)))."))
        value = getproperty(container, name)
        slots = active_parameter_entry!(labels, values, bgc, (path..., name), value, selection, active_index)
        push!(entries, name => slots)
    end

    return (; entries...)
end

function active_parameter_entry!(labels, values, bgc, path::Tuple, value, selected::Bool, active_index)
    selected || throw(ArgumentError("Boolean active parameter selections must be true."))
    validate_runtime_active_parameter(path)
    value isa Number || throw(ArgumentError(
        "Scalar active selector for $(path_label(path)) is not supported because the stored parameter is not scalar."
    ))
    return push_active_value!(labels, values, active_index, path_label(path), value)
end

function active_parameter_entry!(labels, values, bgc, path::Tuple, value, selection::NamedTuple, active_index)
    validate_runtime_active_parameter(path)
    return active_parameter_map!(labels, values, bgc, path, value, selection, active_index)
end

function active_parameter_entry!(labels, values, bgc, path::Tuple, value, selection::Tuple, active_index)
    validate_runtime_active_parameter(path)
    isempty(selection) && return ()

    if all(item -> item isa Symbol, selection)
        return vector_active_parameter_entry!(labels, values, bgc, path, value, selection, active_index)
    elseif all(is_pair_selection, selection)
        return matrix_active_parameter_entry!(labels, values, bgc, path, value, selection, active_index)
    end

    throw(ArgumentError(
        "Active parameter tuple selections must contain tracer symbols, such as (:P1, :P2), " *
        "or row-column tracer pairs, such as ((:Z1, :P1), (:Z1, :P2))."
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
        index = plankton_parameter_index(bgc, tracer)
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
        indices = interaction_parameter_indices(bgc, row, column)
        push_active_slot!(entries, labels, values, active_index, "$(path_label(path))[$row, $column]", value[indices...], indices)
    end

    return Tuple(entries)
end

const CONSTRUCTOR_DERIVED_ACTIVE_PARAMETERS = (
    :palatability_matrix => ":interactions.palatability",
    :assimilation_matrix => ":interactions.assimilation",
    :specificity => ":interactions.palatability",
    :protection => ":interactions.palatability",
    :optimum_predator_prey_ratio => ":interactions.palatability",
    :assimilation_efficiency => ":interactions.assimilation",
)

function validate_runtime_active_parameter(path::Tuple)
    length(path) == 1 || return nothing

    name = only(path)
    for (derived, runtime_path) in CONSTRUCTOR_DERIVED_ACTIVE_PARAMETERS
        name === derived || continue
        throw(ArgumentError(
            "Active parameter :$name is not currently supported because it is used to derive " *
            "the runtime parameter $runtime_path during model construction. Select $runtime_path instead."
        ))
    end

    return nothing
end

function plankton_parameter_index(bgc, tracer::Symbol)
    tracer_names = Tuple(keys(bgc.tracer_functions))
    pos = findfirst(==(tracer), tracer_names)
    pos === nothing && throw(ArgumentError("Unknown tracer :$tracer."))

    plankton_base = bgc.tracers.idx.plankton_base
    plankton_base == 0 && throw(ArgumentError("Model has no plankton parameter axis."))

    parameter_index = pos - plankton_base + 1
    parameter_index >= 1 || throw(ArgumentError("Tracer :$tracer is not a plankton tracer."))
    return parameter_index
end

function interaction_parameter_indices(bgc, consumer::Symbol, prey::Symbol)
    hasproperty(bgc.parameters, :interactions) || throw(ArgumentError("Model has no interaction matrix axes."))

    interactions = bgc.parameters.interactions
    consumer_global = plankton_parameter_index(bgc, consumer)
    prey_global = plankton_parameter_index(bgc, prey)

    consumer_index = interactions.global_to_consumer[consumer_global]
    prey_index = interactions.global_to_prey[prey_global]

    consumer_index > 0 || throw(ArgumentError("Tracer :$consumer is not on the consumer axis."))
    prey_index > 0 || throw(ArgumentError("Tracer :$prey is not on the prey axis."))

    return (consumer_index, prey_index)
end
