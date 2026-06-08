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

"""Internal active view of the runtime interaction matrices."""
struct ActiveInteractions{B,P,M}
    base::B
    p::P
    map::M
end

@inline Base.propertynames(ai::ActiveInteractions) = propertynames(ai.base)

@inline function Base.getproperty(ai::ActiveInteractions, name::Symbol)
    name === :base && return getfield(ai, :base)
    name === :p && return getfield(ai, :p)
    name === :map && return getfield(ai, :map)

    base = getfield(ai, :base)
    map = getfield(ai, :map)

    if name === :palatability && hasproperty(map, :palatability_matrix)
        return ActiveParameterArray(base.palatability, getfield(ai, :p), map.palatability_matrix)
    elseif name === :assimilation && hasproperty(map, :assimilation_matrix)
        return ActiveParameterArray(base.assimilation, getfield(ai, :p), map.assimilation_matrix)
    end

    return getproperty(base, name)
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

    if name === :interactions
        interactions = getproperty(base, :interactions)
        if hasproperty(map, :palatability_matrix) || hasproperty(map, :assimilation_matrix)
            return ActiveInteractions(interactions, p, map)
        end
        return interactions
    end

    if hasproperty(map, name)
        selector = getproperty(map, name)
        selector isa Integer && return p[selector]
        return ActiveParameterArray(getproperty(base, name), p, selector)
    end

    return getproperty(base, name)
end

"""Internal wrapper for evaluating a static BGC with an alternate parameter object."""
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

@inline Adapt.adapt_structure(to, ai::ActiveInteractions) =
    ActiveInteractions(Adapt.adapt(to, ai.base), Adapt.adapt(to, ai.p), Adapt.adapt(to, ai.map))

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

"""Active parameter layout, labels, and reference values for a BGC.

Returned by [`active_parameters`](@ref). `layout` is the internal map from
structured parameter names to positions in a flat vector, `labels` names those
entries, and `values` stores the corresponding values from the BGC used to
create the set.
"""
struct ActiveParameterSet{L,V}
    layout::L
    labels::Tuple{Vararg{String}}
    values::V
end

Base.length(active::ActiveParameterSet) = length(active.values)
Base.isempty(active::ActiveParameterSet) = isempty(active.values)

active_parameter_layout(active::ActiveParameterSet) = active.layout
active_parameter_layout(::Nothing) = (;)

"""Return selected BGC parameters as a labelled flat-vector parameter set.

The keyword arguments name active parameter leaves without exposing vector
indices. Scalar parameters are selected with `true`, vector parameters with
tracer symbols, and supported interaction matrices with consumer-prey pairs.

```julia
active = active_parameters(
    bgc;
    maximum_growth_rate = (:P1, :P2),
    detritus_remineralization = true,
    palatability_matrix = ((:Z1, :P1), (:Z1, :P2)),
)

θ = copy(active.values)
parameterized(bgc, θ; active_parameters = active)
```
"""
function active_parameters(bgc; kwargs...)
    active_index = Ref(1)
    layout_entries = Pair{Symbol, Any}[]
    labels = String[]
    values = Any[]

    for (parameter_name, selection) in pairs((; kwargs...))
        selector = active_parameter_entry!(labels, values, bgc, parameter_name, selection, active_index)
        push!(layout_entries, parameter_name => selector)
    end

    active_values = isempty(values) ? Float64[] : collect(promote(values...))
    return ActiveParameterSet((; layout_entries...), Tuple(labels), active_values)
end

"""Return a BGC wrapper whose selected parameters are read from `p`."""
function parameterized(bgc, p; active_parameters=nothing)
    layout = active_parameter_layout(active_parameters)
    map = active_parameter_map(bgc, layout)
    parameters = ActiveParameters(bgc.parameters, p, map)
    return ParameterizedBGC(bgc, parameters)
end

function active_parameter_entry!(labels, values, bgc, parameter_name::Symbol, selected::Bool, active_index)
    selected || throw(ArgumentError("Boolean active parameter selections must be true."))
    validate_runtime_active_parameter(bgc, parameter_name)

    parameter_value = getproperty(bgc.parameters, parameter_name)
    parameter_value isa Number || throw(ArgumentError(
        "Scalar active selector for :$parameter_name is not supported because the stored parameter is not scalar. " *
        "Use tracer symbols for vector parameters, or consumer-prey pairs for supported interaction matrices."
    ))

    push!(labels, String(parameter_name))
    push!(values, parameter_value)
    return next_active_index!(active_index)
end

function active_parameter_entry!(labels, values, bgc, parameter_name::Symbol, selection::Tuple, active_index)
    isempty(selection) && return (;)

    if all(item -> item isa Symbol, selection)
        return vector_active_parameter_entry!(labels, values, bgc, parameter_name, selection, active_index)
    elseif all(is_pair_selection, selection)
        return matrix_active_parameter_entry!(labels, values, bgc, parameter_name, selection, active_index)
    end

    throw(ArgumentError(
        "Active parameter tuple selections must contain tracer symbols, such as (:P1, :P2), " *
        "or consumer-prey pairs, such as ((:Z1, :P1), (:Z1, :P2))."
    ))
end

function next_active_index!(active_index::Base.RefValue{Int})
    index = active_index[]
    active_index[] += 1
    return index
end

is_pair_selection(item) = item isa Tuple && length(item) == 2 && item[1] isa Symbol && item[2] isa Symbol

function vector_active_parameter_entry!(labels, values, bgc, parameter_name::Symbol, selection::Tuple, active_index)
    validate_runtime_active_parameter(bgc, parameter_name)
    is_interaction_matrix_parameter(parameter_name) && throw(ArgumentError(
        "Vector active selector for :$parameter_name is not supported. " *
        "Use consumer-prey pairs, for example ((:Z1, :P1),)."
    ))

    parameter_value = getproperty(bgc.parameters, parameter_name)
    parameter_value isa AbstractVector || throw(ArgumentError(
        "Vector active selector for :$parameter_name is not supported because the stored parameter is not a vector."
    ))

    entries = Pair{Symbol, Int}[]
    for tracer in selection
        index = next_active_index!(active_index)
        push!(entries, tracer => index)
        push!(labels, "$(parameter_name).$(tracer)")
        push!(values, parameter_value[plankton_parameter_index(bgc, tracer)])
    end

    return (; entries...)
end

function matrix_active_parameter_entry!(labels, values, bgc, parameter_name::Symbol, selection::Tuple, active_index)
    matrix = active_matrix_parameter_value(bgc, parameter_name)
    consumers = Symbol[]
    consumer_entries = Pair{Symbol, Vector{Pair{Symbol, Int}}}[]

    for (consumer, prey) in selection
        index = next_active_index!(active_index)
        pos = findfirst(==(consumer), consumers)
        if pos === nothing
            push!(consumers, consumer)
            push!(consumer_entries, consumer => Pair{Symbol, Int}[])
            pos = length(consumer_entries)
        end

        push!(consumer_entries[pos].second, prey => index)
        push!(labels, "$(parameter_name)[$consumer, $prey]")
        push!(values, matrix[interaction_parameter_indices(bgc, consumer, prey)...])
    end

    return (; (consumer => (; prey_entries...) for (consumer, prey_entries) in consumer_entries)...)
end

function active_matrix_parameter_value(bgc, parameter_name::Symbol)
    parameter_name === :palatability_matrix && return bgc.parameters.interactions.palatability
    parameter_name === :assimilation_matrix && return bgc.parameters.interactions.assimilation
    throw(ArgumentError(
        "Matrix active selector for :$parameter_name is not supported. " *
        "Currently supported matrix parameters are :palatability_matrix and :assimilation_matrix."
    ))
end

function active_parameter_map(bgc, selectors::NamedTuple)
    entries = Pair{Symbol, Any}[]
    for (parameter_name, selector) in pairs(selectors)
        push!(entries, parameter_name => active_parameter_slots(bgc, parameter_name, selector))
    end
    return (; entries...)
end

const CONSTRUCTOR_DERIVED_ACTIVE_PARAMETERS = (
    :specificity,
    :protection,
    :optimum_predator_prey_ratio,
    :assimilation_efficiency,
)

function validate_runtime_active_parameter(bgc, parameter_name::Symbol)
    if parameter_name in CONSTRUCTOR_DERIVED_ACTIVE_PARAMETERS
        throw(ArgumentError(
            "Active parameter :$parameter_name is not currently supported because it is used to derive " *
            "runtime interaction matrices during model construction. Select entries of :palatability_matrix " *
            "or :assimilation_matrix as active parameters instead."
        ))
    end

    hasproperty(bgc.parameters, parameter_name) || throw(ArgumentError("Unknown active parameter :$parameter_name."))
    return nothing
end

function active_parameter_index(active_index)
    active_index isa Integer || throw(ArgumentError("Active parameter indices must be integers."))
    index = Int(active_index)
    index > 0 || throw(ArgumentError("Active parameter indices must be positive."))
    return index
end

function active_parameter_slots(bgc, parameter_name::Symbol, selector::Integer)
    validate_runtime_active_parameter(bgc, parameter_name)
    value = getproperty(bgc.parameters, parameter_name)
    value isa Number || throw(ArgumentError(
        "Scalar active selector for :$parameter_name is not supported because the stored parameter is not scalar. " *
        "Use a NamedTuple selector for vector parameters, or nested NamedTuples for supported interaction matrices."
    ))
    return active_parameter_index(selector)
end

function active_parameter_slots(bgc, parameter_name::Symbol, selector::NamedTuple)
    validate_runtime_active_parameter(bgc, parameter_name)

    values_tuple = values(selector)
    isempty(values_tuple) && return ()

    nested = first(values_tuple) isa NamedTuple
    for value in values_tuple
        if nested != (value isa NamedTuple)
            throw(ArgumentError("Active selector for :$parameter_name mixes vector and matrix syntax."))
        end
    end

    return nested ? matrix_active_parameter_slots(bgc, parameter_name, selector) : vector_active_parameter_slots(bgc, parameter_name, selector)
end

function vector_active_parameter_slots(bgc, parameter_name::Symbol, selector::NamedTuple)
    is_interaction_matrix_parameter(parameter_name) && throw(ArgumentError(
        "Vector active selector for :$parameter_name is not supported. " *
        "Use nested NamedTuples, for example (; Z1 = (; P1 = 1))."
    ))

    value = getproperty(bgc.parameters, parameter_name)
    value isa AbstractVector || throw(ArgumentError(
        "Vector active selector for :$parameter_name is not supported because the stored parameter is not a vector."
    ))

    entries = []
    for (tracer, active_index) in pairs(selector)
        push!(entries, (; indices=(plankton_parameter_index(bgc, tracer),), active_index=active_parameter_index(active_index)))
    end
    return Tuple(entries)
end

function matrix_active_parameter_slots(bgc, parameter_name::Symbol, selector::NamedTuple)
    is_interaction_matrix_parameter(parameter_name) || throw(
        ArgumentError(
            "Matrix active selector for :$parameter_name is not supported. " *
            "Currently supported matrix parameters are :palatability_matrix and :assimilation_matrix."
        )
    )

    entries = []
    for (consumer, prey_selector) in pairs(selector)
        prey_selector isa NamedTuple || throw(ArgumentError("Matrix selector for :$parameter_name requires nested NamedTuples."))
        for (prey, active_index) in pairs(prey_selector)
            push!(entries, (; indices=interaction_parameter_indices(bgc, consumer, prey), active_index=active_parameter_index(active_index)))
        end
    end
    return Tuple(entries)
end

@inline is_interaction_matrix_parameter(parameter_name::Symbol) =
    parameter_name === :palatability_matrix || parameter_name === :assimilation_matrix

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

