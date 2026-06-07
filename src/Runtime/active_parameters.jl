"""Internal vector view for one parameter array."""
struct ActiveParameterVector{B,P,S}
    base::B
    p::P
    slots::S
end

@inline Base.length(v::ActiveParameterVector) = length(v.base)
@inline Base.axes(v::ActiveParameterVector) = axes(v.base)
@inline Base.eachindex(v::ActiveParameterVector) = eachindex(v.base)
@inline Base.IndexStyle(::Type{<:ActiveParameterVector}) = IndexLinear()

@inline function Base.eltype(::Type{<:ActiveParameterVector{B,P}}) where {B,P}
    return promote_type(eltype(B), eltype(P))
end

@inline function Base.getindex(v::ActiveParameterVector, i::Int)
    @inbounds for slot in v.slots
        slot.parameter_index == i && return v.p[slot.active_index]
    end
    return v.base[i]
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
    map = getfield(ap, :map)

    if hasproperty(map, name)
        return ActiveParameterVector(getproperty(base, name), getfield(ap, :p), getproperty(map, name))
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

@inline function (bgc_p::ParameterizedBGC)(::Val{tracer}, args...) where {tracer}
    f = getfield(bgc_p.bgc.tracer_functions, tracer)
    return f(bgc_p, args...)
end

"""Return a BGC wrapper whose selected parameters are read from `p`.

`active_parameters` maps parameter names to tracer names and active-vector
indices. For example:

```julia
parameterized(bgc, p; active_parameters = (; maximum_growth_rate = (; P1 = 1)))
```
"""
function parameterized(bgc, p; active_parameters=(;))
    map = active_parameter_map(bgc, active_parameters)
    parameters = ActiveParameters(bgc.parameters, p, map)
    return ParameterizedBGC(bgc, parameters)
end

function active_parameter_map(bgc, selectors::NamedTuple)
    entries = Pair{Symbol,Any}[]
    for (parameter_name, selector) in pairs(selectors)
        push!(entries, parameter_name => active_parameter_slots(bgc, selector))
    end
    return (; entries...)
end

function active_parameter_slots(bgc, selector::NamedTuple)
    entries = map(pairs(selector)) do pair
        tracer, active_index = pair
        (; parameter_index=plankton_parameter_index(bgc, tracer), active_index=Int(active_index))
    end
    return Tuple(entries)
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

function with_parameters end

function with_active_parameters(bgc, p; kwargs...)
    return parameterized(bgc, p; active_parameters=(; kwargs...))
end
