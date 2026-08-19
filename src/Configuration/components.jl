"""Population state described by intrinsic component properties.

Ecological function is supplied by process participation. `currency` identifies
what conserved material the state represents, while `size_structure` and
`sinking` describe intrinsic realization properties.
"""
struct Population{C,S,V}
    currency::C
    size_structure::S
    sinking::V
end

function Population(; currency, size_structure=nothing, sinking=nothing)
    isnothing(currency) && throw(ArgumentError("Population currency must be specified."))
    return Population(currency, size_structure, sinking)
end

"""Material-pool state described by intrinsic component properties."""
struct Pool{C,S,V}
    currency::C
    size_structure::S
    sinking::V
end

function Pool(currency; size_structure=nothing, sinking=nothing)
    isnothing(currency) && throw(ArgumentError("Pool currency must be specified."))
    return Pool(currency, size_structure, sinking)
end

Pool(; currency, size_structure=nothing, sinking=nothing) =
    Pool(currency; size_structure, sinking)

"""Return the conserved-material currency represented by `component`."""
@inline currency(component::Union{Population,Pool}) = component.currency

"""Return the intrinsic size-structure specification for `component`."""
@inline size_structure(component::Union{Population,Pool}) = component.size_structure

"""Return the intrinsic sinking configuration for `component`."""
@inline sinking(component::Union{Population,Pool}) = component.sinking

"""Setup-time realization of logical components into concrete tracer identities.

`component_tracers` maps each logical component to its realized tracer tuple,
`component_indices` maps each logical component into `tracer_order`, and
`component_diameters` records class diameters for structured populations.
"""
struct ComponentLayout{TR,CT,CI,CD}
    tracer_order::TR
    component_tracers::CT
    component_indices::CI
    component_diameters::CD
end

function _realized_component(component::Population, name::Symbol, ::Type{T}) where {T<:Real}
    structure = size_structure(component)
    if isnothing(structure)
        return ((name,), nothing)
    end

    normalized = normalize_diameters(structure)
    specification = normalized.specification
    n = normalized.n
    isnothing(n) && throw(ArgumentError("component :$name size structure must define `n`."))
    n isa Integer && n > 0 ||
        throw(ArgumentError("component :$name size structure must define a positive `n`."))

    diameters = Tuple(param_compute_diameters(T, n, specification))
    tracers = ntuple(i -> Symbol(string(name), "_", i), n)
    return tracers, diameters
end

function _realized_component(component::Pool, name::Symbol, ::Type{T}) where {T<:Real}
    isnothing(size_structure(component)) || throw(
        ArgumentError(
            "component :$name is a structured Pool; the current component layout realizes scalar Pool state.",
        ),
    )
    return ((name,), nothing)
end

"""Realize a named component collection into a deterministic tracer layout.

Logical components are visited in `NamedTuple` key order. Scalar components use
their logical identity directly; structured populations expand to `name_1`,
`name_2`, ... according to their normalized size structure.
"""
function realize_components(components::NamedTuple; scalar_type::Type{T}=Float64) where {T<:Real}
    names = keys(components)
    n_components = length(names)

    tracer_values = Vector{Any}(undef, n_components)
    index_values = Vector{Any}(undef, n_components)
    diameter_values = Vector{Any}(undef, n_components)
    tracer_order = Symbol[]

    for i in 1:n_components
        name = names[i]
        component = getproperty(components, name)
        component isa Union{Population,Pool} || throw(
            ArgumentError(
                "component :$name must be a Population or Pool; got $(typeof(component)).",
            ),
        )

        tracers, diameters = _realized_component(component, name, T)
        conflicting = filter(tracer -> tracer in tracer_order, tracers)
        isempty(conflicting) || throw(
            ArgumentError("component :$name realizes duplicate tracer identities $(conflicting)."),
        )

        first_index = length(tracer_order) + 1
        append!(tracer_order, tracers)
        last_index = length(tracer_order)

        tracer_values[i] = tracers
        index_values[i] = Tuple(first_index:last_index)
        diameter_values[i] = diameters
    end

    component_tracers = NamedTuple{names}(Tuple(tracer_values))
    component_indices = NamedTuple{names}(Tuple(index_values))
    component_diameters = NamedTuple{names}(Tuple(diameter_values))

    return ComponentLayout(
        Tuple(tracer_order), component_tracers, component_indices, component_diameters
    )
end
