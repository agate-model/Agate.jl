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
`component_diameters` records class diameters for any structured component.
"""
struct ComponentLayout{TR,CT,CI,CD}
    tracer_order::TR
    component_tracers::CT
    component_indices::CI
    component_diameters::CD
end

function _realized_component(
    component::Union{Population,Pool}, name::Symbol, ::Type{T}
) where {T<:Real}
    structure = size_structure(component)
    isnothing(structure) && return ((name,), nothing)

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

function _component_value(layout::ComponentLayout, field::Symbol, component::Symbol)
    values = getfield(layout, field)
    hasproperty(values, component) || throw(ArgumentError("unknown component :$component"))
    return getproperty(values, component)
end

"""Return the concrete tracer identities realized by one logical component."""
component_tracers(layout::ComponentLayout, component::Symbol) =
    _component_value(layout, :component_tracers, component)

"""Return the positions of one logical component in `layout.tracer_order`."""
component_indices(layout::ComponentLayout, component::Symbol) =
    _component_value(layout, :component_indices, component)

"""Return class diameters for one component, or `nothing` for scalar state."""
component_diameters(layout::ComponentLayout, component::Symbol) =
    _component_value(layout, :component_diameters, component)

"""Return the number of concrete classes realized by one logical component."""
component_class_count(layout::ComponentLayout, component::Symbol) =
    length(component_tracers(layout, component))

"""Resolve intrinsic component sinking metadata onto concrete tracer classes.

A tuple supplies one value per realized class. Any other non-`nothing` value is
broadcast to every class of the component. The returned `NamedTuple` is keyed by
concrete tracer identity and remains setup-time structural metadata.
"""
function realize_component_sinking(components::NamedTuple, layout::ComponentLayout)
    entries = Pair{Symbol,Any}[]
    for (name, component) in pairs(components)
        value = sinking(component)
        isnothing(value) && continue
        tracers = component_tracers(layout, name)
        values = if value isa Tuple
            length(value) == length(tracers) || throw(
                ArgumentError(
                    "component :$name sinking tuple must have one value per realized class"
                ),
            )
            value
        else
            ntuple(_ -> value, length(tracers))
        end
        append!(entries, Pair.(tracers, values))
    end
    return (; entries...)
end

"""Realize a named component collection into a deterministic tracer layout.

Logical components are visited in `NamedTuple` key order. Scalar components use
their logical identity directly; any component with a size structure expands to
`name_1`, `name_2`, ... according to its normalized size structure.
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

    component_tracers_values = NamedTuple{names}(Tuple(tracer_values))
    component_indices_values = NamedTuple{names}(Tuple(index_values))
    component_diameters_values = NamedTuple{names}(Tuple(diameter_values))

    return ComponentLayout(
        Tuple(tracer_order),
        component_tracers_values,
        component_indices_values,
        component_diameters_values,
    )
end

function _population_names(components::NamedTuple)
    return Tuple(name for name in keys(components) if getproperty(components, name) isa Population)
end

function _pool_names(components::NamedTuple)
    return Tuple(name for name in keys(components) if getproperty(components, name) isa Pool)
end

function _validate_population_groups(
    components::NamedTuple, population_groups::NamedTuple, context::CommunityContext
)
    population_names = _population_names(components)
    Set(keys(population_groups)) == Set(population_names) || throw(
        ArgumentError(
            "population_groups must map exactly the logical population components $(population_names)"
        ),
    )

    assigned_groups = Symbol[]
    for population in population_names
        groups = getproperty(population_groups, population)
        groups isa Tuple || throw(
            ArgumentError("population component :$population subgroup realization must be a tuple"),
        )
        isempty(groups) && throw(
            ArgumentError("population component :$population must realize at least one subgroup"),
        )
        for group in groups
            group isa Symbol || throw(
                ArgumentError("population component :$population subgroup identities must be Symbols"),
            )
            group in assigned_groups && throw(
                ArgumentError("community subgroup :$group is assigned to more than one population component"),
            )
            haskey(context.group_indices, group) || throw(
                ArgumentError("population component :$population references unknown subgroup :$group"),
            )
            push!(assigned_groups, group)
        end
    end
    Set(assigned_groups) == Set(keys(context.group_indices)) || throw(
        ArgumentError("every runtime community subgroup must realize exactly one population component"),
    )
    return population_names
end

"""Realize logical components over named population subgroups.

Population components are realized from `CommunityContext` subgroup classes.
Pool components are realized directly from their intrinsic component structure,
including size-structured pools. Pool tracers precede population-community
tracers in the runtime layout, preserving the established scalar-pool ordering.
"""
function realize_component_groups(
    components::NamedTuple,
    population_groups::NamedTuple,
    context::CommunityContext,
)
    names = keys(components)
    _validate_population_groups(components, population_groups, context)
    pool_names = _pool_names(components)

    pool_components = NamedTuple{pool_names}(
        Tuple(getproperty(components, name) for name in pool_names)
    )
    pool_layout = realize_components(pool_components; scalar_type=context.scalar_type)
    conflicting = filter(tracer -> tracer in context.plankton_symbols, pool_layout.tracer_order)
    isempty(conflicting) || throw(
        ArgumentError(
            "pool tracer names conflict with realized population tracers: $(unique(conflicting))"
        ),
    )
    tracer_order = (pool_layout.tracer_order..., Tuple(context.plankton_symbols)...)

    tracer_values = ntuple(length(names)) do i
        name = names[i]
        component = getproperty(components, name)
        component isa Pool && return component_tracers(pool_layout, name)

        groups = getproperty(population_groups, name)
        tracers = ()
        for group in groups
            indices = context.group_indices[group]
            tracers = (tracers..., Tuple(context.plankton_symbols[indices])...)
        end
        return tracers
    end
    component_tracers_values = NamedTuple{names}(tracer_values)

    index_values = ntuple(length(names)) do i
        tracers = tracer_values[i]
        indices = Tuple(findfirst(==(tracer), tracer_order) for tracer in tracers)
        any(isnothing, indices) && throw(
            ArgumentError("component :$(names[i]) realizes a tracer absent from the runtime layout"),
        )
        return Tuple(Int(index) for index in indices)
    end
    component_indices_values = NamedTuple{names}(index_values)

    diameter_values = ntuple(length(names)) do i
        name = names[i]
        component = getproperty(components, name)
        component isa Pool && return component_diameters(pool_layout, name)

        tracers = tracer_values[i]
        indices = Tuple(findfirst(==(tracer), context.plankton_symbols) for tracer in tracers)
        any(isnothing, indices) && throw(
            ArgumentError("component :$name realizes a tracer absent from the population layout"),
        )
        return Tuple(context.diameters[Int(index)] for index in indices)
    end
    component_diameters_values = NamedTuple{names}(diameter_values)

    return ComponentLayout(
        tracer_order,
        component_tracers_values,
        component_indices_values,
        component_diameters_values,
    )
end
