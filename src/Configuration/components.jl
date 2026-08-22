"""Population state described by intrinsic component properties.

Ecological function is supplied by process participation. `states` maps each
prognostic state identity to the conserved-material currency it represents,
while `size_structure` describes the shared ecological class realization.
A one-state population may be authored with the `currency=` convenience.
"""
struct Population{ST<:NamedTuple,S}
    states::ST
    size_structure::S

    function Population(states::ST, size_structure::S) where {ST<:NamedTuple,S}
        isempty(states) && throw(ArgumentError("Population must define at least one state."))
        all(value -> value isa Symbol, values(states)) || throw(
            ArgumentError("Population state currencies must be Symbols."),
        )
        return new{ST,S}(states, size_structure)
    end
end

function Population(; currency=nothing, states=nothing, size_structure=nothing)
    isnothing(currency) != isnothing(states) || throw(
        ArgumentError("Population requires exactly one of `currency` or `states`."),
    )
    resolved_states = if isnothing(states)
        currency isa Symbol || throw(
            ArgumentError("Population `currency` must be a Symbol; use `states` for named state mappings."),
        )
        NamedTuple{(currency,)}((currency,))
    else
        states
    end
    resolved_states isa NamedTuple || throw(
        ArgumentError("Population states must be a NamedTuple mapping state names to currencies."),
    )
    return Population(resolved_states, size_structure)
end

"""Material-pool state described by intrinsic component properties."""
struct Pool{C,S}
    currency::C
    size_structure::S
end

function Pool(currency; size_structure=nothing)
    isnothing(currency) && throw(ArgumentError("Pool currency must be specified."))
    return Pool(currency, size_structure)
end

Pool(; currency, size_structure=nothing) = Pool(currency; size_structure)


"""Reference one named prognostic state carried by a logical population."""
struct PopulationStateRef
    population::Symbol
    state::Symbol
end

"""Construct a setup-time reference to one population prognostic state."""
@inline population_state(population::Symbol, state::Symbol) = PopulationStateRef(population, state)

"""Return the prognostic state-to-currency mapping for a population."""
@inline states(component::Population) = component.states

"""Return the currency represented by one named population state."""
function state_currency(component::Population, state::Symbol)
    hasproperty(component.states, state) || throw(
        ArgumentError("Population has no state :$state."),
    )
    return getproperty(component.states, state)
end

"""Return the conserved-material currency represented by a single-state component.

Multi-state populations have no unique component currency; use
`state_currency(population, state)` instead.
"""
@inline currency(component::Pool) = component.currency
function currency(component::Population)
    length(component.states) == 1 || throw(
        ArgumentError(
            "multi-state Population has no unique currency; use state_currency(population, state)",
        ),
    )
    return only(values(component.states))
end

"""Return the intrinsic size-structure specification for `component`."""
@inline size_structure(component::Union{Population,Pool}) = component.size_structure

"""Setup-time realization of logical components into ecological classes and tracers.

`component_classes` maps each logical component to its ecological class identities.
For pools these are also the physical tracer identities. Population
`component_state_tracers` map each prognostic state onto one tracer per ecological
class; `component_tracers` is their flattened physical tracer tuple.
`component_indices` index flattened physical tracers in `tracer_order`; state-specific
positions are derived from the state tracer mapping. `component_diameters` records one
diameter per ecological class.
"""
struct ComponentLayout{TR,CC,CST,CT,CI,CD}
    tracer_order::TR
    component_classes::CC
    component_state_tracers::CST
    component_tracers::CT
    component_indices::CI
    component_diameters::CD
end

function _realized_classes(
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
    classes = ntuple(i -> Symbol(string(name), "_", i), n)
    return classes, diameters
end

@inline function _population_state_tracer(class::Symbol, state::Symbol, nstates::Int)
    return nstates == 1 ? class : Symbol(string(class), "_", state)
end

function _population_state_tracers(component::Population, classes::Tuple)
    state_names = keys(states(component))
    nstates = length(state_names)
    values = ntuple(length(state_names)) do i
        state = state_names[i]
        Tuple(_population_state_tracer(class, state, nstates) for class in classes)
    end
    return NamedTuple{state_names}(values)
end

function _flatten_population_tracers(
    component::Population, classes::Tuple, state_tracers::NamedTuple
)
    state_names = keys(states(component))
    return Tuple(
        getproperty(state_tracers, state)[class_index]
        for class_index in eachindex(classes) for state in state_names
    )
end

function _realized_component(component::Population, name::Symbol, ::Type{T}) where {T<:Real}
    classes, diameters = _realized_classes(component, name, T)
    state_tracers = _population_state_tracers(component, classes)
    tracers = _flatten_population_tracers(component, classes, state_tracers)
    return (; classes, state_tracers, tracers, diameters)
end

function _realized_component(component::Pool, name::Symbol, ::Type{T}) where {T<:Real}
    classes, diameters = _realized_classes(component, name, T)
    return (; classes, state_tracers=nothing, tracers=classes, diameters)
end

function _component_value(layout::ComponentLayout, field::Symbol, component::Symbol)
    values = getfield(layout, field)
    hasproperty(values, component) || throw(ArgumentError("unknown component :$component"))
    return getproperty(values, component)
end

"""Return ecological class identities realized by one logical component."""
component_classes(layout::ComponentLayout, component::Symbol) =
    _component_value(layout, :component_classes, component)

"""Return the state-to-tracer mapping for a population, or `nothing` for a pool."""
component_state_tracers(layout::ComponentLayout, component::Symbol) =
    _component_value(layout, :component_state_tracers, component)

"""Return the flattened concrete tracer identities realized by one logical component."""
component_tracers(layout::ComponentLayout, component::Symbol) =
    _component_value(layout, :component_tracers, component)

"""Return the positions of one logical component's physical tracers in `tracer_order`."""
component_indices(layout::ComponentLayout, component::Symbol) =
    _component_value(layout, :component_indices, component)

"""Return concrete tracers for one population state."""
function state_tracers(layout::ComponentLayout, component::Symbol, state::Symbol)
    mapping = component_state_tracers(layout, component)
    mapping isa NamedTuple || throw(
        ArgumentError("component :$component does not expose population states"),
    )
    hasproperty(mapping, state) || throw(
        ArgumentError("component :$component has no realized state :$state"),
    )
    return getproperty(mapping, state)
end

state_tracers(layout::ComponentLayout, reference::PopulationStateRef) =
    state_tracers(layout, reference.population, reference.state)

"""Return one concrete tracer for a population state and local ecological class ordinal."""
function state_tracer(
    layout::ComponentLayout, reference::PopulationStateRef, class_ordinal::Integer
)
    tracers = state_tracers(layout, reference)
    1 <= class_ordinal <= length(tracers) || throw(
        ArgumentError(
            "class ordinal $class_ordinal is out of bounds for population :$(reference.population) state :$(reference.state)",
        ),
    )
    return tracers[Int(class_ordinal)]
end

"""Return concrete tracer positions for one population state."""
function state_indices(layout::ComponentLayout, component::Symbol, state::Symbol)
    tracers = state_tracers(layout, component, state)
    indices = Tuple(findfirst(==(tracer), layout.tracer_order) for tracer in tracers)
    any(isnothing, indices) && throw(
        ArgumentError("component :$component state :$state is absent from the runtime layout"),
    )
    return Tuple(Int(index) for index in indices)
end

state_indices(layout::ComponentLayout, reference::PopulationStateRef) =
    state_indices(layout, reference.population, reference.state)

"""Return class diameters for one component, or `nothing` for scalar state."""
component_diameters(layout::ComponentLayout, component::Symbol) =
    _component_value(layout, :component_diameters, component)

"""Return the number of ecological classes realized by one logical component."""
component_class_count(layout::ComponentLayout, component::Symbol) =
    length(component_classes(layout, component))

function _physical_indices(tracers::Tuple, tracer_order::Vector{Symbol}, component::Symbol)
    indices = Tuple(findfirst(==(tracer), tracer_order) for tracer in tracers)
    any(isnothing, indices) && throw(
        ArgumentError("component :$component realizes a tracer absent from the runtime layout"),
    )
    return Tuple(Int(index) for index in indices)
end


function _check_component_identity_conflicts(
    name::Symbol,
    classes::Tuple,
    tracers::Tuple,
    prior_classes::Vector{Symbol},
    prior_tracers::Vector{Symbol},
)
    class_conflicts = filter(class -> class in prior_classes || class in prior_tracers, classes)
    tracer_conflicts = filter(tracer -> tracer in prior_tracers || tracer in prior_classes, tracers)
    conflicts = unique((class_conflicts..., tracer_conflicts...))
    isempty(conflicts) || throw(
        ArgumentError("component :$name realizes duplicate class/tracer identities $(conflicts)."),
    )
    return nothing
end

"""Realize a named component collection into deterministic class and tracer layouts.

Logical components are visited in `NamedTuple` key order. Population size structure
creates ecological classes once; each population state then realizes one physical
tracer per class. Single-state populations preserve the established class tracer
names (`P_1`, `P_2`, ...). Multi-state populations use class-qualified state names
such as `P_1_carbon`. Pool realization is unchanged.
"""
function realize_components(components::NamedTuple; scalar_type::Type{T}=Float64) where {T<:Real}
    names = keys(components)
    n_components = length(names)

    class_values = Vector{Any}(undef, n_components)
    state_tracer_values = Vector{Any}(undef, n_components)
    tracer_values = Vector{Any}(undef, n_components)
    index_values = Vector{Any}(undef, n_components)
    diameter_values = Vector{Any}(undef, n_components)
    class_order = Symbol[]
    tracer_order = Symbol[]

    for i in 1:n_components
        name = names[i]
        component = getproperty(components, name)
        component isa Union{Population,Pool} || throw(
            ArgumentError(
                "component :$name must be a Population or Pool; got $(typeof(component)).",
            ),
        )

        realized = _realized_component(component, name, T)
        _check_component_identity_conflicts(
            name, realized.classes, realized.tracers, class_order, tracer_order
        )

        append!(class_order, realized.classes)
        first_index = length(tracer_order) + 1
        append!(tracer_order, realized.tracers)
        last_index = length(tracer_order)

        class_values[i] = realized.classes
        state_tracer_values[i] = realized.state_tracers
        tracer_values[i] = realized.tracers
        index_values[i] = Tuple(first_index:last_index)
        diameter_values[i] = realized.diameters
    end

    return ComponentLayout(
        Tuple(tracer_order),
        NamedTuple{names}(Tuple(class_values)),
        NamedTuple{names}(Tuple(state_tracer_values)),
        NamedTuple{names}(Tuple(tracer_values)),
        NamedTuple{names}(Tuple(index_values)),
        NamedTuple{names}(Tuple(diameter_values)),
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

function _population_group_owner(population_groups::NamedTuple)
    owners = Dict{Symbol,Symbol}()
    for population in keys(population_groups), group in getproperty(population_groups, population)
        owners[group] = population
    end
    return owners
end

function _population_classes(
    population::Symbol, population_groups::NamedTuple, context::CommunityContext
)
    classes = ()
    for group in getproperty(population_groups, population)
        classes = (classes..., Tuple(context.class_symbols[context.group_indices[group]])...)
    end
    return classes
end

"""Realize logical components over named population subgroups.

Population ecological classes come from `CommunityContext`; physical state tracers
are then realized from those classes. Pool components are realized directly from
their intrinsic structure. Pool tracers precede population state tracers, while
population class ordering follows the parsed community ordering.
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
    class_conflicts = filter(class -> class in pool_layout.tracer_order, context.class_symbols)
    isempty(class_conflicts) || throw(
        ArgumentError(
            "pool tracer names conflict with realized population classes: $(unique(class_conflicts))"
        ),
    )

    population_names = _population_names(components)
    population_classes = NamedTuple{population_names}(
        Tuple(_population_classes(name, population_groups, context) for name in population_names)
    )
    population_state_tracers = NamedTuple{population_names}(
        Tuple(
            _population_state_tracers(
                getproperty(components, name), getproperty(population_classes, name)
            ) for name in population_names
        )
    )

    group_owner = _population_group_owner(population_groups)
    population_tracer_order = Symbol[]
    for class_index in eachindex(context.class_symbols)
        population = group_owner[context.group_symbols[class_index]]
        component = getproperty(components, population)
        state_names = keys(states(component))
        class = context.class_symbols[class_index]
        nstates = length(state_names)
        for state in state_names
            push!(population_tracer_order, _population_state_tracer(class, state, nstates))
        end
    end
    length(unique(population_tracer_order)) == length(population_tracer_order) || throw(
        ArgumentError("population states realize duplicate tracer identities"),
    )
    physical_conflicts = filter(tracer -> tracer in pool_layout.tracer_order, population_tracer_order)
    isempty(physical_conflicts) || throw(
        ArgumentError(
            "pool tracer names conflict with realized population state tracers: $(unique(physical_conflicts))"
        ),
    )

    tracer_order = Symbol[pool_layout.tracer_order..., population_tracer_order...]

    class_values = ntuple(length(names)) do i
        name = names[i]
        component = getproperty(components, name)
        component isa Pool && return component_classes(pool_layout, name)
        return getproperty(population_classes, name)
    end
    component_classes_values = NamedTuple{names}(class_values)

    state_tracer_values = ntuple(length(names)) do i
        name = names[i]
        component = getproperty(components, name)
        component isa Pool && return nothing
        return getproperty(population_state_tracers, name)
    end
    component_state_tracers_values = NamedTuple{names}(state_tracer_values)

    tracer_values = ntuple(length(names)) do i
        name = names[i]
        component = getproperty(components, name)
        component isa Pool && return component_tracers(pool_layout, name)
        return _flatten_population_tracers(
            component,
            getproperty(population_classes, name),
            getproperty(population_state_tracers, name),
        )
    end
    component_tracers_values = NamedTuple{names}(tracer_values)

    index_values = ntuple(length(names)) do i
        _physical_indices(tracer_values[i], tracer_order, names[i])
    end
    component_indices_values = NamedTuple{names}(index_values)

    diameter_values = ntuple(length(names)) do i
        name = names[i]
        component = getproperty(components, name)
        component isa Pool && return component_diameters(pool_layout, name)

        classes = getproperty(population_classes, name)
        indices = Tuple(findfirst(==(class), context.class_symbols) for class in classes)
        any(isnothing, indices) && throw(
            ArgumentError("component :$name realizes a class absent from the population layout"),
        )
        return Tuple(context.diameters[Int(index)] for index in indices)
    end
    component_diameters_values = NamedTuple{names}(diameter_values)

    return ComponentLayout(
        Tuple(tracer_order),
        component_classes_values,
        component_state_tracers_values,
        component_tracers_values,
        component_indices_values,
        component_diameters_values,
    )
end
