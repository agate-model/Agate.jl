"""Population state described by intrinsic component properties.

Ecological function is supplied by process participation. `states` maps each
prognostic state identity to the conserved-material currency it represents,
while `size_structure` describes the shared ecological class realization.
A one-state population is authored with its currency as the positional argument.
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

function Population(currency; size_structure=nothing)
    currency isa Symbol || throw(ArgumentError("Population currency must be a Symbol."))
    return Population(NamedTuple{(currency,)}((currency,)), size_structure)
end

function Population(; states, size_structure=nothing)
    states isa NamedTuple || throw(
        ArgumentError("Population states must be a NamedTuple mapping state names to currencies."),
    )
    return Population(states, size_structure)
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

"""Return the conserved-material currency represented by a single-state component."""
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

"""One authoritative setup-time realization of components, ecological classes, and inputs.

`ModelLayout` owns physical tracer positions, logical component classes, population-state
tracers and indices, population subgroup membership, class diameters, interaction axes,
and auxiliary input positions. It is constructed once and consumed by parameter planning,
compilation, manifests, and host-side metadata.
"""
struct ModelLayout{T<:Real,TR,TI,AF,AI,CC,CCI,CST,CSI,CT,CI,CD,PG,CS,CP,GS,GL,CL,GI,D,CO,PR}
    scalar_type::Type{T}
    tracer_order::TR
    tracer_indices::TI
    auxiliary_fields::AF
    auxiliary_indices::AI
    component_classes::CC
    component_class_indices::CCI
    component_state_tracers::CST
    component_state_indices::CSI
    component_tracers::CT
    component_indices::CI
    component_diameters::CD
    population_groups::PG
    class_symbols::CS
    class_populations::CP
    group_symbols::GS
    group_local_indices::GL
    component_local_indices::CL
    group_indices::GI
    diameters::D
    consumer_indices::CO
    prey_indices::PR
end

_population_names(components::NamedTuple) = Tuple(
    name for name in keys(components) if getproperty(components, name) isa Population
)
_pool_names(components::NamedTuple) = Tuple(
    name for name in keys(components) if getproperty(components, name) isa Pool
)

@inline _unspecified_diameter(::Type{T}) where {T<:AbstractFloat} = T(NaN)
@inline _unspecified_diameter(::Type{T}) where {T<:Real} = zero(T)

@inline function _population_state_tracer(class::Symbol, state::Symbol, nstates::Int)
    return nstates == 1 ? class : Symbol(string(class), "_", state)
end

function _component_value(layout::ModelLayout, field::Symbol, component::Symbol)
    values = getfield(layout, field)
    hasproperty(values, component) || throw(ArgumentError("unknown component :$component"))
    return getproperty(values, component)
end

"""Return ecological class identities realized by one logical component."""
component_classes(layout::ModelLayout, component::Symbol) =
    _component_value(layout, :component_classes, component)

"""Return global ecological-class positions for a population, or `nothing` for a pool."""
component_class_indices(layout::ModelLayout, component::Symbol) =
    _component_value(layout, :component_class_indices, component)

"""Return the state-to-tracer mapping for a population, or `nothing` for a pool."""
component_state_tracers(layout::ModelLayout, component::Symbol) =
    _component_value(layout, :component_state_tracers, component)

"""Return the state-to-physical-index mapping for a population, or `nothing` for a pool."""
component_state_indices(layout::ModelLayout, component::Symbol) =
    _component_value(layout, :component_state_indices, component)

"""Return flattened concrete tracer identities realized by one logical component."""
component_tracers(layout::ModelLayout, component::Symbol) =
    _component_value(layout, :component_tracers, component)

"""Return physical tracer positions realized by one logical component."""
component_indices(layout::ModelLayout, component::Symbol) =
    _component_value(layout, :component_indices, component)

"""Return concrete tracers for one population state."""
function state_tracers(layout::ModelLayout, component::Symbol, state::Symbol)
    mapping = component_state_tracers(layout, component)
    mapping isa NamedTuple || throw(
        ArgumentError("component :$component does not expose population states"),
    )
    hasproperty(mapping, state) || throw(
        ArgumentError("component :$component has no realized state :$state"),
    )
    return getproperty(mapping, state)
end

state_tracers(layout::ModelLayout, reference::PopulationStateRef) =
    state_tracers(layout, reference.population, reference.state)

"""Return one concrete tracer for a population state and local ecological class ordinal."""
function state_tracer(
    layout::ModelLayout, reference::PopulationStateRef, class_ordinal::Integer
)
    tracers = state_tracers(layout, reference)
    1 <= class_ordinal <= length(tracers) || throw(
        ArgumentError(
            "class ordinal $class_ordinal is out of bounds for population :$(reference.population) state :$(reference.state)",
        ),
    )
    return tracers[Int(class_ordinal)]
end

"""Return precomputed physical tracer positions for one population state."""
function state_indices(layout::ModelLayout, component::Symbol, state::Symbol)
    mapping = component_state_indices(layout, component)
    mapping isa NamedTuple || throw(
        ArgumentError("component :$component does not expose population states"),
    )
    hasproperty(mapping, state) || throw(
        ArgumentError("component :$component has no realized state :$state"),
    )
    return getproperty(mapping, state)
end

state_indices(layout::ModelLayout, reference::PopulationStateRef) =
    state_indices(layout, reference.population, reference.state)

"""Return class diameters for one component, or `nothing` for a scalar pool."""
component_diameters(layout::ModelLayout, component::Symbol) =
    _component_value(layout, :component_diameters, component)

"""Return the number of ecological classes realized by one logical component."""
component_class_count(layout::ModelLayout, component::Symbol) =
    length(component_classes(layout, component))

function _canonical_population_realization(
    components::NamedTuple,
    population_groups,
    group_diameters;
    intrinsic::Bool,
)
    population_names = _population_names(components)

    if intrinsic
        groups = NamedTuple{population_names}(Tuple((name,) for name in population_names))
        diameters = NamedTuple{population_names}(ntuple(length(population_names)) do i
            name = population_names[i]
            structure = size_structure(getproperty(components, name))
            isnothing(structure) && return nothing
            return normalize_diameters(
                structure; path="component :$name size_structure"
            ).specification
        end)
        return groups, diameters
    end

    population_groups isa NamedTuple || throw(
        ArgumentError("population_groups must be a NamedTuple"),
    )
    group_diameters isa NamedTuple || throw(
        ArgumentError("group_diameters must be a NamedTuple"),
    )
    Set(keys(population_groups)) == Set(population_names) || throw(
        ArgumentError(
            "population_groups must map exactly the logical population components $(population_names)",
        ),
    )

    assigned = Symbol[]
    group_values = ntuple(length(population_names)) do i
        population = population_names[i]
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
            group in assigned && throw(
                ArgumentError("population subgroup :$group is assigned more than once"),
            )
            push!(assigned, group)
        end
        groups
    end
    Set(assigned) == Set(keys(group_diameters)) || throw(
        ArgumentError("population_groups and group_diameters must contain the same subgroups"),
    )

    canonical_groups = NamedTuple{population_names}(group_values)
    diameter_names = keys(group_diameters)
    canonical_diameters = NamedTuple{diameter_names}(ntuple(length(diameter_names)) do i
        group = diameter_names[i]
        normalize_diameters(
            getproperty(group_diameters, group); path="population group :$group diameters"
        ).specification
    end)
    return canonical_groups, canonical_diameters
end

function _role_indices(role, role_name::Symbol, group_indices::NamedTuple, nclasses::Int)
    role === nothing && return Tuple(1:nclasses)
    (role isa Tuple || role isa AbstractVector{Symbol}) || throw(
        ArgumentError("$role_name roles must be a collection of population subgroup Symbols"),
    )
    requested = Set{Symbol}(role)
    for group in requested
        hasproperty(group_indices, group) || throw(
            ArgumentError("unknown population subgroup :$group in $role_name roles"),
        )
    end
    indices = Int[]
    for group in keys(group_indices)
        group in requested || continue
        append!(indices, getproperty(group_indices, group))
    end
    return Tuple(indices)
end

function _validate_auxiliary_fields(auxiliary_fields::Tuple, tracer_order::Tuple)
    all(field -> field isa Symbol, auxiliary_fields) || throw(
        ArgumentError("auxiliary_fields entries must be Symbols"),
    )
    length(unique(auxiliary_fields)) == length(auxiliary_fields) || throw(
        ArgumentError("auxiliary_fields contains duplicate entries"),
    )
    conflicts = Tuple(field for field in auxiliary_fields if field in tracer_order)
    isempty(conflicts) || throw(
        ArgumentError("auxiliary fields conflict with tracer names: $(collect(conflicts))"),
    )
    return nothing
end

function _realize_pool(component::Pool, name::Symbol, ::Type{T}) where {T<:Real}
    structure = size_structure(component)
    isnothing(structure) && return ((name,), nothing)
    normalized = normalize_diameters(structure; path="component :$name size_structure")
    classes = ntuple(i -> Symbol(string(name), "_", i), normalized.n)
    diameters = Tuple(param_compute_diameters(T, normalized.specification))
    return classes, diameters
end

function _check_new_identities!(
    name::Symbol,
    classes,
    tracers,
    prior_classes::Set{Symbol},
    prior_tracers::Set{Symbol},
)
    conflicts = Symbol[]
    for class in classes
        (class in prior_classes || class in prior_tracers) && push!(conflicts, class)
    end
    for tracer in tracers
        (tracer in prior_classes || tracer in prior_tracers) && push!(conflicts, tracer)
    end
    isempty(conflicts) || throw(
        ArgumentError("component/group :$name realizes duplicate class/tracer identities $(unique(conflicts))"),
    )
    union!(prior_classes, classes)
    union!(prior_tracers, tracers)
    return nothing
end

"""Realize authored components and optional registered-family subgroup inputs into one `ModelLayout`.

When subgroup inputs are omitted, each population's intrinsic `size_structure` is used.
Registered families supply `population_groups` and `group_diameters`; both paths then share
the same append/indexing pass and downstream representation.
"""
function realize_model_layout(
    components::NamedTuple;
    scalar_type::Type{T}=Float64,
    population_groups=nothing,
    group_diameters=nothing,
    interaction_roles=nothing,
    auxiliary_fields::Tuple=(),
) where {T<:Real}
    all(component -> component isa Union{Population,Pool}, values(components)) || throw(
        ArgumentError("all model components must be Population or Pool values"),
    )
    intrinsic = isnothing(population_groups) && isnothing(group_diameters)
    xor(isnothing(population_groups), isnothing(group_diameters)) && throw(
        ArgumentError("population_groups and group_diameters must be supplied together"),
    )
    population_groups, group_diameters = _canonical_population_realization(
        components, population_groups, group_diameters; intrinsic
    )

    component_names = keys(components)
    population_names = _population_names(components)
    pool_names = _pool_names(components)

    classes_by_component = Dict{Symbol,Vector{Symbol}}(name => Symbol[] for name in component_names)
    class_indices_by_component = Dict{Symbol,Any}(name => nothing for name in component_names)
    tracers_by_component = Dict{Symbol,Vector{Symbol}}(name => Symbol[] for name in component_names)
    indices_by_component = Dict{Symbol,Vector{Int}}(name => Int[] for name in component_names)
    diameters_by_component = Dict{Symbol,Any}(name => nothing for name in component_names)
    state_tracers_by_component = Dict{Symbol,Any}(name => nothing for name in component_names)
    state_indices_by_component = Dict{Symbol,Any}(name => nothing for name in component_names)

    state_tracer_vectors = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()
    state_index_vectors = Dict{Symbol,Dict{Symbol,Vector{Int}}}()
    ecological_indices_by_population = Dict{Symbol,Vector{Int}}(
        name => Int[] for name in population_names
    )
    component_local_counts = Dict{Symbol,Int}(name => 0 for name in population_names)

    tracer_order = Symbol[]
    seen_classes = Set{Symbol}()
    seen_tracers = Set{Symbol}()

    # Physical pools are always appended before population states.
    for name in pool_names
        component = getproperty(components, name)
        classes, diameters = _realize_pool(component, name, T)
        _check_new_identities!(name, classes, classes, seen_classes, seen_tracers)
        append!(classes_by_component[name], classes)
        diameters_by_component[name] = diameters
        for tracer in classes
            push!(tracer_order, tracer)
            push!(tracers_by_component[name], tracer)
            push!(indices_by_component[name], length(tracer_order))
        end
    end

    for population in population_names
        state_names = keys(states(getproperty(components, population)))
        state_tracer_vectors[population] = Dict(state => Symbol[] for state in state_names)
        state_index_vectors[population] = Dict(state => Int[] for state in state_names)
        groups = getproperty(population_groups, population)
        has_diameters = any(group -> getproperty(group_diameters, group) !== nothing, groups)
        diameters_by_component[population] = has_diameters ? T[] : nothing
    end

    group_owner = Dict{Symbol,Symbol}()
    for population in population_names, group in getproperty(population_groups, population)
        group_owner[group] = population
    end

    class_symbols = Symbol[]
    class_populations = Symbol[]
    group_symbols = Symbol[]
    group_local_indices = Int[]
    component_local_indices = Int[]
    diameters = T[]
    group_index_values = Vector{Any}(undef, length(keys(group_diameters)))
    group_names = keys(group_diameters)

    for (group_position, group) in pairs(group_names)
        population = group_owner[group]
        component = getproperty(components, population)
        specification = getproperty(group_diameters, group)
        group_diameters_realized = if specification === nothing
            T[_unspecified_diameter(T)]
        else
            param_compute_diameters(T, specification)
        end
        nclasses = length(group_diameters_realized)
        classes = if intrinsic && specification === nothing && group === population
            (population,)
        else
            ntuple(i -> Symbol(string(group), "_", i), nclasses)
        end
        state_names = keys(states(component))
        nstates = length(state_names)
        group_global_indices = Int[]

        for group_local in eachindex(classes)
            class = classes[group_local]
            physical_tracers = Tuple(
                _population_state_tracer(class, state, nstates) for state in state_names
            )
            _check_new_identities!(
                group, (class,), physical_tracers, seen_classes, seen_tracers
            )

            push!(class_symbols, class)
            push!(class_populations, population)
            push!(group_symbols, group)
            push!(group_local_indices, group_local)
            component_local_counts[population] += 1
            push!(component_local_indices, component_local_counts[population])
            push!(diameters, group_diameters_realized[group_local])
            global_class_index = length(class_symbols)
            push!(group_global_indices, global_class_index)
            push!(ecological_indices_by_population[population], global_class_index)
            push!(classes_by_component[population], class)
            population_diameters = diameters_by_component[population]
            population_diameters === nothing ||
                push!(population_diameters, group_diameters_realized[group_local])

            for (state_position, state) in pairs(state_names)
                tracer = physical_tracers[state_position]
                push!(tracer_order, tracer)
                physical_index = length(tracer_order)
                push!(tracers_by_component[population], tracer)
                push!(indices_by_component[population], physical_index)
                push!(state_tracer_vectors[population][state], tracer)
                push!(state_index_vectors[population][state], physical_index)
            end
        end
        group_index_values[group_position] = Tuple(group_global_indices)
    end

    group_indices = NamedTuple{group_names}(Tuple(group_index_values))
    roles = isnothing(interaction_roles) ? (consumers=nothing, prey=nothing) : interaction_roles
    hasproperty(roles, :consumers) && hasproperty(roles, :prey) || throw(
        ArgumentError("interaction_roles must define :consumers and :prey"),
    )
    consumer_indices = _role_indices(
        roles.consumers, :consumers, group_indices, length(class_symbols)
    )
    prey_indices = _role_indices(roles.prey, :prey, group_indices, length(class_symbols))

    for population in population_names
        component = getproperty(components, population)
        state_names = keys(states(component))
        state_tracers_by_component[population] = NamedTuple{state_names}(
            Tuple(Tuple(state_tracer_vectors[population][state]) for state in state_names)
        )
        state_indices_by_component[population] = NamedTuple{state_names}(
            Tuple(Tuple(state_index_vectors[population][state]) for state in state_names)
        )
        class_indices_by_component[population] = Tuple(
            ecological_indices_by_population[population]
        )
        population_diameters = diameters_by_component[population]
        population_diameters === nothing ||
            (diameters_by_component[population] = Tuple(population_diameters))
    end

    component_classes_values = NamedTuple{component_names}(
        Tuple(Tuple(classes_by_component[name]) for name in component_names)
    )
    component_class_indices_values = NamedTuple{component_names}(
        Tuple(class_indices_by_component[name] for name in component_names)
    )
    component_state_tracers_values = NamedTuple{component_names}(
        Tuple(state_tracers_by_component[name] for name in component_names)
    )
    component_state_indices_values = NamedTuple{component_names}(
        Tuple(state_indices_by_component[name] for name in component_names)
    )
    component_tracers_values = NamedTuple{component_names}(
        Tuple(Tuple(tracers_by_component[name]) for name in component_names)
    )
    component_indices_values = NamedTuple{component_names}(
        Tuple(Tuple(indices_by_component[name]) for name in component_names)
    )
    component_diameters_values = NamedTuple{component_names}(
        Tuple(diameters_by_component[name] for name in component_names)
    )

    tracer_order_tuple = Tuple(tracer_order)
    _validate_auxiliary_fields(auxiliary_fields, tracer_order_tuple)
    tracer_indices = NamedTuple{tracer_order_tuple}(Tuple(eachindex(tracer_order_tuple)))
    auxiliary_indices = NamedTuple{auxiliary_fields}(
        ntuple(i -> length(tracer_order_tuple) + i, length(auxiliary_fields))
    )

    return ModelLayout(
        T,
        tracer_order_tuple,
        tracer_indices,
        auxiliary_fields,
        auxiliary_indices,
        component_classes_values,
        component_class_indices_values,
        component_state_tracers_values,
        component_state_indices_values,
        component_tracers_values,
        component_indices_values,
        component_diameters_values,
        population_groups,
        Tuple(class_symbols),
        Tuple(class_populations),
        Tuple(group_symbols),
        Tuple(group_local_indices),
        Tuple(component_local_indices),
        group_indices,
        Tuple(diameters),
        consumer_indices,
        prey_indices,
    )
end

"""Convenience realization for direct component definitions."""
realize_components(components::NamedTuple; kwargs...) = realize_model_layout(components; kwargs...)

"""Return compact host-side metadata derived from one authoritative `ModelLayout`."""
function model_metadata(layout::ModelLayout; interaction_axes=nothing, parameter_axes=(;))
    group_names = keys(layout.group_indices)
    group_classes = NamedTuple{group_names}(ntuple(length(group_names)) do i
        indices = getproperty(layout.group_indices, group_names[i])
        Tuple(layout.class_symbols[index] for index in indices)
    end)
    population_tracer_set = Set{Symbol}()
    for name in keys(layout.component_state_tracers)
        getproperty(layout.component_state_tracers, name) isa NamedTuple || continue
        union!(population_tracer_set, getproperty(layout.component_tracers, name))
    end
    population_tracers = Tuple(
        tracer for tracer in layout.tracer_order if tracer in population_tracer_set
    )
    return (;
        group_classes,
        class_symbols=layout.class_symbols,
        population_tracers,
        plankton_diameters=layout.diameters,
        interaction_axes,
        parameter_axes,
    )
end
