"""One authoritative setup-time realization of components, ecological classes, and inputs.

`ModelLayout` owns physical tracer positions, logical component classes, population-state
tracers, class diameters, interaction axes, and final tracer/auxiliary input positions. It is
constructed once and consumed by parameter planning, compilation, manifests, and host-side
metadata.
"""
struct ModelLayout{T<:Real,TR,II,CC,CST,CT,CD,CS,GI,D,CO,PR}
    scalar_type::Type{T}
    tracer_order::TR
    input_indices::II
    component_classes::CC
    component_state_tracers::CST
    component_tracers::CT
    component_diameters::CD
    class_symbols::CS
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

"""Return the state-to-tracer mapping for a population, or `nothing` for a pool."""
component_state_tracers(layout::ModelLayout, component::Symbol) =
    _component_value(layout, :component_state_tracers, component)

"""Return flattened concrete tracer identities realized by one logical component."""
component_tracers(layout::ModelLayout, component::Symbol) =
    _component_value(layout, :component_tracers, component)

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

"""Return class diameters for one component, or `nothing` when no diameter structure is defined."""
component_diameters(layout::ModelLayout, component::Symbol) =
    _component_value(layout, :component_diameters, component)

function canonicalize_population_realization(
    components::NamedTuple,
    population_groups=nothing,
    group_diameters=nothing,
)
    population_names = _population_names(components)

    if isnothing(population_groups) && isnothing(group_diameters)
        groups = NamedTuple{population_names}(Tuple((name,) for name in population_names))
        diameters = NamedTuple{population_names}(ntuple(length(population_names)) do i
            name = population_names[i]
            structure = size_structure(getproperty(components, name))
            isnothing(structure) && return nothing
            return canonicalize_diameters(
                structure; path="component :$name size_structure"
            ).specification
        end)
        return groups, diameters
    end

    xor(isnothing(population_groups), isnothing(group_diameters)) && throw(
        ArgumentError("population_groups and group_diameters must be supplied together"),
    )

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
        canonicalize_diameters(
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
    canonical = canonicalize_diameters(structure; path="component :$name size_structure")
    classes = ntuple(i -> Symbol(string(name), "_", i), canonical.n)
    diameters = Tuple(realize_diameters(T, canonical.specification))
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

"""Realize canonical population/group inputs into one `ModelLayout`."""
function realize_model_layout(
    components::NamedTuple,
    population_groups::NamedTuple,
    group_diameters::NamedTuple;
    scalar_type::Type{T}=Float64,
    interaction_roles=nothing,
    auxiliary_fields::Tuple=(),
) where {T<:Real}
    all(component -> component isa Union{Population,Pool}, values(components)) || throw(
        ArgumentError("all model components must be Population or Pool values"),
    )

    component_names = keys(components)
    population_names = _population_names(components)
    pool_names = _pool_names(components)

    classes_by_component = Dict{Symbol,Vector{Symbol}}(name => Symbol[] for name in component_names)
    tracers_by_component = Dict{Symbol,Vector{Symbol}}(name => Symbol[] for name in component_names)
    diameters_by_component = Dict{Symbol,Any}(name => nothing for name in component_names)
    state_tracers_by_component = Dict{Symbol,Any}(name => nothing for name in component_names)

    state_tracer_vectors = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()

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
        end
    end

    for population in population_names
        state_names = states(getproperty(components, population))
        state_tracer_vectors[population] = Dict(state => Symbol[] for state in state_names)
        groups = getproperty(population_groups, population)
        has_diameters = any(group -> getproperty(group_diameters, group) !== nothing, groups)
        diameters_by_component[population] = has_diameters ? T[] : nothing
    end

    group_owner = Dict{Symbol,Symbol}()
    for population in population_names, group in getproperty(population_groups, population)
        group_owner[group] = population
    end

    class_symbols = Symbol[]
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
            realize_diameters(T, specification)
        end
        nclasses = length(group_diameters_realized)
        classes = specification === nothing ?
            (group,) : ntuple(i -> Symbol(string(group), "_", i), nclasses)
        state_names = states(component)
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
            push!(diameters, group_diameters_realized[group_local])
            global_class_index = length(class_symbols)
            push!(group_global_indices, global_class_index)
            push!(classes_by_component[population], class)
            population_diameters = diameters_by_component[population]
            population_diameters === nothing ||
                push!(population_diameters, group_diameters_realized[group_local])

            for (state_position, state) in pairs(state_names)
                tracer = physical_tracers[state_position]
                push!(tracer_order, tracer)
                push!(tracers_by_component[population], tracer)
                push!(state_tracer_vectors[population][state], tracer)
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
        state_names = states(component)
        state_tracers_by_component[population] = NamedTuple{state_names}(
            Tuple(Tuple(state_tracer_vectors[population][state]) for state in state_names)
        )
        population_diameters = diameters_by_component[population]
        population_diameters === nothing ||
            (diameters_by_component[population] = Tuple(population_diameters))
    end

    component_classes_values = NamedTuple{component_names}(
        Tuple(Tuple(classes_by_component[name]) for name in component_names)
    )
    component_state_tracers_values = NamedTuple{component_names}(
        Tuple(state_tracers_by_component[name] for name in component_names)
    )
    component_tracers_values = NamedTuple{component_names}(
        Tuple(Tuple(tracers_by_component[name]) for name in component_names)
    )
    component_diameters_values = NamedTuple{component_names}(
        Tuple(diameters_by_component[name] for name in component_names)
    )

    tracer_order_tuple = Tuple(tracer_order)
    _validate_auxiliary_fields(auxiliary_fields, tracer_order_tuple)
    input_names = (tracer_order_tuple..., auxiliary_fields...)
    input_indices = NamedTuple{input_names}(Tuple(eachindex(input_names)))

    return ModelLayout(
        T,
        tracer_order_tuple,
        input_indices,
        component_classes_values,
        component_state_tracers_values,
        component_tracers_values,
        component_diameters_values,
        Tuple(class_symbols),
        group_indices,
        Tuple(diameters),
        consumer_indices,
        prey_indices,
    )
end

"""Canonicalize authored population realization inputs and realize one `ModelLayout`."""
function realize_model_layout(
    components::NamedTuple;
    scalar_type::Type{T}=Float64,
    population_groups=nothing,
    group_diameters=nothing,
    interaction_roles=nothing,
    auxiliary_fields::Tuple=(),
) where {T<:Real}
    population_groups, group_diameters = canonicalize_population_realization(
        components, population_groups, group_diameters
    )
    return realize_model_layout(
        components,
        population_groups,
        group_diameters;
        scalar_type=T,
        interaction_roles,
        auxiliary_fields,
    )
end

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
