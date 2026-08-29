"""One authoritative setup-time realization of components, realized entities, and inputs.

`ModelLayout` owns physical tracer positions, logical component entities, plankton-state
tracers, entity diameters, interaction axes, and final tracer/auxiliary input positions. It is
constructed once and consumed by parameter planning, compilation, manifests, and host-side
metadata.
"""
struct ModelLayout{T<:Real,TR,II,CC,CST,CT,CD,CS,GI,D,CO,PR}
    scalar_type::Type{T}
    tracer_order::TR
    input_indices::II
    component_entities::CC
    component_state_tracers::CST
    component_tracers::CT
    component_diameters::CD
    size_classes::CS
    pft_indices::GI
    size_class_diameters::D
    consumer_indices::CO
    prey_indices::PR
end

_plankton_names(components::NamedTuple) = Tuple(
    name for name in keys(components) if getproperty(components, name) isa Plankton
)
_pool_names(components::NamedTuple) = Tuple(
    name for name in keys(components) if getproperty(components, name) isa Pool
)

@inline _unspecified_diameter(::Type{T}) where {T<:AbstractFloat} = T(NaN)
@inline _unspecified_diameter(::Type{T}) where {T<:Real} = zero(T)

@inline function _plankton_state_tracer(size_class::Symbol, state::Symbol, nstates::Int)
    return nstates == 1 ? size_class : Symbol(string(size_class), "_", state)
end

function _component_value(layout::ModelLayout, field::Symbol, component::Symbol)
    values = getfield(layout, field)
    hasproperty(values, component) || throw(ArgumentError("unknown component :$component"))
    return getproperty(values, component)
end

"""Return realized entity identities for one logical component."""
component_entities(layout::ModelLayout, component::Symbol) =
    _component_value(layout, :component_entities, component)

"""Return the state-to-tracer mapping for a plankton, or `nothing` for a pool."""
component_state_tracers(layout::ModelLayout, component::Symbol) =
    _component_value(layout, :component_state_tracers, component)

"""Return flattened concrete tracer identities realized by one logical component."""
component_tracers(layout::ModelLayout, component::Symbol) =
    _component_value(layout, :component_tracers, component)

"""Return concrete tracers for one plankton state."""
function state_tracers(layout::ModelLayout, component::Symbol, state::Symbol)
    mapping = component_state_tracers(layout, component)
    mapping isa NamedTuple || throw(
        ArgumentError("component :$component does not expose plankton states"),
    )
    hasproperty(mapping, state) || throw(
        ArgumentError("component :$component has no realized state :$state"),
    )
    return getproperty(mapping, state)
end

state_tracers(layout::ModelLayout, reference::PlanktonStateRef) =
    state_tracers(layout, reference.plankton, reference.state)

"""Return one concrete tracer for a plankton state and local realized-entity index."""
function state_tracer(
    layout::ModelLayout, reference::PlanktonStateRef, entity_index::Integer
)
    tracers = state_tracers(layout, reference)
    1 <= entity_index <= length(tracers) || throw(
        ArgumentError(
            "entity index $entity_index is out of bounds for plankton :$(reference.plankton) state :$(reference.state)",
        ),
    )
    return tracers[Int(entity_index)]
end

"""Return entity diameters for one component, or `nothing` when no diameter structure is defined."""
component_diameters(layout::ModelLayout, component::Symbol) =
    _component_value(layout, :component_diameters, component)

function canonicalize_plankton_realization(
    components::NamedTuple,
    plankton_pfts=nothing,
    pft_size_structures=nothing,
)
    plankton_names = _plankton_names(components)

    if isnothing(plankton_pfts) && isnothing(pft_size_structures)
        pfts = NamedTuple{plankton_names}(Tuple((name,) for name in plankton_names))
        size_structures = NamedTuple{plankton_names}(ntuple(length(plankton_names)) do i
            name = plankton_names[i]
            structure = size_structure(getproperty(components, name))
            isnothing(structure) && return nothing
            return canonicalize_diameters(
                structure; path="component :$name size_structure"
            ).specification
        end)
        return pfts, size_structures
    end

    xor(isnothing(plankton_pfts), isnothing(pft_size_structures)) && throw(
        ArgumentError("plankton_pfts and pft_size_structures must be supplied together"),
    )

    plankton_pfts isa NamedTuple || throw(
        ArgumentError("plankton_pfts must be a NamedTuple"),
    )
    pft_size_structures isa NamedTuple || throw(
        ArgumentError("pft_size_structures must be a NamedTuple"),
    )
    Set(keys(plankton_pfts)) == Set(plankton_names) || throw(
        ArgumentError(
            "plankton_pfts must map exactly the logical plankton components $(plankton_names)",
        ),
    )

    assigned = Symbol[]
    pft_values = ntuple(length(plankton_names)) do i
        plankton = plankton_names[i]
        pfts = getproperty(plankton_pfts, plankton)
        pfts isa Tuple || throw(
            ArgumentError("plankton component :$plankton PFT realization must be a tuple"),
        )
        isempty(pfts) && throw(
            ArgumentError("plankton component :$plankton must realize at least one PFT"),
        )
        for pft in pfts
            pft isa Symbol || throw(
                ArgumentError("plankton component :$plankton PFT identities must be Symbols"),
            )
            pft in assigned && throw(
                ArgumentError("plankton PFT :$pft is assigned more than once"),
            )
            push!(assigned, pft)
        end
        pfts
    end
    Set(assigned) == Set(keys(pft_size_structures)) || throw(
        ArgumentError("plankton_pfts and pft_size_structures must contain the same PFTs"),
    )

    canonical_pfts = NamedTuple{plankton_names}(pft_values)
    pft_names = Tuple(assigned)
    canonical_size_structures = NamedTuple{pft_names}(ntuple(length(pft_names)) do i
        pft = pft_names[i]
        canonicalize_diameters(
            getproperty(pft_size_structures, pft); path="plankton PFT :$pft size_structure"
        ).specification
    end)
    return canonical_pfts, canonical_size_structures
end

function _role_indices(role, role_name::Symbol, pft_indices::NamedTuple, nsizeclasses::Int)
    role === nothing && return Tuple(1:nsizeclasses)
    (role isa Tuple || role isa AbstractVector{Symbol}) || throw(
        ArgumentError("$role_name roles must be a collection of plankton PFT Symbols"),
    )
    requested = Set{Symbol}(role)
    for pft in requested
        hasproperty(pft_indices, pft) || throw(
            ArgumentError("unknown plankton PFT :$pft in $role_name roles"),
        )
    end
    indices = Int[]
    for pft in keys(pft_indices)
        pft in requested || continue
        append!(indices, getproperty(pft_indices, pft))
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

function _check_new_identities!(
    name::Symbol,
    entities,
    tracers,
    prior_entities::Set{Symbol},
    prior_tracers::Set{Symbol},
)
    conflicts = Symbol[]
    for entity in entities
        (entity in prior_entities || entity in prior_tracers) && push!(conflicts, entity)
    end
    for tracer in tracers
        (tracer in prior_entities || tracer in prior_tracers) && push!(conflicts, tracer)
    end
    isempty(conflicts) || throw(
        ArgumentError("component/PFT :$name realizes duplicate entity/tracer identities $(unique(conflicts))"),
    )
    union!(prior_entities, entities)
    union!(prior_tracers, tracers)
    return nothing
end

"""Realize canonical plankton/PFT inputs into one `ModelLayout`."""
function realize_model_layout(
    components::NamedTuple,
    plankton_pfts::NamedTuple,
    pft_size_structures::NamedTuple;
    scalar_type::Type{T}=Float64,
    interaction_roles=nothing,
    auxiliary_fields::Tuple=(),
) where {T<:Real}
    all(component -> component isa Union{Plankton,Pool}, values(components)) || throw(
        ArgumentError("all model components must be Plankton or Pool values"),
    )

    component_names = keys(components)
    plankton_names = _plankton_names(components)
    pool_names = _pool_names(components)

    entities_by_component = Dict{Symbol,Vector{Symbol}}(name => Symbol[] for name in component_names)
    tracers_by_component = Dict{Symbol,Vector{Symbol}}(name => Symbol[] for name in component_names)
    diameters_by_component = Dict{Symbol,Any}(name => nothing for name in component_names)
    state_tracers_by_component = Dict{Symbol,Any}(name => nothing for name in component_names)

    state_tracer_vectors = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()

    tracer_order = Symbol[]
    seen_entities = Set{Symbol}()
    seen_tracers = Set{Symbol}()

    # Physical pools are always appended before plankton states.
    for name in pool_names
        entities = (name,)
        _check_new_identities!(name, entities, entities, seen_entities, seen_tracers)
        push!(entities_by_component[name], name)
        push!(tracer_order, name)
        push!(tracers_by_component[name], name)
    end

    for plankton in plankton_names
        state_names = states(getproperty(components, plankton))
        state_tracer_vectors[plankton] = Dict(state => Symbol[] for state in state_names)
        pfts = getproperty(plankton_pfts, plankton)
        has_diameters = any(pft -> getproperty(pft_size_structures, pft) !== nothing, pfts)
        diameters_by_component[plankton] = has_diameters ? T[] : nothing
    end

    pft_owner = Dict{Symbol,Symbol}()
    for plankton in plankton_names, pft in getproperty(plankton_pfts, plankton)
        pft_owner[pft] = plankton
    end

    size_classes = Symbol[]
    size_class_diameters = T[]
    pft_names = Tuple(
        pft for plankton in plankton_names for pft in getproperty(plankton_pfts, plankton)
    )
    pft_index_values = Vector{Any}(undef, length(pft_names))

    for (pft_position, pft) in pairs(pft_names)
        plankton = pft_owner[pft]
        component = getproperty(components, plankton)
        specification = getproperty(pft_size_structures, pft)
        realized_diameters = if specification === nothing
            T[_unspecified_diameter(T)]
        else
            realize_diameters(T, specification)
        end
        nsizeclasses = length(realized_diameters)
        realized_size_classes = specification === nothing ?
            (pft,) : ntuple(i -> Symbol(string(pft), "_", i), nsizeclasses)
        state_names = states(component)
        nstates = length(state_names)
        pft_global_indices = Int[]

        for pft_local in eachindex(realized_size_classes)
            size_class = realized_size_classes[pft_local]
            physical_tracers = Tuple(
                _plankton_state_tracer(size_class, state, nstates) for state in state_names
            )
            _check_new_identities!(
                pft, (size_class,), physical_tracers, seen_entities, seen_tracers
            )

            push!(size_classes, size_class)
            push!(size_class_diameters, realized_diameters[pft_local])
            global_size_class_index = length(size_classes)
            push!(pft_global_indices, global_size_class_index)
            push!(entities_by_component[plankton], size_class)
            plankton_diameters = diameters_by_component[plankton]
            plankton_diameters === nothing ||
                push!(plankton_diameters, realized_diameters[pft_local])

            for (state_position, state) in pairs(state_names)
                tracer = physical_tracers[state_position]
                push!(tracer_order, tracer)
                push!(tracers_by_component[plankton], tracer)
                push!(state_tracer_vectors[plankton][state], tracer)
            end
        end
        pft_index_values[pft_position] = Tuple(pft_global_indices)
    end

    pft_indices = NamedTuple{pft_names}(Tuple(pft_index_values))
    roles = isnothing(interaction_roles) ? (consumers=nothing, prey=nothing) : interaction_roles
    hasproperty(roles, :consumers) && hasproperty(roles, :prey) || throw(
        ArgumentError("interaction_roles must define :consumers and :prey"),
    )
    consumer_indices = _role_indices(
        roles.consumers, :consumers, pft_indices, length(size_classes)
    )
    prey_indices = _role_indices(roles.prey, :prey, pft_indices, length(size_classes))

    for plankton in plankton_names
        component = getproperty(components, plankton)
        state_names = states(component)
        state_tracers_by_component[plankton] = NamedTuple{state_names}(
            Tuple(Tuple(state_tracer_vectors[plankton][state]) for state in state_names)
        )
        plankton_diameters = diameters_by_component[plankton]
        plankton_diameters === nothing ||
            (diameters_by_component[plankton] = Tuple(plankton_diameters))
    end

    component_entities_values = NamedTuple{component_names}(
        Tuple(Tuple(entities_by_component[name]) for name in component_names)
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
        component_entities_values,
        component_state_tracers_values,
        component_tracers_values,
        component_diameters_values,
        Tuple(size_classes),
        pft_indices,
        Tuple(size_class_diameters),
        consumer_indices,
        prey_indices,
    )
end

"""Canonicalize authored plankton realization inputs and realize one `ModelLayout`."""
function realize_model_layout(
    components::NamedTuple;
    scalar_type::Type{T}=Float64,
    plankton_pfts=nothing,
    pft_size_structures=nothing,
    interaction_roles=nothing,
    auxiliary_fields::Tuple=(),
) where {T<:Real}
    plankton_pfts, pft_size_structures = canonicalize_plankton_realization(
        components, plankton_pfts, pft_size_structures
    )
    return realize_model_layout(
        components,
        plankton_pfts,
        pft_size_structures;
        scalar_type=T,
        interaction_roles,
        auxiliary_fields,
    )
end

"""Return compact host-side metadata derived from one authoritative `ModelLayout`."""
function model_metadata(layout::ModelLayout; interaction_axes=nothing, parameter_axes=(;))
    pft_names = keys(layout.pft_indices)
    pft_size_classes = NamedTuple{pft_names}(ntuple(length(pft_names)) do i
        indices = getproperty(layout.pft_indices, pft_names[i])
        Tuple(layout.size_classes[index] for index in indices)
    end)
    plankton_tracer_set = Set{Symbol}()
    for name in keys(layout.component_state_tracers)
        getproperty(layout.component_state_tracers, name) isa NamedTuple || continue
        union!(plankton_tracer_set, getproperty(layout.component_tracers, name))
    end
    plankton_tracers = Tuple(
        tracer for tracer in layout.tracer_order if tracer in plankton_tracer_set
    )
    return (;
        pft_size_classes,
        size_classes=layout.size_classes,
        plankton_tracers,
        plankton_diameters=layout.size_class_diameters,
        interaction_axes,
        parameter_axes,
    )
end
