using ..Parameters: ConstructionParameter
using ..Runtime: active_parameters

# The numerical map retains only the fixed numeric baseline, derivations and
# index maps. Structural validation and compilation happen on the host.
struct TraitParameterMap{T,Numerics,Slots,Gather}
    numerics::Numerics
    slots::Slots
    gather::Gather
end

@inline _trait_update(tau, baseline, index::Int) = tau[index]
@inline _trait_update(tau, baseline, indices::Tuple) = [
    index == 0 ? baseline[i] : tau[index] for (i, index) in enumerate(indices)
]

@inline function _trait_updates(tau, baseline, slots::NamedTuple{Names}) where {Names}
    # Select the fixed fields before mapping over values. A runtime Symbol
    # lookup on a mixed scalar/array NamedTuple introduces boxed unions in AD.
    selected = NamedTuple{Names}(baseline)
    return NamedTuple{Names}(map(
        (base, indices) -> _trait_update(tau, base, indices), values(selected), values(slots)
    ))
end

@inline function _realized_trait_values(f::TraitParameterMap, tau)
    updates = _trait_updates(tau, f.numerics.baseline, f.slots)
    return _rematerialize_traits(f.numerics, updates)
end

@inline function (f::TraitParameterMap{T})(tau::AbstractVector{T}) where {T}
    return _gather_runtime_vector(T, _realized_trait_values(f, tau), f.gather)
end

"""Selected continuous construction traits and their dependent runtime parameters.

Returned by [`active_traits`](@ref). `labels` and `values` describe the flat trait
vector; `active_parameters` describes the corresponding runtime parameter vector.
"""
struct ActiveTraitSet{BGC,Active,Values,Labels,Selection,Overrides,Map,Realization}
    bgc::BGC
    active_parameters::Active
    values::Values
    labels::Labels
    trait_selection::Selection
    explicit_override_keys::Overrides
    parameter_map::Map
    realization::Realization
end

function _trait_slots(prepared, selection)
    T = prepared.layout.scalar_type
    T <: AbstractFloat || throw(ArgumentError("trait realization requires floating-point parameters"))
    isempty(selection) && throw(ArgumentError("select at least one construction trait"))
    tau = T[]
    labels = String[]
    slots = map(keys(selection), values(selection)) do name, selector
        parameter = planned_parameter(prepared.plan, name)
        parameter.definition isa ConstructionParameter || throw(ArgumentError(
            "parameter :$name is not a ConstructionParameter; use Runtime.active_parameters for runtime parameters"
        ))
        parameter.definition.default isa DerivedDefault && throw(ArgumentError(
            "construction trait :$name is itself derived; select its independent inputs"
        ))
        baseline = getproperty(prepared.values, name)
        _validate_trait_input(baseline, name)
        if parameter.rank == 0
            selector === true || throw(ArgumentError("select scalar trait :$name with true"))
            push!(tau, baseline)
            push!(labels, string(name))
            return length(tau)
        end
        selector isa Tuple && !isempty(selector) && all(x -> x isa Symbol, selector) &&
            allunique(selector) || throw(ArgumentError(
                "select vector trait :$name with a nonempty tuple of distinct SizeClass symbols"
            ))
        axis = only(parameter.storage_labels)
        indices = zeros(Int, length(axis))
        for label in selector
            i = findfirst(==(label), axis)
            isnothing(i) && throw(ArgumentError(
                "unknown SizeClass :$label for trait :$name; expected one of $axis"
            ))
            push!(tau, baseline[i])
            push!(labels, "$name.$label")
            indices[i] = length(tau)
        end
        return Tuple(indices)
    end
    return NamedTuple{keys(selection)}(slots), tau, Tuple(labels)
end

function _trait_runtime_selection(prepared, traits)
    # A declared dependency names an entire parameter, not individual entries.
    # Preserve plan order and select all applicable entries of each live output.
    targets = Tuple(name for name in keys(prepared.derivations) if begin
        parameter = planned_parameter(prepared.plan, name)
        parameter.runtime_bound && any(dep -> dep in traits, parameter.definition.default.deps)
    end)
    for trait in traits
        any(name -> trait in planned_parameter(prepared.plan, name).definition.default.deps, targets) ||
            throw(ArgumentError("construction trait :$trait has no live derived runtime output (check explicit output overrides)"))
    end
    return NamedTuple{targets}(map(targets) do name
        parameter = planned_parameter(prepared.plan, name)
        if parameter.rank == 0
            typeof(getproperty(prepared.values, name)) === prepared.layout.scalar_type || throw(ArgumentError(
                "derived runtime parameter :$name must be a continuous numeric value"))
            return true
        end
        parameter.rank == 1 && return only(parameter.storage_labels)
        return parameter.applicability
    end)
end

function _active_trait_set(bgc, prepared, selection)
    slots, tau, labels = _trait_slots(prepared, selection)
    runtime_selection = _trait_runtime_selection(prepared, keys(selection))
    active = active_parameters(bgc; runtime_selection...)
    numerical = _prepare_trait_realization(prepared, keys(selection), active)
    T = prepared.layout.scalar_type
    f = TraitParameterMap{T,typeof(numerical.numerics),typeof(slots),typeof(numerical.gather)}(
        numerical.numerics, slots, numerical.gather
    )
    return ActiveTraitSet(bgc, active, tau, labels, selection,
                          prepared.explicit_override_keys, f, prepared)
end

function _check_trait_architecture(grid, arch)
    (isnothing(arch) || arch isa CPU) && (isnothing(grid) || architecture(grid) isa CPU) ||
        throw(ArgumentError("active_traits requires a CPU model for construction-trait sensitivity"))
    return nothing
end

"""
    active_traits(definition::ModelDefinition; trait_selection, kwargs...)
    active_traits(family::AbstractModelFamily; trait_selection, kwargs...)
    active_traits(family::AbstractModelFamily, recipe::ModelRecipe; trait_selection, kwargs...)

Select continuous construction traits for sensitivity analysis in a fixed CPU model.
Construction keywords match `construct` for the corresponding source. Select scalar
traits with `true` and vector traits with tuples of realized SizeClass symbols, e.g.
`trait_selection=(specificity=(:Z_1,),)`. Trait labels and `values` follow selection
order; `active_parameters` contains the union of all applicable runtime entries of
live derived outputs in parameter-plan order. Explicit derived-output overrides stay fixed.

The returned [`ActiveTraitSet`](@ref) exposes `bgc`, `labels`, `values`,
`trait_selection`, `explicit_override_keys`, `active_parameters`, and `parameter_map`.
Use [`realize_trait_parameters`](@ref) for checked primal evaluation. `parameter_map(tau)` is the numerical map for AD at a checked trait
vector: it expects the same length and floating-point element type as `values` and
returns the same runtime vector, without repeating host-side validation. Treat the
trait set metadata and base BGC as fixed; copy `values` before changing the trait vector.
"""
function active_traits(
    definition::ModelDefinition; trait_selection::NamedTuple,
    plankton_pfts=nothing, parameter_overrides::NamedTuple=(;),
    sinking_tracers=nothing, open_bottom::Bool=true, grid=nothing, arch=nothing,
    scalar_type=nothing,
)
    _check_trait_architecture(grid, arch)
    bgc, _, prepared = _construct_process_definition(definition;
        plankton_pfts, parameter_overrides, sinking_tracers, open_bottom, grid, arch, scalar_type)
    return _active_trait_set(bgc, prepared, trait_selection)
end

function active_traits(
    family::AbstractModelFamily; trait_selection::NamedTuple, plankton_pfts::NamedTuple,
    parameter_overrides::NamedTuple=(;), setting_overrides::NamedTuple=(;),
    sinking_tracers=nothing, open_bottom::Bool=true, grid=nothing, arch=nothing,
    scalar_type=nothing,
)
    _check_trait_architecture(grid, arch)
    realization = (; plankton_pfts, parameter_overrides, setting_overrides, sinking_tracers, open_bottom)
    bgc, _, prepared = _construct_family_model(family, realization; grid, arch, scalar_type)
    return _active_trait_set(bgc, prepared, trait_selection)
end

function active_traits(
    family::AbstractModelFamily, recipe::ModelRecipe; trait_selection::NamedTuple,
    grid=nothing, arch=nothing, scalar_type=nothing,
)
    _check_trait_architecture(grid, arch)
    bgc, _, prepared = _construct_recipe(family, recipe; grid, arch, scalar_type)
    return _active_trait_set(bgc, prepared, trait_selection)
end

"""
    realize_trait_parameters(traits, tau) -> p

Validate a flat trait vector and realize its dependent runtime parameters in
`traits.active_parameters` order. Pass `p` and that active set to `Runtime.parameterized`
or `Runtime.ode_problem`, using `traits.bgc`. Each call checks the realized scientific
parameter domains and coupled constraints, and leaves the baseline unchanged.
"""
function realize_trait_parameters(traits::ActiveTraitSet, tau::AbstractVector)
    realization = traits.realization
    T = realization.layout.scalar_type
    length(tau) == length(traits.values) && eltype(tau) === T && all(isfinite, tau) ||
        throw(ArgumentError("expected $(length(traits.values)) finite trait values with eltype $T"))
    Base.require_one_based_indexing(tau)
    f = traits.parameter_map
    resolved = _realized_trait_values(f, tau)
    checked = _checked_parameter_values(realization.canonical, realization.layout, realization.plan,
        resolved, f.numerics.derivations, realization.constraints)
    p = _gather_runtime_vector(T, checked, f.gather)
    all(isfinite, p) || throw(ArgumentError("realized runtime parameter vector must be finite"))
    return p
end
