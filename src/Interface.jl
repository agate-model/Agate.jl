"""Model/Factory interface for defaults and group-level hooks.

Agate's public model constructors (`NiPiZD.construct`, `DARWIN.construct`, …) are
implemented as thin wrappers around factories.

Factories provide explicit defaults for:
- community structure (groups, size structure)
- tracer dynamics builders

This module also defines optional *group-level hook points* intended for
OceanBioME-style extensibility.

GPU notes
---------
Hook dispatch uses `Val{:Group}` so downstream code can specialize without any
runtime `Symbol` lookups inside kernels.

These hooks are optional; if unused they impose no kernel cost.
"""

module Interface

using ..Utils: AbstractBGCFactory, ParameterSpec, parameter_directory

export ParameterSpec
export parameter_directory
export default_plankton_dynamics
export default_community
export default_biogeochem_dynamics

export sinking_velocity
export grazing_kernel

"""Default plankton dynamics for a factory.

Returns a `NamedTuple` mapping group symbols (e.g. `:Z`, `:P`) to dynamics builder
functions.
"""
function default_plankton_dynamics(::AbstractBGCFactory)
    throw(
        ArgumentError(
            "No method `default_plankton_dynamics(factory)` is defined for this factory."
        ),
    )
end

"""Default plankton community structure for a factory.

Returns a `NamedTuple` mapping group symbols to group specifications.

This is structural information only (group symbols, diameter specifications,
PFT specifications, etc.). Numeric parameter defaults are sourced from the
factory's default parameter generator (`Constructor.default_parameters`).
"""
function default_community(::AbstractBGCFactory)
    throw(
        ArgumentError("No method `default_community(factory)` is defined for this factory.")
    )
end

"""Default non-plankton tracer dynamics for a factory.

Returns a `NamedTuple` mapping tracer symbols (e.g. `:N`, `:DIC`) to dynamics
builder functions.
"""
function default_biogeochem_dynamics(::AbstractBGCFactory)
    throw(
        ArgumentError(
            "No method `default_biogeochem_dynamics(factory)` is defined for this factory."
        ),
    )
end

# -----------------------------------------------------------------------------
# Optional group-level hooks
# -----------------------------------------------------------------------------

"""Return a sinking velocity for a plankton group/class.

    sinking_velocity(model_or_factory, ::Val{G}, class_id::Int, ctx, params) -> Union{Nothing,Number}

- `G` is a group symbol (e.g. `:P`, `:Z`).
- `class_id` is the within-group ordinal (1..n_group).
- Return `nothing` to indicate no sinking.

This hook is intended for host-side construction of sinking-velocity fields.
"""
@inline sinking_velocity(::Any, ::Val, ::Int, ctx, params) = nothing

"""Group-pair hook for selecting a grazing/predation kernel.

    grazing_kernel(model_or_factory, ::Val{PredG}, ::Val{PreyG}, predator_class::Int, prey_class::Int, ctx, params)

Implementations should return a callable `k(prey, predator)` that can be used in
kernel-callable code.

`predator_class` and `prey_class` are within-group ordinals.

The default throws to make missing implementations explicit when you opt in.
"""
function grazing_kernel(::Any, ::Val, ::Val, ::Int, ::Int, ctx, params)
    throw(
        ArgumentError(
            "No method `grazing_kernel(model, ::Val{PredG}, ::Val{PreyG}, ...)` is defined for this model/factory.",
        ),
    )
end

end # module Interface
