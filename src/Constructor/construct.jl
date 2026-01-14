using OceanBioME: BoxModelGrid, setup_velocity_fields

import Adapt

using ..Utils:
    AbstractBGCFactory,
    define_tracer_functions,
    normalize_interactions,
    parse_community,
    validate_plankton_inputs

using ..Models:
    default_plankton_dynamics,
    default_biogeochem_dynamics

using ..Library.Equations: Equation, expr, requirements, req, merge_requirements
using ..Library.Allometry: allometric_palatability_unimodal_protection

using ..Parameters: resolve_runtime_parameters

# Bring the `Agate.Models` module into scope for qualified calls like
# `Models.default_parameter_args` without relying on the parent module name.
import ..Models

# -----------------------------------------------------------------------------
# Args helper
# -----------------------------------------------------------------------------

"""A lightweight construction-time bundle for parameter resolution.

`ParameterRegistryArgs` is CPU-only. It holds:

- `community`: the plankton community specification (size structure + PFT overrides)
- `params`: keyword overrides for parameters declared in the registry
- `backend`: optional `Adapt` backend (e.g. `CUDA.CuArray`) to move runtime arrays

The resolved runtime parameter bundle remains GPU-safe (plain scalars and arrays).
"""
struct ParameterRegistryArgs{C,O,B,R}
    community::C
    params::O
    backend::B
    registry::R
end

"""Create a default parameter argument bundle.

This is the recommended way to supply parameter overrides without duplicating defaults outside the registry.
"""
function default_parameter_args(factory::AbstractBGCFactory;
    FT::Type{<:AbstractFloat},
    community = Models.default_parameter_args(factory; FT=FT),
    params::NamedTuple=NamedTuple(),
    backend=nothing,
    registry=nothing,
)
    return ParameterRegistryArgs(community, params, backend, registry)
end

# -----------------------------------------------------------------------------
# Public constructor
# -----------------------------------------------------------------------------

"""
    construct(factory::AbstractBGCFactory; kw...) -> Type

Compile a concrete biogeochemistry type from a factory and optional overrides.

The return value is a *type*; instantiate with `bgc = bgc_type()`.

Key keyword arguments
---------------------

- `plankton_dynamics`: group dynamics builders (must return `Equation`)
- `biogeochem_dynamics`: non-plankton tracer dynamics builders (must return `Equation`)
- `parameter_args`: a `ParameterRegistryArgs` (from `default_parameter_args`) or a raw community `NamedTuple`
- `interactions`: optional interaction matrix overrides
- `sinking_tracers`: optional sinking velocity specification

Registry-owned defaults
-----------------------

All parameter defaults are sourced from the model's parameter registry. This
constructor only supplies structural defaults (community sizes/diameters) via
`default_parameter_args(factory; FT=FT)`.
"""
function construct(
    factory::AbstractBGCFactory;
    FT::Type{<:AbstractFloat},

    plankton_dynamics=default_plankton_dynamics(factory),
    biogeochem_dynamics=default_biogeochem_dynamics(factory),

    # Parameter inputs (community + overrides + backend)
    parameter_args=default_parameter_args(factory; FT=FT),

    # Optional interaction matrices (NamedTuple/InteractionMatrices/Function)
    interactions=nothing,

    # Optional sinking velocities
    sinking_tracers=nothing,
    grid=BoxModelGrid(),
    open_bottom::Bool=true,

    # Default palatability rule used by registry defaults
    palatability_fn=allometric_palatability_unimodal_protection,

    # Convenience explicit matrix overrides (take precedence over `interactions`)
    palatability_matrix=nothing,
    assimilation_efficiency_matrix=nothing,
)
    return _construct_impl(
        factory,
        FT,
        plankton_dynamics,
        biogeochem_dynamics,
        parameter_args,
        interactions,
        sinking_tracers,
        grid,
        open_bottom,
        palatability_fn,
        palatability_matrix,
        assimilation_efficiency_matrix,
    )
end

# -----------------------------------------------------------------------------
# Internal implementation
# -----------------------------------------------------------------------------

@inline function _as_registry_args(factory::AbstractBGCFactory, parameter_args, FT)
    if parameter_args isa ParameterRegistryArgs
        return parameter_args
    elseif parameter_args isa NamedTuple
        return ParameterRegistryArgs(parameter_args, NamedTuple(), nothing, nothing)
    else
        throw(ArgumentError("parameter_args must be a ParameterRegistryArgs or a community NamedTuple, got $(typeof(parameter_args))"))
    end
end

function _construct_impl(
    factory::AbstractBGCFactory,
    FT::Type{<:AbstractFloat},
    plankton_dynamics,
    biogeochem_dynamics,
    parameter_args,
    interactions,
    sinking_tracers,
    grid,
    open_bottom::Bool,
    palatability_fn,
    palatability_matrix,
    assimilation_efficiency_matrix,
)

    # Normalize parameter_args and validate structural inputs.
    argpack = _as_registry_args(factory, parameter_args, FT)
    community = argpack.community

    validate_plankton_inputs(plankton_dynamics, community)
    biogeochem_dynamics isa NamedTuple || throw(ArgumentError("biogeochem_dynamics must be a NamedTuple"))

    # Parse community and normalize any explicit interaction matrices.
    ctx = parse_community(
        FT,
        community;
        plankton_dynamics=plankton_dynamics,
        biogeochem_dynamics=biogeochem_dynamics,
    )

    explicit_mats = normalize_interactions(factory, FT, ctx, interactions)

    # ---------------------------------------------------------------------
    # Build tracer expressions and collect parameter requirements.
    # ---------------------------------------------------------------------

    plankton_syms = ctx.plankton_symbols

    tracer_names = Symbol[collect(keys(biogeochem_dynamics))...]
    append!(tracer_names, plankton_syms)

    tracer_exprs = Expr[]
    merged = req()

    for (_, f) in pairs(biogeochem_dynamics)
        eq = f(plankton_syms)
        eq isa Equation || throw(ArgumentError("biogeochem dynamics $(nameof(f)) must return Equation"))
        push!(tracer_exprs, expr(eq))
        merged = merge_requirements(merged, requirements(eq))
    end

    for idx in eachindex(plankton_syms)
        g = ctx.group_symbols[idx]
        f = getfield(plankton_dynamics, g)
        tr = plankton_syms[idx]

        eq = f(plankton_syms, tr, idx)
        eq isa Equation || throw(ArgumentError("plankton dynamics $(nameof(f)) must return Equation"))
        push!(tracer_exprs, expr(eq))

        merged = merge_requirements(merged, requirements(eq))
    end

    tracers = NamedTuple{Tuple(tracer_names)}(Tuple(tracer_exprs))

    # ---------------------------------------------------------------------
    # Parameter resolution (CPU) -> optional backend adaptation (CPU->GPU).
    # ---------------------------------------------------------------------

    # Apply explicit matrix overrides with precedence:
    # keyword args > `interactions` > registry defaults
    overrides = argpack.params

    # Merge explicit interactions matrices into overrides.
    overrides = merge(overrides, explicit_mats)

    if palatability_matrix !== nothing
        overrides = merge(overrides, (; palatability_matrix = palatability_matrix))
    end
    if assimilation_efficiency_matrix !== nothing
        overrides = merge(overrides, (; assimilation_efficiency_matrix = assimilation_efficiency_matrix))
    end

    params = resolve_runtime_parameters(
        factory,
        ctx,
        merged;
        FT=FT,
        overrides=overrides,
        palatability_fn=palatability_fn,
        registry=argpack.registry,
    )

	# Optional backend adaptation (CPU->GPU) at the final construction stage.
	if argpack.backend != nothing
		params = Adapt.adapt(argpack.backend, params)
	end

    # Optional sinking velocities.
    if isnothing(sinking_tracers)
        return define_tracer_functions(params, tracers)
    end

    sinking_velocities = setup_velocity_fields(sinking_tracers, grid, open_bottom)
    return define_tracer_functions(params, tracers; sinking_velocities=sinking_velocities)
end
