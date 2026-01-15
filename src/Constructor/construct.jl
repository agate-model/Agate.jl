using OceanBioME: BoxModelGrid, setup_velocity_fields

using Adapt: adapt

using Oceananigans.Architectures: architecture, device, CPU, GPU



using Agate.Utils:    AbstractBGCFactory,
    normalize_interactions,
    parse_community,
    validate_plankton_inputs

using Agate.Models:    default_plankton_dynamics,
    default_biogeochem_dynamics

# Bring the `Agate.Models` module into scope for qualified calls like
# `Models.default_community` without relying on the parent module name.
import Agate.Models
# For qualified calls inside registry update helpers.
import Agate.Parameters
using Agate.Equations: Equation, expr, requirements, req, merge_requirements
using Agate.Equations: declare_parameter_vars!
using Agate.Library.Allometry: allometric_palatability_unimodal_protection
import Agate.ParamVars
using Agate.Parameters:    resolve_runtime_parameters,
    parameter_registry

"""Apply `interactions` overrides by updating/extending the parameter registry.

`interactions` is intended for interaction matrices and related knobs. Any
provided matrix key must end in `_matrix`. Unknown keys error unless they are
new `_matrix` keys, in which case a matrix `ParamSpec` is appended.
"""
function _apply_interactions_to_registry(ctx, registry, overrides::NamedTuple)
    isempty(overrides) && return registry

    update_pairs = Pair{Symbol,Any}[]
    new_specs = Parameters.ParamSpec[]
    known = Set(spec.name for spec in registry.specs)

    for (k, v) in pairs(overrides)
        k_str = String(k)
        is_matrix_val = v isa AbstractMatrix
        if is_matrix_val && !endswith(k_str, "_matrix")
            throw(ArgumentError("interactions: matrix key :$k must end with `_matrix`"))
        end

        # Basic shape validation for interaction matrices.
        if is_matrix_val && endswith(k_str, "_matrix")
            n = ctx.n_total
            (size(v, 1) == n && size(v, 2) == n) || throw(ArgumentError("interactions: :$k must be size ($n,$n)"))
        end

        if k in known
            push!(update_pairs, k => v)
        else
            # Allow extending the registry with new matrices via `interactions`.
            if endswith(k_str, "_matrix")
                push!(new_specs, Parameters.ParamSpec(k, "Interaction matrix provided via `interactions`.", v; missing_policy=:zero_silent, value_kind=:real))
                push!(known, k)
            else
                throw(ArgumentError("interactions: unknown parameter key :$k (not present in registry)."))
            end
        end
    end

    reg = registry
    if !isempty(new_specs)
        reg = Parameters.extend_registry(reg, new_specs...)
    end
    if !isempty(update_pairs)
        nt = (; update_pairs...)
        reg = Parameters.update_registry(reg; nt...)
    end

    return reg
end

"""
    construct(factory::AbstractBGCFactory; kw...) -> bgc

Construct and compile a concrete biogeochemistry *instance* from a factory and
optional overrides.

Design principles
-----------------
- Structural defaults (plankton community size structure) are provided by
  `Models.default_community(factory; FT=...)`.
- Parameter defaults are provided by `Parameters.parameter_registry(factory)`.
- User overrides flow through the registry (no separate `params` keyword).
- The returned instance is `Adapt.jl`-compatible (CPU <-> GPU).

Key keyword arguments
---------------------
- `FT::Type{<:AbstractFloat}=Float64`: floating point type for runtime values.
- `community`: plankton community structure (size classes, diameters, PFT specs).
- `registry`: parameter registry (defaults/specs), typically updated/extended by the user.
- `interactions`: optional NamedTuple or function `(ctx)->NamedTuple` providing
  interaction-related parameter overrides (e.g. matrices). Unknown keys should error.
- `arch=CPU()`: `CPU()` or `GPU()` used to set architecture.
"""
function construct(
    factory::AbstractBGCFactory;
    FT::Type{<:AbstractFloat}=Float64,

    plankton_dynamics=default_plankton_dynamics(factory),
    biogeochem_dynamics=default_biogeochem_dynamics(factory),

    community = Models.default_community(factory; FT=FT),
    registry = parameter_registry(factory),

    arch = CPU(),
    interactions::Union{Nothing,NamedTuple,Function}=nothing,

    sinking_tracers=nothing,
    grid=BoxModelGrid(),
    open_bottom::Bool=true,

    palatability_fn=allometric_palatability_unimodal_protection,
)
    validate_plankton_inputs(plankton_dynamics, community)
    biogeochem_dynamics isa NamedTuple || throw(ArgumentError("biogeochem_dynamics must be a NamedTuple"))

    # Parse community and normalize interaction overrides.
    ctx = parse_community(
        FT,
        community;
        plankton_dynamics=plankton_dynamics,
        biogeochem_dynamics=biogeochem_dynamics,
    )

    overrides = normalize_interactions(factory, FT, ctx, interactions)
    final_registry = _apply_interactions_to_registry(ctx, registry, overrides)

    # Ensure canonical parameter placeholders exist for the final registry.
    # Equation authoring (Library + user extensions) should reference parameters as
    # `PV.<name>` where `const PV = Agate.ParamVars`.
    declare_parameter_vars!(ParamVars, map(s -> s.name, final_registry.specs); export_vars=false)

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
    # Parameter resolution -> architecture adaptation (CPU/GPU).
    # ---------------------------------------------------------------------

    params = resolve_runtime_parameters(
        factory,
        ctx,
        merged;
        FT=FT,
        palatability_fn=palatability_fn,
        registry=final_registry,
    )

    # Optional sinking velocities.
    if isnothing(sinking_tracers)
        bgc_type = define_tracer_functions(params, tracers)
        # NOTE: bgc_type is created via eval in `define_tracer_functions`, which can
        # trigger Julia's world-age restriction if we call it immediately.
        bgc = Base.invokelatest(bgc_type, params)
    else
        sinking_velocities = setup_velocity_fields(sinking_tracers, grid, open_bottom)
        bgc_type = define_tracer_functions(params, tracers; sinking_velocities=sinking_velocities)
        bgc = Base.invokelatest(bgc_type, params, sinking_velocities)
    end

    # move any arrays inside `bgc` onto the requested architecture.
    bgc = adapt(device(arch), bgc)

    return bgc
end
