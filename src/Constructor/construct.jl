using OceanBioME: BoxModelGrid, setup_velocity_fields

using Adapt: adapt

using Oceananigans.Architectures: architecture, device, CPU, GPU

using Agate.Utils:
    AbstractBGCFactory,
    normalize_interactions,
    parse_community,
    validate_plankton_inputs

using Agate.Models:
    default_plankton_dynamics,
    default_biogeochem_dynamics

# Bring the `Agate.Models` module into scope for qualified calls like
# `Models.default_community` without relying on the parent module name.
import Agate.Models
# For qualified calls inside registry update helpers.
import Agate.Parameters
using Agate.Equations: Equation, expr, requirements, req, merge_requirements
using Agate.Equations: ParamVar

# Local construction-time parameter placeholder namespace.
#
# Equation builders accept a first argument `PV` and reference placeholders as `PV.<name>`.
# We construct `PV` as a `NamedTuple` of `ParamVar{:<name>}` values so equation authors can
# keep the terse `PV.foo[i]` syntax without relying on a global module namespace.
@inline function parameter_vars(names::AbstractVector{Symbol})
    vals = map(names) do name
        # Create ParamVar{name} at runtime.
        Core.apply_type(ParamVar, name)()
    end
    NT = NamedTuple{Tuple(names)}
    return NT(Tuple(vals))
end

using Agate.Parameters: resolve_runtime_parameters, parameter_registry

"""Apply `interactions` overrides by updating the parameter registry.

`interactions` is intended primarily for interaction matrices and other parameter overrides.

Strictness rules
----------------
- Keys must already exist in the registry. Unknown keys throw to catch typos early.
- To add new parameters, extend the registry explicitly with `extend_registry(...)`.

For any matrix-shaped parameter, concrete matrix values are validated to be size
`(ctx.n_total, ctx.n_total)`.
"""
function _apply_interactions_to_registry(
    ctx,
    registry,
    overrides::NamedTuple,
)
    isempty(overrides) && return registry

    # Validate keys and eager matrix sizes (for concrete matrices).
    n = ctx.n_total
    for (k, v) in pairs(overrides)
        spec = Parameters.lookup(registry, k)
        spec === nothing && throw(ArgumentError(
            "interactions: unknown parameter key :$k. Add it explicitly with extend_registry(...) or fix the typo.",
        ))

        if v isa AbstractMatrix
            spec.shape === :matrix || throw(ArgumentError(
                "interactions: :$k is a $(spec.shape) parameter, but a matrix value was provided.",
            ))
            (size(v, 1) == n && size(v, 2) == n) || throw(ArgumentError(
                "interactions: :$k must be size ($n,$n).",
            ))
        end
    end

    return Parameters.update_registry(registry; overrides...)
end

"""
    construct(factory::AbstractBGCFactory; kw...) -> bgc

Construct and compile a concrete biogeochemistry *instance* from a factory and
optional overrides.

Design principles
-----------------
- Structural defaults (plankton community size structure) are provided by
  `Models.default_community(factory)`.
- Parameter defaults are provided by `Parameters.parameter_registry(factory)`.
- User overrides flow through the registry (no separate `params` keyword).
- The returned instance is `Adapt.jl`-compatible (CPU <-> GPU).

Key keyword arguments
---------------------
- `grid=nothing`: optional grid used for sinking-velocity fields and for choosing
  the floating point type when interfacing with Oceananigans / OceanBioME.
  Precision is determined by `eltype(grid)`. When `grid` is not provided,
  Agate constructs a `Float64` instance.
- `arch=nothing`: `CPU()` or `GPU()`; when omitted and `grid` is provided, defaults
  to `architecture(grid)`.
- `community`: plankton community structure (size classes, diameters, PFT specs).
- `registry`: parameter registry (defaults/specs), typically updated/extended by the user.
- `interactions`: optional `NamedTuple` or function `(ctx)->NamedTuple` providing
  interaction-related parameter overrides (e.g. matrices).
"""
function construct(
    factory::AbstractBGCFactory;
    plankton_dynamics=default_plankton_dynamics(factory),
    biogeochem_dynamics=default_biogeochem_dynamics(factory),

    community=Models.default_community(factory),
    registry=parameter_registry(factory),

    arch=nothing,
    interactions::Union{Nothing,NamedTuple,Function}=nothing,

    sinking_tracers=nothing,
    grid=nothing,
    open_bottom::Bool=true,
)

    # ---------------------------------------------------------------------
    # Precision and architecture selection.
    #
    # When `grid` is provided, its element type is the source-of-truth for `FT`.
    # This mirrors OceanBioME / Oceananigans, where model precision is determined
    # by the grid / model configuration rather than an independent BGC keyword.
    # ---------------------------------------------------------------------

    # If sinking is requested and no grid was supplied, fall back to a BoxModelGrid.
    if isnothing(grid) && !isnothing(sinking_tracers)
        grid = BoxModelGrid()
    end

    if !isnothing(grid)
        FT = eltype(grid)
        arch_grid = architecture(grid)
        if isnothing(arch)
            arch = arch_grid
        elseif typeof(arch) !== typeof(arch_grid)
            throw(ArgumentError(
                "arch=$arch does not match architecture(grid)=$arch_grid. Architecture is determined by the grid; either omit arch or construct a grid for $arch.",
            ))
        end
    else
        FT = Float64
        isnothing(arch) && (arch = CPU())
    end

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

    # Construction-time parameter placeholder namespace.
    PV = parameter_vars(map(s -> s.name, final_registry.specs))

    # ---------------------------------------------------------------------
    # Build tracer expressions and collect parameter requirements.
    # ---------------------------------------------------------------------

    plankton_syms = ctx.plankton_symbols

    tracer_names = Symbol[collect(keys(biogeochem_dynamics))...]
    append!(tracer_names, plankton_syms)

    tracer_exprs = Expr[]
    merged = req()

    for (_, f) in pairs(biogeochem_dynamics)
        eq = f(PV, plankton_syms)
        eq isa Equation || throw(ArgumentError("biogeochem dynamics $(nameof(f)) must return Equation"))
        push!(tracer_exprs, expr(eq))
        merged = merge_requirements(merged, requirements(eq))
    end

    for idx in eachindex(plankton_syms)
        g = ctx.group_symbols[idx]
        f = getfield(plankton_dynamics, g)
        tr = plankton_syms[idx]

        eq = f(PV, plankton_syms, tr, idx)
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

    # Move any arrays inside `bgc` onto the requested architecture.
    bgc = adapt(device(arch), bgc)

    return bgc
end
