using OceanBioME: BoxModelGrid, setup_velocity_fields

using Adapt: adapt

# We use Oceananigans' architecture abstraction (CPU/GPU + backend-specific array types).
# Importing from a submodule does not bind the `Oceananigans` name, and we also
# reference `Oceananigans.Architectures` directly below.
import Oceananigans

using Oceananigans.Architectures: architecture, CPU, GPU

using ..Utils:
    AbstractBGCFactory,
    normalize_interactions,
    parse_community,
    validate_plankton_inputs

using ..FactoryInterface:
    default_plankton_dynamics,
    default_biogeochem_dynamics,
    default_community
# For qualified calls inside registry update helpers.
import ..Parameters
using ..Equations: Equation, expr, requirements, req, merge_requirements
using ..Equations: ParamVar

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

using ..Parameters: resolve_runtime_parameters, parameter_registry

"""Move `x` to the requested Oceananigans architecture.

This is a construction-time helper that uses `Oceananigans.Architectures.array_type(arch)`
to pick a storage type (e.g. `Array` on CPU, `CUDA.CuArray` on CUDA GPU) and then
`Adapt.adapt` to recursively move arrays inside `x`.
"""
function _on_architecture(arch, x)
    arch === nothing && return x
    arch isa CPU && return x

    # We deliberately avoid `Oceananigans.Architectures.on_architecture` here.
    # Oceananigans only defines `on_architecture` for its own types, and the
    # generic fallback for unknown structs can be a no-op. Instead, we always
    # Adapt directly to the architecture's preferred array storage type.
    # IMPORTANT: `construct` generates the biogeochemistry type (and its
    # `Adapt.@adapt_structure` methods) at runtime via `eval`. Without
    # `invokelatest` here, Julia's world-age restriction can cause `Adapt.adapt`
    # to miss those freshly-defined methods and silently behave like a no-op.
    return Base.invokelatest(Adapt.adapt, _architecture_array_type(arch), x)
end

"""Return the preferred array storage type for `arch`.

This delegates to Oceananigans' public `array_type(arch)` API.
"""
function _architecture_array_type(arch)
    arch isa CPU && return Array
    arch isa GPU && return Oceananigans.Architectures.array_type(arch)
    return Array
end

"""Apply `interactions` overrides by updating the parameter registry.

`interactions` is intended primarily for interaction matrices and other parameter overrides.

Strictness rules
----------------
- Keys must already exist in the registry. Unknown keys throw to catch typos early.

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
            "interactions: unknown parameter key :$k. Fix the typo or add the parameter to the model's parameter registry.",
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
  `default_community(factory)`.
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

    community=default_community(factory),
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

    # Parameter resolution is independent of the concrete factory type.
    # Everything needed to resolve scalars/vectors/matrices is encoded in:
    #   - `ctx` (community + trait metadata)
    #   - `final_registry` (declared parameter providers)
    #   - `merged` (equation requirements)
    params = resolve_runtime_parameters(ctx, final_registry, merged, FT)

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
    bgc = _on_architecture(arch, bgc)

    return bgc
end
