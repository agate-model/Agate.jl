using OceanBioME: BoxModelGrid, setup_velocity_fields

using Adapt: adapt

# We use Oceananigans' architecture abstraction (CPU/GPU + backend-specific array types).
# Importing from a submodule does not bind the `Oceananigans` name, and we also
# reference `Oceananigans.Architectures` directly below.
import Oceananigans

using Oceananigans.Architectures: architecture, CPU, GPU

using ..Utils:
    AbstractBGCFactory,
    parameter_spec,
    normalize_interactions,
    resolve_derived_matrices,
    finalize_interaction_parameters,
    parse_community,
    validate_plankton_inputs

using ..FactoryInterface:
    default_plankton_dynamics,
    default_biogeochem_dynamics,
    default_community

using ..Functors: CompiledEquation, requirements, req, merge_requirements

"""Return factory-specific default runtime parameters.

Factories can extend this method to compute default parameters from the parsed
community (`ctx`) and the target float type (`FT`).

The returned `NamedTuple` is expected to contain *resolved* values (scalars,
vectors, matrices), not higher-level providers.
"""
default_parameters(::AbstractBGCFactory, ctx, ::Type{FT}) where {FT} = (;)

"""Move `x` to the requested Oceananigans architecture."""
function _on_architecture(arch, x)
    arch === nothing && return x
    arch isa CPU && return x

    return Base.invokelatest(Adapt.adapt, _architecture_array_type(arch), x)
end

"""Return the preferred array storage type for `arch`."""
function _architecture_array_type(arch)
    arch isa CPU && return Array
    arch isa GPU && return Oceananigans.Architectures.array_type(arch)
    return Array
end

@inline function _required_keys(r)
    return (r.scalars..., r.vectors..., r.matrices...)
end

function _validate_parameter_directory(factory::AbstractBGCFactory, r)
    # Ensure every required parameter has explicit metadata.
    for k in _required_keys(r)
        spec = parameter_spec(factory, k)
        spec === nothing && throw(ArgumentError(
            "Factory $(typeof(factory)) is missing a ParameterSpec for required parameter :$k. " *
            "Add it to parameter_directory(::$(typeof(factory))).",
        ))
    end

    # Ensure the directory is internally consistent with compiled-equation requirements.
    for k in r.scalars
        spec = parameter_spec(factory, k)
        spec.shape === :scalar || throw(ArgumentError(
            "parameter_directory(::$(typeof(factory))) declares :$k as shape $(spec.shape), but compiled equations require a scalar.",
        ))
    end

    for k in r.vectors
        spec = parameter_spec(factory, k)
        spec.shape === :vector || throw(ArgumentError(
            "parameter_directory(::$(typeof(factory))) declares :$k as shape $(spec.shape), but compiled equations require a vector.",
        ))
    end

    for k in r.matrices
        spec = parameter_spec(factory, k)
        spec.shape === :matrix || throw(ArgumentError(
            "parameter_directory(::$(typeof(factory))) declares :$k as shape $(spec.shape), but compiled equations require a matrix.",
        ))
    end

    return nothing
end

function _validate_parameter_shapes(ctx, params::NamedTuple, r)
    n = ctx.n_total

    for k in r.vectors
        v = getproperty(params, k)
        length(v) == n || throw(ArgumentError(
            "parameter :$k must have length $n (got $(length(v))).",
        ))
    end

    for k in r.matrices
        m = getproperty(params, k)
        (size(m, 1) == n && size(m, 2) == n) || throw(ArgumentError(
            "parameter :$k must have size ($n,$n) (got $(size(m))).",
        ))
    end

    return nothing
end

function _validate_override_keys(where_, overrides::NamedTuple, required::Tuple, factory::AbstractBGCFactory)
    isempty(overrides) && return nothing

    required_set = Set(required)
    for k in keys(overrides)
        k in required_set && continue

        # Optional overrides are permitted only for keys declared in the
        # factory's parameter directory.
        parameter_spec(factory, k) === nothing && throw(ArgumentError(
            "$(where_): unknown parameter key :$k.",
        ))
    end

    return nothing
end

@inline _contains_missing(x) = x === missing

function _contains_missing(x::NamedTuple)
    for v in values(x)
        _contains_missing(v) && return true
    end
    return false
end

function _contains_missing(x::AbstractArray)
    return any(ismissing, x)
end

function _reject_missing_values(params::NamedTuple)
    for (k, v) in pairs(params)
        _contains_missing(v) && throw(ArgumentError(
            "parameter :$k contains `missing`; all required parameters must be explicitly defined.",
        ))
    end
    return nothing
end

"""
    construct_factory(factory::AbstractBGCFactory; kw...) -> bgc

Construct and compile a concrete biogeochemistry *instance* from a factory and
optional overrides.

Key keyword arguments
---------------------
- `grid=nothing`: optional grid used for sinking-velocity fields and for choosing
  the floating point type when interfacing with Oceananigans / OceanBioME.
  Precision is determined by `eltype(grid)`. When `grid` is not provided,
  Agate constructs a `Float64` instance.
- `arch=nothing`: `CPU()` or `GPU()`; when omitted and `grid` is provided, defaults
  to `architecture(grid)`.
- `community`: plankton community structure (size classes, diameters, PFT specs).
- `parameters`: `NamedTuple` of fully-resolved parameter overrides.
- `interactions`: optional `NamedTuple` of interaction parameter overrides (often matrices such as `:palatability_matrix` and `:assimilation_matrix`).
  Values may be concrete objects or provider functions callable as `f(ctx)`.
  For matrix parameters, overrides may be full `(n_total, n_total)` matrices. A group-block `(n_groups, n_groups)` matrix may be supplied and expanded during construction; when the parameter declares role-aware axes, wrap the block matrix as `GroupBlockMatrix(B)` to avoid ambiguity. When axes are declared, rectangular consumer-by-prey matrices sized to those axes (for example `(n_consumer, n_prey)`) are also accepted, as are axis-local group-block matrices.
"""
function construct_factory(
    factory::AbstractBGCFactory;
    plankton_dynamics = default_plankton_dynamics(factory),
    biogeochem_dynamics = default_biogeochem_dynamics(factory),
    community = default_community(factory),
    parameters::NamedTuple = (;),
    interactions::Union{Nothing,NamedTuple} = nothing,
    arch = nothing,
    sinking_tracers = nothing,
    grid = nothing,
    open_bottom::Bool = true,
)

    # If sinking is requested and no grid was supplied, fall back to a BoxModelGrid.
    if isnothing(grid) && !isnothing(sinking_tracers)
        grid = BoxModelGrid()
    end

    # Precision + architecture selection.
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

    # Parse community.
    ctx = parse_community(
        factory,
        FT,
        community;
        plankton_dynamics = plankton_dynamics,
        biogeochem_dynamics = biogeochem_dynamics,
    )

    # ---------------------------------------------------------------------
    # Build tracer expressions and collect parameter requirements.
    # ---------------------------------------------------------------------

    plankton_syms = ctx.plankton_symbols

    tracer_names = Symbol[collect(keys(biogeochem_dynamics))...]
    append!(tracer_names, plankton_syms)

    tracer_defs = Any[]
    merged = req()

    for (_, f) in pairs(biogeochem_dynamics)
        tr = f(plankton_syms)
        (tr isa CompiledEquation) || throw(ArgumentError(
            "biogeochem dynamics $(nameof(f)) must return CompiledEquation",
        ))
        push!(tracer_defs, tr)
        merged = merge_requirements(merged, requirements(tr))
    end

    for idx in eachindex(plankton_syms)
        g = ctx.group_symbols[idx]
        f = getfield(plankton_dynamics, g)
        trsym = plankton_syms[idx]

        tr = f(plankton_syms, trsym, idx)
        (tr isa CompiledEquation) || throw(ArgumentError(
            "plankton dynamics $(nameof(f)) must return CompiledEquation",
        ))
        push!(tracer_defs, tr)

        merged = merge_requirements(merged, requirements(tr))
    end

    tracers = NamedTuple{Tuple(tracer_names)}(Tuple(tracer_defs))

    # ---------------------------------------------------------------------
    # Parameters
    # ---------------------------------------------------------------------

    _validate_parameter_directory(factory, merged)

    required = _required_keys(merged)

    overrides = normalize_interactions(factory, ctx, interactions)

    _validate_override_keys("parameters", parameters, required, factory)
    _validate_override_keys("interactions", overrides, required, factory)

    defaults = default_parameters(factory, ctx, FT)

    # Merge precedence: defaults < parameters < interactions
    params_full = merge(defaults, parameters, overrides)

    # Recompute any derived matrices affected by explicit trait overrides.
    explicit_override_keys = (keys(parameters)..., keys(overrides)...)
    params_full = resolve_derived_matrices(factory, ctx, params_full, explicit_override_keys)

    # Ensure the required keys are present.
    missing = Symbol[]
    for k in required
        hasproperty(params_full, k) || push!(missing, k)
    end
    isempty(missing) || throw(ArgumentError(
        "missing required parameters: $(join(string.(missing), ", "))",
    ))

    # Finalize any role-aware interaction matrices.
    params_full = finalize_interaction_parameters(factory, ctx, params_full)

    # Slice down to the required keys (plus selected internal helpers) for better type stability.
    internal = hasproperty(params_full, :interactions) ? (:interactions,) : ()
    all_keys = (required..., internal...)
    params_req = NamedTuple{all_keys}(Tuple(getproperty(params_full, k) for k in all_keys))

    _reject_missing_values(params_req)

    _validate_parameter_shapes(ctx, params_req, merged)

    # ---------------------------------------------------------------------
    # Compile + instantiate.
    # ---------------------------------------------------------------------

    if isnothing(sinking_tracers)
        bgc_factory = define_tracer_functions(params_req, tracers)
        bgc = bgc_factory(params_req)
    else
        sinking_velocities = setup_velocity_fields(sinking_tracers, grid, open_bottom)
        bgc_factory = define_tracer_functions(params_req, tracers; sinking_velocities = sinking_velocities)
        bgc = bgc_factory(params_req, sinking_velocities)
    end

    # Move any arrays inside `bgc` onto the requested architecture.
    bgc = _on_architecture(arch, bgc)

    return bgc
end
