using OceanBioME: BoxModelGrid, setup_velocity_fields

using Adapt: adapt

# We use Oceananigans' architecture abstraction (CPU/GPU + backend-specific array types).
# Importing from a submodule does not bind the `Oceananigans` name, and we also
# reference `Oceananigans.Architectures` directly below.
import Oceananigans

using Oceananigans.Architectures: architecture, CPU, GPU

using ..Factories: AbstractBGCFactory

using ..Utils:
    parameter_definitions,
    ConstDefault,
    NoDefault,
    FillDefault,
    DiameterIndexedVectorDefault,
    parameter_spec,
    axis_indices,
    normalize_interaction_overrides,
    resolve_derived_matrices,
    finalize_interaction_parameters,
    parse_community,
    build_tracer_index,
    validate_community_inputs

using ..Interface: default_plankton_dynamics, default_biogeochem_dynamics, default_community

using ..Equations: CompiledEquation, requirements, EquationRequirements, merge_requirements

using ..Library.Allometry: resolve_diameter_indexed_vector

"""A factory-supplied callable that builds a non-plankton tracer equation.

Biogeochemical dynamics builders are stored in `biogeochem_dynamics` and are called
once per tracer symbol during construction.

Expected signature
------------------
`builder() -> CompiledEquation`.
"""
const BiogeochemDynamicsBuilder = Function

"""A factory-supplied callable that builds a plankton tracer equation.

Plankton dynamics builders are stored in `plankton_dynamics` under their group
symbol (e.g. `P`, `Z`) and are called once per plankton class.

Expected signature
------------------
`builder(global_index::Int) -> CompiledEquation`

Arguments
---------
- `global_index`: global plankton class index (ordered as in `community_context.plankton_symbols`)
"""
const PlanktonDynamicsBuilder = Function

"""Evaluate `parameter_definitions(factory)` to produce baseline parameter defaults.

Defaults are evaluated on the host during model construction. The returned
`NamedTuple` maps parameter keys to concrete values (scalars, vectors, matrices).

Parameters declared with `NoDefault()` are omitted; they are expected to be
provided by user overrides or derived-matrix providers later in construction.
"""
function build_parameter_defaults(
    factory::AbstractBGCFactory, community_context, ::Type{FT}
) where {FT}
    defs = parameter_definitions(factory)
    isempty(defs) && return (;)

    keys_ = map(d -> d.spec.name, defs)
    length(unique(keys_)) == length(keys_) || throw(
        ArgumentError(
            "parameter_definitions(::$(typeof(factory))) contains duplicate keys.",
        ),
    )

    pairs = Pair{Symbol,Any}[]
    for def in defs
        spec = def.spec
        provider = def.default
        provider isa NoDefault && continue
        value = _evaluate_default(provider, spec, factory, community_context, FT)
        push!(pairs, spec.name => value)
    end

    return (; pairs...)
end

@inline function _evaluate_default(
    provider::ConstDefault, spec, ::AbstractBGCFactory, ::Any, ::Type{FT}
) where {FT}
    spec.shape === :scalar || throw(
        ArgumentError("ConstDefault can only be used for scalar parameters (:$((spec.name)))."),
    )
    v = provider.value
    return v isa Bool ? v : FT(v)
end

@inline function _evaluate_default(
    provider::FillDefault, spec, ::AbstractBGCFactory, community_context, ::Type{FT}
) where {FT}
    spec.shape in (:vector, :matrix) || throw(
        ArgumentError("FillDefault can only be used for vector or matrix parameters (:$((spec.name)))."),
    )

    v = provider.value
    v = v isa Bool ? v : FT(v)

    if spec.shape === :vector
        return fill(v, community_context.n_total)
    end

    if spec.axes === nothing
        n = community_context.n_total
        return fill(v, n, n)
    end

    row_axis, col_axis = spec.axes
    nr = length(axis_indices(community_context, row_axis))
    nc = length(axis_indices(community_context, col_axis))
    return fill(v, nr, nc)
end

@inline function _evaluate_default(
    provider::DiameterIndexedVectorDefault,
    spec,
    ::AbstractBGCFactory,
    community_context,
    ::Type{FT},
) where {FT}
    spec.shape === :vector || throw(
        ArgumentError(
            "DiameterIndexedVectorDefault can only be used for vector parameters (:$((spec.name))).",
        ),
    )

    indices = getproperty(community_context, provider.indices_field)
    default = FT(provider.default)
    return resolve_diameter_indexed_vector(
        FT,
        community_context.diameters,
        indices,
        provider.value;
        default=default,
    )
end

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
        spec === nothing && throw(
            ArgumentError(
                "Factory $(typeof(factory)) is missing a ParameterSpec for required parameter :$k. " *
                "Add it to parameter_definitions(::$(typeof(factory))).",
            ),
        )
    end

    # Ensure the directory is internally consistent with compiled-equation requirements.
    for k in r.scalars
        spec = parameter_spec(factory, k)
        spec.shape === :scalar || throw(
            ArgumentError(
                "parameter_definitions(::$(typeof(factory))) declares :$k as shape $(spec.shape), but compiled equations require a scalar.",
            ),
        )
    end

    for k in r.vectors
        spec = parameter_spec(factory, k)
        spec.shape === :vector || throw(
            ArgumentError(
                "parameter_definitions(::$(typeof(factory))) declares :$k as shape $(spec.shape), but compiled equations require a vector.",
            ),
        )
    end

    for k in r.matrices
        spec = parameter_spec(factory, k)
        spec.shape === :matrix || throw(
            ArgumentError(
                "parameter_definitions(::$(typeof(factory))) declares :$k as shape $(spec.shape), but compiled equations require a matrix.",
            ),
        )
    end

    return nothing
end

function _validate_parameter_shapes(
    factory::AbstractBGCFactory, context, params::NamedTuple, r
)
    n = context.n_total

    for k in r.vectors
        v = getproperty(params, k)
        length(v) == n ||
            throw(ArgumentError("parameter :$k must have length $n (got $(length(v)))."))
    end

    for k in r.matrices
        m = getproperty(params, k)
        spec = parameter_spec(factory, k)
        spec === nothing && throw(
            ArgumentError(
                "Factory $(typeof(factory)) is missing a ParameterSpec for required parameter :$k.",
            ),
        )

        if spec.axes === nothing
            (size(m, 1) == n && size(m, 2) == n) || throw(
                ArgumentError("parameter :$k must have size ($n,$n) (got $(size(m))).")
            )
        else
            row_axis, col_axis = spec.axes
            nr = length(axis_indices(context, row_axis))
            nc = length(axis_indices(context, col_axis))
            (size(m, 1) == nr && size(m, 2) == nc) || throw(
                ArgumentError(
                    "parameter :$k must have size ($nr,$nc) for axes $(spec.axes) (got $(size(m))).",
                ),
            )
        end
    end

    return nothing
end

function _validate_override_keys(
    where_, overrides::NamedTuple, required::Tuple, factory::AbstractBGCFactory
)
    isempty(overrides) && return nothing

    required_set = Set(required)
    for k in keys(overrides)
        k in required_set && continue

        # Optional overrides are permitted only for keys declared in the
        # factory's parameter directory.
        parameter_spec(factory, k) === nothing &&
            throw(ArgumentError("$(where_): unknown parameter key :$k."))
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
        _contains_missing(v) && throw(
            ArgumentError(
                "parameter :$k contains `missing`; all required parameters must be explicitly defined.",
            ),
        )
    end
    return nothing
end

"""
    construct_factory(factory::AbstractBGCFactory; kw...) -> bgc

Construct and compile a concrete biogeochemistry *instance* from a factory and optional parameter overrides.

Key keyword arguments
---------------------
- `grid=nothing`: optional grid used for sinking-velocity fields and for choosing
  the floating point type when interfacing with Oceananigans / OceanBioME.
  Precision is determined by `eltype(grid)`. When `grid` is not provided,
  Agate constructs a `Float64` instance.
- `arch=nothing`: `CPU()` or `GPU()`; when omitted and `grid` is provided, defaults
  to `architecture(grid)`.
- `community`: plankton community structure (size classes, diameters, PFT specs).
- `parameters`: `NamedTuple` of user-supplied parameter overrides.
- `interaction_roles=nothing`: optional role membership for consumer/prey axes. Provide as a `NamedTuple` with fields `consumers` and `prey`, each either `nothing` (all classes), a collection of group `Symbol`s, an index vector, or a boolean mask. When omitted, both roles default to `nothing` (all classes).
- `default_parameter_roles=nothing`: optional role membership used **only** when generating default parameter vectors (e.g. producer- vs consumer-like trait defaults). Provide as `(; producers=..., consumers=...)` using the same formats as `interaction_roles`. When omitted, defaults to matching the interaction axes: `(; producers=interaction_roles.prey, consumers=interaction_roles.consumers)`.
- `interaction_overrides`: optional `NamedTuple` of interaction-parameter overrides (often matrices such as `:palatability_matrix` and `:assimilation_matrix`).
  Values may be concrete objects or provider functions callable as `f(community_context)`.
  For matrix parameters, overrides may be full `(n_total, n_total)` matrices. A group-block `(n_groups, n_groups)` matrix may be supplied and expanded during construction; when the parameter declares role-aware axes, wrap the block matrix as `GroupBlockMatrix(B)` to avoid ambiguity. When axes are declared, rectangular consumer-by-prey matrices sized to those axes (for example `(n_consumer, n_prey)`) are also accepted, as are axis-local group-block matrices.
"""

# ---------------------------------------------------------------------
# Auxiliary field validation
# ---------------------------------------------------------------------

function _validate_auxiliary_fields(auxiliary_fields::Tuple, tracer_names::Tuple)
    isempty(auxiliary_fields) && return nothing

    seen = Set{Symbol}()
    for s in auxiliary_fields
        s isa Symbol || throw(
            ArgumentError("auxiliary_fields entries must be Symbols, got $(typeof(s))")
        )
        (s ∉ seen) || throw(ArgumentError("auxiliary_fields contains duplicate entry :$s"))
        push!(seen, s)
        (s ∉ tracer_names) || throw(
            ArgumentError("auxiliary field :$s conflicts with an existing tracer name")
        )
    end

    return nothing
end

"""
    construct_factory(factory::AbstractBGCFactory; kwargs...) -> bgc

Generic constructor used by model-specific entrypoints like `NiPiZD.construct` and `DARWIN.construct`.

This function assembles a biogeochemistry instance in four conceptual steps:

1. **Parse community structure** into a `CommunityContext` (group ordering, diameter classes, role axes).
2. **Build baseline parameters** by evaluating `parameter_definitions(factory)` into concrete numeric values
   (scalars, vectors, matrices) with element type `FT`.
3. **Apply overrides** (`parameters` and interaction overrides), then compute or recompute any derived interaction
   matrices declared by `derived_matrix_specs(factory)`.
4. **Finalize interactions** into a canonical consumer-by-prey representation and adapt the instance to the
   requested architecture.

Precision is determined by the grid when `grid` is provided: `FT = eltype(grid)` (the grid is the single source of
truth for floating-point type). GPU compatibility is provided by an architecture-wide `Adapt.adapt` pass at the end of
construction; default providers run on the host and should allocate host arrays.

Keyword arguments
-----------------
- `plankton_dynamics`, `biogeochem_dynamics`: dynamics builders (defaults come from the factory).
- `community`: community specification (defaults come from the factory).
- `parameters`: a `NamedTuple` of parameter overrides.
- `interaction_overrides`: optional `NamedTuple` of interaction matrix overrides (model-specific constructors may also
  expose these as convenience keywords).
- `grid`, `arch`: optionally choose precision and architecture; when `grid` is provided, architecture is determined by
  `architecture(grid)`.

The returned object stores the fully **resolved** parameter set at `bgc.parameters`.
"""

function construct_factory(
    factory::AbstractBGCFactory;
    plankton_dynamics=default_plankton_dynamics(factory),
    biogeochem_dynamics=default_biogeochem_dynamics(factory),
    community=default_community(factory),
    parameters::NamedTuple=(;),
    interaction_overrides::Union{Nothing,NamedTuple}=nothing,
    interaction_roles=nothing,
    default_parameter_roles=nothing,
    auxiliary_fields::Tuple=(:PAR,),
    arch=nothing,
    sinking_tracers=nothing,
    grid=nothing,
    open_bottom::Bool=true,
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
            throw(
                ArgumentError(
                    "arch=$arch does not match architecture(grid)=$arch_grid. Architecture is determined by the grid; either omit arch or construct a grid for $arch.",
                ),
            )
        end
    else
        FT = Float64
        isnothing(arch) && (arch = CPU())
    end

    validate_community_inputs(plankton_dynamics, community)
    biogeochem_dynamics isa NamedTuple ||
        throw(ArgumentError("biogeochem_dynamics must be a NamedTuple"))

    # Parse community.
    community_context = parse_community(
        factory,
        FT,
        community;
        plankton_dynamics=plankton_dynamics,
        biogeochem_dynamics=biogeochem_dynamics,
        interaction_roles=interaction_roles,
        default_parameter_roles=default_parameter_roles,
    )

    # ---------------------------------------------------------------------
    # Build tracer expressions and collect parameter requirements.
    # ---------------------------------------------------------------------

    plankton_syms = community_context.plankton_symbols

    # Keep tracer names as a tuple so downstream NamedTuple construction preserves concrete types.
    tracer_names = (keys(biogeochem_dynamics)..., Tuple(plankton_syms)...)

    # Accumulate compiled tracer definitions in a tuple to avoid `Any` erasure.
    tracer_defs = ()
    merged = EquationRequirements()

    for (_, f) in pairs(biogeochem_dynamics)
        tr = f()
        (tr isa CompiledEquation) || throw(
            ArgumentError("biogeochem dynamics $(nameof(f)) must return CompiledEquation"),
        )
        tracer_defs = (tracer_defs..., tr)
        merged = merge_requirements(merged, requirements(tr))
    end

    for idx in eachindex(plankton_syms)
        g = community_context.group_symbols[idx]
        f = getfield(plankton_dynamics, g)
        tr = f(idx)
        (tr isa CompiledEquation) || throw(
            ArgumentError("plankton dynamics $(nameof(f)) must return CompiledEquation")
        )
        tracer_defs = (tracer_defs..., tr)

        merged = merge_requirements(merged, requirements(tr))
    end

    tracers = NamedTuple{tracer_names}(tracer_defs)

    # ---------------------------------------------------------------------
    # Parameters
    # ---------------------------------------------------------------------

    _validate_parameter_directory(factory, merged)

    required = _required_keys(merged)

    interaction_parameter_overrides = normalize_interaction_overrides(
        factory, community_context, interaction_overrides
    )

    _validate_override_keys("parameters", parameters, required, factory)
    _validate_override_keys("interaction_overrides", interaction_parameter_overrides, required, factory)

    parameter_defaults = build_parameter_defaults(factory, community_context, FT)

    # Merge precedence: parameter_defaults < parameters < interaction_overrides
    merged_parameters = merge(parameter_defaults, parameters, interaction_parameter_overrides)

    # Recompute any derived matrices affected by explicit trait overrides.
    explicit_override_keys = (keys(parameters)..., keys(interaction_parameter_overrides)...)
    merged_parameters = resolve_derived_matrices(
        factory, community_context, merged_parameters, explicit_override_keys
    )

    # Ensure the required keys are present.
    missing = Symbol[]
    for k in required
        hasproperty(merged_parameters, k) || push!(missing, k)
    end
    isempty(missing) ||
        throw(ArgumentError("missing required parameters: $(join(string.(missing), ", "))"))

    # Finalize any role-aware interaction matrices.
    merged_parameters = finalize_interaction_parameters(
        factory, community_context, merged_parameters
    )

    # Slice down to the required keys (plus selected internal helpers) for better type stability.
    internal = hasproperty(merged_parameters, :interactions) ? (:interactions,) : ()
    all_keys = (required..., internal...)
    resolved_parameters = NamedTuple{all_keys}(
        Tuple(getproperty(merged_parameters, k) for k in all_keys)
    )

    _reject_missing_values(resolved_parameters)

    _validate_parameter_shapes(factory, community_context, resolved_parameters, merged)

    # ---------------------------------------------------------------------
    # Compile + instantiate.
    # ---------------------------------------------------------------------

    # Auxiliary fields are appended to the Oceananigans kernel argument list.
    # The symbol tuple lives only on the host side; the compiled model stores
    # a GPU-safe integer `TracerIndex`.
    #
    # NOTE: `auxiliary_fields` is a host-side tuple of `Symbol`s only. The compiled
    # model stores a GPU-safe integer `TracerIndex`, so no runtime Symbol work
    # occurs inside kernels.
    _validate_auxiliary_fields(auxiliary_fields, tracer_names)
    tracer_index = build_tracer_index(
        community_context,
        tracer_names,
        auxiliary_fields;
        n_biogeochem_tracers=length(keys(biogeochem_dynamics)),
    )

    if isnothing(sinking_tracers)
        bgc_factory = define_tracer_functions(
            resolved_parameters,
            tracers;
            auxiliary_fields=auxiliary_fields,
            tracer_index=tracer_index,
        )
        bgc = bgc_factory(resolved_parameters)
    else
        sinking_velocities = setup_velocity_fields(sinking_tracers, grid, open_bottom)
        bgc_factory = define_tracer_functions(
            resolved_parameters,
            tracers;
            auxiliary_fields=auxiliary_fields,
            tracer_index=tracer_index,
            sinking_velocities=sinking_velocities,
        )
        bgc = bgc_factory(resolved_parameters, sinking_velocities)
    end

    # Move any arrays inside `bgc` onto the requested architecture.
    bgc = _on_architecture(arch, bgc)

    return bgc
end