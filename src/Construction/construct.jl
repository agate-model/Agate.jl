using OceanBioME: BoxModelGrid, setup_velocity_fields

using Adapt: adapt

import Oceananigans

using Oceananigans.Architectures: architecture, CPU, GPU

using ..Factories:
    AbstractBGCFactory,
    parameter_definitions,
    ConstDefault,
    NoDefault,
    FillDefault,
    DiameterIndexedVectorDefault,
    parameter_spec,
    default_plankton_dynamics,
    default_biogeochem_dynamics,
    default_community

using ..Configuration:
    axis_indices,
    normalize_interaction_overrides,
    resolve_derived_matrices,
    finalize_interaction_parameters,
    parse_community,
    validate_community_inputs

using ..Runtime: build_tracer_index

using ..Equations: CompiledEquation

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
            "parameter_definitions(::$(typeof(factory))) contains duplicate keys."
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
        ArgumentError(
            "ConstDefault can only be used for scalar parameters (:$((spec.name)))."
        ),
    )
    v = provider.value
    return v isa Bool ? v : FT(v)
end

@inline function _evaluate_default(
    provider::FillDefault, spec, ::AbstractBGCFactory, community_context, ::Type{FT}
) where {FT}
    spec.shape in (:vector, :matrix) || throw(
        ArgumentError(
            "FillDefault can only be used for vector or matrix parameters (:$((spec.name))).",
        ),
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
        FT, community_context.diameters, indices, provider.value; default=default
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

const _RESERVED_PARAMETER_KEYS = (:x, :y, :z, :t)

function _validate_parameter_directory(factory::AbstractBGCFactory)
    defs = parameter_definitions(factory)
    isempty(defs) && return ()

    keys_ = map(d -> d.spec.name, defs)
    length(unique(keys_)) == length(keys_) || throw(
        ArgumentError(
            "parameter_definitions(::$(typeof(factory))) contains duplicate keys."
        ),
    )

    for k in keys_
        (k in _RESERVED_PARAMETER_KEYS) && throw(
            ArgumentError(
                "parameter_definitions(::$(typeof(factory))) declares reserved parameter key :$k.",
            ),
        )
    end

    for def in defs
        spec = def.spec
        spec.shape in (:scalar, :vector, :matrix) || throw(
            ArgumentError(
                "parameter_definitions(::$(typeof(factory))) declares :$(spec.name) with invalid shape $(spec.shape).",
            ),
        )

        if spec.shape === :matrix
            if spec.axes !== nothing
                (spec.axes isa Tuple && length(spec.axes) == 2) || throw(
                    ArgumentError(
                        "parameter :$(spec.name) axes must be a 2-tuple of Symbols (got $(typeof(spec.axes))).",
                    ),
                )
                row_axis, col_axis = spec.axes
                (row_axis isa Symbol && col_axis isa Symbol) || throw(
                    ArgumentError(
                        "parameter :$(spec.name) axes must be Symbols (got $(spec.axes)).",
                    ),
                )
            end
        else
            spec.axes === nothing || throw(
                ArgumentError(
                    "parameter :$(spec.name) has axes=$(spec.axes) but is not a matrix."
                ),
            )
        end
    end

    return Tuple(keys_)
end

function _validate_parameter_shapes(
    factory::AbstractBGCFactory, context, params::NamedTuple, required::Tuple
)
    n = context.n_total

    for k in required
        spec = parameter_spec(factory, k)
        spec === nothing && throw(
            ArgumentError(
                "Factory $(typeof(factory)) is missing a ParameterSpec for parameter :$k.",
            ),
        )

        if spec.shape === :vector
            v = getproperty(params, k)
            length(v) == n || throw(
                ArgumentError("parameter :$k must have length $n (got $(length(v))).")
            )
        elseif spec.shape === :matrix
            m = getproperty(params, k)

            if spec.axes === nothing
                (size(m, 1) == n && size(m, 2) == n) || throw(
                    ArgumentError("parameter :$k must have size ($n,$n) (got $(size(m)))."),
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

Construct a concrete biogeochemistry model from `factory` and optional
configuration overrides.

Construction proceeds in four stages:

1. Parse the community into a `CommunityContext`.
2. Evaluate `parameter_definitions(factory)` into concrete defaults.
3. Apply user overrides and resolve any derived interaction matrices.
4. Finalize interaction parameters, compile tracer functions, and adapt the
   result to the requested architecture.

Keyword arguments
-----------------
- `plankton_dynamics`, `biogeochem_dynamics`: dynamics builders.
- `community`: plankton community specification.
- `parameters`: `NamedTuple` of parameter overrides.
- `interaction_overrides`: `NamedTuple` of explicit interaction-matrix
  overrides.
- `interaction_roles`: optional `NamedTuple` with `consumers` and `prey`
  membership for interaction axes.
- `default_parameter_roles`: optional `NamedTuple` with `producers` and
  `consumers` membership used only when generating default parameter vectors.
- `auxiliary_fields`: auxiliary values appended to the tracer argument list.
- `grid`, `arch`: optional precision and architecture inputs.
- `sinking_tracers`, `open_bottom`: sinking-velocity configuration.

The returned object stores the fully resolved parameter set in
`bgc.parameters`.
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
    if isnothing(grid) && !isnothing(sinking_tracers)
        grid = BoxModelGrid()
    end
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
    community_context = parse_community(
        factory,
        FT,
        community;
        plankton_dynamics=plankton_dynamics,
        biogeochem_dynamics=biogeochem_dynamics,
        interaction_roles=interaction_roles,
        default_parameter_roles=default_parameter_roles,
    )

    plankton_syms = community_context.plankton_symbols
    tracer_names = (keys(biogeochem_dynamics)..., Tuple(plankton_syms)...)
    tracer_defs = ()

    for (_, f) in pairs(biogeochem_dynamics)
        tr = f()
        (tr isa CompiledEquation) || throw(
            ArgumentError("biogeochem dynamics $(nameof(f)) must return CompiledEquation"),
        )
        tracer_defs = (tracer_defs..., tr)
    end

    for idx in eachindex(plankton_syms)
        g = community_context.group_symbols[idx]
        f = getfield(plankton_dynamics, g)
        tr = f(idx)
        (tr isa CompiledEquation) || throw(
            ArgumentError("plankton dynamics $(nameof(f)) must return CompiledEquation")
        )
        tracer_defs = (tracer_defs..., tr)
    end

    tracers = NamedTuple{tracer_names}(tracer_defs)

    required = _validate_parameter_directory(factory)

    interaction_parameter_overrides = normalize_interaction_overrides(
        factory, community_context, interaction_overrides
    )

    _validate_override_keys("parameters", parameters, required, factory)
    _validate_override_keys(
        "interaction_overrides", interaction_parameter_overrides, required, factory
    )

    parameter_defaults = build_parameter_defaults(factory, community_context, FT)
    merged_parameters = merge(
        parameter_defaults, parameters, interaction_parameter_overrides
    )
    explicit_override_keys = (keys(parameters)..., keys(interaction_parameter_overrides)...)
    merged_parameters = resolve_derived_matrices(
        factory, community_context, merged_parameters, explicit_override_keys
    )
    missing = Symbol[]
    for k in required
        hasproperty(merged_parameters, k) || push!(missing, k)
    end
    isempty(missing) ||
        throw(ArgumentError("missing required parameters: $(join(string.(missing), ", "))"))
    merged_parameters = finalize_interaction_parameters(
        factory, community_context, merged_parameters
    )
    internal = hasproperty(merged_parameters, :interactions) ? (:interactions,) : ()
    all_keys = (required..., internal...)
    resolved_parameters = NamedTuple{all_keys}(
        Tuple(getproperty(merged_parameters, k) for k in all_keys)
    )

    _reject_missing_values(resolved_parameters)

    _validate_parameter_shapes(factory, community_context, resolved_parameters, required)
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
    bgc = _on_architecture(arch, bgc)

    return bgc
end
