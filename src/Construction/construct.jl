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
    DiameterIndexedMaterialization,
    parameter_directory,
    parameter_spec,
    default_plankton_dynamics,
    default_biogeochem_dynamics,
    default_community

using ..Configuration:
    axis_indices,
    normalize_interaction_overrides,
    matrix_definitions,
    resolve_derived_matrices,
    finalize_interaction_parameters,
    parse_community,
    parameter_role_indices,
    validate_community_inputs

using ..Runtime: build_tracer_index

using ..Equations: CompiledEquation

using ..Library.Allometry: AbstractParamDef, resolve_diameter_indexed_vector

"""Evaluate `parameter_definitions(factory)` to produce baseline parameter defaults.

Defaults are evaluated on the host during model construction. The returned
`NamedTuple` maps parameter keys to concrete values (scalars, vectors, matrices).

Parameters declared with `NoDefault()` are omitted; they are expected to be
provided by user overrides or derived-matrix providers later in construction.
"""
function build_parameter_defaults(
    factory::AbstractBGCFactory, community_context, ::Type{T}
) where {T<:Real}
    defs = parameter_definitions(factory)
    isempty(defs) && return (;)

    pairs = Pair{Symbol,Any}[]
    for def in defs
        spec = def.spec
        provider = def.default
        provider isa NoDefault && continue
        value = evaluate_default(provider, spec, factory, community_context, T)
        push!(pairs, spec.name => value)
    end

    return (; pairs...)
end

@inline function evaluate_default(
    provider::ConstDefault, spec, ::AbstractBGCFactory, ::Any, ::Type{T}
) where {T<:Real}
    spec.shape === :scalar || throw(
        ArgumentError(
            "ConstDefault can only be used for scalar parameters (:$((spec.name)))."
        ),
    )
    v = provider.value
    return v isa Bool ? v : T(v)
end

@inline function evaluate_default(
    provider::FillDefault, spec, ::AbstractBGCFactory, community_context, ::Type{T}
) where {T<:Real}
    spec.shape in (:vector, :matrix) || throw(
        ArgumentError(
            "FillDefault can only be used for vector or matrix parameters (:$((spec.name))).",
        ),
    )

    v = provider.value
    v = v isa Bool ? v : T(v)

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

@inline function evaluate_default(
    provider::DiameterIndexedVectorDefault,
    spec,
    ::AbstractBGCFactory,
    community_context,
    ::Type{T},
) where {T<:Real}
    spec.shape === :vector || throw(
        ArgumentError(
            "DiameterIndexedVectorDefault can only be used for vector parameters (:$((spec.name))).",
        ),
    )

    indices = parameter_role_indices(community_context, provider.role)
    default = T(provider.default)
    return resolve_diameter_indexed_vector(
        T, community_context.diameters, indices, provider.value; default=default
    )
end

"""Move `x` to the requested Oceananigans architecture."""
function on_architecture(arch, x)
    arch === nothing && return x
    arch isa CPU && return x

    return Base.invokelatest(Adapt.adapt, architecture_array_type(arch), x)
end

"""Return the preferred array storage type for `arch`."""
function architecture_array_type(arch)
    arch isa CPU && return Array
    arch isa GPU && return Oceananigans.Architectures.array_type(arch)
    return Array
end

const RESERVED_PARAMETER_KEYS = (:x, :y, :z, :t)

function validate_parameter_directory(factory::AbstractBGCFactory)
    defs = parameter_definitions(factory)
    isempty(defs) && return ()

    keys_ = map(d -> d.spec.name, defs)
    length(unique(keys_)) == length(keys_) || throw(
        ArgumentError(
            "parameter_definitions(::$(typeof(factory))) contains duplicate keys."
        ),
    )

    for k in keys_
        (k in RESERVED_PARAMETER_KEYS) && throw(
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

        if spec.shape === :vector
            if spec.axes !== nothing
                spec.axes isa Symbol || throw(
                    ArgumentError(
                        "parameter :$(spec.name) vector axis must be a Symbol (got $(typeof(spec.axes))).",
                    ),
                )
            end
        elseif spec.shape === :matrix
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
                    "parameter :$(spec.name) has axes=$(spec.axes) but is not vector or matrix."
                ),
            )
        end

        if spec.materialization isa DiameterIndexedMaterialization
            spec.shape === :vector || throw(
                ArgumentError(
                    "parameter :$(spec.name) declares diameter-indexed materialization but is not vector-shaped."
                ),
            )
        end
    end

    return Tuple(keys_)
end

function validate_parameter_shapes(
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


function parameter_axis_names(context, axis::Symbol, parameter_name::Symbol)
    axis === :plankton && return context.plankton_symbols
    throw(ArgumentError("parameter :$parameter_name has unsupported vector axis :$axis."))
end

function expand_named_vector_override(
    spec, default_value, user_value::NamedTuple, context, ::Type{T}
) where {T<:Real}
    spec.axes === nothing && throw(
        ArgumentError(
            "parameter :$(spec.name) does not support NamedTuple overrides because it has no named vector axis."
        ),
    )

    names = parameter_axis_names(context, spec.axes, spec.name)
    length(default_value) == length(names) || throw(
        ArgumentError(
            "parameter :$(spec.name) default length $(length(default_value)) does not match axis :$(spec.axes) length $(length(names))."
        ),
    )

    expanded = copy(default_value)
    for (key, value) in pairs(user_value)
        idx = findfirst(==(key), names)
        if idx === nothing
            expected = join(string.(names), ", ")
            throw(
                ArgumentError(
                    "Unknown key `$(key)` for parameter `$(spec.name)`. Expected one of: $(expected)."
                ),
            )
        end
        expanded[idx] = value isa Bool ? value : T(value)
    end
    return expanded
end

function materialize_parameter_law_override(
    factory::AbstractBGCFactory, context, key::Symbol, value::AbstractParamDef, ::Type{T}
) where {T<:Real}
    spec = parameter_spec(factory, key)
    materialization = spec === nothing ? nothing : spec.materialization

    materialization isa DiameterIndexedMaterialization || throw(
        ArgumentError(
            "parameter :$key only supports parameter-law overrides for parameters with declared diameter-indexed vector materialization."
        ),
    )

    indices = if isnothing(materialization.role)
        eachindex(context.diameters)
    else
        parameter_role_indices(context, materialization.role)
    end
    fill_value = T(materialization.fill_value)
    return resolve_diameter_indexed_vector(
        T, context.diameters, indices, value; default=fill_value
    )
end

function materialize_parameter_value(spec, value, ::Type{T}) where {T<:Real}
    if spec.shape === :scalar
        return value isa Bool ? value : T(value)
    elseif spec.shape === :vector
        value isa AbstractVector || return value
        eltype(value) === T && return value
        out = similar(value, T, axes(value))
        copyto!(out, value)
        return out
    elseif spec.shape === :matrix
        value isa AbstractMatrix || return value
        eltype(value) === T && return value
        out = similar(value, T, axes(value))
        copyto!(out, value)
        return out
    end

    return value
end

function materialize_parameter_overrides(
    factory::AbstractBGCFactory,
    context,
    defaults::NamedTuple,
    overrides::NamedTuple,
    ::Type{T},
) where {T<:Real}
    isempty(overrides) && return overrides

    entries = Pair{Symbol,Any}[]
    for (key, value) in Base.pairs(overrides)
        spec = parameter_spec(factory, key)
        spec === nothing && begin
            push!(entries, key => value)
            continue
        end

        if value isa AbstractParamDef
            push!(
                entries,
                key => materialize_parameter_law_override(
                    factory, context, key, value, T
                ),
            )
        elseif value isa NamedTuple
            spec.shape === :vector || throw(
                ArgumentError(
                    "parameter :$key does not support NamedTuple overrides because it is $(spec.shape)-shaped."
                ),
            )
            hasproperty(defaults, key) || throw(
                ArgumentError(
                    "parameter :$key does not support partial NamedTuple overrides because it has no direct default value."
                ),
            )
            push!(
                entries,
                key => expand_named_vector_override(
                    spec, getproperty(defaults, key), value, context, T
                ),
            )
        else
            push!(entries, key => materialize_parameter_value(spec, value, T))
        end
    end

    return (; entries...)
end

function validate_override_keys(
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

@inline contains_missing(x) = x === missing

function contains_missing(x::NamedTuple)
    for v in values(x)
        contains_missing(v) && return true
    end
    return false
end

function contains_missing(x::AbstractArray)
    return any(ismissing, x)
end

function reject_missing_values(params::NamedTuple)
    for (k, v) in pairs(params)
        contains_missing(v) && throw(
            ArgumentError(
                "parameter :$k contains `missing`; all required parameters must be explicitly defined.",
            ),
        )
    end
    return nothing
end

function validate_auxiliary_fields(auxiliary_fields::Tuple, tracer_names::Tuple)
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


"""Resolve the scalar type from an explicit choice, the grid, or `Float64`."""
function resolve_construction_scalar_type(grid, scalar_type)
    if scalar_type !== nothing
        scalar_type isa Type || throw(
            ArgumentError(
                "scalar_type must be a concrete subtype of Real; got $(scalar_type)"
            ),
        )
        scalar_type <: Real || throw(
            ArgumentError(
                "scalar_type must be a concrete subtype of Real; got $(scalar_type)"
            ),
        )
        isconcretetype(scalar_type) ||
            throw(ArgumentError("scalar_type must be concrete; got $(scalar_type)"))
        return scalar_type
    end

    grid !== nothing && return eltype(grid)
    return Float64
end

convert_sinking_tracers(::Type{T}, ::Nothing) where {T<:Real} = nothing
function convert_sinking_tracers(::Type{T}, sinking_tracers::NamedTuple) where {T<:Real}
    return NamedTuple{keys(sinking_tracers)}(
        Tuple(convert(T, velocity) for velocity in values(sinking_tracers))
    )
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
- `ecological_roles`: optional model-defined ecological group identities retained in manifest state.
- `interaction_roles`: optional `NamedTuple` with `consumers` and `prey`
  membership for interaction axes.
- `parameter_roles`: optional `NamedTuple` of named parameter-applicability roles.
- `auxiliary_fields`: auxiliary values appended to the tracer argument list.
- `grid`, `arch`: optional grid and architecture inputs.
- `scalar_type`: explicit runtime scalar type; when provided it takes precedence over
  `eltype(grid)`. When omitted, construction uses `eltype(grid)` or `Float64` if no
  grid is supplied.
- `sinking_tracers`, `open_bottom`: sinking-velocity configuration.

The returned object stores the fully resolved parameter set in
`bgc.parameters`.
"""
function construct_factory(factory::AbstractBGCFactory; kwargs...)
    bgc, _ = _construct_factory(factory; build_manifest=false, kwargs...)
    return bgc
end

function _recipe_realization_inputs(factory::AbstractBGCFactory, recipe::ModelRecipe)
    realization = deepcopy(recipe)
    return (;
        plankton_dynamics=default_plankton_dynamics(
            factory, realization.community, realization.ecological_roles
        ),
        biogeochem_dynamics=default_biogeochem_dynamics(factory),
        community=realization.community,
        parameters=realization.parameter_overrides,
        interaction_overrides=realization.interaction_overrides,
        ecological_roles=realization.ecological_roles,
        interaction_roles=realization.interaction_roles,
        parameter_roles=realization.parameter_roles,
        auxiliary_fields=realization.auxiliary_fields,
        sinking_tracers=realization.sinking_tracers,
        open_bottom=realization.open_bottom,
    )
end

"""Realize a canonical `ModelRecipe` in the supplied execution environment."""
function construct_factory(
    recipe::ModelRecipe; grid=nothing, arch=nothing, scalar_type=nothing
)
    factory = replay_factory(recipe)
    inputs = _recipe_realization_inputs(factory, recipe)
    return construct_factory(factory; inputs..., grid, arch, scalar_type)
end

"""Construct a factory-defined model and return it with its resolved `ModelManifest`."""
function construct_factory_plus_manifest(factory::AbstractBGCFactory; kwargs...)
    return _construct_factory(factory; build_manifest=true, kwargs...)
end

"""Realize a canonical `ModelRecipe` and return its resolved `ModelManifest`."""
function construct_factory_plus_manifest(
    recipe::ModelRecipe; grid=nothing, arch=nothing, scalar_type=nothing
)
    factory = replay_factory(recipe)
    inputs = _recipe_realization_inputs(factory, recipe)
    return construct_factory_plus_manifest(factory; inputs..., grid, arch, scalar_type)
end

function _construct_factory(
    factory::AbstractBGCFactory;
    plankton_dynamics=default_plankton_dynamics(factory),
    biogeochem_dynamics=default_biogeochem_dynamics(factory),
    community=default_community(factory),
    parameters::NamedTuple=(;),
    interaction_overrides::NamedTuple=(;),
    ecological_roles::NamedTuple=(;),
    interaction_roles=nothing,
    parameter_roles=nothing,
    auxiliary_fields::Tuple=(:PAR,),
    arch=nothing,
    sinking_tracers=nothing,
    grid=nothing,
    scalar_type=nothing,
    open_bottom::Bool=true,
    build_manifest::Bool=false,
)
    if isnothing(grid) && !isnothing(sinking_tracers)
        grid = BoxModelGrid()
    end
    T = resolve_construction_scalar_type(grid, scalar_type)
    sinking_tracers = convert_sinking_tracers(T, sinking_tracers)

    if !isnothing(grid)
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
        isnothing(arch) && (arch = CPU())
    end

    validate_community_inputs(plankton_dynamics, community)
    biogeochem_dynamics isa NamedTuple ||
        throw(ArgumentError("biogeochem_dynamics must be a NamedTuple"))
    community_context = parse_community(
        T,
        community;
        biogeochem_tracers=keys(biogeochem_dynamics),
        interaction_roles=interaction_roles,
        parameter_roles=parameter_roles,
    )

    plankton_syms = community_context.plankton_symbols
    tracer_names = (keys(biogeochem_dynamics)..., Tuple(plankton_syms)...)
    tracer_defs = ()

    for (k, f) in pairs(biogeochem_dynamics)
        tr = f()
        (tr isa CompiledEquation) || throw(
            ArgumentError("biogeochem dynamics :$k must return CompiledEquation"),
        )
        tracer_defs = (tracer_defs..., tr)
    end

    for idx in eachindex(plankton_syms)
        g = community_context.group_symbols[idx]
        f = getfield(plankton_dynamics, g)
        tr = f(idx)
        (tr isa CompiledEquation) || throw(
            ArgumentError("plankton dynamics :$g must return CompiledEquation")
        )
        tracer_defs = (tracer_defs..., tr)
    end

    tracers = NamedTuple{tracer_names}(tracer_defs)

    required = validate_parameter_directory(factory)

    interaction_parameter_overrides = normalize_interaction_overrides(
        factory, community_context, interaction_overrides
    )

    validate_override_keys("parameters", parameters, required, factory)
    validate_override_keys(
        "interaction_overrides", interaction_parameter_overrides, required, factory
    )

    parameter_defaults = build_parameter_defaults(factory, community_context, T)
    parameter_overrides = materialize_parameter_overrides(
        factory, community_context, parameter_defaults, parameters, T
    )

    merged_parameters = merge(
        parameter_defaults, parameter_overrides, interaction_parameter_overrides
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

    reject_missing_values(resolved_parameters)

    validate_parameter_shapes(factory, community_context, resolved_parameters, required)
    validate_auxiliary_fields(auxiliary_fields, tracer_names)
    tracer_index = build_tracer_index(
        community_context,
        tracer_names,
        auxiliary_fields;
        n_biogeochem_tracers=length(keys(biogeochem_dynamics)),
    )

    plankton_diameter_metadata = Tuple(community_context.diameters)

    if isnothing(sinking_tracers)
        bgc_factory = define_tracer_functions(
            resolved_parameters,
            tracers;
            auxiliary_fields=auxiliary_fields,
            tracer_index=tracer_index,
        )
        bgc = bgc_factory(
            resolved_parameters; plankton_diameters=plankton_diameter_metadata
        )
    else
        sinking_velocities = setup_velocity_fields(sinking_tracers, grid, open_bottom)
        bgc_factory = define_tracer_functions(
            resolved_parameters,
            tracers;
            auxiliary_fields=auxiliary_fields,
            tracer_index=tracer_index,
            sinking_velocities=sinking_velocities,
        )
        bgc = bgc_factory(
            resolved_parameters,
            sinking_velocities;
            plankton_diameters=plankton_diameter_metadata,
        )
    end
    manifest = if build_manifest
        capture_model_manifest(
            factory,
            resolved_parameters,
            community_context;
            tracer_order=tracer_names,
            auxiliary_fields,
            ecological_roles,
            explicit_override_keys,
            sinking_tracers,
            open_bottom,
            scalar_type=T,
        )
    else
        nothing
    end

    bgc = on_architecture(arch, bgc)

    return bgc, manifest
end
