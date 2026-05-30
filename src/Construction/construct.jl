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

using ..Library.Allometry: AllometricParam, resolve_diameter_indexed_vector

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
    factory::AbstractBGCFactory, community_context, ::Type{T}
) where {T<:Real}
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

    indices = getproperty(community_context, provider.indices_field)
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

function parameter_definition(factory::AbstractBGCFactory, key::Symbol)
    for def in parameter_definitions(factory)
        def.spec.name === key && return def
    end
    return nothing
end

function materialize_allometric_parameter_override(
    factory::AbstractBGCFactory, context, key::Symbol, value::AllometricParam, ::Type{T}
) where {T<:Real}
    def = parameter_definition(factory, key)
    provider = def === nothing ? nothing : def.default

    provider isa DiameterIndexedVectorDefault || throw(
        ArgumentError(
            "parameter :$key only supports AllometricParam overrides for diameter-indexed vector parameters."
        ),
    )

    indices = getproperty(context, provider.indices_field)
    default = T(provider.default)
    return resolve_diameter_indexed_vector(
        T, context.diameters, indices, value; default=default
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

        if value isa AllometricParam
            push!(
                entries,
                key => materialize_allometric_parameter_override(
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


function manifest_axis(axis)
    axis === nothing && return nothing
    axis isa Symbol && return string(axis)
    axis isa Tuple && return Any[string(a) for a in axis]
    return string(axis)
end

function manifest_parameter_value(x, name=nothing)
    if x === nothing || x isa Bool || x isa Integer || x isa AbstractString
        return x
    elseif x isa AbstractFloat
        return isfinite(x) ? x : string(x)
    elseif x isa Symbol
        return string(x)
    elseif x isa NamedTuple
        return Dict{String,Any}(
            string(k) => manifest_parameter_value(v, k) for (k, v) in pairs(x)
        )
    elseif x isa AbstractDict
        return Dict{String,Any}(
            string(k) => manifest_parameter_value(v, k) for (k, v) in pairs(x)
        )
    elseif x isa AbstractMatrix
        return Any[
            Any[manifest_parameter_value(x[i, j], name) for j in axes(x, 2)]
            for i in axes(x, 1)
        ]
    elseif x isa AbstractVector || x isa Tuple
        return Any[manifest_parameter_value(v, name) for v in x]
    else
        label = isnothing(name) ? "" : " $(repr(name))"
        throw(
            ArgumentError(
                "Cannot serialize manifest parameter$(label) of type $(typeof(x))."
            )
        )
    end
end

function parameter_record_manifest(spec, value)
    return Dict{String,Any}(
        "shape" => string(spec.shape),
        "axes" => manifest_axis(spec.axes),
        "doc" => spec.doc,
        "value" => manifest_parameter_value(value),
    )
end

function parameter_manifest(factory::AbstractBGCFactory, params::NamedTuple, required::Tuple)
    records = Dict{String,Any}()
    for key in required
        spec = parameter_spec(factory, key)
        spec === nothing && throw(
            ArgumentError(
                "Factory $(typeof(factory)) is missing a ParameterSpec for parameter :$key."
            ),
        )
        records[string(key)] = parameter_record_manifest(spec, getproperty(params, key))
    end
    return records
end

function parameter_value_manifest(params::NamedTuple, required::Tuple)
    return Dict{String,Any}(
        string(key) => manifest_parameter_value(getproperty(params, key)) for key in required
    )
end

function plankton_diameter_group_manifest(context)
    return Dict{String,Any}(
        string(group) => Any[manifest_parameter_value(context.diameters[i]) for i in indices]
        for (group, indices) in context.group_indices
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
- `interaction_roles`: optional `NamedTuple` with `consumers` and `prey`
  membership for interaction axes.
- `default_parameter_roles`: optional `NamedTuple` with `producers` and
  `consumers` membership used only when generating default parameter vectors.
- `auxiliary_fields`: auxiliary values appended to the tracer argument list.
- `grid`, `arch`: optional grid and architecture inputs.
- `scalar_type`: explicit runtime scalar type; when omitted, construction uses `eltype(grid)` or `Float64` if no grid is supplied.
- `sinking_tracers`, `open_bottom`: sinking-velocity configuration.

The returned object stores the fully resolved parameter set in
`bgc.parameters`.
"""

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

function construct_factory(factory::AbstractBGCFactory; kwargs...)
    bgc, _ = _construct_factory_with_context(factory; build_context=false, kwargs...)
    return bgc
end

function construct_factory_with_context(factory::AbstractBGCFactory; kwargs...)
    return _construct_factory_with_context(factory; build_context=true, kwargs...)
end

function _construct_factory_with_context(
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
    scalar_type=nothing,
    open_bottom::Bool=true,
    build_context::Bool=true,
)
    if isnothing(grid) && !isnothing(sinking_tracers)
        grid = BoxModelGrid()
    end
    T = resolve_construction_scalar_type(grid, scalar_type)

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
        factory,
        T,
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
    bgc = on_architecture(arch, bgc)

    construction_context = build_context ? (
        tracers=Any[string(name) for name in tracer_names],
        auxiliary_fields=Any[string(name) for name in auxiliary_fields],
        parameters=parameter_manifest(factory, resolved_parameters, required),
        parameter_values=parameter_value_manifest(resolved_parameters, required),
        plankton_diameters=Any[manifest_parameter_value(d) for d in plankton_diameter_metadata],
        plankton_diameters_by_group=plankton_diameter_group_manifest(community_context),
        scalar_type=string(T),
        architecture=string(typeof(arch)),
        has_sinking_velocities=!isnothing(sinking_tracers),
        sinking=Dict{String,Any}(
            "enabled" => !isnothing(sinking_tracers),
            "tracers" => manifest_parameter_value(sinking_tracers),
            "open_bottom" => open_bottom,
        ),
    ) : nothing

    return bgc, construction_context
end
