using OceanBioME: BoxModelGrid, setup_velocity_fields

using Adapt: adapt

import Oceananigans

using Oceananigans.Architectures: architecture, CPU, GPU

using ..Factories:
    AbstractBGCFactory,
    parameter_definitions,
    ConstDefault,
    DerivedDefault,
    NoDefault,
    FillDefault,
    DiameterIndexedVectorDefault,
    DiameterIndexedMaterialization,
    derive_default,
    parameter_directory,
    parameter_spec,
    default_plankton_dynamics,
    default_biogeochem_dynamics,
    default_community

using ..Configuration:
    axis_indices,
    normalize_interaction_overrides,
    finalize_interaction_parameters,
    parse_community,
    parameter_role_indices,
    validate_community_inputs,
    Pool,
    realize_components,
    component_tracers,
    realize_component_groups

using ..Runtime: build_tracer_index

using ..Processes:
    ModelDefinition, normalize_model, driver_identities, participants, process_kind,
    resolve_parameter_applicability

using ..Compilation: compile_model_tendencies

using ..Equations: CompiledEquation

using ..Library.Allometry: AbstractParamDef, resolve_diameter_indexed_vector

"""Evaluate `parameter_definitions(factory)` to produce baseline parameter defaults.

Defaults are evaluated on the host during model construction. The returned
`NamedTuple` maps parameter keys to concrete values (scalars, vectors, matrices).

Parameters declared with `NoDefault()` or `DerivedDefault(...)` are omitted from
this first pass. Derived defaults are resolved after user overrides are materialized.
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
        (provider isa NoDefault || provider isa DerivedDefault) && continue
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

    key_set = Set(keys_)
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

        provider = def.default
        if provider isa DerivedDefault
            for dep in provider.deps
                dep in key_set || throw(
                    ArgumentError(
                        "derived default :$(spec.name) depends on undeclared parameter :$dep.",
                    ),
                )
            end
        end
    end

    return Tuple(keys_)
end

function derived_default_order(factory::AbstractBGCFactory)
    derived = Any[
        definition for definition in parameter_definitions(factory) if
        definition.default isa DerivedDefault
    ]
    isempty(derived) && return ()

    derived_keys = Tuple(definition.spec.name for definition in derived)
    resolved_keys = Set{Symbol}()
    ordered = Any[]
    pending = derived

    while !isempty(pending)
        remaining = Any[]
        progressed = false

        for definition in pending
            dependencies = Tuple(
                dep for dep in definition.default.deps if dep in derived_keys
            )
            if all(dep -> dep in resolved_keys, dependencies)
                push!(ordered, definition)
                push!(resolved_keys, definition.spec.name)
                progressed = true
            else
                push!(remaining, definition)
            end
        end

        progressed || throw(
            ArgumentError(
                "derived parameter defaults contain a dependency cycle among: " *
                join(string.(Tuple(definition.spec.name for definition in remaining)), ", "),
            ),
        )
        pending = remaining
    end

    return Tuple(ordered)
end

function validate_derived_parameter_result(spec, context, value)
    T = context.scalar_type

    if spec.shape === :scalar
        value isa Bool && return nothing
        value isa Number || throw(
            ArgumentError(
                "derived default :$(spec.name) must be scalar; got $(typeof(value)).",
            ),
        )
        typeof(value) === T || throw(
            ArgumentError(
                "derived default :$(spec.name) must have type $(T); got $(typeof(value)). No implicit casting is performed.",
            ),
        )
        return nothing
    elseif spec.shape === :vector
        value isa AbstractVector || throw(
            ArgumentError(
                "derived default :$(spec.name) must be a vector; got $(typeof(value)).",
            ),
        )
        eltype(value) === T || throw(
            ArgumentError(
                "derived default :$(spec.name) must have eltype $(T); got eltype $(eltype(value)). No implicit casting is performed.",
            ),
        )
        length(value) == context.n_total || throw(
            ArgumentError(
                "derived default :$(spec.name) must have length $(context.n_total); got $(length(value)).",
            ),
        )
        return nothing
    end

    value isa AbstractMatrix || throw(
        ArgumentError(
            "derived default :$(spec.name) must be a matrix; got $(typeof(value)).",
        ),
    )
    eltype(value) === T || throw(
        ArgumentError(
            "derived default :$(spec.name) must have eltype $(T); got eltype $(eltype(value)). No implicit casting is performed.",
        ),
    )

    if spec.axes === nothing
        expected = (context.n_total, context.n_total)
    else
        row_axis, col_axis = spec.axes
        expected = (length(axis_indices(context, row_axis)), length(axis_indices(context, col_axis)))
    end
    size(value) == expected || throw(
        ArgumentError(
            "derived default :$(spec.name) must have size $expected; got $(size(value)).",
        ),
    )
    return nothing
end

"""Resolve all `DerivedDefault` parameter definitions in deterministic dependency order.

An explicit override of a derived parameter wins. Otherwise a derived value is computed
when missing and recomputed whenever a dependency was explicitly overridden or was itself
recomputed. Dependencies between derived defaults are topologically ordered and cycles are
rejected during setup.
"""
function resolve_parameter_defaults(
    factory::AbstractBGCFactory,
    context,
    params::NamedTuple,
    explicit_override_keys::Tuple{Vararg{Symbol}},
)
    ordered = derived_default_order(factory)
    isempty(ordered) && return params

    override_set = Set{Symbol}(explicit_override_keys)
    changed = Set{Symbol}(explicit_override_keys)
    resolved = params

    for definition in ordered
        key = definition.spec.name
        provider = definition.default

        if key in override_set && hasproperty(resolved, key)
            continue
        end

        missing_deps = Tuple(dep for dep in provider.deps if !hasproperty(resolved, dep))
        isempty(missing_deps) || throw(
            ArgumentError(
                "derived default :$key is missing dependencies: " *
                join(string.(missing_deps), ", "),
            ),
        )

        needs_compute = !hasproperty(resolved, key)
        needs_recompute = any(dep -> dep in changed, provider.deps)
        if needs_compute || needs_recompute
            value = derive_default(provider.deriver, factory, context, resolved)
            validate_derived_parameter_result(definition.spec, context, value)
            resolved = merge(resolved, NamedTuple{(key,)}((value,)))
            push!(changed, key)
        end
    end

    return resolved
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


function _groups_for_components(population_groups::NamedTuple, components::Tuple)
    groups = Symbol[]
    for component in components
        hasproperty(population_groups, component) || continue
        for group in getproperty(population_groups, component)
            group in groups || push!(groups, group)
        end
    end
    return Tuple(groups)
end

function _process_interaction_roles(definition, population_groups::NamedTuple)
    consumer_components = Symbol[]
    resource_components = Symbol[]
    for named in values(definition.processes)
        process_kind(named) === :grazing || continue
        process_participants = participants(named)
        hasproperty(process_participants, :consumer) &&
            append!(consumer_components, process_participants.consumer)
        hasproperty(process_participants, :resource) &&
            append!(resource_components, process_participants.resource)
    end
    consumers = _groups_for_components(population_groups, Tuple(unique(consumer_components)))
    resources = _groups_for_components(population_groups, Tuple(unique(resource_components)))
    isempty(consumers) && isempty(resources) && return nothing
    return (consumers=consumers, prey=resources)
end

function _process_manifest_roles(definition, population_groups::NamedTuple)
    growth_components = Symbol[]
    consumer_components = Symbol[]
    for named in values(definition.processes)
        kind = process_kind(named)
        process_participants = participants(named)
        if kind === :growth && hasproperty(process_participants, :population)
            append!(growth_components, process_participants.population)
        elseif kind === :grazing && hasproperty(process_participants, :consumer)
            append!(consumer_components, process_participants.consumer)
        end
    end
    return (
        phytoplankton=_groups_for_components(population_groups, Tuple(unique(growth_components))),
        zooplankton=_groups_for_components(population_groups, Tuple(unique(consumer_components))),
    )
end

function _process_parameter_indices(definition, layout, context, parameter::Symbol)
    selected = Set{Symbol}()
    for applicability in resolve_parameter_applicability(definition, layout)
        applicability.binding.parameter === parameter || continue
        for axis in applicability.axis_tracers, tracer in axis
            tracer in context.plankton_symbols && push!(selected, tracer)
        end
    end
    return [i for (i, tracer) in pairs(context.plankton_symbols) if tracer in selected]
end

function evaluate_process_default(
    provider::ConstDefault, spec, ::AbstractBGCFactory, definition, layout, context, ::Type{T}
) where {T<:Real}
    spec.shape === :scalar || throw(
        ArgumentError("ConstDefault can only be used for scalar parameters (:$(spec.name)).")
    )
    value = provider.value
    return value isa Bool ? value : T(value)
end

function evaluate_process_default(
    provider::FillDefault, spec, factory::AbstractBGCFactory, definition, layout, context, ::Type{T}
) where {T<:Real}
    value = provider.value
    value = value isa Bool ? value : T(value)
    if spec.shape === :vector
        result = fill(zero(value), context.n_total)
        indices = _process_parameter_indices(definition, layout, context, spec.name)
        isempty(indices) && (indices = collect(eachindex(result)))
        result[indices] .= value
        return result
    elseif spec.shape === :matrix
        return evaluate_default(provider, spec, factory, context, T)
    end
    throw(ArgumentError("FillDefault requires vector or matrix parameter storage."))
end

function evaluate_process_default(
    provider::DiameterIndexedVectorDefault,
    spec,
    ::AbstractBGCFactory,
    definition,
    layout,
    context,
    ::Type{T},
) where {T<:Real}
    spec.shape === :vector || throw(
        ArgumentError("DiameterIndexedVectorDefault requires vector parameter storage.")
    )
    indices = _process_parameter_indices(definition, layout, context, spec.name)
    default = T(provider.default)
    return resolve_diameter_indexed_vector(
        T, context.diameters, indices, provider.value; default
    )
end

function build_process_parameter_defaults(
    factory::AbstractBGCFactory, definition, layout, context, ::Type{T}
) where {T<:Real}
    entries = Pair{Symbol,Any}[]
    for parameter_definition in parameter_definitions(factory)
        provider = parameter_definition.default
        (provider isa NoDefault || provider isa DerivedDefault) && continue
        spec = parameter_definition.spec
        value = evaluate_process_default(provider, spec, factory, definition, layout, context, T)
        push!(entries, spec.name => value)
    end
    return (; entries...)
end

function materialize_process_parameter_law_override(
    context, definition, layout, spec, value::AbstractParamDef, ::Type{T}
) where {T<:Real}
    materialization = spec.materialization
    materialization isa DiameterIndexedMaterialization || throw(
        ArgumentError(
            "parameter :$(spec.name) only supports parameter-law overrides for parameters with declared diameter-indexed vector materialization."
        ),
    )
    spec.shape === :vector || throw(
        ArgumentError("parameter :$(spec.name) diameter-indexed materialization requires vector storage")
    )
    indices = _process_parameter_indices(definition, layout, context, spec.name)
    return resolve_diameter_indexed_vector(
        T, context.diameters, indices, value; default=T(materialization.fill_value)
    )
end

function materialize_process_parameter_overrides(
    factory::AbstractBGCFactory,
    context,
    definition,
    layout,
    defaults::NamedTuple,
    overrides::NamedTuple,
    ::Type{T},
) where {T<:Real}
    isempty(overrides) && return overrides
    entries = Pair{Symbol,Any}[]
    for (key, value) in Base.pairs(overrides)
        spec = parameter_spec(factory, key)
        if spec === nothing
            push!(entries, key => value)
        elseif value isa AbstractParamDef
            push!(entries, key => materialize_process_parameter_law_override(
                context, definition, layout, spec, value, T
            ))
        elseif value isa NamedTuple
            spec.shape === :vector || throw(
                ArgumentError("parameter :$key does not support NamedTuple overrides because it is $(spec.shape)-shaped.")
            )
            hasproperty(defaults, key) || throw(
                ArgumentError("parameter :$key has no direct default for partial overrides.")
            )
            push!(entries, key => expand_named_vector_override(
                spec, getproperty(defaults, key), value, context, T
            ))
        else
            push!(entries, key => materialize_parameter_value(spec, value, T))
        end
    end
    return (; entries...)
end

"""
    construct_factory(factory::AbstractBGCFactory; kwargs...) -> bgc

Construct a concrete biogeochemistry model from `factory` and optional
configuration overrides.

Construction proceeds in four stages:

1. Parse the community into a `CommunityContext`.
2. Evaluate `parameter_definitions(factory)` into concrete defaults.
3. Apply user overrides and resolve derived parameter defaults.
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


function _construct_process_factory(
    factory::AbstractBGCFactory,
    recipe::ProcessModelRecipe;
    grid=nothing,
    arch=nothing,
    scalar_type=nothing,
    build_manifest::Bool=false,
)
    if isnothing(grid) && !isnothing(recipe.sinking_tracers)
        grid = BoxModelGrid()
    end
    T = resolve_construction_scalar_type(grid, scalar_type)
    sinking_tracers = convert_sinking_tracers(T, recipe.sinking_tracers)

    if !isnothing(grid)
        arch_grid = architecture(grid)
        if isnothing(arch)
            arch = arch_grid
        elseif typeof(arch) !== typeof(arch_grid)
            throw(ArgumentError(
                "arch=$arch does not match architecture(grid)=$arch_grid. Architecture is determined by the grid."
            ))
        end
    else
        isnothing(arch) && (arch = CPU())
    end

    definition = normalize_model(ModelDefinition(;
        components=recipe.components,
        processes=recipe.processes,
        parameters=parameter_definitions(factory),
    ))
    placeholder_dynamics = NamedTuple{keys(recipe.community)}(
        ntuple(_ -> identity, length(recipe.community))
    )
    validate_community_inputs(placeholder_dynamics, recipe.community)
    interaction_roles = _process_interaction_roles(definition, recipe.population_groups)
    pool_names = Tuple(
        name for name in keys(recipe.components) if getproperty(recipe.components, name) isa Pool
    )
    pool_components = NamedTuple{pool_names}(
        Tuple(getproperty(recipe.components, name) for name in pool_names)
    )
    pool_layout = realize_components(pool_components; scalar_type=T)
    community_context = parse_community(
        T,
        recipe.community;
        biogeochem_tracers=pool_layout.tracer_order,
        interaction_roles,
        parameter_roles=NamedTuple(),
    )
    layout = realize_component_groups(
        recipe.components, recipe.population_groups, community_context
    )
    tracer_names = layout.tracer_order
    auxiliary_fields = driver_identities(definition)

    required = validate_parameter_directory(factory)
    interaction_parameter_overrides = normalize_interaction_overrides(
        factory, community_context, deepcopy(recipe.interaction_overrides)
    )
    validate_override_keys("parameters", recipe.parameter_overrides, required, factory)
    validate_override_keys(
        "interaction_overrides", interaction_parameter_overrides, required, factory
    )

    parameter_defaults = build_process_parameter_defaults(
        factory, definition, layout, community_context, T
    )
    parameter_overrides = materialize_process_parameter_overrides(
        factory,
        community_context,
        definition,
        layout,
        parameter_defaults,
        recipe.parameter_overrides,
        T,
    )
    merged_parameters = merge(
        parameter_defaults, parameter_overrides, interaction_parameter_overrides
    )
    explicit_override_keys = (
        keys(recipe.parameter_overrides)..., keys(interaction_parameter_overrides)...
    )
    merged_parameters = resolve_parameter_defaults(
        factory, community_context, merged_parameters, explicit_override_keys
    )
    missing = Symbol[k for k in required if !hasproperty(merged_parameters, k)]
    isempty(missing) || throw(
        ArgumentError("missing required parameters: $(join(string.(missing), ", "))")
    )
    merged_parameters = finalize_interaction_parameters(
        factory, community_context, merged_parameters
    )
    internal = hasproperty(merged_parameters, :interactions) ? (:interactions,) : ()
    all_keys = (required..., internal...)
    resolved_parameters = NamedTuple{all_keys}(
        Tuple(getproperty(merged_parameters, key) for key in all_keys)
    )
    reject_missing_values(resolved_parameters)
    validate_parameter_shapes(factory, community_context, resolved_parameters, required)
    validate_auxiliary_fields(auxiliary_fields, tracer_names)

    tracers = compile_model_tendencies(
        definition, layout, community_context; target_order=tracer_names
    )
    tracer_index = build_tracer_index(
        community_context,
        tracer_names,
        auxiliary_fields;
        n_biogeochem_tracers=sum(
            length(component_tracers(layout, name)) for name in pool_names
        ),
    )
    plankton_diameter_metadata = Tuple(community_context.diameters)

    bgc = if isnothing(sinking_tracers)
        bgc_factory = define_tracer_functions(
            resolved_parameters,
            tracers;
            auxiliary_fields,
            tracer_index,
        )
        bgc_factory(resolved_parameters; plankton_diameters=plankton_diameter_metadata)
    else
        sinking_velocities = setup_velocity_fields(
            sinking_tracers, grid, recipe.open_bottom
        )
        bgc_factory = define_tracer_functions(
            resolved_parameters,
            tracers;
            auxiliary_fields,
            tracer_index,
            sinking_velocities,
        )
        bgc_factory(
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
            ecological_roles=_process_manifest_roles(definition, recipe.population_groups),
            explicit_override_keys,
            sinking_tracers,
            open_bottom=recipe.open_bottom,
            scalar_type=T,
        )
    else
        nothing
    end

    return on_architecture(arch, bgc), manifest
end

"""Realize a v3 component/process recipe in the supplied execution environment."""
function construct_factory(
    recipe::ProcessModelRecipe; grid=nothing, arch=nothing, scalar_type=nothing
)
    factory = replay_factory(recipe)
    bgc, _ = _construct_process_factory(
        factory, recipe; grid, arch, scalar_type, build_manifest=false
    )
    return bgc
end

"""Realize a v3 component/process recipe and return its resolved manifest."""
function construct_factory_plus_manifest(
    recipe::ProcessModelRecipe; grid=nothing, arch=nothing, scalar_type=nothing
)
    factory = replay_factory(recipe)
    return _construct_process_factory(
        factory, recipe; grid, arch, scalar_type, build_manifest=true
    )
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
    merged_parameters = resolve_parameter_defaults(
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
