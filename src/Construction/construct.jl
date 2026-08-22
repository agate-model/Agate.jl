using OceanBioME: BoxModelGrid, setup_velocity_fields

using Adapt: adapt

import Oceananigans

using Oceananigans.Architectures: architecture, CPU, GPU

using ..ModelFamilies: AbstractModelFamily

using ..Parameters:
    Parameter,
    ConstantDefault,
    DerivedDefault,
    NoDefault,
    DiameterIndexedVectorDefault,
    derive_default,
    parameter_spec

import ..Parameters: parameter_definitions

using ..Configuration:
    axis_indices,
    interaction_axis_metadata,
    parse_community,
    validate_community,
    Population,
    Pool,
    size_structure,
    realize_components,
    component_classes,
    component_tracers,
    _realize_component_groups

using ..Runtime: build_tracer_index

using ..Processes:
    ModelDefinition, normalize_model, driver_identities, participants, formulation,
    uses_living_interactions, resolve_parameter_applicability

using ..Compilation: compile_model_tendencies

using ..Library.Allometry: AbstractParamDef, resolve_diameter_indexed_vector

struct _DefinitionParameterSource{P}
    definitions::P
end

parameter_definitions(source::_DefinitionParameterSource) = source.definitions

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

function validate_parameter_directory(source)
    definitions = parameter_definitions(source)
    definitions isa NamedTuple || throw(
        ArgumentError("parameter_definitions(::$(typeof(source))) must return a NamedTuple"),
    )
    isempty(definitions) && return ()
    all(parameter -> parameter isa Parameter, values(definitions)) || throw(
        ArgumentError("parameter_definitions(::$(typeof(source))) must contain only Parameter values"),
    )

    parameter_keys = Tuple(keys(definitions))
    key_set = Set(parameter_keys)
    for name in parameter_keys
        (name in RESERVED_PARAMETER_KEYS) && throw(
            ArgumentError(
                "parameter_definitions(::$(typeof(source))) declares reserved parameter key :$name.",
            ),
        )
    end

    for (name, parameter) in pairs(definitions)
        spec = parameter.spec
        spec.shape in (:scalar, :vector, :matrix) || throw(
            ArgumentError(
                "parameter_definitions(::$(typeof(source))) declares :$name with invalid shape $(spec.shape).",
            ),
        )

        if spec.shape === :vector
            if spec.axes !== nothing
                spec.axes isa Symbol || throw(
                    ArgumentError(
                        "parameter :$name vector axis must be a Symbol (got $(typeof(spec.axes))).",
                    ),
                )
            end
        elseif spec.shape === :matrix
            if spec.axes !== nothing
                (spec.axes isa Tuple && length(spec.axes) == 2) || throw(
                    ArgumentError(
                        "parameter :$name axes must be a 2-tuple of Symbols (got $(typeof(spec.axes))).",
                    ),
                )
                row_axis, col_axis = spec.axes
                (row_axis isa Symbol && col_axis isa Symbol) || throw(
                    ArgumentError(
                        "parameter :$name axes must be Symbols (got $(spec.axes)).",
                    ),
                )
            end
        else
            spec.axes === nothing || throw(
                ArgumentError(
                    "parameter :$name has axes=$(spec.axes) but is not vector or matrix."
                ),
            )
        end

        provider = parameter.default
        if provider isa DerivedDefault
            for dep in provider.deps
                dep in key_set || throw(
                    ArgumentError(
                        "derived default :$name depends on undeclared parameter :$dep.",
                    ),
                )
            end
        end
    end

    return parameter_keys
end

function derived_default_order(source)
    derived = Pair{Symbol,Any}[
        name => parameter for (name, parameter) in pairs(parameter_definitions(source)) if
        parameter.default isa DerivedDefault
    ]
    isempty(derived) && return ()

    derived_keys = Tuple(first(entry) for entry in derived)
    resolved_keys = Set{Symbol}()
    ordered = Pair{Symbol,Any}[]
    pending = derived

    while !isempty(pending)
        remaining = Pair{Symbol,Any}[]
        progressed = false

        for (name, parameter) in pending
            dependencies = Tuple(
                dep for dep in parameter.default.deps if dep in derived_keys
            )
            if all(dep -> dep in resolved_keys, dependencies)
                push!(ordered, name => parameter)
                push!(resolved_keys, name)
                progressed = true
            else
                push!(remaining, name => parameter)
            end
        end

        progressed || throw(
            ArgumentError(
                "derived parameter defaults contain a dependency cycle among: " *
                join(string.(Tuple(first(entry) for entry in remaining)), ", "),
            ),
        )
        pending = remaining
    end

    return Tuple(ordered)
end

function _local_parameter_axis_classes(definition, layout, parameter::Symbol, shape::Symbol)
    matches = Tuple(
        applicability.axis_classes
        for applicability in resolve_parameter_applicability(definition, layout)
        if applicability.binding.parameter === parameter
    )
    isempty(matches) && throw(
        ArgumentError(
            "parameter :$parameter has local storage but no resolved process applicability",
        ),
    )

    expected_axes = shape === :vector ? 1 : shape === :matrix ? 2 : 0
    all(axes -> length(axes) == expected_axes, matches) || throw(
        ArgumentError(
            "parameter :$parameter local storage does not match its declared $shape shape",
        ),
    )
    first_axes = first(matches)
    all(==(first_axes), matches) || throw(
        ArgumentError(
            "parameter :$parameter supplies incompatible process-local axes; declare explicit storage axes",
        ),
    )
    return first_axes
end

function _parameter_storage_shape(definition, layout, context, name::Symbol, spec)
    spec.shape === :scalar && return ()

    if spec.axes === nothing
        (definition === nothing || layout === nothing) && throw(
            ArgumentError(
                "parameter :$name local storage requires a normalized model definition and component layout",
            ),
        )
        local_axes = _local_parameter_axis_classes(
            definition, layout, name, spec.shape
        )
        return Tuple(length(axis) for axis in local_axes)
    elseif spec.shape === :vector
        n = spec.axes === :plankton ? context.n_total : length(axis_indices(context, spec.axes))
        return (n,)
    end

    row_axis, col_axis = spec.axes
    return (length(axis_indices(context, row_axis)), length(axis_indices(context, col_axis)))
end

function validate_derived_parameter_result(
    name::Symbol, spec, context, value; definition=nothing, layout=nothing
)
    T = context.scalar_type

    if spec.shape === :scalar
        value isa Bool && return nothing
        value isa Number || throw(
            ArgumentError(
                "derived default :$name must be scalar; got $(typeof(value)).",
            ),
        )
        typeof(value) === T || throw(
            ArgumentError(
                "derived default :$name must have type $(T); got $(typeof(value)). No implicit casting is performed.",
            ),
        )
        return nothing
    elseif spec.shape === :vector
        value isa AbstractVector || throw(
            ArgumentError(
                "derived default :$name must be a vector; got $(typeof(value)).",
            ),
        )
        eltype(value) === T || throw(
            ArgumentError(
                "derived default :$name must have eltype $(T); got eltype $(eltype(value)). No implicit casting is performed.",
            ),
        )
        expected = only(_parameter_storage_shape(definition, layout, context, name, spec))
        length(value) == expected || throw(
            ArgumentError(
                "derived default :$name must have length $expected; got $(length(value)).",
            ),
        )
        return nothing
    end

    value isa AbstractMatrix || throw(
        ArgumentError(
            "derived default :$name must be a matrix; got $(typeof(value)).",
        ),
    )
    eltype(value) === T || throw(
        ArgumentError(
            "derived default :$name must have eltype $(T); got eltype $(eltype(value)). No implicit casting is performed.",
        ),
    )

    expected = _parameter_storage_shape(definition, layout, context, name, spec)
    size(value) == expected || throw(
        ArgumentError(
            "derived default :$name must have size $expected; got $(size(value)).",
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
    source,
    context,
    params::NamedTuple,
    explicit_override_keys::Tuple{Vararg{Symbol}};
    derivation_owner=source,
    normalized_definition=nothing,
    layout=nothing,
)
    ordered = derived_default_order(source)
    isempty(ordered) && return params

    override_set = Set{Symbol}(explicit_override_keys)
    changed = Set{Symbol}(explicit_override_keys)
    resolved = params

    for (key, parameter) in ordered
        provider = parameter.default

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
            value = derive_default(provider.deriver, derivation_owner, context, resolved)
            validate_derived_parameter_result(
                key,
                parameter.spec,
                context,
                value;
                definition=normalized_definition,
                layout=layout,
            )
            resolved = merge(resolved, NamedTuple{(key,)}((value,)))
            push!(changed, key)
        end
    end

    return resolved
end

function validate_parameter_shapes(
    source, definition, layout, context, params::NamedTuple, required::Tuple
)
    for k in required
        spec = parameter_spec(source, k)
        spec === nothing && throw(
            ArgumentError(
                "Parameter source $(typeof(source)) is missing a ParameterSpec for parameter :$k.",
            ),
        )

        expected = _parameter_storage_shape(definition, layout, context, k, spec)
        if spec.shape === :vector
            v = getproperty(params, k)
            v isa AbstractVector || throw(
                ArgumentError("parameter :$k must be a vector; got $(typeof(v))."),
            )
            length(v) == only(expected) || throw(
                ArgumentError(
                    "parameter :$k must have length $(only(expected)) (got $(length(v))).",
                ),
            )
        elseif spec.shape === :matrix
            m = getproperty(params, k)
            m isa AbstractMatrix || throw(
                ArgumentError("parameter :$k must be a matrix; got $(typeof(m))."),
            )
            size(m) == expected || throw(
                ArgumentError(
                    "parameter :$k must have size $expected (got $(size(m))).",
                ),
            )
        end
    end

    return nothing
end


function parameter_axis_names(context, axis::Symbol, parameter_name::Symbol)
    axis === :plankton && return context.class_symbols
    throw(ArgumentError("parameter :$parameter_name has unsupported vector axis :$axis."))
end

function expand_named_vector_override(
    name::Symbol, spec, default_value, user_value::NamedTuple, context, ::Type{T}
) where {T<:Real}
    spec.axes === nothing && throw(
        ArgumentError(
            "parameter :$name does not support NamedTuple overrides because it has no named vector axis."
        ),
    )

    names = parameter_axis_names(context, spec.axes, name)
    length(default_value) == length(names) || throw(
        ArgumentError(
            "parameter :$name default length $(length(default_value)) does not match axis :$(spec.axes) length $(length(names))."
        ),
    )

    expanded = copy(default_value)
    for (key, value) in pairs(user_value)
        idx = findfirst(==(key), names)
        if idx === nothing
            expected = join(string.(names), ", ")
            throw(
                ArgumentError(
                    "Unknown key `$(key)` for parameter `$name`. Expected one of: $(expected)."
                ),
            )
        end
        expanded[idx] = value isa Bool ? value : T(value)
    end
    return expanded
end

function materialize_parameter_value(spec, value, ::Type{T}) where {T<:Real}
    if spec.shape === :scalar
        return value isa Bool ? value : T(value)
    elseif spec.shape === :vector
        value isa AbstractVector || return value
        eltype(value) === T && return copy(value)
        out = similar(value, T, axes(value))
        copyto!(out, value)
        return out
    elseif spec.shape === :matrix
        value isa AbstractMatrix || return value
        eltype(value) === T && return copy(value)
        out = similar(value, T, axes(value))
        copyto!(out, value)
        return out
    end

    return value
end

function validate_override_keys(where_, overrides::NamedTuple, required::Tuple, source)
    isempty(overrides) && return nothing

    required_set = Set(required)
    for k in keys(overrides)
        k in required_set && continue
        parameter_spec(source, k) === nothing &&
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
        uses_living_interactions(formulation(named.process)) || continue
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

function _process_parameter_indices(definition, layout, context, parameter::Symbol)
    selected = Set{Symbol}()
    has_applicability = false
    for applicability in resolve_parameter_applicability(definition, layout)
        applicability.binding.parameter === parameter || continue
        has_applicability = true
        for axis in applicability.axis_classes, class in axis
            class in context.class_symbols && push!(selected, class)
        end
    end
    has_applicability || return collect(eachindex(context.class_symbols))
    return [i for (i, class) in pairs(context.class_symbols) if class in selected]
end

function evaluate_process_default(
    provider::ConstantDefault, name::Symbol, spec, source, definition, layout, context, ::Type{T}
) where {T<:Real}
    value = provider.value
    value = value isa Bool ? value : T(value)
    spec.shape === :scalar && return value

    expected = _parameter_storage_shape(definition, layout, context, name, spec)
    if spec.shape === :vector
        if spec.axes === :plankton
            result = fill(zero(value), only(expected))
            indices = _process_parameter_indices(definition, layout, context, name)
            result[indices] .= value
            return result
        end
        return fill(value, only(expected))
    elseif spec.shape === :matrix
        return fill(value, expected...)
    end
    throw(ArgumentError("unsupported parameter shape $(spec.shape) for :$name"))
end

function evaluate_process_default(
    provider::DiameterIndexedVectorDefault,
    name::Symbol,
    spec,
    source,
    definition,
    layout,
    context,
    ::Type{T},
) where {T<:Real}
    spec.shape === :vector || throw(
        ArgumentError("DiameterIndexedVectorDefault requires vector parameter storage.")
    )
    indices = _process_parameter_indices(definition, layout, context, name)
    default = T(provider.default)
    return resolve_diameter_indexed_vector(
        T, context.diameters, indices, provider.value; default
    )
end

function build_process_parameter_defaults(source, definition, layout, context, ::Type{T}) where {T<:Real}
    entries = Pair{Symbol,Any}[]
    for (name, parameter) in pairs(parameter_definitions(source))
        provider = parameter.default
        (provider isa NoDefault || provider isa DerivedDefault) && continue
        spec = parameter.spec
        value = evaluate_process_default(provider, name, spec, source, definition, layout, context, T)
        push!(entries, name => value)
    end
    return (; entries...)
end

function materialize_process_parameter_law_override(
    context, definition, layout, name::Symbol, parameter, value::AbstractParamDef, ::Type{T}
) where {T<:Real}
    spec = parameter.spec
    provider = parameter.default
    provider isa DiameterIndexedVectorDefault || throw(
        ArgumentError(
            "parameter :$name only supports parameter-law overrides with a diameter-indexed vector default provider (DiameterIndexedVectorDefault)."
        ),
    )
    spec.shape === :vector || throw(
        ArgumentError("parameter :$name diameter-indexed override requires vector storage")
    )
    indices = _process_parameter_indices(definition, layout, context, name)
    return resolve_diameter_indexed_vector(
        T, context.diameters, indices, value; default=T(provider.default)
    )
end

function materialize_process_parameter_overrides(
    source,
    context,
    definition,
    layout,
    defaults::NamedTuple,
    overrides::NamedTuple,
    ::Type{T},
) where {T<:Real}
    isempty(overrides) && return overrides
    entries = Pair{Symbol,Any}[]
    definitions = parameter_definitions(source)
    for (key, value) in Base.pairs(overrides)
        parameter = hasproperty(definitions, key) ? getproperty(definitions, key) : nothing
        if parameter === nothing
            push!(entries, key => value)
            continue
        end
        spec = parameter.spec
        if value isa AbstractParamDef
            push!(entries, key => materialize_process_parameter_law_override(
                context, definition, layout, key, parameter, value, T
            ))
        elseif value isa NamedTuple
            spec.shape === :vector || throw(
                ArgumentError("parameter :$key does not support NamedTuple overrides because it is $(spec.shape)-shaped.")
            )
            hasproperty(defaults, key) || throw(
                ArgumentError("parameter :$key has no direct default for partial overrides.")
            )
            push!(entries, key => expand_named_vector_override(
                key, spec, getproperty(defaults, key), value, context, T
            ))
        else
            push!(entries, key => materialize_parameter_value(spec, value, T))
        end
    end
    return (; entries...)
end

function _intrinsic_population_groups(components::NamedTuple)
    names = Tuple(
        name for name in keys(components) if getproperty(components, name) isa Population
    )
    return NamedTuple{names}(Tuple((name,) for name in names))
end

@inline _unspecified_diameter(::Type{T}) where {T<:AbstractFloat} = T(NaN)
@inline _unspecified_diameter(::Type{T}) where {T<:Real} = zero(T)

function _intrinsic_population_community(
    components::NamedTuple,
    population_groups::NamedTuple,
    layout,
    ::Type{T},
) where {T<:Real}
    names = keys(population_groups)
    specs = ntuple(length(names)) do i
        name = names[i]
        component = getproperty(components, name)
        structure = size_structure(component)
        diameters = isnothing(structure) ? T[_unspecified_diameter(T)] : structure
        class_symbols = component_classes(layout, name)
        return (; diameters, pft=PFTSpecification(), class_symbols)
    end
    return NamedTuple{names}(specs)
end

function _realize_process_definition(
    definition,
    ::Type{T};
    population_groups=nothing,
    community=nothing,
) where {T<:Real}
    intrinsic = isnothing(population_groups) && isnothing(community)
    xor(isnothing(population_groups), isnothing(community)) && throw(
        ArgumentError("population_groups and community must be supplied together."),
    )

    groups = intrinsic ? _intrinsic_population_groups(definition.components) : population_groups
    pool_names = Tuple(
        name for name in keys(definition.components) if
        getproperty(definition.components, name) isa Pool
    )

    layout = nothing
    pool_tracers = ()
    community_input = community

    if intrinsic
        population_names = keys(groups)
        realization_names = (pool_names..., population_names...)
        realization_components = NamedTuple{realization_names}(
            Tuple(getproperty(definition.components, name) for name in realization_names)
        )
        layout = realize_components(realization_components; scalar_type=T)
        pool_tracers = Tuple(
            tracer for name in pool_names for tracer in component_tracers(layout, name)
        )
        community_input = _intrinsic_population_community(
            definition.components, groups, layout, T
        )
    else
        pool_components = NamedTuple{pool_names}(
            Tuple(getproperty(definition.components, name) for name in pool_names)
        )
        pool_layout = realize_components(pool_components; scalar_type=T)
        pool_tracers = pool_layout.tracer_order
    end

    validate_community(community_input)
    interaction_roles = _process_interaction_roles(definition, groups)
    context = parse_community(
        T,
        community_input;
        biogeochem_tracers=Tuple(pool_tracers),
        interaction_roles,
    )

    intrinsic || (layout = _realize_component_groups(definition.components, groups, context, pool_layout))
    return (; layout, context, population_groups=groups, pool_names)
end

function _construct_process_definition(
    definition::ModelDefinition;
    population_groups=nothing,
    community=nothing,
    parameter_overrides::NamedTuple=(;),
    sinking_tracers=nothing,
    open_bottom::Bool=true,
    grid=nothing,
    arch=nothing,
    scalar_type=nothing,
    build_manifest::Bool=false,
    derivation_owner=nothing,
    manifest_family=nothing,
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
            throw(ArgumentError(
                "arch=$arch does not match architecture(grid)=$arch_grid. Architecture is determined by the grid."
            ))
        end
    else
        isnothing(arch) && (arch = CPU())
    end

    normalized = normalize_model(definition)
    definitions = isnothing(normalized.parameters) ? (;) : normalized.parameters
    parameter_source = _DefinitionParameterSource(definitions)
    isnothing(derivation_owner) && (derivation_owner = definition)
    isnothing(definition.parameters) && !isempty(normalized.parameter_requirements) && throw(
        ArgumentError(
            "construct(definition) requires ModelDefinition.parameters to provide the process parameter requirements."
        ),
    )
    realization = _realize_process_definition(
        normalized, T; population_groups, community
    )
    layout = realization.layout
    community_context = realization.context
    tracer_names = layout.tracer_order
    auxiliary_fields = driver_identities(normalized)

    required = validate_parameter_directory(parameter_source)
    validate_override_keys(
        "parameters", parameter_overrides, required, parameter_source
    )

    parameter_defaults = build_process_parameter_defaults(
        parameter_source, normalized, layout, community_context, T
    )
    materialized_overrides = materialize_process_parameter_overrides(
        parameter_source,
        community_context,
        normalized,
        layout,
        parameter_defaults,
        parameter_overrides,
        T,
    )
    merged_parameters = merge(parameter_defaults, materialized_overrides)
    explicit_override_keys = Tuple(keys(parameter_overrides))
    merged_parameters = resolve_parameter_defaults(
        parameter_source,
        community_context,
        merged_parameters,
        explicit_override_keys;
        derivation_owner,
        normalized_definition=normalized,
        layout,
    )
    missing = Symbol[k for k in required if !hasproperty(merged_parameters, k)]
    isempty(missing) || throw(
        ArgumentError("missing required parameters: $(join(string.(missing), ", "))")
    )
    resolved_parameters = NamedTuple{required}(
        Tuple(getproperty(merged_parameters, key) for key in required)
    )
    reject_missing_values(resolved_parameters)
    validate_parameter_shapes(
        parameter_source, normalized, layout, community_context, resolved_parameters, required
    )
    validate_auxiliary_fields(auxiliary_fields, tracer_names)

    tracers = compile_model_tendencies(
        normalized, layout, community_context; target_order=tracer_names
    )
    tracer_index = build_tracer_index(
        community_context,
        tracer_names,
        auxiliary_fields;
        n_biogeochem_tracers=sum(
            (length(component_tracers(layout, name)) for name in realization.pool_names);
            init=0,
        ),
    )
    plankton_diameter_metadata = Tuple(community_context.diameters)
    interaction_axes = interaction_axis_metadata(parameter_source, community_context)

    bgc = if isnothing(sinking_tracers)
        bgc_factory = define_tracer_functions(
            resolved_parameters,
            tracers;
            auxiliary_fields,
            tracer_index,
        )
        bgc_factory(
            resolved_parameters;
            plankton_diameters=plankton_diameter_metadata,
            interaction_axes,
        )
    else
        sinking_velocities = setup_velocity_fields(sinking_tracers, grid, open_bottom)
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
            interaction_axes,
        )
    end

    manifest = if build_manifest
        isnothing(manifest_family) && throw(
            ArgumentError("a registered model family is required to capture a replay manifest")
        )
        capture_model_manifest(
            manifest_family,
            resolved_parameters,
            community_context;
            tracer_order=tracer_names,
            auxiliary_fields,
            explicit_override_keys,
            sinking_tracers,
            open_bottom,
            scalar_type=T,
        )
    else
        nothing
    end

    return on_architecture(arch, bgc), manifest
end

function _construct_registered_model(
    family::AbstractModelFamily,
    recipe::ProcessModelRecipe;
    grid=nothing,
    arch=nothing,
    scalar_type=nothing,
    build_manifest::Bool=false,
)
    definition = ModelDefinition(;
        components=recipe.components,
        processes=recipe.processes,
        parameters=parameter_definitions(family),
    )
    return _construct_process_definition(
        definition;
        population_groups=recipe.population_groups,
        community=recipe.community,
        parameter_overrides=recipe.parameter_overrides,
        sinking_tracers=recipe.sinking_tracers,
        open_bottom=recipe.open_bottom,
        grid,
        arch,
        scalar_type,
        build_manifest,
        derivation_owner=family,
        manifest_family=family,
    )
end

"""
    construct(definition::ModelDefinition; kwargs...) -> bgc

Construct a model directly from authored components, named processes, and parameter
definitions. Population and pool size structures are realized from the definition,
process participation determines interaction axes and required auxiliary drivers, and
runtime tracer equations are compiled during setup.

`parameter_overrides` supplies concrete parameter values over the defaults declared in
`definition.parameters`, including explicit axis-sized interaction matrices. Runtime grid,
architecture, and scalar precision remain execution choices rather than part of the
scientific definition.
"""
function construct(
    definition::ModelDefinition;
    parameter_overrides::NamedTuple=(;),
    sinking_tracers=nothing,
    open_bottom::Bool=true,
    grid=nothing,
    arch=nothing,
    scalar_type=nothing,
)
    bgc, _ = _construct_process_definition(
        definition;
        parameter_overrides,
        sinking_tracers,
        open_bottom,
        grid,
        arch,
        scalar_type,
    )
    return bgc
end


"""Realize a component/process recipe in the supplied execution environment."""
function construct(
    recipe::ProcessModelRecipe; grid=nothing, arch=nothing, scalar_type=nothing
)
    family = replay_family(recipe)
    bgc, _ = _construct_registered_model(
        family, recipe; grid, arch, scalar_type, build_manifest=false
    )
    return bgc
end

"""Realize a component/process recipe and return its resolved manifest."""
function construct_plus_manifest(
    recipe::ProcessModelRecipe; grid=nothing, arch=nothing, scalar_type=nothing
)
    family = replay_family(recipe)
    return _construct_registered_model(
        family, recipe; grid, arch, scalar_type, build_manifest=true
    )
end
