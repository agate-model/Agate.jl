using OceanBioME: BoxModelGrid, setup_velocity_fields

using Adapt: adapt

import Oceananigans

using Oceananigans.Architectures: architecture, CPU, GPU

using ..ModelFamilies: AbstractModelFamily

using ..Parameters:
    ConstantDefault,
    DerivedDefault,
    NoDefault,
    DiameterIndexedVectorDefault,
    derive_default

using ..Configuration:
    canonicalize_population_realization, interaction_axis_metadata, realize_model_layout,
    model_metadata


using ..Processes:
    ModelDefinition, canonicalize_model, driver_identities, participants, formulation,
    uses_living_interactions, build_parameter_plan, planned_parameter,
    runtime_parameter_values, parameter_plan_metadata, validate_realized_science

using ..Compilation: CompileContext, compile_model_tendencies

using ..Library.Allometry: AbstractParamDef, resolve_diameter_indexed_vector

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

function validate_parameter_value(parameter, value, ::Type{T}; derived::Bool=false) where {T<:Real}
    rank = parameter.rank
    name = parameter.name
    shape = parameter.storage_shape

    if rank == 0
        value isa Bool && return nothing
        value isa Number || throw(ArgumentError(
            "parameter :$name must be scalar; got $(typeof(value)).",
        ))
        if derived && typeof(value) !== T
            throw(ArgumentError(
                "derived default :$name must have type $(T); got $(typeof(value)). No implicit casting is performed.",
            ))
        end
        return nothing
    elseif rank == 1
        value isa AbstractVector || throw(ArgumentError(
            "parameter :$name must be a vector; got $(typeof(value)).",
        ))
        length(value) == only(shape) || throw(ArgumentError(
            "parameter :$name must have length $(only(shape)) (got $(length(value))).",
        ))
        if derived && eltype(value) !== T
            throw(ArgumentError(
                "derived default :$name must have eltype $(T); got eltype $(eltype(value)). No implicit casting is performed.",
            ))
        end
        return nothing
    elseif rank == 2
        value isa AbstractMatrix || throw(ArgumentError(
            "parameter :$name must be a matrix; got $(typeof(value)).",
        ))
        size(value) == shape || throw(ArgumentError(
            "parameter :$name must have size $shape (got $(size(value))).",
        ))
        if derived && eltype(value) !== T
            throw(ArgumentError(
                "derived default :$name must have eltype $(T); got eltype $(eltype(value)). No implicit casting is performed.",
            ))
        end
        return nothing
    end

    throw(ArgumentError("parameter :$name has unsupported rank $rank"))
end

function validate_parameter_storage(plan, values::NamedTuple, ::Type{T}) where {T<:Real}
    for name in keys(plan.parameters)
        validate_parameter_value(getproperty(plan.parameters, name), getproperty(values, name), T)
    end
    return nothing
end

function parameter_axis_names(parameter)
    parameter.rank == 1 || throw(ArgumentError(
        "parameter :$(parameter.name) does not have one vector storage axis",
    ))
    return only(parameter.storage_labels)
end

function expand_named_vector_override(
    parameter, default_value, user_value::NamedTuple, ::Type{T}
) where {T<:Real}
    name = parameter.name
    names = parameter_axis_names(parameter)
    length(default_value) == length(names) || throw(ArgumentError(
        "parameter :$name default length $(length(default_value)) does not match its planned axis length $(length(names)).",
    ))

    expanded = copy(default_value)
    for (key, value) in pairs(user_value)
        idx = findfirst(==(key), names)
        if idx === nothing
            expected = join(string.(names), ", ")
            throw(ArgumentError(
                "Unknown key `$(key)` for parameter `$name`. Expected one of: $(expected).",
            ))
        end
        expanded[idx] = value isa Bool ? value : T(value)
    end
    return expanded
end

function materialize_parameter_value(parameter, value, ::Type{T}) where {T<:Real}
    rank = parameter.rank
    if rank == 0
        return value isa Bool ? value : T(value)
    elseif rank == 1
        value isa AbstractVector || return value
        eltype(value) === T && return copy(value)
        out = similar(value, T, axes(value))
        copyto!(out, value)
        return out
    elseif rank == 2
        value isa AbstractMatrix || return value
        eltype(value) === T && return copy(value)
        out = similar(value, T, axes(value))
        copyto!(out, value)
        return out
    end
    throw(ArgumentError("parameter :$(parameter.name) has unsupported rank $rank"))
end

function validate_override_keys(plan, overrides::NamedTuple)
    for key in keys(overrides)
        hasproperty(plan.parameters, key) || throw(
            ArgumentError("parameters: unknown parameter key :$key."),
        )
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

function evaluate_process_default(
    provider::ConstantDefault, parameter, ::Type{T}
) where {T<:Real}
    value = provider.value
    value = value isa Bool ? value : T(value)
    rank = parameter.rank
    rank == 0 && return value

    expected = parameter.storage_shape
    if rank == 1
        return fill(value, only(expected))
    elseif rank == 2
        return fill(value, expected...)
    end
    throw(ArgumentError("parameter :$(parameter.name) has unsupported rank $rank"))
end

function evaluate_process_default(
    provider::DiameterIndexedVectorDefault, parameter, ::Type{T}
) where {T<:Real}
    parameter.rank == 1 || throw(
        ArgumentError("DiameterIndexedVectorDefault requires vector parameter storage."),
    )
    diameters = parameter.storage_diameters
    diameters === nothing && throw(ArgumentError(
        "parameter :$(parameter.name) has no realized diameter axis for DiameterIndexedVectorDefault",
    ))
    default = T(provider.default)
    return resolve_diameter_indexed_vector(
        T, diameters, Tuple(eachindex(diameters)), provider.value; default
    )
end

function build_process_parameter_defaults(plan, ::Type{T}) where {T<:Real}
    entries = Pair{Symbol,Any}[]
    for (name, parameter) in pairs(plan.parameters)
        provider = parameter.definition.default
        (provider isa NoDefault || provider isa DerivedDefault) && continue
        push!(entries, name => evaluate_process_default(provider, parameter, T))
    end
    return (; entries...)
end

function materialize_process_parameter_law_override(
    parameter, value::AbstractParamDef, ::Type{T}
) where {T<:Real}
    provider = parameter.definition.default
    provider isa DiameterIndexedVectorDefault || throw(ArgumentError(
        "parameter :$(parameter.name) only supports parameter-law overrides with a diameter-indexed vector default provider (DiameterIndexedVectorDefault).",
    ))
    parameter.rank == 1 || throw(ArgumentError(
        "parameter :$(parameter.name) diameter-indexed override requires vector storage",
    ))
    diameters = parameter.storage_diameters
    diameters === nothing && throw(ArgumentError(
        "parameter :$(parameter.name) has no realized diameter axis for a diameter-indexed override",
    ))
    return resolve_diameter_indexed_vector(
        T,
        diameters,
        Tuple(eachindex(diameters)),
        value;
        default=T(provider.default),
    )
end

function materialize_process_parameter_overrides(
    plan,
    defaults::NamedTuple,
    overrides::NamedTuple,
    ::Type{T},
) where {T<:Real}
    isempty(overrides) && return overrides
    entries = Pair{Symbol,Any}[]
    for (key, value) in pairs(overrides)
        parameter = planned_parameter(plan, key)
        if value isa AbstractParamDef
            push!(entries, key => materialize_process_parameter_law_override(
                parameter, value, T
            ))
        elseif value isa NamedTuple
            parameter.rank == 1 || throw(ArgumentError(
                "parameter :$key does not support NamedTuple overrides because it is not vector-valued.",
            ))
            hasproperty(defaults, key) || throw(ArgumentError(
                "parameter :$key has no direct default for partial overrides.",
            ))
            push!(entries, key => expand_named_vector_override(
                parameter, getproperty(defaults, key), value, T
            ))
        else
            push!(entries, key => materialize_parameter_value(parameter, value, T))
        end
    end
    return (; entries...)
end

"""Resolve one-level `DerivedDefault` values from already-materialized parameters."""
function resolve_parameter_defaults(
    plan,
    layout,
    params::NamedTuple;
    derivation_owner,
)
    resolved = params
    T = layout.scalar_type

    for (key, parameter) in pairs(plan.parameters)
        provider = parameter.definition.default
        provider isa DerivedDefault || continue
        hasproperty(resolved, key) && continue

        missing_deps = Tuple(dep for dep in provider.deps if !hasproperty(params, dep))
        isempty(missing_deps) || throw(ArgumentError(
            "derived default :$key is missing dependencies: " * join(string.(missing_deps), ", "),
        ))
        dependencies = NamedTuple{provider.deps}(
            Tuple(getproperty(params, dep) for dep in provider.deps)
        )

        value = derive_default(provider.deriver, derivation_owner, layout, dependencies)
        validate_parameter_value(parameter, value, T; derived=true)
        resolved = merge(resolved, NamedTuple{(key,)}((value,)))
    end
    return resolved
end

function _realize_process_definition(
    definition,
    ::Type{T};
    population_groups=nothing,
    group_diameters=nothing,
    auxiliary_fields::Tuple=(),
) where {T<:Real}
    population_groups, group_diameters = canonicalize_population_realization(
        definition.components, population_groups, group_diameters
    )
    interaction_roles = _process_interaction_roles(definition, population_groups)
    return realize_model_layout(
        definition.components,
        population_groups,
        group_diameters;
        scalar_type=T,
        interaction_roles,
        auxiliary_fields,
    )
end

function _construct_process_definition(
    definition::ModelDefinition;
    population_groups=nothing,
    group_diameters=nothing,
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

    canonical = canonicalize_model(definition)
    isnothing(derivation_owner) && (derivation_owner = definition)
    isnothing(definition.parameters) && !isempty(canonical.parameter_bindings) && throw(
        ArgumentError(
            "construct(definition) requires ModelDefinition.parameters for the declared process parameter slots."
        ),
    )
    auxiliary_fields = driver_identities(canonical)
    layout = _realize_process_definition(
        canonical, T; population_groups, group_diameters, auxiliary_fields
    )
    tracer_names = layout.tracer_order
    parameter_plan = build_parameter_plan(canonical, layout)
    required = Tuple(keys(parameter_plan.parameters))
    validate_override_keys(parameter_plan, parameter_overrides)

    parameter_defaults = build_process_parameter_defaults(parameter_plan, T)
    materialized_overrides = materialize_process_parameter_overrides(
        parameter_plan, parameter_defaults, parameter_overrides, T
    )
    explicit_override_keys = Tuple(keys(parameter_overrides))
    resolved = resolve_parameter_defaults(
        parameter_plan,
        layout,
        merge(parameter_defaults, materialized_overrides);
        derivation_owner,
    )
    missing = Symbol[key for key in required if !hasproperty(resolved, key)]
    isempty(missing) || throw(
        ArgumentError("missing required parameters: $(join(string.(missing), ", "))")
    )
    resolved_parameters = NamedTuple{required}(
        Tuple(getproperty(resolved, key) for key in required)
    )
    reject_missing_values(resolved_parameters)
    validate_parameter_storage(parameter_plan, resolved_parameters, T)
    validate_realized_science(canonical, layout, parameter_plan, resolved_parameters)

    runtime_parameters = runtime_parameter_values(parameter_plan, resolved_parameters)
    compile_context = CompileContext(canonical, layout, parameter_plan)
    equations = compile_model_tendencies(compile_context; target_order=tracer_names)
    interaction_axes = interaction_axis_metadata(parameter_plan, layout)
    metadata = model_metadata(
        layout; interaction_axes, parameter_axes=parameter_plan_metadata(parameter_plan)
    )
    sinking_velocities = isnothing(sinking_tracers) ? nothing :
        setup_velocity_fields(sinking_tracers, grid, open_bottom)
    bgc = AgateBGC(
        runtime_parameters, equations, auxiliary_fields, sinking_velocities, metadata
    )

    manifest = if build_manifest
        isnothing(manifest_family) && throw(
            ArgumentError("a registered model family is required to capture a replay manifest")
        )
        capture_model_manifest(
            manifest_family,
            resolved_parameters,
            layout,
            parameter_plan;
            interaction_axes,
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
    realization::NamedTuple;
    grid=nothing,
    arch=nothing,
    scalar_type=nothing,
    build_manifest::Bool=false,
)
    return _construct_process_definition(
        ModelDefinition(family);
        realization...,
        grid,
        arch,
        scalar_type,
        build_manifest,
        derivation_owner=family,
        manifest_family=family,
    )
end

function _construct_family(inputs::NamedTuple)
    bgc, _ = _construct_registered_model(
        inputs.family, _family_realization(inputs); inputs.execution...
    )
    return bgc
end

function _construct_recipe(
    recipe::ProcessModelRecipe;
    grid=nothing,
    arch=nothing,
    scalar_type=nothing,
    build_manifest::Bool=false,
)
    return _construct_registered_model(
        replay_family(recipe),
        _family_realization(recipe);
        grid,
        arch,
        scalar_type,
        build_manifest,
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


"""Replay a versioned family recipe in the supplied execution environment."""
function construct(
    recipe::ProcessModelRecipe; grid=nothing, arch=nothing, scalar_type=nothing
)
    bgc, _ = _construct_recipe(recipe; grid, arch, scalar_type)
    return bgc
end

"""Replay a versioned family recipe and return its resolved manifest."""
function construct_plus_manifest(
    recipe::ProcessModelRecipe; grid=nothing, arch=nothing, scalar_type=nothing
)
    return _construct_recipe(recipe; grid, arch, scalar_type, build_manifest=true)
end
