using OceanBioME: BoxModelGrid, setup_velocity_fields

using Adapt: adapt

import Oceananigans

using Oceananigans.Architectures: architecture, CPU, GPU

using ..ModelFamilies: AbstractModelFamily

using ..Components:
    canonicalize_plankton_realization, realize_model_layout, model_metadata

using ..Processes:
    ModelDefinition, canonicalize_model, driver_identities, build_parameter_plan,
    runtime_parameter_values, parameter_plan_metadata, validate_realized_parameters

using ..Compilation: CompileContext, compile_model_tendencies

"""Move `x` to the requested Oceananigans architecture."""
function on_architecture(arch, x)
    arch === nothing && return x
    arch isa CPU && return x

    return adapt(architecture_array_type(arch), x)
end

function on_architecture(
    arch, bgc::AgateBGC{Parameters,Equations,SinkingVelocities,Metadata,AuxiliaryFields}
) where {Parameters,Equations,SinkingVelocities,Metadata,AuxiliaryFields}
    arch === nothing && return bgc
    arch isa CPU && return bgc

    to = architecture_array_type(arch)
    # Move runtime storage while retaining metadata on the host-side object.
    return AgateBGC(
        adapt(to, bgc.parameters),
        adapt(to, bgc.equations),
        AuxiliaryFields,
        adapt(to, bgc.sinking_velocities),
        bgc.metadata,
    )
end

"""Return the preferred array storage type for `arch`."""
function architecture_array_type(arch)
    arch isa CPU && return Array
    arch isa GPU && return Oceananigans.Architectures.array_type(arch)
    return Array
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
    names = Tuple(sort!(collect(keys(sinking_tracers)); by=String))
    return NamedTuple{names}(
        Tuple(convert(T, getproperty(sinking_tracers, name)) for name in names)
    )
end

function validate_sinking_tracers(sinking_tracers, layout)
    isnothing(sinking_tracers) && return nothing
    sinking_tracers isa NamedTuple || throw(ArgumentError(
        "sinking_tracers must be a NamedTuple mapping realized tracer names to nonnegative finite speeds."
    ))

    realized = Set(layout.tracer_order)
    for (tracer, speed) in pairs(sinking_tracers)
        tracer in realized || throw(ArgumentError(
            "Unknown sinking tracer :$tracer; expected a realized tracer in $(layout.tracer_order)."
        ))
        speed isa Real && !(speed isa Bool) || throw(ArgumentError(
            "Sinking speed for :$tracer must be a real number; got $(repr(speed))."
        ))
        isfinite(speed) || throw(ArgumentError(
            "Sinking speed for :$tracer must be finite; got $(repr(speed))."
        ))
        speed >= zero(speed) || throw(ArgumentError(
            "Sinking speed for :$tracer must be nonnegative; got $(repr(speed))."
        ))
    end
    return nothing
end

"""Normalize named-family `(n=0,)` shorthand to one implicit SizeClass (`nothing`)."""
@inline function normalize_pft_size_structure(specification)
    specification isa NamedTuple || return specification
    keys(specification) == (:n,) || return specification
    n = specification.n
    return n isa Integer && !(n isa Bool) && n == 0 ? nothing : specification
end


function _realize_process_definition(
    definition,
    ::Type{T};
    plankton_pfts=nothing,
    auxiliary_fields::Tuple=(),
) where {T<:Real}
    plankton_pfts = canonicalize_plankton_realization(definition.components, plankton_pfts)
    return realize_model_layout(
        definition.components, plankton_pfts; scalar_type=T, auxiliary_fields
    )
end

function _construct_process_definition(
    definition::ModelDefinition;
    plankton_pfts=nothing,
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
        canonical, T; plankton_pfts, auxiliary_fields
    )
    validate_sinking_tracers(sinking_tracers, layout)
    sinking_tracers = convert_sinking_tracers(T, sinking_tracers)
    validate_sinking_tracers(sinking_tracers, layout)
    tracer_names = layout.tracer_order
    parameter_plan = build_parameter_plan(canonical, layout)
    required = Tuple(keys(parameter_plan.parameters))
    validate_override_keys(parameter_plan, parameter_overrides)

    parameter_defaults = materialize_parameter_defaults(
        parameter_plan, T, parameter_overrides
    )
    materialized_overrides = materialize_parameter_overrides(
        parameter_plan, parameter_defaults, parameter_overrides, T
    )
    explicit_override_keys = Tuple(keys(parameter_overrides))
    resolved = resolve_derived_parameter_defaults(
        parameter_plan,
        layout,
        merge(parameter_defaults, materialized_overrides);
        derivation_owner,
    )
    missing_names = Symbol[key for key in required if !hasproperty(resolved, key)]
    isempty(missing_names) || throw(
        ArgumentError("missing required parameters: $(join(string.(missing_names), ", "))")
    )
    resolved_parameters = NamedTuple{required}(
        Tuple(getproperty(resolved, key) for key in required)
    )
    reject_missing_parameter_values(resolved_parameters)
    validate_parameter_storage(parameter_plan, resolved_parameters, T)
    validate_realized_parameters(canonical, resolved_parameters)

    runtime_parameters = runtime_parameter_values(parameter_plan, resolved_parameters)
    compile_context = CompileContext(canonical, layout, parameter_plan)
    equations = compile_model_tendencies(compile_context; target_order=tracer_names)
    metadata = model_metadata(
        layout; parameter_axes=parameter_plan_metadata(canonical, parameter_plan)
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

function _construct_recipe(
    recipe::ModelRecipe;
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
    construct(family::AbstractModelFamily;
              plankton_pfts, parameter_overrides=(;),
              sinking_tracers=nothing, open_bottom=true, grid=nothing,
              arch=nothing, scalar_type=nothing) -> bgc

Construct a registered model family from its resolved family realization. This is the
supported construction seam for external family packages after their own user-facing
constructor syntax has been translated into the nested `plankton_pfts` mapping and
parameter overrides. Runtime grid, architecture, and scalar precision remain execution
choices.
"""
function construct(
    family::AbstractModelFamily;
    plankton_pfts::NamedTuple,
    parameter_overrides::NamedTuple=(;),
    sinking_tracers=nothing,
    open_bottom::Bool=true,
    grid=nothing,
    arch=nothing,
    scalar_type=nothing,
)
    realization = (; plankton_pfts, parameter_overrides, sinking_tracers, open_bottom)
    bgc, _ = _construct_registered_model(family, realization; grid, arch, scalar_type)
    return bgc
end

"""
    construct(definition::ModelDefinition; kwargs...) -> bgc

Construct a model directly from authored components, named processes, and parameter
definitions. `plankton_pfts` optionally replaces each logical plankton component's intrinsic
size structure with the same named PFT realization vocabulary used by registered families.
Process participation determines interaction axes and required auxiliary drivers, and runtime
tracer equations are compiled during setup.

`parameter_overrides` supplies concrete parameter values over the defaults declared in
`definition.parameters`, including explicit axis-sized interaction matrices. Runtime grid,
architecture, and scalar precision remain execution choices rather than part of the
scientific definition.
"""
function construct(
    definition::ModelDefinition;
    plankton_pfts=nothing,
    parameter_overrides::NamedTuple=(;),
    sinking_tracers=nothing,
    open_bottom::Bool=true,
    grid=nothing,
    arch=nothing,
    scalar_type=nothing,
)
    bgc, _ = _construct_process_definition(
        definition;
        plankton_pfts,
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
    recipe::ModelRecipe; grid=nothing, arch=nothing, scalar_type=nothing
)
    bgc, _ = _construct_recipe(recipe; grid, arch, scalar_type)
    return bgc
end

"""Replay a versioned family recipe and return its resolved manifest."""
function construct_plus_manifest(
    recipe::ModelRecipe; grid=nothing, arch=nothing, scalar_type=nothing
)
    return _construct_recipe(recipe; grid, arch, scalar_type, build_manifest=true)
end
