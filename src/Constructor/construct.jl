using OceanBioME: BoxModelGrid, setup_velocity_fields

using Agate.Utils:
    AbstractBGCFactory,
    build_tracer_expressions,
    create_parameters,
    define_tracer_functions,
    normalize_interactions,
    parse_community,
    required_parameters,
    validate_plankton_inputs

using Agate.Models:
    default_plankton_dynamics,
    default_plankton_args,
    default_biogeochem_dynamics,
    default_biogeochem_args

"""Gather required parameter symbols from registered dynamics."""
function _required_from_dynamics(plankton_dynamics::NamedTuple, biogeochem_dynamics::NamedTuple)
    required = Symbol[]
    for f in values(plankton_dynamics)
        append!(required, required_parameters(f))
    end
    for f in values(biogeochem_dynamics)
        append!(required, required_parameters(f))
    end
    return unique(required)
end

function _validate_interaction_matrices(mats::NamedTuple, required_mats::Vector{Symbol}, n::Int)
    issues = String[]

    provided = collect(keys(mats))
    missing = setdiff(required_mats, provided)
    !isempty(missing) && push!(issues, "missing required interaction matrices: $(missing)")

    for (k, M) in pairs(mats)
        if !(M isa AbstractMatrix)
            push!(issues, "interaction $(k) must be a matrix, got $(typeof(M))")
            continue
        end
        if size(M, 1) != n || size(M, 2) != n
            push!(issues, "interaction $(k) must have size ($(n), $(n)), got $(size(M))")
        end
    end

    if !isempty(issues)
        throw(ArgumentError(join(unique(issues), "\n")))
    end

    return nothing
end

"""
    construct(factory::AbstractBGCFactory; kwargs...) -> Type

Compile a concrete biogeochemistry type from a factory and optional overrides.

The return value is a *type*; instantiate with `bgc = bgc_type()`.
"""

function construct(
    factory::AbstractBGCFactory;
    FT::Type{<:AbstractFloat}=Float64,
    plankton_dynamics=default_plankton_dynamics(factory),
    plankton_args=default_plankton_args(factory, FT),
    biogeochem_dynamics=default_biogeochem_dynamics(factory),
    biogeochem_args=default_biogeochem_args(factory, FT),
    interactions=nothing,
    sinking_tracers=nothing,
    grid=BoxModelGrid(),
    open_bottom::Bool=true,
)
    return _construct_impl(factory, FT, plankton_dynamics, plankton_args, biogeochem_dynamics, biogeochem_args,
                           interactions, sinking_tracers, grid, open_bottom)
end

function _construct_impl(
    factory::AbstractBGCFactory,
    FT::Type{<:AbstractFloat},
    plankton_dynamics,
    plankton_args,
    biogeochem_dynamics,
    biogeochem_args,
    interactions,
    sinking_tracers,
    grid,
    open_bottom::Bool,
)
    Base.@nospecialize factory plankton_dynamics plankton_args biogeochem_dynamics biogeochem_args interactions sinking_tracers grid

    # 1) Fill defaults and validate input shapes.
    validate_plankton_inputs(plankton_dynamics, plankton_args)

    biogeochem_dynamics isa NamedTuple || throw(ArgumentError("biogeochem_dynamics must be a NamedTuple"))

    # 2) Parse community and normalize interactions.
    ctx = parse_community(
        FT,
        plankton_args;
        plankton_dynamics=plankton_dynamics,
        biogeochem_dynamics=biogeochem_dynamics,
    )
    mats = normalize_interactions(factory, FT, ctx, interactions)

    # 3) Validate interactions vs required names and sizes.
    required = _required_from_dynamics(plankton_dynamics, biogeochem_dynamics)
    required_mats = Symbol[s for s in required if endswith(String(s), "_matrix")]
    _validate_interaction_matrices(mats, required_mats, ctx.n_total)

    # 4) Build runtime parameters.
    params = create_parameters(factory, FT, ctx, biogeochem_args; interactions_matrices=mats)

    # 5) Build tracer expressions.
    tracers = build_tracer_expressions(plankton_dynamics, biogeochem_dynamics, ctx)

    # 6) Optional sinking velocities.
    if isnothing(sinking_tracers)
        return define_tracer_functions(params, tracers)
    end

    sinking_velocities = setup_velocity_fields(sinking_tracers, grid, open_bottom)
    return define_tracer_functions(params, tracers; sinking_velocities=sinking_velocities)
end
