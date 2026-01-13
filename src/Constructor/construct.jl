using OceanBioME: BoxModelGrid, setup_velocity_fields

using Agate.Utils:
    AbstractBGCFactory,
    define_tracer_functions,
    normalize_interactions,
    parse_community,
    validate_plankton_inputs,
    param_cast_matrix

using Agate.Utils.Specifications: BiogeochemistrySpecification, cast_spec, ModelSpecification, pft_get, pft_has

using Agate.Models:
    default_plankton_dynamics,
    default_plankton_args,
    default_biogeochem_dynamics,
    default_biogeochem_args

using Agate.Library.Equations: Equation, expr, requirements, req, merge_requirements
using Agate.Library.Allometry:
    allometric_palatability_unimodal_protection,
    build_palatability_matrix,
    build_assimilation_matrix,
    resolve_param

# -----------------------------------------------------------------------------
# Public constructor
# -----------------------------------------------------------------------------

"""
    construct(factory::AbstractBGCFactory; kwargs...) -> Type

Compile a concrete biogeochemistry type from a factory and optional overrides.

The return value is a *type*; instantiate with `bgc = bgc_type()`.

Dynamics builders **must** return `Agate.Library.Equations.Equation`.
"""
function construct(
    factory::AbstractBGCFactory;
    FT::Type{<:AbstractFloat}=Float64,

    # Interaction matrices (default: built from PFT traits when required by equations)
    palat_matrix=nothing,
    assim_matrix=nothing,

    # Default palatability rule used by the library matrix builder.
    palatability_fn=allometric_palatability_unimodal_protection,

    plankton_dynamics=default_plankton_dynamics(factory),
    plankton_args=default_plankton_args(factory, FT),
    biogeochem_dynamics=default_biogeochem_dynamics(factory),
    biogeochem_args=default_biogeochem_args(factory, FT),

    # Legacy/advanced: allow providing extra interaction matrices.
    interactions=nothing,

    sinking_tracers=nothing,
    grid=BoxModelGrid(),
    open_bottom::Bool=true,
)
    return _construct_impl(
        factory,
        FT,
        palatability_fn,
        palat_matrix,
        assim_matrix,
        plankton_dynamics,
        plankton_args,
        biogeochem_dynamics,
        biogeochem_args,
        interactions,
        sinking_tracers,
        grid,
        open_bottom,
    )
end

# -----------------------------------------------------------------------------
# Internal implementation
# -----------------------------------------------------------------------------

@inline function _push_unique!(v::Vector{Symbol}, items)
    for s in items
        if s ∉ v
            push!(v, s)
        end
    end
    return v
end

function _construct_impl(
    factory::AbstractBGCFactory,
    FT::Type{<:AbstractFloat},
    palatability_fn,
    palat_matrix,
    assim_matrix,
    plankton_dynamics,
    plankton_args,
    biogeochem_dynamics,
    biogeochem_args,
    interactions,
    sinking_tracers,
    grid,
    open_bottom::Bool,
)

    # 1) Fill defaults and validate input shapes.
    validate_plankton_inputs(plankton_dynamics, plankton_args)
    biogeochem_dynamics isa NamedTuple || throw(ArgumentError("biogeochem_dynamics must be a NamedTuple"))

    # 2) Parse community and normalize any explicit interaction matrices.
    ctx = parse_community(
        FT,
        plankton_args;
        plankton_dynamics=plankton_dynamics,
        biogeochem_dynamics=biogeochem_dynamics,
    )
    explicit_mats = normalize_interactions(factory, FT, ctx, interactions)

    # 3) Build equations and collect requirements.
    plankton_syms = ctx.plankton_symbols

    tracer_names = Symbol[collect(keys(biogeochem_dynamics))...]
    append!(tracer_names, plankton_syms)

    tracer_exprs = Expr[]
    merged = req()

    # Per-group parameter requirements (for typo protection).
    group_required = Dict{Symbol, Vector{Symbol}}()

    # Biogeochemical tracers first.
    for (_, f) in pairs(biogeochem_dynamics)
        eq = f(plankton_syms)
        eq isa Equation || throw(ArgumentError("biogeochem dynamics $(nameof(f)) must return Equation"))
        push!(tracer_exprs, expr(eq))
        merged = merge_requirements(merged, requirements(eq))
    end

    # Plankton tracers in parsed community order.
    for idx in eachindex(plankton_syms)
        g = ctx.group_symbols[idx]
        f = getfield(plankton_dynamics, g)
        tr = plankton_syms[idx]

        eq = f(plankton_syms, tr, idx)
        eq isa Equation || throw(ArgumentError("plankton dynamics $(nameof(f)) must return Equation"))
        push!(tracer_exprs, expr(eq))

        r = requirements(eq)
        merged = merge_requirements(merged, r)

        # Track which group-owned parameters are required by this group's own dynamics.
        reqvec = get!(group_required, g) do
            Symbol[]
        end
        _push_unique!(reqvec, r.group_params)
    end

    tracers = NamedTuple{Tuple(tracer_names)}(Tuple(tracer_exprs))

    # 4) Validate required BGC scalars and cast the biogeochemistry spec.
    spec = biogeochem_args isa BiogeochemistrySpecification ? biogeochem_args : BiogeochemistrySpecification(biogeochem_args)
    spec = cast_spec(FT, spec)

    for s in merged.scalars
        hasproperty(spec.data, s) || throw(ArgumentError("Missing required biogeochemical scalar: $(s)"))
    end

    # 5) Validate group-owned required parameters (missing vs explicit `nothing`).
    for (g, keys) in group_required
        rep = findfirst(i -> ctx.group_symbols[i] == g, eachindex(ctx.group_symbols))
        rep === nothing && continue
        pft = ctx.pfts[rep]
        for k in keys
            pft_has(pft, k) || throw(ArgumentError("Group $(g) is missing required parameter: $(k)"))
        end
    end

    # 6) Allocate and fill parameter containers.
    #
    # - For group parameters: missing => error (validated above); explicit `nothing` => zeros.
    # - For community parameters: missing or `nothing` => zeros.

    data_pairs = Pair{Symbol, Any}[]

    # Canonical always-present vector.
    diameters = ctx.diameters
    push!(data_pairs, :diameters => diameters)

    # Biogeochemical scalars (always included; validated if referenced by equations).
    for k in keys(spec.data)
        push!(data_pairs, k => getproperty(spec.data, k))
    end

    # Vector parameters derived entirely from equation requirements.
    vector_keys = Symbol[]
    _push_unique!(vector_keys, merged.group_params)
    _push_unique!(vector_keys, merged.community_params)

    for k in vector_keys
        v = zeros(FT, ctx.n_total)
        @inbounds for i in 1:ctx.n_total
            pft = ctx.pfts[i]
            if pft_has(pft, k)
                val = pft_get(pft, k)
                if val !== nothing
                    v[i] = resolve_param(FT, val, diameters[i])
                end
            end
        end
        push!(data_pairs, k => v)
    end

    # Matrices derived from equation requirements.
    pft_data = map(p -> p.data, ctx.pfts)

    # Apply explicit overrides with precedence:
    # keyword args > `interactions` NamedTuple
    pal_override = palat_matrix
    pal_override === nothing && (pal_override = hasproperty(explicit_mats, :palatability_matrix) ? getproperty(explicit_mats, :palatability_matrix) : nothing)

    assim_override = assim_matrix
    assim_override === nothing && (assim_override = hasproperty(explicit_mats, :assimilation_efficiency_matrix) ? getproperty(explicit_mats, :assimilation_efficiency_matrix) : nothing)

    for m in merged.matrices
        if m == :palatability_matrix
            M = build_palatability_matrix(FT, pft_data, diameters; overrides=pal_override, palatability_fn=palatability_fn)
            push!(data_pairs, m => M)
        elseif m == :assimilation_efficiency_matrix
            M = build_assimilation_matrix(FT, pft_data, diameters; overrides=assim_override)
            push!(data_pairs, m => M)
        else
            hasproperty(explicit_mats, m) || throw(ArgumentError("Missing required interaction matrix: $(m)"))
            push!(data_pairs, m => param_cast_matrix(FT, getproperty(explicit_mats, m)))
        end
    end

    params = ModelSpecification((; data_pairs...))

    # 7) Optional sinking velocities.
    if isnothing(sinking_tracers)
        return define_tracer_functions(params, tracers)
    end

    sinking_velocities = setup_velocity_fields(sinking_tracers, grid, open_bottom)
    return define_tracer_functions(params, tracers; sinking_velocities=sinking_velocities)
end
