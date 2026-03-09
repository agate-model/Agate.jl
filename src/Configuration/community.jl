"""Build a plankton community from a base community spec.

This helper updates structural fields (`n` and `diameters`) for one or
more plankton groups, leaving all other group fields intact.

The ordering of groups in the returned `NamedTuple` is the ordering of
`base` (i.e. `keys(base)`).

Keywords
--------
- `n=(; ...)`: optional per-group size-class counts, keyed by group symbol.
- `diameters=(; ...)`: optional per-group diameter specifications, keyed by group symbol.

Examples
--------
```julia
community = build_plankton_community(base;
    n=(Z=10, P=20),
    diameters=(Z=zoo_diameters, P=phyto_diameters),
)
```
"""
function build_plankton_community(
    base::NamedTuple; n::NamedTuple=NamedTuple(), diameters::NamedTuple=NamedTuple()
)
    base_keys = keys(base)

    # Explicit errors for unknown update keys.
    for k in keys(n)
        k in base_keys || throw(ArgumentError("n: unknown group symbol $k"))
    end
    for k in keys(diameters)
        k in base_keys || throw(ArgumentError("diameters: unknown group symbol $k"))
    end

    values_ = ntuple(length(base_keys)) do i
        g = base_keys[i]
        spec = getfield(base, g)

        # Only touch structural fields; preserve everything else.
        new_d = if hasproperty(diameters, g)
            getfield(diameters, g)
        else
            getproperty(spec, :diameters)
        end

        # `n` is optional when diameters are given explicitly as a vector/list.
        # Prefer an explicit `n` override, otherwise infer from explicit diameters,
        # otherwise fall back to the base spec.
        new_n = if hasproperty(n, g)
            getfield(n, g)
        elseif new_d isa AbstractVector
            length(new_d)
        elseif hasproperty(spec, :n)
            getproperty(spec, :n)
        else
            throw(
                ArgumentError(
                    "group $g: missing `n` and diameters are not an explicit vector; provide `n` explicitly",
                ),
            )
        end

        return (; spec..., n=new_n, diameters=new_d)
    end

    return NamedTuple{base_keys}(values_)
end

# -----------------------------------------------------------------------------
# Diameter specifications
# -----------------------------------------------------------------------------

"""Abstract supertype for diameter specifications."""
abstract type AbstractDiameterSpecification end

"""A diameter specification defined by an explicit list of diameters."""
struct DiameterListSpecification{T,VT<:AbstractVector{T}} <: AbstractDiameterSpecification
    diameters::VT
end

"""A diameter specification defined by a range and a splitting method."""
struct DiameterRangeSpecification{T} <: AbstractDiameterSpecification
    min_diameter::T
    max_diameter::T
    splitting::Symbol
end

##############################################################################
# Construction context
##############################################################################

"""Context used during model construction (community parsing, parameter resolution, interaction normalization, derived matrices).

`CommunityContext` exists at **construction time** and is used by factory/variant code to compute concrete parameter defaults. User-facing interaction overrides are data-only and are validated separately. It is distinct from `TendencyContext`, which is used at **tendency time** inside Oceananigans kernels.
"""
struct CommunityContext{FT<:AbstractFloat,VT<:AbstractVector{FT}}
    FT::Type{FT}
    n_total::Int
    diameters::VT
    pfts::Vector{PFTSpecification}
    plankton_symbols::Vector{Symbol}
    group_symbols::Vector{Symbol}
    group_local_index::Vector{Int}
    group_indices::Dict{Symbol,Vector{Int}}

    # Interaction-role axes (used for matrix shapes + tendency kernels)
    consumer_indices::Vector{Int}
    prey_indices::Vector{Int}

    # Default-parameter membership (used only when generating default parameter vectors)
    default_producer_indices::Vector{Int}
    default_consumer_indices::Vector{Int}

    plankton_dynamics::NamedTuple
    biogeochem_dynamics::NamedTuple
end

# -----------------------------------------------------------------------------
# Community parsing
# -----------------------------------------------------------------------------

"""Return a diameter specification for an explicit diameter list."""
diameter_specification(diameters::AbstractVector) = DiameterListSpecification(diameters)

"""Return a diameter specification defined by (min, max, splitting)."""
diameter_specification(spec::Tuple{Any,Any,Symbol}) =
    DiameterRangeSpecification(spec[1], spec[2], spec[3])

"""Return the diameter specification when one is already provided."""
diameter_specification(spec::AbstractDiameterSpecification) = spec

"""Validate `plankton_dynamics` and `community` inputs.

Throws a single `ArgumentError` listing all issues.
"""
function validate_community_inputs(plankton_dynamics, community)
    issues = String[]

    if !(plankton_dynamics isa NamedTuple)
        push!(issues, "plankton_dynamics must be a NamedTuple")
    end
    if !(community isa NamedTuple)
        push!(issues, "community must be a NamedTuple")
    end

    if !isempty(issues)
        throw(ArgumentError(join(issues, "\n")))
    end

    dyn_keys = collect(keys(plankton_dynamics))
    arg_keys = collect(keys(community))

    missing = setdiff(dyn_keys, arg_keys)
    extra = setdiff(arg_keys, dyn_keys)
    !isempty(missing) && push!(issues, "community is missing groups: $(missing)")
    !isempty(extra) && push!(issues, "community has extra groups: $(extra)")

    for k in arg_keys
        if !haskey(community, k)
            continue
        end
        spec = getfield(community, k)

        if !hasproperty(spec, :diameters)
            push!(issues, "group $(k): missing required field `diameters`")
        else
            d = getproperty(spec, :diameters)
            if !(
                d isa AbstractVector ||
                d isa AbstractDiameterSpecification ||
                (d isa Tuple && length(d) == 3)
            )
                push!(issues, "group $(k): invalid `diameters` specification")
            end

            needs_n = !(d isa AbstractVector)

            # For non-explicit diameter specifications (range/splitting or pre-built specs), `n` is required
            # and must be a positive integer. Without this check the downstream community parser will throw
            # a `MethodError` (e.g., when `n === nothing`) instead of a user-facing `ArgumentError`.
            if needs_n
                if !hasproperty(spec, :n)
                    push!(
                        issues,
                        "group $(k): missing required field `n` for non-explicit diameters",
                    )
                else
                    n = getproperty(spec, :n)
                    if !(n isa Integer) || n < 1
                        push!(
                            issues,
                            "group $(k): `n` must be a positive integer for non-explicit diameters",
                        )
                    end
                end
            end

            if d isa AbstractVector && hasproperty(spec, :n)
                n = getproperty(spec, :n)
                if n != length(d)
                    push!(
                        issues,
                        "group $(k): `n` ($(n)) does not match length(diameters) ($(length(d)))",
                    )
                end
            end
        end

        if !hasproperty(spec, :pft)
            push!(issues, "group $(k): missing required field `pft`")
        else
            pft = getproperty(spec, :pft)
            ok = pft isa PFTSpecification || pft isa NamedTuple
            ok || push!(issues, "group $(k): `pft` must be PFTSpecification or NamedTuple")
        end
    end

    if !isempty(issues)
        throw(ArgumentError(join(issues, "\n")))
    end

    return nothing
end

"""Parse and flatten a plankton community into a `CommunityContext`.

Role axes (consumer/prey membership) are defined by `roles`, a `NamedTuple` with fields:
- `consumers`: `nothing` (all groups/classes) or an iterable of group `Symbol`s, indices, or a boolean mask
- `prey`: `nothing` (all groups/classes) or an iterable of group `Symbol`s, indices, or a boolean mask

When `roles` is omitted, both roles default to `nothing` (all classes).
"""
function parse_community(
    factory::AbstractBGCFactory,
    ::Type{FT},
    community::NamedTuple;
    plankton_dynamics::NamedTuple=NamedTuple(),
    biogeochem_dynamics::NamedTuple=NamedTuple(),
    interaction_roles=nothing,
    default_parameter_roles=nothing,
) where {FT<:AbstractFloat}
    if !isnothing(interaction_roles)
        (
            hasproperty(interaction_roles, :consumers) &&
            hasproperty(interaction_roles, :prey)
        ) || throw(
            ArgumentError(
                "interaction_roles must have fields :consumers and :prey (each may be `nothing`, group Symbols, indices, or boolean masks).",
            ),
        )
    end

    if !isnothing(default_parameter_roles)
        (
            hasproperty(default_parameter_roles, :producers) &&
            hasproperty(default_parameter_roles, :consumers)
        ) || throw(
            ArgumentError(
                "default_parameter_roles must have fields :producers and :consumers (each may be `nothing`, group Symbols, indices, or boolean masks).",
            ),
        )
    end
    # Canonical group ordering is the (explicit, stable) ordering of `community`.
    # This makes ordering decisions visible to the caller and avoids hidden
    # factory-specific group ordering.
    group_order = keys(community)

    # `validate_community_inputs(plankton_dynamics, community)` is responsible for
    # ensuring the dynamics and community group sets match.
    plankton_symbols = Symbol[]
    group_of = Symbol[]
    local_idx = Int[]
    pfts = PFTSpecification[]
    diameters = FT[]

    for g in group_order
        spec = getfield(community, g)
        dspec = diameter_specification(getproperty(spec, :diameters))
        n = if dspec isa DiameterListSpecification
            length(dspec.diameters)
        else
            getproperty(spec, :n)
        end
        ds = param_compute_diameters(FT, n, dspec)
        pft_raw = getproperty(spec, :pft)
        pft = pft_raw isa PFTSpecification ? pft_raw : PFTSpecification(pft_raw)

        for i in 1:n
            push!(plankton_symbols, Symbol(string(g), i))
            push!(group_of, g)
            push!(local_idx, i)
            push!(pfts, pft)
            push!(diameters, ds[i])
        end
    end

    n_total = length(plankton_symbols)

    group_indices = Dict{Symbol,Vector{Int}}()
    for (i, g) in enumerate(group_of)
        push!(get!(group_indices, g, Int[]), i)
    end

    interaction_roles_resolved =
        isnothing(interaction_roles) ? (consumers=nothing, prey=nothing) : interaction_roles
    hasproperty(interaction_roles_resolved, :consumers) ||
        throw(ArgumentError("interaction_roles must define :consumers"))
    hasproperty(interaction_roles_resolved, :prey) ||
        throw(ArgumentError("interaction_roles must define :prey"))

    # Default-parameter membership is controlled separately from interaction roles.
    # When `default_parameter_roles` is omitted, we match the interaction axes:
    # - producers := interaction_roles.prey
    # - consumers := interaction_roles.consumers
    default_parameter_roles_resolved = if isnothing(default_parameter_roles)
        (
            producers=getproperty(interaction_roles_resolved, :prey),
            consumers=getproperty(interaction_roles_resolved, :consumers),
        )
    else
        default_parameter_roles
    end
    hasproperty(default_parameter_roles_resolved, :producers) ||
        throw(ArgumentError("default_parameter_roles must define :producers"))
    hasproperty(default_parameter_roles_resolved, :consumers) ||
        throw(ArgumentError("default_parameter_roles must define :consumers"))

    function _indices_for_role(role, role_name::Symbol)
        if role === nothing
            return collect(1:n_total)
        elseif role isa AbstractVector{Bool}
            length(role) == n_total ||
                throw(ArgumentError("$role_name mask must have length $n_total"))
            return findall(role)
        elseif role isa AbstractVector{Int}
            idx = collect(role)
            all(1 .<= idx .<= n_total) ||
                throw(ArgumentError("$role_name indices must be in 1:$n_total"))
            return idx
        elseif role isa Tuple || role isa AbstractVector{Symbol}
            # Treat symbol collections as *membership* (not ordering): preserve the
            # canonical group ordering from `keys(community)` to avoid
            # surprising axis re-ordering when users pass e.g. `(:P, :Z)`.
            requested = Set{Symbol}(role)

            # Validate upfront to give a clean error message.
            for g in requested
                haskey(group_indices, g) ||
                    throw(ArgumentError("Unknown group symbol $g in $role_name roles"))
            end

            idx = Int[]
            for g in group_order
                g ∈ requested || continue
                append!(idx, group_indices[g])
            end
            return idx
        else
            throw(
                ArgumentError(
                    "$role_name roles must be nothing, a Bool mask, an Int index vector, or a collection of group Symbols",
                ),
            )
        end
    end

    consumer_indices = _indices_for_role(
        getproperty(interaction_roles_resolved, :consumers), :consumers
    )
    prey_indices = _indices_for_role(getproperty(interaction_roles_resolved, :prey), :prey)

    default_producer_indices = _indices_for_role(
        getproperty(default_parameter_roles_resolved, :producers), :default_producers
    )
    default_consumer_indices = _indices_for_role(
        getproperty(default_parameter_roles_resolved, :consumers), :default_consumers
    )

    community_context = CommunityContext{FT,typeof(diameters)}(
        FT,
        n_total,
        diameters,
        pfts,
        plankton_symbols,
        group_of,
        local_idx,
        group_indices,
        consumer_indices,
        prey_indices,
        default_producer_indices,
        default_consumer_indices,
        plankton_dynamics,
        biogeochem_dynamics,
    )

    return community_context
end

##############################################################################
# Parameter Utils
###############################################################################

@inline function param_check_length(name::Symbol, expected::Int, got::Int)
    if expected != got
        throw(ArgumentError("$(name) must have length $(expected) but has length $(got)"))
    end
    return nothing
end

function param_compute_diameters(
    ::Type{FT}, n::Int, spec::DiameterRangeSpecification
) where {FT<:AbstractFloat}
    min_d = FT(spec.min_diameter)
    max_d = FT(spec.max_diameter)

    if n == 1
        return FT[min_d]
    end

    diameters = Vector{FT}(undef, n)

    if spec.splitting === :log_splitting
        log_min = log(min_d)
        log_max = log(max_d)
        step = (log_max - log_min) / FT(n - 1)
        @inbounds for i in 1:n
            diameters[i] = exp(log_min + FT(i - 1) * step)
        end
    elseif spec.splitting === :linear_splitting
        step = (max_d - min_d) / FT(n - 1)
        @inbounds for i in 1:n
            diameters[i] = min_d + FT(i - 1) * step
        end
    else
        throw(ArgumentError("Unsupported splitting method: $(spec.splitting)"))
    end

    return diameters
end

function param_compute_diameters(
    ::Type{FT}, n::Int, spec::DiameterListSpecification
) where {FT<:AbstractFloat}
    param_check_length(:diameters, n, length(spec.diameters))
    diameters = Vector{FT}(undef, n)
    @inbounds for i in 1:n
        diameters[i] = FT(spec.diameters[i])
    end
    return diameters
end
