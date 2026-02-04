"""
Utilities for configuring and inspecting Agate models.

Contains:
- factory interfaces (`AbstractBGCFactory`)
- community parsing helpers (`parse_community`, diameter specifications)
- lightweight box-model diagnostics (`box_model_budget`, `box_model_mass_balance`)
"""
module Utils

include("Specifications.jl")

using .Specifications: PFTSpecification

import Oceananigans: time_step!

export AbstractDiameterSpecification
export DiameterListSpecification
export DiameterRangeSpecification

# Model-agnostic construction/runtime containers
export AbstractBGCFactory

# Parameter metadata
export ParameterSpec
export parameter_directory
export parameter_spec

# Interactions API
export CommunityContext
export axis_indices
export normalize_interaction_overrides
export GroupBlockMatrix
export InteractionBlocks
export roles_from_groups
export interaction_blocks
export forbid_link!
export set_block!
export scale_block!
export InteractionMatrices

export TendencyContext
export TracerValues
export tendency_views
export MatrixProvider
export derived_matrix_specs
export resolve_derived_matrices

# Class references (group + ordinal)
export ClassRef
export class
export resolve_class
export class_count

# GPU-safe tracer access
export TracerIndex
export Tracers
export build_tracer_index

"""Wrapper to explicitly mark a matrix as a *group-block* interaction override.

For interaction parameters declared as matrices, `normalize_interaction_overrides` accepts several
shapes (full `(n_total, n_total)`, axis-sized matrices when axes are declared, etc.).

When an interaction parameter declares role-aware axes (for example `axes = (:consumer, :prey)`),
plain `(n_groups, n_groups)` matrices are ambiguous: they could be interpreted as an axis-sized
matrix or as a block matrix over all groups.

Wrap a matrix in `GroupBlockMatrix(B)` to force it to be interpreted as a group-block matrix
over *all* groups and expanded to `(n_total, n_total)` via `expand_group_block_matrix`.
"""
struct GroupBlockMatrix{T<:AbstractMatrix}
    B::T
end
export expand_group_block_matrix

# Community parsing
export validate_community_inputs
export parse_community

export param_check_length
export box_model_mass_balance
export box_model_budget
export param_compute_diameters

# -----------------------------------------------------------------------------
# Model-agnostic factories
# -----------------------------------------------------------------------------

"""Abstract supertype for biogeochemical model factories."""
abstract type AbstractBGCFactory end

"""Model families provide implementations via the `Interface` module.

Factories are intentionally minimal: they declare defaults for community structure and
dynamics builders. Group ordering is inferred from the *explicit* ordering of the
`community::NamedTuple` passed to `Constructor.construct_factory`.
"""

include("ParameterDirectory.jl")

include("TendencyContext.jl")

include("ClassRefs.jl")
include("TracerAccessors.jl")

include("InteractionBlocks.jl")

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

"""Context passed to constructor-time provider functions (e.g. `palatability_matrix = (ctx) -> ...`).

`CommunityContext` exists at **construction time** (community parsing, parameter
resolution, interaction normalization, derived matrices). It is distinct from
`TendencyContext`, which is used at **tendency time** inside Oceananigans kernels.
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

    # Parameter-default classification (used by `default_parameters`)
    producer_param_indices::Vector{Int}
    consumer_param_indices::Vector{Int}

    plankton_dynamics::NamedTuple
    biogeochem_dynamics::NamedTuple
end
include("InteractionMatrices.jl")
include("DerivedMatrices.jl")

"""Normalize `interaction_overrides` into a `NamedTuple` of parameter overrides.

`interaction_overrides` may be:

- `nothing` (no overrides)
- a `NamedTuple` of updates

Within the `NamedTuple`, values may be concrete objects (e.g. matrices) *or*
provider functions. Provider functions are evaluated against a `CommunityContext`
and must be callable as:

- `f(ctx)`

For matrix parameters, users may pass either:

- a full `(n_total, n_total)` matrix,
- a group-block `(n_groups, n_groups)` matrix which will be expanded (and, for
  role-aware parameters, sliced to the declared axes) (see `expand_group_block_matrix`). When the parameter
  declares role-aware axes, a plain `(n_groups, n_groups)` matrix is ambiguous;
  wrap it in `GroupBlockMatrix(B)` to force group-block expansion,
- or, when the parameter directory declares role-aware axes,
  a rectangular matrix sized to those axes (for example `(n_consumer, n_prey)`).

For consumer-by-prey matrices with axes `(:consumer, :prey)`, users may also pass
an `InteractionBlocks` object created by `interaction_blocks(roles; init=...)`.

When a rectangular matrix is provided for a role-aware parameter (one that
declares `axes` in `parameter_directory(factory)`), it is kept rectangular.
During construction, `finalize_interaction_parameters` will populate the
canonical `parameters.interactions` container and keep the matrix parameter
itself in its axis-local rectangular form.

Shape validation is driven by `parameter_directory(factory)`.
Final key validation and full parameter shape checks occur during model
construction.
"""
function normalize_interaction_overrides(
    factory::AbstractBGCFactory,
    ctx::CommunityContext{FT},
    interaction_overrides::Union{Nothing,NamedTuple},
) where {FT}
    interaction_overrides === nothing && return (;)

    resolved = ()
    for (key, value) in pairs(interaction_overrides)
        spec = parameter_spec(factory, key)
        spec !== nothing || throw(
            ArgumentError(
                "interaction override '$key' is missing a ParameterSpec in parameter_directory(::$(typeof(factory))).",
            ),
        )

        resolved_value =
            value isa Function ? _call_interaction_provider(value, ctx, key) : value
        resolved_value = _normalize_interaction_value(ctx, spec, key, resolved_value)
        _validate_interaction_override(ctx, spec, key, resolved_value)
        resolved = (resolved..., key => resolved_value)
    end

    return (; resolved...)
end

@inline function _call_interaction_provider(f::Function, ctx::CommunityContext, key::Symbol)
    applicable(f, ctx) && return f(ctx)

    throw(ArgumentError("interaction override '$key' provider must be callable as f(ctx)"))
end

"""Return the global plankton indices for an interaction axis.

Axes may be:

- `:consumer` (role-defined consumer axis)
- `:prey` (role-defined prey axis)
- any existing group `Symbol` present in `ctx.group_indices`
"""
@inline axis_indices(ctx::CommunityContext, axis::Symbol) =
    if axis === :consumer
        ctx.consumer_indices
    elseif axis === :prey
        ctx.prey_indices
    elseif haskey(ctx.group_indices, axis)
        ctx.group_indices[axis]
    else
        throw(
            ArgumentError(
                "Unknown interaction axis '$axis'. Valid axes are :consumer, :prey, or an existing group symbol.",
            ),
        )
    end

@inline function _normalize_interaction_value(
    ctx::CommunityContext, spec::ParameterSpec, key::Symbol, value
)
    spec.shape === :matrix || return value
    n_total = ctx.n_total

    groups = unique(ctx.group_symbols)
    ng = length(groups)

    # Axes-aware matrices are normalized into their axis-local rectangular form.
    if spec.axes !== nothing
        row_axis, col_axis = spec.axes
        row_indices = axis_indices(ctx, row_axis)
        col_indices = axis_indices(ctx, col_axis)

        nr = length(row_indices)
        nc = length(col_indices)

        # Convenience wrapper: explicit group order for a consumer-by-prey block.
        if value isa InteractionBlocks
            spec.axes == (:consumer, :prey) || throw(
                ArgumentError(
                    "InteractionBlocks can only be used for axes (:consumer, :prey); got axes $(spec.axes) for '$key'.",
                ),
            )

            row_groups_axis = unique(ctx.group_symbols[row_indices])
            col_groups_axis = unique(ctx.group_symbols[col_indices])

            Set(row_groups_axis) == Set(value.consumer_groups) || throw(
                ArgumentError(
                    "InteractionBlocks consumer_groups $(value.consumer_groups) do not match consumer axis groups $(Tuple(row_groups_axis)) for '$key'.",
                ),
            )
            Set(col_groups_axis) == Set(value.prey_groups) || throw(
                ArgumentError(
                    "InteractionBlocks prey_groups $(value.prey_groups) do not match prey axis groups $(Tuple(col_groups_axis)) for '$key'.",
                ),
            )

            size(value.B) == (length(value.consumer_groups), length(value.prey_groups)) ||
                throw(
                    ArgumentError(
                        "InteractionBlocks.B must have size $(length(value.consumer_groups))x$(length(value.prey_groups)); got size $(size(value.B)) for '$key'.",
                    ),
                )

            return _expand_group_block_axis_matrix(
                ctx,
                value.B,
                row_indices,
                col_indices,
                value.consumer_groups,
                value.prey_groups,
            )
        end

        # Explicit wrapper forces group-block expansion (then slice to axes).
        if value isa GroupBlockMatrix
            full = expand_group_block_matrix(ctx, value.B)
            return full[row_indices, col_indices]
        end

        value isa AbstractMatrix || return value

        if size(value) == (nr, nc)
            return _convert_matrix_eltype(ctx, value)
        end

        # Axis-local group-block matrix (groups present on each axis, in order of appearance).
        row_groups = unique(ctx.group_symbols[row_indices])
        col_groups = unique(ctx.group_symbols[col_indices])
        ngr = length(row_groups)
        ngc = length(col_groups)

        if size(value) == (ngr, ngc)
            return _expand_group_block_axis_matrix(
                ctx, value, row_indices, col_indices, row_groups, col_groups
            )
        end

        # Full-square override: slice to axes.
        if size(value) == (n_total, n_total)
            return _convert_matrix_eltype(ctx, value[row_indices, col_indices])
        end

        # Ambiguous: a plain (n_groups, n_groups) matrix could be intended as a group-block matrix.
        if size(value) == (ng, ng)
            throw(
                ArgumentError(
                    "interaction override '$key' is ambiguous for axes $(spec.axes): a $(ng)x$(ng) matrix could be either an axis-sized matrix or a group-block matrix. " *
                    "Wrap it as GroupBlockMatrix(B) to force group-block expansion, or pass a $(nr)x$(nc) axis matrix.",
                ),
            )
        end

        throw(
            ArgumentError(
                "interaction override '$key' must be a $(nr)x$(nc) axis matrix, a $(ngr)x$(ngc) axis block, or a $(n_total)x$(n_total) full matrix; got size $(size(value))",
            ),
        )
    end

    # Non-axes matrices normalize to full square matrices.
    if value isa GroupBlockMatrix
        return expand_group_block_matrix(ctx, value.B)
    end

    value isa AbstractMatrix || return value

    value isa AbstractMatrix || return value

    size(value) == (n_total, n_total) && return _convert_matrix_eltype(ctx, value)

    # Group-block override over *all* groups.
    size(value) == (ng, ng) && return expand_group_block_matrix(ctx, value)

    return value
end

"""Expand a group-block matrix to a full `(n_total, n_total)` matrix.

The group ordering is the order of first appearance in `ctx.group_symbols`.

For a block matrix `B` sized `(n_groups, n_groups)`, the expanded matrix `M` is
defined as:

`M[predator, prey] = B[group(predator), group(prey)]`.
"""
function expand_group_block_matrix(ctx::CommunityContext, B::AbstractMatrix)
    groups = unique(ctx.group_symbols)
    ng = length(groups)
    size(B) == (ng, ng) || throw(
        ArgumentError("Expected a $(ng)x$(ng) group-block matrix (got size $(size(B)))."),
    )

    group_to_index = Dict{Symbol,Int}(g => i for (i, g) in pairs(groups))
    group_idx = map(g -> group_to_index[g], ctx.group_symbols)

    return _convert_matrix_eltype(ctx, B[group_idx, group_idx])
end

@inline function _convert_matrix_eltype(ctx::CommunityContext, A::AbstractMatrix)
    eltype(A) === ctx.FT && return A
    R = similar(A, ctx.FT, size(A)...)
    copyto!(R, A)
    return R
end

@inline function _expand_group_block_axis_matrix(
    ctx::CommunityContext,
    B::AbstractMatrix,
    row_indices::AbstractVector{<:Integer},
    col_indices::AbstractVector{<:Integer},
    row_groups,
    col_groups,
)
    row_map = Dict{Symbol,Int}(g => i for (i, g) in pairs(row_groups))
    col_map = Dict{Symbol,Int}(g => i for (i, g) in pairs(col_groups))

    nr = length(row_indices)
    nc = length(col_indices)
    R = similar(B, ctx.FT, nr, nc)

    for (ii, gi) in enumerate(row_indices)
        rg = ctx.group_symbols[gi]
        ri = row_map[rg]
        for (jj, gj) in enumerate(col_indices)
            cg = ctx.group_symbols[gj]
            cj = col_map[cg]
            @inbounds R[ii, jj] = B[ri, cj]
        end
    end

    return R
end

@inline function _validate_interaction_override(
    ctx::CommunityContext, spec::ParameterSpec, key::Symbol, value
)
    if spec.shape === :matrix
        value isa AbstractMatrix || throw(
            ArgumentError(
                "interaction override '$key' must be a matrix; got $(typeof(value))"
            ),
        )

        if spec.axes === nothing
            n_total = ctx.n_total
            size(value) == (n_total, n_total) || throw(
                ArgumentError(
                    "interaction override '$key' must be a $(n_total)x$(n_total) matrix after normalization; got size $(size(value))",
                ),
            )
        else
            row_axis, col_axis = spec.axes
            row_indices = axis_indices(ctx, row_axis)
            col_indices = axis_indices(ctx, col_axis)
            nr = length(row_indices)
            nc = length(col_indices)
            size(value) == (nr, nc) || throw(
                ArgumentError(
                    "interaction override '$key' must be a $(nr)x$(nc) matrix for axes $(spec.axes) after normalization; got size $(size(value))",
                ),
            )
        end

    elseif spec.shape === :vector
        value isa AbstractVector || throw(
            ArgumentError(
                "interaction override '$key' must be a vector; got $(typeof(value))"
            ),
        )

        n_total = ctx.n_total
        length(value) == n_total || throw(
            ArgumentError(
                "interaction override '$key' must have length $n_total (got $(length(value)))",
            ),
        )
    end

    return nothing
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
    roles=nothing,
    parameter_groups=nothing,
) where {FT<:AbstractFloat}
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

    roles_resolved = isnothing(roles) ? (consumers=nothing, prey=nothing) : roles
    hasproperty(roles_resolved, :consumers) ||
        throw(ArgumentError("roles must define :consumers"))
    hasproperty(roles_resolved, :prey) || throw(ArgumentError("roles must define :prey"))

    # Parameter-default group membership is controlled separately from interaction roles.
    # When `parameter_groups` is omitted, we default to matching the role axes:
    # - producers := roles.prey
    # - consumers := roles.consumers
    parameter_groups_resolved = if isnothing(parameter_groups)
        (
            producers=getproperty(roles_resolved, :prey),
            consumers=getproperty(roles_resolved, :consumers),
        )
    else
        parameter_groups
    end
    hasproperty(parameter_groups_resolved, :producers) ||
        throw(ArgumentError("parameter_groups must define :producers"))
    hasproperty(parameter_groups_resolved, :consumers) ||
        throw(ArgumentError("parameter_groups must define :consumers"))

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
        getproperty(roles_resolved, :consumers), :consumers
    )
    prey_indices = _indices_for_role(getproperty(roles_resolved, :prey), :prey)

    producer_param_indices = _indices_for_role(
        getproperty(parameter_groups_resolved, :producers), :producers
    )
    consumer_param_indices = _indices_for_role(
        getproperty(parameter_groups_resolved, :consumers), :parameter_consumers
    )

    ctx = CommunityContext{FT,typeof(diameters)}(
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
        producer_param_indices,
        consumer_param_indices,
        plankton_dynamics,
        biogeochem_dynamics,
    )

    return ctx
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

# -----------------------------------------------------------------------------
# Box model mass balance utilities
# -----------------------------------------------------------------------------

"""
    box_model_budget(box_model, terms; location=(1, 1, 1))

Compute a weighted sum of tracer values in `box_model` at the given grid `location`.

`terms` may be either:
- an `AbstractVector` of pairs `tracer_symbol => weight`, e.g. `[:N => 1, :P1 => 1]`, or
- a `NamedTuple` of weights, e.g. `(N=1, P1=1)`.

The function assumes tracer fields are accessible as `box_model.fields.<tracer>` and store
their data in `.data`.

This is intended as a small helper for mass/budget diagnostics and tests.
"""
function box_model_budget(box_model, terms; location::NTuple{3,Int}=(1, 1, 1))
    i, j, k = location
    pairs_iter = terms isa NamedTuple ? pairs(terms) : terms

    s = 0.0
    for (tracer, weight) in pairs_iter
        fld = getproperty(box_model.fields, tracer)
        s += float(weight) * float(fld.data[i, j, k])
    end
    return s
end

"""
    box_model_mass_balance(box_model, budgets; dt, nsteps, location=(1, 1, 1))

Advance `box_model` forward `nsteps` with timestep `dt` and return a NamedTuple:

`(initial=..., final=..., drift=..., relative_drift=...)`

where each of `initial`, `final`, `drift`, and `relative_drift` is a NamedTuple with the
same keys as `budgets`.

`budgets` is a NamedTuple mapping budget names to a `terms` specification accepted by
`box_model_budget`.

This function does not depend on `Test` and can be used for lightweight model diagnostics.
"""
function box_model_mass_balance(
    box_model, budgets::NamedTuple; dt, nsteps::Integer, location::NTuple{3,Int}=(1, 1, 1)
)
    initial = map(terms -> box_model_budget(box_model, terms; location), budgets)

    for _ in 1:nsteps
        time_step!(box_model, dt)
    end

    final = map(terms -> box_model_budget(box_model, terms; location), budgets)
    drift = map((a, b) -> b - a, initial, final)
    relative_drift = map((a, b) -> a == 0 ? (b - a) : (b - a) / a, initial, final)

    return (; initial, final, drift, relative_drift)
end

end # module
