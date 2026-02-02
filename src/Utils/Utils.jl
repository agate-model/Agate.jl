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
export default_roles

# Parameter metadata
export ParameterSpec
export parameter_directory
export parameter_spec

# Interactions API
export InteractionContext
export axis_indices
export normalize_interactions
export GroupBlockMatrix
export InteractionMatrices
export sum_over

export KernelBundle
export MatrixFn
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

For interaction parameters declared as matrices, `normalize_interactions` accepts several
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
export validate_plankton_inputs
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

"""Return default role membership for a factory.

Factories may override this to define which groups participate as consumers and prey
in interaction matrices and predation sums.

Return a `NamedTuple` with fields:
- `consumers`: `nothing` (all groups) or an iterable of group `Symbol`s
- `prey`: `nothing` (all groups) or an iterable of group `Symbol`s

These roles define the interaction axes. Overlap is allowed.
"""
default_roles(::AbstractBGCFactory) = (consumers=nothing, prey=nothing)

"""Sum `f(i)` for `i` in `itr`, starting from `init`.

Designed for use with Julia's `do`-block syntax, e.g.

```julia
sum_over(n, zero(T)) do i
    ...
end
```

`itr` may be an integer `n` (interpreted as `1:n`), a range, or `eachindex(x)`.
"""
@inline function sum_over(f, itr, init)
    acc = init
    @inbounds for i in itr
        acc += f(i)
    end
    return acc
end

@inline sum_over(f, n::Int, init) = sum_over(f, Base.OneTo(n), init)

include("ParameterDirectory.jl")

include("KernelBundle.jl")

include("ClassRefs.jl")
include("TracerAccessors.jl")

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

"""Context passed to constructor callbacks (e.g. `interactions(ctx)`)."""
struct InteractionContext{FT<:AbstractFloat,VT<:AbstractVector{FT}}
    FT::Type{FT}
    n_total::Int
    diameters::VT
    pfts::Vector{PFTSpecification}
    plankton_symbols::Vector{Symbol}
    group_symbols::Vector{Symbol}
    group_local_index::Vector{Int}
    group_indices::Dict{Symbol,Vector{Int}}
    consumer_indices::Vector{Int}
    prey_indices::Vector{Int}
    plankton_dynamics::NamedTuple
    biogeochem_dynamics::NamedTuple
end
include("InteractionMatrices.jl")
include("DerivedMatrices.jl")

"""Normalize `interactions` into a `NamedTuple` of parameter overrides.

`interactions` may be:

- `nothing` (no overrides)
- a `NamedTuple` of updates

Within the `NamedTuple`, values may be concrete objects (e.g. matrices) *or*
provider functions. Provider functions are evaluated against an `InteractionContext`
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

When a rectangular matrix is provided for a role-aware parameter (one that
declares `axes` in `parameter_directory(factory)`), it is kept rectangular.
During construction, `finalize_interaction_parameters` will populate the
canonical `parameters.interactions` container and keep the matrix parameter
itself in its axis-local rectangular form.

Shape validation is driven by `parameter_directory(factory)`.
Final key validation and full parameter shape checks occur during model
construction.
"""
function normalize_interactions(
    factory::AbstractBGCFactory,
    ctx::InteractionContext{FT},
    interactions::Union{Nothing,NamedTuple},
) where {FT}
    interactions === nothing && return (;)

    resolved = Pair{Symbol,Any}[]
    for (key, value) in pairs(interactions)
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
        push!(resolved, key => resolved_value)
    end

    return (; resolved...)
end

@inline function _call_interaction_provider(
    f::Function, ctx::InteractionContext, key::Symbol
)
    applicable(f, ctx) && return f(ctx)

    throw(ArgumentError("interaction override '$key' provider must be callable as f(ctx)"))
end

"""Return the global plankton indices for an interaction axis.

Axes may be:

- `:consumer` (role-defined consumer axis)
- `:prey` (role-defined prey axis)
- any existing group `Symbol` present in `ctx.group_indices`
"""
@inline axis_indices(ctx::InteractionContext, axis::Symbol) = if axis === :consumer
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
    ctx::InteractionContext, spec::ParameterSpec, key::Symbol, value
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

        # Explicit wrapper forces group-block expansion (then slice to axes).
        if value isa GroupBlockMatrix
            full = expand_group_block_matrix(ctx, value.B)
            return full[row_indices, col_indices]
        end

        value isa AbstractMatrix || return value

        if size(value) == (nr, nc)
            return value
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
            return value[row_indices, col_indices]
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

    size(value) == (n_total, n_total) && return value

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
function expand_group_block_matrix(ctx::InteractionContext, B::AbstractMatrix)
    groups = unique(ctx.group_symbols)
    ng = length(groups)
    size(B) == (ng, ng) || throw(
        ArgumentError("Expected a $(ng)x$(ng) group-block matrix (got size $(size(B)))."),
    )

    group_to_index = Dict{Symbol,Int}(g => i for (i, g) in pairs(groups))
    group_idx = map(g -> group_to_index[g], ctx.group_symbols)

    return B[group_idx, group_idx]
end

@inline function _expand_group_block_axis_matrix(
    ctx::InteractionContext,
    B::AbstractMatrix,
    row_indices::AbstractVector{Int},
    col_indices::AbstractVector{Int},
    row_groups::AbstractVector{Symbol},
    col_groups::AbstractVector{Symbol},
)
    row_map = Dict{Symbol,Int}(g => i for (i, g) in pairs(row_groups))
    col_map = Dict{Symbol,Int}(g => i for (i, g) in pairs(col_groups))

    nr = length(row_indices)
    nc = length(col_indices)
    R = similar(B, nr, nc)

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
    ctx::InteractionContext, spec::ParameterSpec, key::Symbol, value
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
# Plankton community parsing
# -----------------------------------------------------------------------------

"""Return a diameter specification for an explicit diameter list."""
diameter_specification(diameters::AbstractVector) = DiameterListSpecification(diameters)

"""Return a diameter specification defined by (min, max, splitting)."""
diameter_specification(spec::Tuple{Any,Any,Symbol}) =
    DiameterRangeSpecification(spec[1], spec[2], spec[3])

"""Return the diameter specification when one is already provided."""
diameter_specification(spec::AbstractDiameterSpecification) = spec

"""Validate `plankton_dynamics` and `plankton_args` inputs.

Throws a single `ArgumentError` listing all issues.
"""
function validate_plankton_inputs(plankton_dynamics, plankton_args)
    issues = String[]

    if !(plankton_dynamics isa NamedTuple)
        push!(issues, "plankton_dynamics must be a NamedTuple")
    end
    if !(plankton_args isa NamedTuple)
        push!(issues, "plankton_args must be a NamedTuple")
    end

    if !isempty(issues)
        throw(ArgumentError(join(issues, "\n")))
    end

    dyn_keys = collect(keys(plankton_dynamics))
    arg_keys = collect(keys(plankton_args))

    missing = setdiff(dyn_keys, arg_keys)
    extra = setdiff(arg_keys, dyn_keys)
    !isempty(missing) && push!(issues, "plankton_args is missing groups: $(missing)")
    !isempty(extra) && push!(issues, "plankton_args has extra groups: $(extra)")

    for k in arg_keys
        if !haskey(plankton_args, k)
            continue
        end
        spec = getfield(plankton_args, k)

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

        if !(hasproperty(spec, :args) || hasproperty(spec, :pft))
            push!(issues, "group $(k): must provide `args` or `pft`")
        else
            pft =
                hasproperty(spec, :pft) ? getproperty(spec, :pft) : getproperty(spec, :args)
            ok = pft isa PFTSpecification || pft isa NamedTuple
            ok || push!(
                issues,
                "group $(k): `args`/`pft` must be PFTSpecification or NamedTuple",
            )
        end
    end

    if !isempty(issues)
        throw(ArgumentError(join(issues, "\n")))
    end

    return nothing
end

"""Parse and flatten a plankton community into a construction context.

Role axes (consumer/prey membership) are defined by `roles`, a `NamedTuple` with fields:
- `consumers`: `nothing` (all groups/classes) or an iterable of group `Symbol`s, indices, or a boolean mask
- `prey`: `nothing` (all groups/classes) or an iterable of group `Symbol`s, indices, or a boolean mask

When `roles` is omitted, `default_roles(factory)` is used.
"""
function parse_community(
    factory::AbstractBGCFactory,
    ::Type{FT},
    plankton_args::NamedTuple;
    plankton_dynamics::NamedTuple=NamedTuple(),
    biogeochem_dynamics::NamedTuple=NamedTuple(),
    roles=nothing,
) where {FT<:AbstractFloat}
    group_symbols = collect(keys(plankton_args))
    plankton_symbols = Symbol[]
    group_of = Symbol[]
    local_idx = Int[]
    pfts = PFTSpecification[]
    diameters = FT[]

    for g in group_symbols
        spec = getfield(plankton_args, g)
        dspec = diameter_specification(getproperty(spec, :diameters))
        n = if dspec isa DiameterListSpecification
            length(dspec.diameters)
        else
            getproperty(spec, :n)
        end
        ds = param_compute_diameters(FT, n, dspec)
        pft_raw =
            hasproperty(spec, :pft) ? getproperty(spec, :pft) : getproperty(spec, :args)
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

    roles_resolved = isnothing(roles) ? default_roles(factory) : roles
    hasproperty(roles_resolved, :consumers) ||
        throw(ArgumentError("roles must define :consumers"))
    hasproperty(roles_resolved, :prey) || throw(ArgumentError("roles must define :prey"))

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
            idx = Int[]
            for g in role
                haskey(group_indices, g) ||
                    throw(ArgumentError("Unknown group symbol $g in $role_name roles"))
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

    consumer_indices = _indices_for_role(getproperty(roles_resolved, :consumers), :consumers)
    prey_indices = _indices_for_role(getproperty(roles_resolved, :prey), :prey)

    ctx = InteractionContext{FT,typeof(diameters)}(
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
