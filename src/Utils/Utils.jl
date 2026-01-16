"""
Utilities for configuring and inspecting Agate models.

Contains:
- factory interfaces (`AbstractBGCFactory`)
- community parsing helpers (`parse_community`, diameter specifications)
- lightweight box-model diagnostics (`box_model_budget`, `box_model_mass_balance`)
"""
module Utils

include("Specifications.jl")

using .Specifications: PFTSpecification, pft_get, pft_has, ModelSpecification

import Oceananigans: time_step!

export AbstractDiameterSpecification
export DiameterListSpecification
export DiameterRangeSpecification

# Model-agnostic construction/runtime containers
export AbstractBGCFactory

# Interactions API
export InteractionContext
export normalize_interactions

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
    plankton_dynamics::NamedTuple
    biogeochem_dynamics::NamedTuple
end

"""Normalize `interactions` into a `NamedTuple` of registry updates.

`interactions` may be:

- `nothing` (no overrides)
- a `NamedTuple` of updates (typically matrices, but may include any parameter keys)
- a function `(ctx) -> NamedTuple` producing updates from the parsed community context

Shape/type validation is performed during parameter resolution based on the active
parameter registry.
"""
function normalize_interactions(
    factory::AbstractBGCFactory,
    ::Type{FT},
    ctx,
    interactions,
) where {FT<:AbstractFloat}
    interactions === nothing && (interactions = default_interactions(factory))
    interactions === nothing && return NamedTuple()

    overrides = if interactions isa NamedTuple
        interactions
    elseif interactions isa Function
        interactions(ctx)
    else
        throw(ArgumentError("interactions must be a NamedTuple or a function `(ctx)->NamedTuple`, got $(typeof(interactions))"))
    end

    overrides isa NamedTuple ||
        throw(ArgumentError("interactions must return a NamedTuple, got $(typeof(overrides))"))

    return overrides
end
"""Fallback: no default interactions."""
default_interactions(::AbstractBGCFactory) = nothing

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
            if !(d isa AbstractVector || d isa AbstractDiameterSpecification || (d isa Tuple && length(d) == 3))
                push!(issues, "group $(k): invalid `diameters` specification")
            end

            needs_n = !(d isa AbstractVector)

            # For non-explicit diameter specifications (range/splitting or pre-built specs), `n` is required
            # and must be a positive integer. Without this check the downstream community parser will throw
            # a `MethodError` (e.g., when `n === nothing`) instead of a user-facing `ArgumentError`.
            if needs_n
                if !hasproperty(spec, :n)
                    push!(issues, "group $(k): missing required field `n` for non-explicit diameters")
                else
                    n = getproperty(spec, :n)
                    if !(n isa Integer) || n < 1
                        push!(issues, "group $(k): `n` must be a positive integer for non-explicit diameters")
                    end
                end
            end

            if d isa AbstractVector && hasproperty(spec, :n)
                n = getproperty(spec, :n)
                if n != length(d)
                    push!(issues, "group $(k): `n` ($(n)) does not match length(diameters) ($(length(d)))")
                end
            end
        end

        if !(hasproperty(spec, :args) || hasproperty(spec, :pft))
            push!(issues, "group $(k): must provide `args` or `pft`")
        else
            pft = hasproperty(spec, :pft) ? getproperty(spec, :pft) : getproperty(spec, :args)
            ok = pft isa PFTSpecification || pft isa NamedTuple
            ok || push!(issues, "group $(k): `args`/`pft` must be PFTSpecification or NamedTuple")
        end
    end

    if !isempty(issues)
        throw(ArgumentError(join(issues, "\n")))
    end

    return nothing
end

"""Parse and flatten a plankton community into a construction context."""
function parse_community(
    ::Type{FT},
    plankton_args::NamedTuple;
    plankton_dynamics::NamedTuple=NamedTuple(),
    biogeochem_dynamics::NamedTuple=NamedTuple(),
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
        n = dspec isa DiameterListSpecification ? length(dspec.diameters) : getproperty(spec, :n)
        ds = param_compute_diameters(FT, n, dspec)
        pft_raw = hasproperty(spec, :pft) ? getproperty(spec, :pft) : getproperty(spec, :args)
        pft = pft_raw isa PFTSpecification ? pft_raw : PFTSpecification(pft_raw)
        # Keep specifications as-authored; numeric literals are converted to FT when resolved
        # into runtime parameters.

        for i in 1:n
            push!(plankton_symbols, Symbol(string(g), i))
            push!(group_of, g)
            push!(local_idx, i)
            push!(pfts, pft)
            push!(diameters, ds[i])
        end
    end

    ctx = InteractionContext{FT,typeof(diameters)}(
        FT,
        length(plankton_symbols),
        diameters,
        pfts,
        plankton_symbols,
        group_of,
        local_idx,
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
    box_model,
    budgets::NamedTuple;
    dt,
    nsteps::Integer,
    location::NTuple{3,Int}=(1, 1, 1),
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
