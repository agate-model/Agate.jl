"""
Utilities to construct Oceananigans biogeochemistry types.

Agate constructs `AbstractContinuousFormBiogeochemistry` subtypes from:
- a runtime parameter struct stored as `bgc.parameters`, and
- tracer tendency expressions (`Expr`) keyed by tracer name.

The generated types are compatible with Oceananigans and OceanBioME.
Runtime parameter structs can be adapted between CPU and GPU with `Adapt.jl`.
"""
module Utils

using Adapt

include("Specifications.jl")

using .Specifications: PFTSpecification, pft_get, pft_has, cast_pft,
    BiogeochemistrySpecification, cast_spec, ModelSpecification

using Agate.Library.Allometry: allometric_scaling_power

using ..Library.ExprUtils: sum_expr

using Agate.Library.Mortality
using Agate.Library.Nutrients
using Agate.Library.Photosynthesis
using Agate.Library.Predation
using Agate.Library.Remineralization

using OceanBioME
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using Oceananigans.Fields: ZeroField

import Oceananigans: time_step!

import Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity,
    required_biogeochemical_auxiliary_fields,
    required_biogeochemical_tracers

export AbstractDiameterSpecification
export DiameterListSpecification
export DiameterRangeSpecification

# Model-agnostic construction/runtime containers
export AbstractBGCFactory

# Option 1 interactions API
export AbstractInteractions
export InteractionMatrices
export InteractionDynamics
export InteractionContext
export normalize_interactions

# Dynamics parameter registry
export @register_dynamics
export required_parameters

# Community parsing
export validate_plankton_inputs
export parse_community
export build_tracer_expressions

# Model construction helpers
export create_parameters

export add_bgc_methods!, create_bgc_struct, define_tracer_functions, expression_check

export param_check_length
export param_check_square_matrix
export param_cast_matrix
export box_model_mass_balance
export box_model_budget
export sum_expr
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

# -----------------------------------------------------------------------------
# Dynamics parameter registry
# -----------------------------------------------------------------------------

const _DYNAMICS_REQUIRED = IdDict{Any, Vector{Symbol}}()
const _UTILS_MODULE = @__MODULE__

function _register_dynamics!(f, symbols)
    req = Symbol[]
    for s in symbols
        push!(req, s)
    end
    _DYNAMICS_REQUIRED[f] = req
    return nothing
end

"""Register required parameter symbols for a dynamics function."""
macro register_dynamics(f, syms)
    return :($(QuoteNode(_UTILS_MODULE))._register_dynamics!($(esc(f)), $(esc(syms))))
end

"""Return required parameter symbols for `f` (empty if not registered)."""
required_parameters(f) = get(_DYNAMICS_REQUIRED, f, Symbol[])

# -----------------------------------------------------------------------------
# Option 1 interactions API
# -----------------------------------------------------------------------------

abstract type AbstractInteractions end

"""Explicit interaction matrices container."""
struct InteractionMatrices{NT<:NamedTuple} <: AbstractInteractions
    matrices::NT
end

InteractionMatrices(; kwargs...) = InteractionMatrices((; kwargs...))

"""Construction-time interaction builder container."""
struct InteractionDynamics{F,A} <: AbstractInteractions
    f::F
    args::A
end

InteractionDynamics(f; args=NamedTuple()) = InteractionDynamics(f, args)

"""Context passed to interaction builders."""
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

"""Normalize `interactions` into a `NamedTuple` of matrices."""
function normalize_interactions(factory::AbstractBGCFactory, ::Type{FT}, ctx, interactions) where {FT<:AbstractFloat}
    interactions === nothing && (interactions = default_interactions(factory))
    interactions === nothing && return NamedTuple()

    mats = if interactions isa InteractionMatrices
        interactions.matrices
    elseif interactions isa InteractionDynamics
        interactions.f(ctx, interactions.args)
    elseif interactions isa NamedTuple
        interactions
    elseif interactions isa Function
        interactions(ctx)
    else
        throw(ArgumentError("Unsupported interactions type: $(typeof(interactions))"))
    end

    # Cast all matrices to FT.
    return map(M -> param_cast_matrix(FT, M), mats)
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
        pft = cast_pft(FT, pft)

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

"""Build tracer tendency expressions with deterministic tracer ordering."""
function build_tracer_expressions(plankton_dynamics::NamedTuple, biogeochem_dynamics::NamedTuple, ctx)
    plankton_syms = ctx.plankton_symbols

    tracer_names = Symbol[collect(keys(biogeochem_dynamics))...]
    append!(tracer_names, plankton_syms)

    tracer_exprs = Expr[]
    for (tr, f) in pairs(biogeochem_dynamics)
        push!(tracer_exprs, f(plankton_syms))
    end

    # Preserve the ordering of `plankton_args` via `ctx`.
    for idx in eachindex(plankton_syms)
        g = ctx.group_symbols[idx]
        f = getfield(plankton_dynamics, g)
        tracer_sym = plankton_syms[idx]
        push!(tracer_exprs, f(plankton_syms, tracer_sym, idx))
    end

    return NamedTuple{Tuple(tracer_names)}(Tuple(tracer_exprs))
end

# -----------------------------------------------------------------------------
# Parameter construction
# -----------------------------------------------------------------------------

const _SAFE_PFT_DEFAULTS = Set([
    :maximum_growth_rate,
    :nutrient_half_saturation,
    :alpha,
    :photosynthetic_slope,
    :chlorophyll_to_carbon_ratio,
    :maximum_predation_rate,
    :holling_half_saturation,
    :quadratic_mortality,
    :half_saturation_DIN,
    :half_saturation_PO4,
])

@inline function _pft_can_supply(pft::PFTSpecification, key::Symbol)
    pft_has(pft, key) && return true
    a = Symbol(string(key), "_a")
    b = Symbol(string(key), "_b")
    return pft_has(pft, a) && pft_has(pft, b)
end

@inline function _pft_value(::Type{FT}, pft::PFTSpecification, diameter, key::Symbol) where {FT<:AbstractFloat}
    if pft_has(pft, key)
        v = pft_get(pft, key)
        return v isa Bool ? v : FT(v)
    end
    a = Symbol(string(key), "_a")
    b = Symbol(string(key), "_b")
    if pft_has(pft, a) && pft_has(pft, b)
        return allometric_scaling_power(FT(pft_get(pft, a)), FT(pft_get(pft, b)), FT(diameter))
    end
    return nothing
end

"""Create model-agnostic runtime parameters from dynamics requirements, interactions, and specs."""
function create_parameters(
    factory::AbstractBGCFactory,
    ::Type{FT},
    ctx::InteractionContext,
    biogeochem_args;
    interactions_matrices::NamedTuple=NamedTuple(),
) where {FT<:AbstractFloat}
    plankton_dynamics = ctx.plankton_dynamics
    biogeochem_dynamics = ctx.biogeochem_dynamics
    spec = biogeochem_args isa BiogeochemistrySpecification ? cast_spec(FT, biogeochem_args) : cast_spec(FT, BiogeochemistrySpecification(biogeochem_args))
    spec_keys = collect(keys(spec.data))

    required = Symbol[]
    for f in values(plankton_dynamics)
        append!(required, required_parameters(f))
    end
    for f in values(biogeochem_dynamics)
        append!(required, required_parameters(f))
    end
    required = unique(required)

    # Split requirements
    matrix_keys = Symbol[s for s in required if endswith(String(s), "_matrix")]
    scalar_keys = Symbol[s for s in required if s in spec_keys]
    vector_keys = Symbol[s for s in required if !(s in matrix_keys) && !(s in scalar_keys)]

    # Always include diameters
    :diameters in vector_keys || push!(vector_keys, :diameters)

    # Validate required pft-provided parameters for plankton dynamics
    issues = String[]
    for (g, f) in pairs(plankton_dynamics)
        req = required_parameters(f)
        for key in req
            if key in scalar_keys || key in matrix_keys || key == :diameters
                continue
            end
            # This group can default some parameters safely.
            g_pft = ctx.pfts[findfirst(==(g), ctx.group_symbols)]
            if !_pft_can_supply(g_pft, key) && !(key in _SAFE_PFT_DEFAULTS)
                push!(issues, "group $(g) dynamics $(nameof(f)) is missing parameter $(key)")
            end
        end
    end
    if !isempty(issues)
        throw(ArgumentError(join(unique(issues), "\n")))
    end

    n = ctx.n_total
    data_pairs = Pair{Symbol,Any}[]

    # Vector parameters
    for key in vector_keys
        if key == :diameters
            push!(data_pairs, :diameters => ctx.diameters)
            continue
        end
        v = zeros(FT, n)
        @inbounds for i in 1:n
            val = _pft_value(FT, ctx.pfts[i], ctx.diameters[i], key)
            if val === nothing
                v[i] = zero(FT)
            else
                v[i] = FT(val)
            end
        end
        push!(data_pairs, key => v)
    end

    # Scalar parameters (include all spec keys, not just required)
    for (k, v) in pairs(spec.data)
        push!(data_pairs, k => v)
    end

    # Interaction matrices
    for (k, M) in pairs(interactions_matrices)
        push!(data_pairs, k => M)
    end

    return ModelSpecification((; data_pairs...))
end

# -----------------------------------------------------------------------------
# Expression validation
# -----------------------------------------------------------------------------

"""Return all Symbols referenced by an expression tree."""
function parse_expression(f_expr)
    symbols = Symbol[]
    expressions = Expr[f_expr]

    for exp in expressions
        for arg in exp.args
            if arg isa Expr
                push!(expressions, arg)
            elseif arg isa Symbol
                push!(symbols, arg)
            end
        end
    end

    return symbols
end

"""
    expression_check(allowed_symbols, f_expr; module_name=Utils)

Validate that all Symbols referenced in `f_expr` are either:
- present in `allowed_symbols`, or
- defined in `module_name`.

Throws `UndefVarError` when an undefined Symbol is found.
"""
function expression_check(allowed_symbols, f_expr; module_name=Utils)
    symbols = parse_expression(f_expr)

    for s in symbols
        if s ∉ allowed_symbols && !isdefined(module_name, s)
            throw(UndefVarError(s))
        end
    end

    return nothing
end

# -----------------------------------------------------------------------------
# Biogeochemistry type construction
# -----------------------------------------------------------------------------

"""
    define_tracer_functions(
        parameters,
        tracers;
        auxiliary_fields=(:PAR,),
        helper_functions=nothing,
        sinking_velocities=nothing,
    )

Create an Oceananigans biogeochemistry model type.

# Arguments
- `parameters`: runtime parameter struct stored as `bgc.parameters`
- `tracers`: `NamedTuple` mapping tracer names (Symbols) to tracer tendency expressions (`Expr`)

# Keywords
- `auxiliary_fields`: tuple of auxiliary field names (Symbols)
- `helper_functions`: optional file path containing helper functions referenced by `tracers`
- `sinking_velocities`: optional `NamedTuple` mapping sinking tracer names to vertical velocity fields
"""
function define_tracer_functions(
    parameters,
    tracers::NamedTuple;
    auxiliary_fields::Tuple=(:PAR,),
    helper_functions=nothing,
    sinking_velocities=nothing,
)
    Base.@nospecialize parameters tracers auxiliary_fields helper_functions sinking_velocities
    model_name = gensym(:AgateBGC)

    bgc_type = create_bgc_struct(
        model_name, parameters; sinking_velocities=sinking_velocities
    )

    add_bgc_methods!(
        bgc_type,
        tracers,
        parameters;
        auxiliary_fields=auxiliary_fields,
        helper_functions=helper_functions,
        sinking_velocities=sinking_velocities,
    )

    return bgc_type
end

"""
    create_bgc_struct(struct_name, parameters; sinking_velocities=nothing)

Create a subtype of `AbstractContinuousFormBiogeochemistry` that stores a runtime
parameter struct (and optional sinking velocities).

The generated type is `Adapt.jl`-compatible.

Returns the wrapper type (a `UnionAll`) so it remains valid after `Adapt.adapt`
changes the concrete type parameters.
"""
function create_bgc_struct(struct_name::Symbol, parameters; sinking_velocities=nothing)
    if isnothing(sinking_velocities)
        type_expr = quote
            Base.@kwdef struct $struct_name{PT} <: AbstractContinuousFormBiogeochemistry
                parameters::PT = $parameters
            end
            $struct_name
        end

        T = eval(type_expr)
        eval(:(Adapt.@adapt_structure $struct_name))
        return T
    end

    type_expr = quote
        Base.@kwdef struct $struct_name{PT,W} <: AbstractContinuousFormBiogeochemistry
            parameters::PT = $parameters
            sinking_velocities::W = $sinking_velocities
        end
        $struct_name
    end

    T = eval(type_expr)
    eval(:(Adapt.@adapt_structure $struct_name))
    return T
end

# -----------------------------------------------------------------------------
# Method attachment 
# -----------------------------------------------------------------------------

@inline function _bgc_wrapper(bgc_type)
    return bgc_type isa UnionAll ? bgc_type : bgc_type.name.wrapper
end

@inline function _parameter_view(parameters)
    # `Agate.Models.create_parameters` returns a wrapper with a `.data` payload.
    # But several tests (and potential user code) pass a plain struct without
    # a `.data` field (e.g. `LVParameters`).
    return Base.hasproperty(parameters, :data) ? getproperty(parameters, :data) : parameters
end

"""
    add_bgc_methods!(
        bgc_type,
        tracers,
        parameters;
        auxiliary_fields=(),
        helper_functions=nothing,
        sinking_velocities=nothing,
    )

Internal implementation: attaches Oceananigans-required methods.

Important:
- Methods are attached to the *wrapper* type so they remain valid after `Adapt.adapt`
  changes type parameters (CPU -> GPU).
- We DO NOT define/override any `OceanBioME.ContinuousBiogeochemistry` trait methods here.
  OceanBioME already provides those; overriding them is version-fragile (e.g., `.model` field
  does not exist in OceanBioME 0.14.0).
"""
function add_bgc_methods!(
    bgc_type,
    tracers::NamedTuple,
    parameters;
    auxiliary_fields::Tuple=(),
    helper_functions=nothing,
    sinking_velocities=nothing,
)
    Base.@nospecialize bgc_type tracers parameters helper_functions sinking_velocities
    if !isnothing(helper_functions)
        include(helper_functions)
    end

    coordinates = (:x, :y, :z, :t)
    tracer_vars = keys(tracers)
    aux_field_vars = auxiliary_fields
    all_state_vars = (coordinates..., tracer_vars..., aux_field_vars...)

    wrapper = _bgc_wrapper(bgc_type)

    # Traits on the inner model type (wrapper dispatch survives Adapt.adapt).
    eval(:(required_biogeochemical_tracers(::$(wrapper)) = $(tracer_vars)))
    eval(:(required_biogeochemical_auxiliary_fields(::$(wrapper)) = $(aux_field_vars)))

    # Bind runtime parameters into local variables to match tracer expression expectations.
    #
    # To reduce compilation pressure (especially for large models like DARWIN), we only
    # materialize *parameters actually referenced* by each tracer expression.
    #
    # Parameters may be either:
    #   - an Agate wrapper with `.data` payload (from `create_parameters`), or
    #   - a plain struct (used in utility tests and potentially user code).
    parameter_keys = collect(propertynames(_parameter_view(parameters)))


    # Define callable tracer tendency methods.
    for (tracer_name, tracer_expression) in pairs(tracers)
        # Determine which parameters this tracer expression references.
        used_symbols = parse_expression(tracer_expression)
        used_params = Symbol[]
        @inbounds for k in parameter_keys
            (k in used_symbols) && push!(used_params, k)
        end

        # Guard against reserved coordinate names appearing as parameters.
        for k in used_params
            if k in coordinates
                throw(ArgumentError("Parameter name $(k) is reserved for coordinates."))
            end
        end

        method_vars = Vector{Expr}(undef, length(used_params))
        @inbounds for (i, key) in enumerate(used_params)
            # NOTE: keep the reference unqualified. This method body is `eval`'d
            # into the `Agate.Utils` module, and the parent module name `Agate`
            # is not guaranteed to be bound in that module. Qualifying with
            # `Agate.Utils` can therefore trigger `UndefVarError: Agate not defined`
            # when the generated tracer method is executed.
            method_vars[i] = :($key = _parameter_view(bgc.parameters).$key)
        end

        allowed_symbols = (all_state_vars..., used_params...)
        expression_check(allowed_symbols, tracer_expression)

        tracer_method = quote
            function (bgc::$(wrapper))(::Val{$(QuoteNode(tracer_name))}, $(all_state_vars...))
                $(method_vars...)
                return $(tracer_expression)
            end
        end

        eval(tracer_method)
    end

    # Optional sinking velocities.

    # Optional sinking velocities.
    if sinking_velocities !== nothing
        sinking_fields = fieldnames(typeof(sinking_velocities))

        # Fallback: no drift unless explicitly provided.
        eval(
            quote
                function biogeochemical_drift_velocity(
                    ::$(wrapper), ::Val{tracer_name}
                ) where {tracer_name}
                    return (u=ZeroField(), v=ZeroField(), w=ZeroField())
                end
            end,
        )

        for tracer_name in sinking_fields
            sinking_method = quote
                function biogeochemical_drift_velocity(
                    bgc::$(wrapper), ::Val{$(QuoteNode(tracer_name))}
                )
                    return (u=ZeroField(), v=ZeroField(), w=bgc.sinking_velocities.$tracer_name)
                end
            end
            eval(sinking_method)
        end
    end

    return bgc_type
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

@inline function param_check_square_matrix(name::Symbol, n::Int, M::AbstractMatrix)
    if size(M, 1) != n || size(M, 2) != n
        throw(ArgumentError("$(name) must have size ($(n), $(n))"))
    end
    return nothing
end

@inline function param_cast_matrix(::Type{FT}, M::AbstractMatrix) where {FT<:AbstractFloat}
    return eltype(M) === FT ? M : FT.(M)
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


# -----------------------------------------------------------------------------
# Namespaces / submodules
# -----------------------------------------------------------------------------

"""Namespace for community parsing and interaction-matrix utilities."""
module EcosystemInteractions
    import ..Utils:
        # diameter specifications / community layout
        AbstractDiameterSpecification,
        DiameterListSpecification,
        DiameterRangeSpecification,
        diameter_specification,
        param_compute_diameters,
        validate_plankton_inputs,
        parse_community,
        # interactions
        AbstractInteractions,
        InteractionMatrices,
        InteractionDynamics,
        InteractionContext,
        normalize_interactions

    export AbstractDiameterSpecification, DiameterListSpecification, DiameterRangeSpecification
    export diameter_specification, param_compute_diameters
    export validate_plankton_inputs, parse_community

    export AbstractInteractions, InteractionMatrices, InteractionDynamics, InteractionContext
    export normalize_interactions
end

"""Namespace for runtime generation of tracer expressions and BGC types."""
module Generator
    import ..Utils:
        @register_dynamics,
        required_parameters,
        build_tracer_expressions,
        create_parameters,
        expression_check,
        create_bgc_struct,
        add_bgc_methods!,
        define_tracer_functions

    export @register_dynamics, required_parameters
    export build_tracer_expressions, create_parameters
    export expression_check
    export create_bgc_struct, add_bgc_methods!, define_tracer_functions
end

"""Namespace for model diagnostics (e.g. box-model mass balance checks)."""
module Diagnostics
    import ..Utils: box_model_budget, box_model_mass_balance
    export box_model_budget, box_model_mass_balance
end


end # module
