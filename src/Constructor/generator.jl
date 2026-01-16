# -----------------------------------------------------------------------------
# Construction-time code generation utilities
# -----------------------------------------------------------------------------

using Adapt

using Agate.Library.Mortality
using Agate.Library.Nutrients
using Agate.Library.Photosynthesis
using Agate.Library.Predation
using Agate.Library.Remineralization
using Agate.Equations: Equation, expr
using OceanBioME
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using Oceananigans.Fields: ZeroField

import Oceananigans: time_step!

import Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity,
    required_biogeochemical_auxiliary_fields,
    required_biogeochemical_tracers

export build_tracer_expressions
export expression_check
export add_bgc_methods!, create_bgc_struct, define_tracer_functions

"""Build tracer tendency expressions with deterministic tracer ordering.

This is a small convenience helper used by the constructor and tests. All
provided dynamics builders must return `Agate.Equations.Equation`.
"""
function build_tracer_expressions(plankton_dynamics::NamedTuple, biogeochem_dynamics::NamedTuple, ctx)
    plankton_syms = ctx.plankton_symbols

    tracer_names = Symbol[collect(keys(biogeochem_dynamics))...]
    append!(tracer_names, plankton_syms)

    tracer_exprs = Expr[]
    for (tr, f) in pairs(biogeochem_dynamics)
        eq = f(plankton_syms)
        eq isa Equation || throw(ArgumentError("biogeochem dynamics $(nameof(f)) must return Equation"))
        push!(tracer_exprs, expr(eq))
    end

    # Preserve the ordering of `plankton_args` via `ctx`.
    for idx in eachindex(plankton_syms)
        g = ctx.group_symbols[idx]
        f = getfield(plankton_dynamics, g)
        tracer_sym = plankton_syms[idx]
        eq = f(plankton_syms, tracer_sym, idx)
        eq isa Equation || throw(ArgumentError("plankton dynamics $(nameof(f)) must return Equation"))
        push!(tracer_exprs, expr(eq))
    end

    return NamedTuple{Tuple(tracer_names)}(Tuple(tracer_exprs))
end

# -----------------------------------------------------------------------------
# Expression validation
# -----------------------------------------------------------------------------

"""Return all `Symbol`s referenced by an expression tree.

Supports `Expr`, `Symbol`, and literal roots (numbers). This is used at
construction time to determine which runtime parameters need to be bound into
generated tracer methods.
"""
function parse_expression(f_expr)
    symbols = Symbol[]
    stack = Any[f_expr]

    while !isempty(stack)
        x = pop!(stack)
        if x isa Expr
            append!(stack, x.args)
        elseif x isa Symbol
            push!(symbols, x)
        end
    end

    return symbols
end

"""
    expression_check(allowed_symbols, f_expr; module_name=Constructor)

Validate that all Symbols referenced in `f_expr` are either:
- present in `allowed_symbols`, or
- defined in `module_name`.

Throws `UndefVarError` when an undefined Symbol is found.
"""
function expression_check(allowed_symbols, f_expr; module_name=Constructor)
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

    bgc_type = create_bgc_struct(model_name, parameters; sinking_velocities=sinking_velocities)

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
            struct $struct_name{PT} <: AbstractContinuousFormBiogeochemistry
                parameters::PT
            end
            $struct_name
        end

        T = eval(type_expr)
        eval(:(Adapt.@adapt_structure $struct_name))
        return T
    end

    type_expr = quote
        struct $struct_name{PT,W} <: AbstractContinuousFormBiogeochemistry
            parameters::PT
            sinking_velocities::W
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
    # `Agate.Models.compute_model_parameters` returns a wrapper with a `.data` payload.
    # But several tests (and potential user code) pass a plain struct without
    # a `.data` field.
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
  OceanBioME already provides those; overriding them is version-fragile.
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
    parameter_keys = collect(propertynames(_parameter_view(parameters)))

    for (tracer_name, tracer_expression) in pairs(tracers)
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
            # Keep the reference unqualified: the generated method body is `eval`'d
            # into `Agate.Constructor`, and relying on parent-module bindings is
            # brittle.
            method_vars[i] = :($key = _parameter_view(bgc.parameters).$key)
        end

        allowed_symbols = (all_state_vars..., used_params...)
        try
            expression_check(allowed_symbols, tracer_expression)
        catch e
            if e isa UndefVarError
                sym = getfield(e, :var)
                throw(ArgumentError(
                    "Tracer :$(tracer_name) expression references undefined symbol :$(sym). " *
                    "Provide it as a parameter, a state variable, or define it in the helper_functions module.",
                ))
            end
            rethrow()
        end

        tracer_method = quote
            function (bgc::$(wrapper))(::Val{$(QuoteNode(tracer_name))}, $(all_state_vars...))
                $(method_vars...)
                return $(tracer_expression)
            end
        end

        try
            eval(tracer_method)
        catch e
            throw(ArgumentError(
                "Failed to define tracer :$(tracer_name) method. Underlying error: " * sprint(showerror, e),
            ))
        end
    end

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
