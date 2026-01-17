# -----------------------------------------------------------------------------
# Construction-time code generation utilities
# -----------------------------------------------------------------------------

using Adapt

using ..Equations: Equation, expr, requirements
using ..Functors: CompiledEquation
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
function build_tracer_expressions(PV, plankton_dynamics::NamedTuple, biogeochem_dynamics::NamedTuple, ctx)
    plankton_syms = ctx.plankton_symbols

    tracer_names = Symbol[collect(keys(biogeochem_dynamics))...]
    append!(tracer_names, plankton_syms)

    tracer_exprs = Expr[]
    for (tr, f) in pairs(biogeochem_dynamics)
        eq = f(PV, plankton_syms)
        eq isa Equation || throw(ArgumentError("biogeochem dynamics $(nameof(f)) must return Equation"))
        push!(tracer_exprs, expr(eq))
    end

    # Preserve the ordering of `plankton_args` via `ctx`.
    for idx in eachindex(plankton_syms)
        g = ctx.group_symbols[idx]
        f = getfield(plankton_dynamics, g)
        tracer_sym = plankton_syms[idx]
        eq = f(PV, plankton_syms, tracer_sym, idx)
        eq isa Equation || throw(ArgumentError("plankton dynamics $(nameof(f)) must return Equation"))
        push!(tracer_exprs, expr(eq))
    end

    return NamedTuple{Tuple(tracer_names)}(Tuple(tracer_exprs))
end

# -----------------------------------------------------------------------------
# Expression validation
# -----------------------------------------------------------------------------

"""Return all `Symbol`s referenced by an expression tree.

Supports `Expr`, `Symbol`, and literal roots (numbers). This is used by
`expression_check` to validate expressions at construction time.
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
        if s ∉ allowed_symbols && !isdefined(module_name, s) && !isdefined(Base, s) && !isdefined(Core, s)
            throw(UndefVarError(s))
        end
    end

    return nothing
end



# -----------------------------------------------------------------------------
# Functor-based biogeochemistry (3A-style driver)
# -----------------------------------------------------------------------------

"""A concrete Oceananigans biogeochemistry model.

This model dispatches tracer tendencies through a stored `NamedTuple` of callables
instead of defining per-tracer methods via runtime `eval`.

This is the 3A-style approach: *store callables, not generated methods*.
"""
struct AgateBGC{PT,TF,AF,SV} <: AbstractContinuousFormBiogeochemistry
    parameters::PT
    tracer_functions::TF
    auxiliary_fields::AF
    sinking_velocities::SV
end

Adapt.@adapt_structure AgateBGC

@inline required_biogeochemical_tracers(bgc::AgateBGC) = keys(bgc.tracer_functions)
@inline required_biogeochemical_auxiliary_fields(bgc::AgateBGC) = bgc.auxiliary_fields

@inline function (bgc::AgateBGC)(::Val{tracer}, args...) where {tracer}
    f = getfield(bgc.tracer_functions, tracer)
    return f(bgc, args...)
end

@inline function biogeochemical_drift_velocity(bgc::AgateBGC, ::Val{tracer}) where {tracer}
    sv = bgc.sinking_velocities
    if sv === nothing
        return (u=ZeroField(), v=ZeroField(), w=ZeroField())
    end

    if hasproperty(sv, tracer)
        return (u=ZeroField(), v=ZeroField(), w=getproperty(sv, tracer))
    end

    return (u=ZeroField(), v=ZeroField(), w=ZeroField())
end

"""A callable factory returned by `define_tracer_functions`.

The factory validates parameters and produces `AgateBGC` instances.
"""
struct AgateBGCFactory{TF,AF,RP,SV}
    tracer_functions::TF
    auxiliary_fields::AF
    required_params::RP
    default_sinking_velocities::SV
end

@inline function _validate_parameters(parameters, required_params)
    isempty(required_params) && return nothing
    keys = propertynames(_parameter_view(parameters))
    for k in required_params
        if k ∉ keys
            throw(ArgumentError(
                "Provided parameters are missing required field :$(k). " *
                "Resolve parameters using `resolve_runtime_parameters` (merged requirements), " *
                "or include :$(k) in your custom parameter struct.",
            ))
        end
    end
    return nothing
end

function (f::AgateBGCFactory)(parameters)
    _validate_parameters(parameters, f.required_params)
    return AgateBGC(parameters, f.tracer_functions, f.auxiliary_fields, f.default_sinking_velocities)
end

function (f::AgateBGCFactory)(parameters, sinking_velocities)
    _validate_parameters(parameters, f.required_params)
    return AgateBGC(parameters, f.tracer_functions, f.auxiliary_fields, sinking_velocities)
end

@inline function _compile_tracer_function(tracer_expression, used_params::Vector{Symbol}, all_state_vars)
    # Arguments: bgc, x, y, z, t, tracers..., aux...
    arg_tuple = Expr(:tuple, :bgc, all_state_vars...)

    assigns = Any[]
    for k in used_params
        push!(assigns, :($k = _parameter_view(bgc.parameters).$k))
    end

    body = Expr(:block, assigns..., Expr(:return, tracer_expression))
    lam = Expr(:->, arg_tuple, body)
    return eval(lam)
end

function _compile_tracer_functions(parameters, tracers::NamedTuple; auxiliary_fields::Tuple=())
    coordinates = (:x, :y, :z, :t)
    tracer_vars = keys(tracers)
    aux_field_vars = auxiliary_fields
    all_state_vars = (coordinates..., tracer_vars..., aux_field_vars...)

    parameter_keys = collect(propertynames(_parameter_view(parameters)))

    compiled = Any[]
    required_params = Symbol[]

    for (tracer_name, tracer_val) in pairs(tracers)
        kind, payload, tracer_req, has_req = _tracer_payload(tracer_val)

        # Fast path: store user-provided callables directly.
        # Note: we must not store `CompiledEquation` wrappers because they carry
        # heap-allocated `Requirements` metadata that is not GPU-friendly.
        if kind === :callable
            if has_req
                used_params = _unique_params_from_requirements(tracer_req)
                for k in used_params
                    if k ∉ parameter_keys
                        throw(ArgumentError(
                            "Tracer :$(tracer_name) declares requirement :$(k), but it is not present in the provided parameters. " *
                            "If this callable does not need parameters, wrap it with `CompiledEquation(f, req())`.",
                        ))
                    end
                end

                # Track requirements so `AgateBGCFactory` can validate user-supplied
                # parameter bundles.
                for k in used_params
                    k in required_params || push!(required_params, k)
                    if k in coordinates
                        throw(ArgumentError("Parameter name $(k) is reserved for coordinates."))
                    end
                end
            end

            push!(compiled, payload)
            continue
        end

        tracer_expression = payload

        used_params = Symbol[]
        if has_req
            used_params = _unique_params_from_requirements(tracer_req)

            # Fail fast if the caller provided an incompatible parameter bundle.
            for k in used_params
                if k ∉ parameter_keys
                    throw(ArgumentError(
                        "Tracer :$(tracer_name) requires parameter :$(k), but it is not present in the provided parameters. " *
                        "Resolve parameters using `resolve_runtime_parameters` (merged requirements), or include :$(k) in your custom parameter struct.",
                    ))
                end
            end
        else
            used_symbols = parse_expression(tracer_expression)
            @inbounds for k in parameter_keys
                (k in used_symbols) && push!(used_params, k)
            end
        end

        # Guard against reserved coordinate names appearing as parameters.
        for k in used_params
            k in required_params || push!(required_params, k)
            if k in coordinates
                throw(ArgumentError("Parameter name $(k) is reserved for coordinates."))
            end
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

        f = _compile_tracer_function(tracer_expression, used_params, all_state_vars)
        push!(compiled, f)
    end

    tracer_functions = NamedTuple{Tuple(tracer_vars)}(Tuple(compiled))
    return tracer_functions, required_params
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

Create a callable Oceananigans biogeochemistry model factory.

The returned object can be called like a type constructor:

    bgc_factory = define_tracer_functions(parameters, tracers)
    bgc = bgc_factory(parameters)


# Arguments
- `parameters`: runtime parameter struct stored as `bgc.parameters`
- `tracers`: `NamedTuple` mapping tracer names (Symbols) to tracer tendencies. Each value may be:
  - `Expr` (legacy; parameters inferred by scanning symbols),
  - `Agate.Equations.Equation` (preferred symbolic constructor API; explicit requirements),
  - `Function` (already-compiled tracer function with full signature),
  - `Agate.Functors.CompiledEquation` (a callable + explicit requirements; only the callable is stored at runtime).

Callable tracer values must accept the same arguments as Oceananigans biogeochemistry kernels:

    f(bgc, x, y, z, t, tracers..., auxiliary_fields...)

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

    if !isnothing(helper_functions)
        include(helper_functions)
    end

    tracer_functions, required_params = _compile_tracer_functions(
        parameters,
        tracers;
        auxiliary_fields=auxiliary_fields,
    )

    return AgateBGCFactory(tracer_functions, auxiliary_fields, required_params, sinking_velocities)
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

@inline function _unique_params_from_requirements(r)
    out = Symbol[]
    for k in r.vectors
        k in out || push!(out, k)
    end
    for k in r.matrices
        k in out || push!(out, k)
    end
    for k in r.scalars
        k in out || push!(out, k)
    end
    return out
end

@inline function _tracer_payload(tr)
    if tr isa Equation
        return :expr, expr(tr), requirements(tr), true
    elseif tr isa Expr
        return :expr, tr, nothing, false
    elseif tr isa CompiledEquation
        return :callable, tr.f, requirements(tr), true
    elseif tr isa Function
        return :callable, tr, nothing, false
    else
        throw(ArgumentError(
            "Tracer map values must be Expr, Agate.Equations.Equation, Agate.Functors.CompiledEquation, or a Function; got $(typeof(tr)).",
        ))
    end
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
    # Preferred path: if tracer dynamics were authored via `Agate.Equations`, then parameter
    # requirements are the source of truth and we do NOT re-parse the expression to infer
    # parameter usage.
    #
    # Backwards-compatible path: if a tracer is provided as a plain `Expr`, we infer used
    # parameters by scanning for parameter names.
    parameter_keys = collect(propertynames(_parameter_view(parameters)))

    for (tracer_name, tracer_val) in pairs(tracers)
        kind, tracer_expression, tracer_req, has_req = _tracer_payload(tracer_val)

        if kind === :callable
            throw(ArgumentError(
                "add_bgc_methods! only supports Expr and Agate.Equations.Equation tracer values. " *
                "For callable tracers, use define_tracer_functions (AgateBGCFactory) instead.",
            ))
        end

        # Determine which parameters to bind into the method body.
        used_params = Symbol[]
        if has_req
            used_params = _unique_params_from_requirements(tracer_req)

            # Fail fast if the caller provided an incompatible parameter bundle.
            for k in used_params
                if k ∉ parameter_keys
                    throw(ArgumentError(
                        "Tracer :$(tracer_name) requires parameter :$(k), but it is not present in the provided parameters. " *
                        "Resolve parameters using `resolve_runtime_parameters` (merged requirements), or include :$(k) in your custom parameter struct.",
                    ))
                end
            end
        else
            used_symbols = parse_expression(tracer_expression)
            @inbounds for k in parameter_keys
                (k in used_symbols) && push!(used_params, k)
            end
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
