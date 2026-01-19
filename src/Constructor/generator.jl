# -----------------------------------------------------------------------------
# Functor-based biogeochemistry driver
# -----------------------------------------------------------------------------

using Adapt

using ..Functors: CompiledEquation, requirements
using OceanBioME
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using Oceananigans.Fields: ZeroField

import Oceananigans: time_step!

import Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity,
    required_biogeochemical_auxiliary_fields,
    required_biogeochemical_tracers

export define_tracer_functions

"""A concrete Oceananigans biogeochemistry model.

Tracer tendencies are dispatched through a stored `NamedTuple` of callables.
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

@inline _parameter_view(parameters) = parameters

@inline function _validate_parameters(parameters, required_params)
    isempty(required_params) && return nothing
    keys = propertynames(_parameter_view(parameters))
    for k in required_params
        if k ∉ keys
            throw(ArgumentError("Provided parameters are missing required field :$(k)."))
        end
    end
    return nothing
end

function (f::AgateBGCFactory)(parameters)
    _validate_parameters(parameters, f.required_params)
    return AgateBGC(
        parameters, f.tracer_functions, f.auxiliary_fields, f.default_sinking_velocities
    )
end

function (f::AgateBGCFactory)(parameters, sinking_velocities)
    _validate_parameters(parameters, f.required_params)
    return AgateBGC(parameters, f.tracer_functions, f.auxiliary_fields, sinking_velocities)
end

@inline function _unique_params_from_requirements(r)
    out = Symbol[]
    for k in r.vectors
        (k in out) || push!(out, k)
    end
    for k in r.matrices
        (k in out) || push!(out, k)
    end
    for k in r.scalars
        (k in out) || push!(out, k)
    end
    return out
end

function _compile_tracer_functions(parameters, tracers::NamedTuple)
    coordinates = (:x, :y, :z, :t)
    parameter_keys = collect(propertynames(_parameter_view(parameters)))

    compiled = Any[]
    required_params = Symbol[]

    for (tracer_name, tracer_val) in pairs(tracers)
        tracer_val isa CompiledEquation || throw(
            ArgumentError(
                "Tracer map values must be Agate.Functors.CompiledEquation; got $(typeof(tracer_val)).",
            ),
        )

        r = requirements(tracer_val)
        used_params = _unique_params_from_requirements(r)

        for k in used_params
            if k in coordinates
                throw(
                    ArgumentError(
                        "Tracer :$(tracer_name) declares reserved parameter name :$(k)."
                    ),
                )
            end
            if k ∉ parameter_keys
                throw(
                    ArgumentError(
                        "Tracer :$(tracer_name) requires parameter :$(k), but it is not present in the provided parameters.",
                    ),
                )
            end
            (k in required_params) || push!(required_params, k)
        end

        # Store only the callable (GPU-friendly).
        push!(compiled, tracer_val.f)
    end

    tracer_vars = keys(tracers)
    tracer_functions = NamedTuple{Tuple(tracer_vars)}(Tuple(compiled))
    return tracer_functions, required_params
end

"""    define_tracer_functions(parameters, tracers; auxiliary_fields=(:PAR,), sinking_velocities=nothing)

Create a callable Oceananigans biogeochemistry model factory.

`tracers` is a `NamedTuple` mapping tracer names to `CompiledEquation` values.

Each compiled equation wraps a callable `f` plus a `Requirements` object describing which model
parameters are accessed by the callable.

Callable tracer values must accept the Oceananigans biogeochemistry kernel signature:

    f(bgc, x, y, z, t, tracers..., auxiliary_fields...)
"""
function define_tracer_functions(
    parameters,
    tracers::NamedTuple;
    auxiliary_fields::Tuple=(:PAR,),
    sinking_velocities=nothing,
)
    Base.@nospecialize parameters tracers auxiliary_fields sinking_velocities

    tracer_functions, required_params = _compile_tracer_functions(parameters, tracers)

    return AgateBGCFactory(
        tracer_functions, auxiliary_fields, required_params, sinking_velocities
    )
end
