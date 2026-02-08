# -----------------------------------------------------------------------------
# Functor-based biogeochemistry driver
# -----------------------------------------------------------------------------

using Adapt

using ..Equations: CompiledEquation
using ..Runtime: Tracers, TracerIndex, build_tracer_index

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

The model stores a small, GPU-safe `Tracers` accessor (`bgc.tracers`) that
converts human-friendly names (e.g. `tr.N`, `tr.PAR`, `tr.plankton`) into
integer indexing into the positional argument list passed to kernels.
"""
struct AgateBGC{PT,TF,TR,SV} <: AbstractContinuousFormBiogeochemistry
    parameters::PT
    tracer_functions::TF
    tracers::TR
    sinking_velocities::SV
end

Adapt.@adapt_structure AgateBGC

@inline required_biogeochemical_tracers(bgc::AgateBGC) = keys(bgc.tracer_functions)

@generated function _auxiliary_fields_from_tracers(::Type{TR}) where {TR}
    # TR is expected to be `Tracers{TracerIndex{TRSYMS,GS,AF,NG}}`.
    TI = TR.parameters[1]
    AF = TI.parameters[3]
    return :($AF)
end

@inline required_biogeochemical_auxiliary_fields(bgc::AgateBGC) =
    _auxiliary_fields_from_tracers(typeof(bgc.tracers))

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
struct AgateBGCFactory{TF,TI,RP,SV}
    tracer_functions::TF
    tracer_index::TI
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
    tr = Tracers(f.tracer_index)
    return AgateBGC(parameters, f.tracer_functions, tr, f.default_sinking_velocities)
end

function (f::AgateBGCFactory)(parameters, sinking_velocities)
    _validate_parameters(parameters, f.required_params)
    tr = Tracers(f.tracer_index)
    return AgateBGC(parameters, f.tracer_functions, tr, sinking_velocities)
end

function _compile_tracer_functions(parameters, tracers::NamedTuple)
    coordinates = (:x, :y, :z, :t)

    # Keep this as a tuple for fast `in` checks without allocations.
    parameter_keys = propertynames(_parameter_view(parameters))

    # Disallow reserved names that are used as Oceananigans kernel arguments.
    for k in parameter_keys
        if k in coordinates
            throw(ArgumentError("Parameters contain reserved field :$(k)."))
        end
    end

    for (_, tracer_val) in pairs(tracers)
        tracer_val isa CompiledEquation || throw(
            ArgumentError(
                "Tracer map values must be Agate.Equations.CompiledEquation; got $(typeof(tracer_val)).",
            ),
        )
    end

    # IMPORTANT: avoid type erasure.
    # Building via `Vector{Any}` (or `Tuple(vec)`) would produce `Any`-typed fields,
    # which triggers dynamic dispatch in kernels and breaks GPU compilation.
    tracer_functions = map(tr -> tr.f, tracers)

    # By construction, the parameter object passed to this factory defines the required
    # fields expected when instantiating `AgateBGC`.
    required_params = parameter_keys

    return tracer_functions, required_params
end

"""    define_tracer_functions(parameters, tracers; auxiliary_fields=(:PAR,), tracer_index=nothing, sinking_velocities=nothing)

Create a callable Oceananigans biogeochemistry model factory.

`tracers` is a `NamedTuple` mapping tracer names to `CompiledEquation` values.

Each compiled equation wraps a callable `f`. Parameter validation is handled upstream by model construction.

Callable tracer values must accept the Oceananigans biogeochemistry kernel signature:

    f(bgc, x, y, z, t, tracers..., auxiliary_fields...)

Notes
-----
- `auxiliary_fields` defines the ordered list of auxiliary values appended to the tracer argument list.
- `tracer_index` controls how `bgc.tracers` maps names to positional indices. If omitted, a scalar-only
  index is generated from `keys(tracers)`.
"""
function define_tracer_functions(
    parameters,
    tracers::NamedTuple;
    auxiliary_fields::Tuple=(:PAR,),
    tracer_index=nothing,
    sinking_velocities=nothing,
)
    tracer_functions, required_params = _compile_tracer_functions(parameters, tracers)

    idx = if tracer_index === nothing
        build_tracer_index(keys(tracers), auxiliary_fields)
    else
        tracer_index
    end

    return AgateBGCFactory(tracer_functions, idx, required_params, sinking_velocities)
end
