using Adapt

using ..Runtime: ParameterizedBGC, evaluate_tendency
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using Oceananigans.Fields: ZeroField


import Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity,
    required_biogeochemical_auxiliary_fields,
    required_biogeochemical_tracers

"""Concrete Oceananigans biogeochemistry built from static compiled equations."""
struct AgateBGC{PT,EQ,SV,MD,AF} <: AbstractContinuousFormBiogeochemistry
    parameters::PT
    equations::EQ
    sinking_velocities::SV
    metadata::MD
end

function AgateBGC(parameters, equations::NamedTuple, auxiliary_fields::Tuple, sinking_velocities, metadata)
    return AgateBGC{
        typeof(parameters), typeof(equations), typeof(sinking_velocities), typeof(metadata), auxiliary_fields,
    }(parameters, equations, sinking_velocities, metadata)
end

@inline function Adapt.adapt_structure(
    to, bgc::AgateBGC{PT,EQ,SV,MD,AF}
) where {PT,EQ,SV,MD,AF}
    # `metadata` is host-side introspection/setup state. Device kernel copies must not
    # carry Symbol-valued metadata because GPU kernel arguments must be bitstypes.
    adapted_metadata = to === identity ? bgc.metadata : nothing
    return AgateBGC(
        Adapt.adapt(to, bgc.parameters),
        Adapt.adapt(to, bgc.equations),
        AF,
        Adapt.adapt(to, bgc.sinking_velocities),
        adapted_metadata,
    )
end

@inline required_biogeochemical_tracers(bgc::AgateBGC) = keys(bgc.equations)

@inline required_biogeochemical_tracers(bgc::ParameterizedBGC) =
    required_biogeochemical_tracers(bgc.bgc)

@inline function required_biogeochemical_tracers(
    ::Type{<:AgateBGC{PT,EQ}}
) where {PT,EQ<:NamedTuple}
    return fieldnames(EQ)
end

@inline required_biogeochemical_auxiliary_fields(::AgateBGC{PT,EQ,SV,MD,AF}) where {PT,EQ,SV,MD,AF} = AF

@inline required_biogeochemical_auxiliary_fields(bgc::ParameterizedBGC) =
    required_biogeochemical_auxiliary_fields(bgc.bgc)

@inline function required_biogeochemical_auxiliary_fields(
    ::Type{<:AgateBGC{PT,EQ,SV,MD,AF}}
) where {PT,EQ,SV,MD,AF}
    return AF
end

@inline function (bgc::AgateBGC)(val_name::Val, args...)
    return evaluate_tendency(bgc, val_name, args...)
end

@inline function biogeochemical_drift_velocity(bgc::ParameterizedBGC, val_name::Val)
    return biogeochemical_drift_velocity(bgc.bgc, val_name)
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
