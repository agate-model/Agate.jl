"""Constructors for the size-structured NiPiZD biogeochemistry model.

This module constructs Oceananigans-compatible `AbstractContinuousFormBiogeochemistry`
subtypes whose tracer tendencies are defined by expressions from
`Agate.Models.NiPiZD.Tracers`.

Agate-owned runtime data is stored in a single concrete parameter container in the
biogeochemistry object: `bgc.parameters`, where `parameters` is a `NiPiZDParameters`
struct created by expanding construction specifications.

Tracer tendency containers are supplied as `NamedTuple`s at the Oceananigans boundary.
"""

module Constructor

using OceanBioME
using OceanBioME: BoxModelGrid, setup_velocity_fields
using Oceananigans.Units

using Agate.Utils: define_tracer_functions

using Agate.Utils:
    AbstractDiameterSpecification, DiameterListSpecification, DiameterRangeSpecification

using Agate.Models.Parameters:
    NiPiZDBiogeochemistrySpecification,
    PhytoPFTParameters,
    ZooPFTParameters,
    create_nipizd_parameters,
    default_phyto_pft_parameters,
    default_zoo_pft_parameters,
    default_nipizd_bgc_specification,
    default_phyto_geider_pft_parameters

using Agate.Models.NiPiZD.Tracers:
    detritus_default,
    nutrient_default,
    nutrient_geider_light,
    phytoplankton_default,
    phytoplankton_geider_light,
    zooplankton_default

export construct
export instantiate
export default_phyto_pft_parameters
export default_zoo_pft_parameters
export default_nipizd_bgc_specification
export default_phyto_geider_pft_parameters

"""Return a diameter specification for an explicit diameter list."""
diameter_specification(diameters::AbstractVector) = DiameterListSpecification(diameters)

"""Return the diameter specification when one is already provided."""
diameter_specification(spec::AbstractDiameterSpecification) = spec

function plankton_symbols(n_P::Int, n_Z::Int)
    names = Vector{Symbol}(undef, n_Z + n_P)

    @inbounds for i in 1:n_Z
        names[i] = Symbol(:Z, i)
    end

    @inbounds for i in 1:n_P
        names[n_Z + i] = Symbol(:P, i)
    end

    return names
end

function build_tracer_expressions(
    n_P::Int,
    n_Z::Int;
    nutrient_dynamics=nutrient_default,
    detritus_dynamics=detritus_default,
    phyto_dynamics=phytoplankton_default,
    zoo_dynamics=zooplankton_default,
)
    plankton_syms = plankton_symbols(n_P, n_Z)

    tracer_names = Symbol[:N, :D]
    append!(tracer_names, [Symbol("Z", i) for i in 1:n_Z])
    append!(tracer_names, [Symbol("P", i) for i in 1:n_P])

    tracer_exprs = Expr[]
    push!(tracer_exprs, nutrient_dynamics(plankton_syms))
    push!(tracer_exprs, detritus_dynamics(plankton_syms))

    for i in 1:n_Z
        sym = Symbol("Z", i)
        push!(tracer_exprs, zoo_dynamics(plankton_syms, sym, i))
    end

    for i in 1:n_P
        sym = Symbol("P", i)
        plankton_idx = n_Z + i
        push!(tracer_exprs, phyto_dynamics(plankton_syms, sym, plankton_idx))
    end

    return NamedTuple{Tuple(tracer_names)}(Tuple(tracer_exprs))
end

"""
    construct(; FT=Float64, n_phyto=2, n_zoo=2, phyto_diameters=..., zoo_diameters=...,
              phyto_pft_parameters=..., zoo_pft_parameters=..., bgc_specification=...,
              nutrient_dynamics=..., detritus_dynamics=..., phyto_dynamics=..., zoo_dynamics=...,
              palatability_matrix=nothing, assimilation_efficiency_matrix=nothing,
              sinking_tracers=nothing, grid=BoxModelGrid(), open_bottom=true)

Construct a NiPiZD biogeochemistry model type.

`FT` sets the floating point type stored in runtime parameter containers. All float
parameters are cast exactly once during construction.

If `sinking_tracers` is provided and `grid` is not provided, the default `BoxModelGrid()`
is used (matching the original Agate behavior).
"""
function construct(;
    FT::Type{<:AbstractFloat}=Float64,
    n_phyto::Int=2,
    n_zoo::Int=2,
    phyto_diameters::Union{AbstractDiameterSpecification,AbstractVector}=DiameterRangeSpecification(
        2, 10, :log_splitting
    ),
    zoo_diameters::Union{AbstractDiameterSpecification,AbstractVector}=DiameterRangeSpecification(
        20, 100, :linear_splitting
    ),
    phyto_pft_parameters=default_phyto_pft_parameters(FT),
    zoo_pft_parameters=default_zoo_pft_parameters(FT),
    bgc_specification=default_nipizd_bgc_specification(FT),
    nutrient_dynamics=nutrient_default,
    detritus_dynamics=detritus_default,
    phyto_dynamics=phytoplankton_default,
    zoo_dynamics=zooplankton_default,
    palatability_matrix=nothing,
    assimilation_efficiency_matrix=nothing,
    sinking_tracers=nothing,
    grid=BoxModelGrid(),
    open_bottom::Bool=true,
)
    phyto_spec = diameter_specification(phyto_diameters)
    zoo_spec = diameter_specification(zoo_diameters)

    parameters = create_nipizd_parameters(
        FT;
        n_phyto=n_phyto,
        n_zoo=n_zoo,
        phyto_diameters=phyto_spec,
        zoo_diameters=zoo_spec,
        phyto_pft_parameters=phyto_pft_parameters,
        zoo_pft_parameters=zoo_pft_parameters,
        bgc_specification=bgc_specification,
        palatability_matrix=palatability_matrix,
        assimilation_efficiency_matrix=assimilation_efficiency_matrix,
    )

    tracers = build_tracer_expressions(
        n_phyto,
        n_zoo;
        nutrient_dynamics=nutrient_dynamics,
        detritus_dynamics=detritus_dynamics,
        phyto_dynamics=phyto_dynamics,
        zoo_dynamics=zoo_dynamics,
    )

    if isnothing(sinking_tracers)
        return define_tracer_functions(parameters, tracers)
    end

    sinking_velocities = setup_velocity_fields(sinking_tracers, grid, open_bottom)
    return define_tracer_functions(
        parameters, tracers; sinking_velocities=sinking_velocities
    )
end

"""
    instantiate(bgc_type; kwargs...)

Instantiate an existing NiPiZD biogeochemistry type.

This method re-expands runtime parameters, enabling overrides of diameter
specifications, PFT parameters, and interaction matrices.

If `bgc_type` defines sinking velocities and `sinking_tracers` is provided, sinking
velocity fields are rebuilt using `grid` (defaulting to `BoxModelGrid()`), matching the
original Agate behavior.
"""
function instantiate(
    bgc_type;
    FT::Type{<:AbstractFloat}=Float64,
    phyto_diameters=nothing,
    zoo_diameters=nothing,
    phyto_pft_parameters=nothing,
    zoo_pft_parameters=nothing,
    bgc_specification=nothing,
    palatability_matrix=nothing,
    assimilation_efficiency_matrix=nothing,
    sinking_tracers=nothing,
    grid=BoxModelGrid(),
    open_bottom::Bool=true,
)
    defaults = bgc_type()
    p0 = defaults.parameters

    n_phyto = p0.n_P
    n_zoo = p0.n_Z

    phyto_spec = if isnothing(phyto_diameters)
        DiameterListSpecification(p0.diameters[(n_zoo + 1):(n_zoo + n_phyto)])
    else
        diameter_specification(phyto_diameters)
    end

    zoo_spec = if isnothing(zoo_diameters)
        DiameterListSpecification(p0.diameters[1:n_zoo])
    else
        diameter_specification(zoo_diameters)
    end

    phyto_pft = if isnothing(phyto_pft_parameters)
        default_phyto_pft_parameters(FT)
    else
        phyto_pft_parameters
    end
    zoo_pft =
        isnothing(zoo_pft_parameters) ? default_zoo_pft_parameters(FT) : zoo_pft_parameters

    bgc_spec = if isnothing(bgc_specification)
        NiPiZDBiogeochemistrySpecification{FT}(
            FT(p0.detritus_remineralization), FT(p0.mortality_export_fraction)
        )
    else
        bgc_specification
    end

    parameters = create_nipizd_parameters(
        FT;
        n_phyto=n_phyto,
        n_zoo=n_zoo,
        phyto_diameters=phyto_spec,
        zoo_diameters=zoo_spec,
        phyto_pft_parameters=phyto_pft,
        zoo_pft_parameters=zoo_pft,
        bgc_specification=bgc_spec,
        palatability_matrix=palatability_matrix,
        assimilation_efficiency_matrix=assimilation_efficiency_matrix,
    )

    if :sinking_velocities in fieldnames(bgc_type)
        sinking_velocities = if isnothing(sinking_tracers)
            defaults.sinking_velocities
        else
            setup_velocity_fields(sinking_tracers, grid, open_bottom)
        end
        return bgc_type(; parameters=parameters, sinking_velocities=sinking_velocities)
    end

    if !isnothing(sinking_tracers)
        throw(ArgumentError("biogeochemistry type does not define sinking velocities"))
    end

    return bgc_type(; parameters=parameters)
end

end # module
