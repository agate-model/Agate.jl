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

using Agate.Models.NiPiZD.Parameters:
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

Construct a size-structured NiPiZD model abstract type.

This constructor builds a size-structured plankton model with two plankton functional types:
phytoplankton (P) and zooplankton (Z), each of which can be specified to have any number of
size classes (`n_phyto` and `n_zoo`). In addition to plankton, the constructor implements
idealized detritus (D) and nutrient (N) cycling by default, although more complex N and D
cycling can also be defined using the `nutrient_dynamics` and `detritus_dynamics` arguments.

During model construction, the size of each plankton determines photosynthetic growth rates,
nutrient half saturation constants, predation rates, and predator-prey assimilation and
palatability values. Alternatively, if manually defined predator-prey assimilation and
palatability values are desired, these can be specified using the `palatability_matrix` and
`assimilation_efficiency_matrix` arguments.

Note that if non-default `*_dynamics` expressions are passed, the relevant parameter
containers (`phyto_pft_parameters`, `zoo_pft_parameters`, and/or `bgc_specification`) also
need to be specified consistently.

The type specification includes a photosynthetic active radiation (PAR) auxiliary field.

`FT` sets the floating point type stored in runtime parameter containers. All float
parameters are cast exactly once during construction.

If `sinking_tracers` is provided and `grid` is not provided, the default `BoxModelGrid()`
is used (matching the original Agate behavior).

# Keywords
- `FT`: floating point type used for runtime parameters (e.g. `Float32`, `Float64`)
- `n_phyto`: number of phytoplankton in the model
- `n_zoo`: number of zooplankton in the model
- `phyto_diameters`: diameter specification (e.g. `AbstractDiameterSpecification`) or a
    list of values to use (as many as the model expects)
- `zoo_diameters`: diameter specification (e.g. `AbstractDiameterSpecification`) or a list
    of values to use (as many as the model expects)
- `nutrient_dynamics`: expression describing how nutrients change over time, see
    `Agate.Models.Tracers`
- `detritus_dynamics`: expression describing how detritus evolves over time, see
    `Agate.Models.Tracers`
- `phyto_dynamics`: expression describing how phytoplankton grow, see `Agate.Models.Tracers`
- `zoo_dynamics`: expression describing how zooplankton grow, see `Agate.Models.Tracers`
- `phyto_pft_parameters`: phytoplankton PFT parameter container; for default values see
    `default_phyto_pft_parameters(FT)`
- `zoo_pft_parameters`: zooplankton PFT parameter container; for default values see
    `default_zoo_pft_parameters(FT)`
- `bgc_specification`: biogeochemistry specification (nutrient/detritus coupling etc.); for
    default values see `default_nipizd_bgc_specification(FT)`
- `palatability_matrix`: optional palatability matrix passed as an Array; if provided,
    then any internally computed palatability values are overridden
- `assimilation_efficiency_matrix`: optional assimilation efficiency matrix passed as an
    Array; if provided, then any internally computed assimilation efficiency values are overridden
- `sinking_tracers`: optional NamedTuple of sinking speeds (passed as positive values) of
   the form (<tracer name expressed as symbol> = <speed>, ...)
- `grid`: optional Oceananigans grid object defining the geometry to build the model on, must
   be passed if `sinking_tracers` is defined, defaults to `BoxModelGrid`
- `open_bottom`: indicates whether the sinking velocity should be smoothly brought to zero
   at the bottom to prevent the tracers leaving the domain, defaults to `true`, which means
   the bottom is open and the tracers leave (i.e., no slowing of velocity to 0 is applied)

# Example
```julia
using Agate.Constructors: NiPiZD

n2p2zd = NiPiZD.construct()
n2p2zd_model_obj = n2p2zd()
```
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

A function to instantiate an object of `bgc_type` returned by `NiPiZD.construct()`.

The type specifies the number of phytoplankton and zooplankton in the model and includes
default parameter values. The instantiate method is used to override the default values
of any of the model parameters or plankton diameters.

This method re-expands runtime parameters, enabling overrides of diameter
specifications, PFT parameters, and interaction matrices.

If `bgc_type` defines sinking velocities and `sinking_tracers` is provided, sinking
velocity fields are rebuilt using `grid` (defaulting to `BoxModelGrid()`), matching the
original Agate behavior.

!!! tip

    Changing the parameter values of an existing NiPiZD model type using `instantiate()` is useful in
    dynamic programming contexts such as `for` loops.


# Arguments
- `bgc_type`: subtype of Oceananigans.Biogeochemistry returned by `NiPiZD.construct()`
   with a specified number of phytoplankton and zooplankton

# Keywords
- `FT`: floating point type used for rebuilt runtime parameters (e.g. `Float32`, `Float64`)
- `phyto_diameters`: diameter specification (e.g. `AbstractDiameterSpecification`) or a list of
    values to use (as many as the model expects). If `nothing`, uses the phytoplankton
    diameters encoded in `bgc_type`.
- `zoo_diameters`: diameter specification (e.g. `AbstractDiameterSpecification`) or a list of
    values to use (as many as the model expects). If `nothing`, uses the zooplankton
    diameters encoded in `bgc_type`.
- `phyto_pft_parameters`: phytoplankton PFT parameter container; if `nothing`, defaults to
    `default_phyto_pft_parameters(FT)`
- `zoo_pft_parameters`: zooplankton PFT parameter container; if `nothing`, defaults to
    `default_zoo_pft_parameters(FT)`
- `bgc_specification`: biogeochemistry specification; if `nothing`, a specification is
    reconstructed from the defaults in `bgc_type()`
- `palatability_matrix`: optional palatability matrix passed as an Array; if provided,
    overrides any internally computed palatability values
- `assimilation_efficiency_matrix`: optional assimilation efficiency matrix passed as an
    Array; if provided, overrides any internally computed assimilation efficiency values
- `sinking_tracers`: optional NamedTuple of sinking speeds (passed as positive values) of
   the form (<tracer name expressed as symbol> = <speed>, ...)
- `grid`: optional Oceananigans grid object defining the geometry to build sinking velocity
   fields on (used when `sinking_tracers` is provided), defaults to `BoxModelGrid`
- `open_bottom`: indicates whether the sinking velocity should be smoothly brought to zero
   at the bottom to prevent the tracers leaving the domain, defaults to `true`, which means
   the bottom is open and the tracers leave (i.e., no slowing of velocity to 0 is applied)

# Example
```julia
using Agate.Constructors: NiPiZD

n2p2zd = NiPiZD.construct()

# change some parameter values
phyto_args = deepcopy(NiPiZD.DEFAULT_PHYTO_ARGS)
phyto_args["allometry"]["maximum_growth_rate"]["a"] = 2
n2p2zd_model_obj = NiPiZD.instantiate(n2p2zd; phyto_args=phyto_args)
```
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
