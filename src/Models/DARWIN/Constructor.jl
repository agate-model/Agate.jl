"""
Module to construct an instance of an Agate.jl-DARWIN model.
"""

module Constructor

using OceanBioME
using OceanBioME: BoxModelGrid, setup_velocity_fields
using Oceananigans.Units

using Agate.Utils:
    define_tracer_functions,
    AbstractDiameterSpecification,
    DiameterListSpecification,
    DiameterRangeSpecification

using Agate.Models.DARWIN.Parameters:
    DarwinBiogeochemistrySpecification,
    create_darwin_parameters,
    default_darwin_phyto_parameters,
    default_darwin_zoo_parameters,
    default_darwin_bgc_specification

using Agate.Models.DARWIN.Tracers:
    DIC_geider_light,
    DIN_geider_light,
    PO4_geider_light,
    DOC_default,
    POC_default,
    DON_default,
    PON_default,
    DOP_default,
    POP_default,
    phytoplankton_growth_two_nutrients_geider_light,
    zooplankton_default

export construct
export instantiate

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
    DIC_dynamics=DIC_geider_light,
    DIN_dynamics=DIN_geider_light,
    PO4_dynamics=PO4_geider_light,
    DOC_dynamics=DOC_default,
    POC_dynamics=POC_default,
    DON_dynamics=DON_default,
    PON_dynamics=PON_default,
    DOP_dynamics=DOP_default,
    POP_dynamics=POP_default,
    phyto_dynamics=phytoplankton_growth_two_nutrients_geider_light,
    zoo_dynamics=zooplankton_default,
)
    plankton_syms = plankton_symbols(n_P, n_Z)

    tracer_names = Symbol[:DIC, :DIN, :PO4, :DOC, :POC, :DON, :PON, :DOP, :POP]
    append!(tracer_names, [Symbol("Z", i) for i in 1:n_Z])
    append!(tracer_names, [Symbol("P", i) for i in 1:n_P])

    tracer_exprs = Expr[]
    push!(tracer_exprs, DIC_dynamics(plankton_syms))
    push!(tracer_exprs, DIN_dynamics(plankton_syms))
    push!(tracer_exprs, PO4_dynamics(plankton_syms))
    push!(tracer_exprs, DOC_dynamics(plankton_syms))
    push!(tracer_exprs, POC_dynamics(plankton_syms))
    push!(tracer_exprs, DON_dynamics(plankton_syms))
    push!(tracer_exprs, PON_dynamics(plankton_syms))
    push!(tracer_exprs, DOP_dynamics(plankton_syms))
    push!(tracer_exprs, POP_dynamics(plankton_syms))

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
    construct(;
        FT::Type{<:AbstractFloat}=Float64,
        n_phyto=2,
        n_zoo=2,
        phyto_diameters=DiameterRangeSpecification(2, 10, :log_splitting),
        zoo_diameters=DiameterRangeSpecification(20, 100, :linear_splitting),
        DIC_dynamics=DIC_geider_light,
        PO4_dynamics=PO4_geider_light,
        DIN_dynamics=DIN_geider_light,
        POC_dynamics=POC_default,
        DOC_dynamics=DOC_default,
        PON_dynamics=PON_default,
        DON_dynamics=DON_default,
        POP_dynamics=POP_default,
        DOP_dynamics=DOP_default,
        phyto_dynamics=phytoplankton_growth_two_nutrients_geider_light,
        zoo_dynamics=zooplankton_default,
        phyto_pft_parameters=default_darwin_phyto_parameters(FT),
        zoo_pft_parameters=default_darwin_zoo_parameters(FT),
        bgc_specification=default_darwin_bgc_specification(FT),
        palatability_matrix=nothing,
        assimilation_efficiency_matrix=nothing,
        sinking_tracers=nothing,
        grid=BoxModelGrid(),
        open_bottom=true,
    )

Construct an Agate.jl-DARWIN model abstract type.

!!! info

    This model is in active development and has not been validated against `MITgcm-DARWIN`.

!!! formulation

    TRACERS:

    ∂t cⱼ = ``Uⱼ``DIC - ``Mⱼ`` + ``Gⱼ`` - ``gⱼ``

    ∂t DIC = ∑(``Uⱼ`` DIC) + ``R``DOC + ``R``POC

    ∂t DIN = ∑(``Uⱼ``DIC * ``Qⱼ``N)  + ``R``DON + ``R``PON

    ∂t PO4 = ∑(``Uⱼ``DIC * ``Qⱼ``P)  + ``R``DOP + ``R``POP

    ∂t DOC = ∑(``Mⱼ``DOC) + ``g``DOC - ``R``DOC

    ∂t DON = ∑(``Mⱼ``DOC * ``Qⱼ``N) + ``g``DON - ``R``DON

    ∂t DOP = ∑(``Mⱼ``DOC * ``Qⱼ``P) + ``g``DOP - ``R``DOP

    ∂t POC = ∑(``Mⱼ``POC) + ``g``POC - ``R``POC

    ∂t PON = ∑(``Mⱼ``POC * ``Qⱼ``N) + ``g``PON - ``R``PON

    ∂t POP = ∑(``Mⱼ``POC * ``Qⱼ``P) + ``g``POP - ``R``POP

    where:
    - ``U`` = uptake
    - ``R`` = remineralization
    - ``M`` = mortality
    - ``g, G`` = grazing losses and gains
    - ``Q`` = plankton elemental ratios

    TRAITS:

    μmax, KR, gmax = a*Volume^b

    palat = η/(1+(``ratio``-``opt``)^2)^σ

    where:
    - μmax = maximum photosynthetic growth
    - KR = nutrient half saturation
    - gmax = maximum predation rate
    - palat = palatability
    - ``ratio`` = predator to prey size ratio (diameter)
    - ``opt`` = predator to prey size optimum (diameter)
    - η = prey protection
    - σ = predator specificity

This constructor builds an Agate.jl-DARWIN model with two plankton functional types:
phytoplankton (P) and zooplankton (Z), each of which can be specified to have any number of
size classes (`n_phyto` and `n_zoo`). In addition to plankton, the constructor implements
idealized dissolved inorganic carbon (DIC), particulate organic matter (POC, POP, PON), dissolved organic
matter (DOC, DOP, DON) and two nutrients (DIN and PO4) cycling by default, although more complex elemental cycling can also be defined using the tracer `*_dynamics`
arguments.

During model construction, the size of each plankton determines photosynthetic growth rates,
nutrient half saturation constants, predation rates, and optionally predator-prey assimilation
and palatability values. Alternatively, if manually defined predator-prey assimilation and
palatability values are desired, these can be defined using the `palatability_matrix` and
`assimilation_efficiency_matrix` arguments.

Note that if non-default `*_dynamics` expressions are passed, the relevant `*_pft_parameters` / `bgc_specification` also
need to be specified.

The type specification includes a photosynthetic active radiation (PAR) auxiliary field.

# Arguments
- `FT`: floating point type used to store parameters in runtime containers
- `n_phyto`: number of phytoplankton to include in the model
- `n_zoo`: number of zooplankton to include in the model
- `phyto_diameters`: diameter specification (e.g., `DiameterRangeSpecification`) from which `n_phyto` diameters
    can be computed or a list of values to use
- `zoo_diameters`: diameter specification (e.g., `DiameterRangeSpecification`) from which `n_zoo` diameters
    can be computed or a list of values to use
- `DIC_dynamics`: expression describing how DIC changes over time, see `Agate.Models.DARWIN.Tracers`
- `PO4_dynamics`: expression describing how PO4 changes over time, see `Agate.Models.DARWIN.Tracers`
- `DIN_dynamics`: expression describing how DIN changes over time, see `Agate.Models.DARWIN.Tracers`
- `POC_dynamics`: expression describing how POC changes over time, see `Agate.Models.DARWIN.Tracers`
- `DOC_dynamics`: expression describing how DOC changes over time, see `Agate.Models.DARWIN.Tracers`
- `PON_dynamics`: expression describing how PON changes over time, see `Agate.Models.DARWIN.Tracers`
- `DON_dynamics`: expression describing how DON changes over time, see `Agate.Models.DARWIN.Tracers`
- `POP_dynamics`: expression describing how POP changes over time, see `Agate.Models.DARWIN.Tracers`
- `DOP_dynamics`: expression describing how DOP changes over time, see `Agate.Models.DARWIN.Tracers`
- `phyto_dynamics`: expression describing how phytoplankton grow, see `Agate.Models.DARWIN.Tracers`
- `zoo_dynamics`: expression describing how zooplankton grow, see `Agate.Models.DARWIN.Tracers`
- `phyto_pft_parameters`: phytoplankton PFT parameters (`DarwinPhytoPFTParameters`), for default values see
    `DARWIN.default_darwin_phyto_parameters(FT)`
- `zoo_pft_parameters`: zooplankton PFT parameters (`DarwinZooPFTParameters`), for default values see
    `DARWIN.default_darwin_zoo_parameters(FT)`
- `bgc_specification`: elemental cycling specification (`DarwinBiogeochemistrySpecification`), for default values see
    `DARWIN.default_darwin_bgc_specification(FT)`
- `palatability_matrix`: optional palatability matrix passed as an Array; if not provided this is computed from
   PFT trait parameters during construction
- `assimilation_efficiency_matrix`: optional assimilation efficiency matrix passed as an Array; if not provided this
   is computed from PFT trait parameters during construction
- `sinking_tracers`: optional NamedTuple of sinking speeds (passed as positive values) of
   the form (<tracer name expressed as symbol> = <speed>, ...)
- `grid`: optional Oceananigans grid object defining the geometry to build the model on, must
   be passed if `sinking_tracers` is defined, defaults to BoxModelGrid
- `open_bottom`: indicates whether the sinking velocity should be smoothly brought to zero
   at the bottom to prevent the tracers leaving the domain, defaults to `true`, which means
   the bottom is open and the tracers leave (i.e., no slowing of velocity to 0 is applied)

# Example
```julia
using Agate.Models: DARWIN

darwin_2p_2z = DARWIN.construct()
darwin_2p_2z_model_obj = darwin_2p_2z()
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
    phyto_pft_parameters=default_darwin_phyto_parameters(FT),
    zoo_pft_parameters=default_darwin_zoo_parameters(FT),
    bgc_specification=default_darwin_bgc_specification(FT),
    DIC_dynamics=DIC_geider_light,
    DIN_dynamics=DIN_geider_light,
    PO4_dynamics=PO4_geider_light,
    DOC_dynamics=DOC_default,
    POC_dynamics=POC_default,
    DON_dynamics=DON_default,
    PON_dynamics=PON_default,
    DOP_dynamics=DOP_default,
    POP_dynamics=POP_default,
    phyto_dynamics=phytoplankton_growth_two_nutrients_geider_light,
    zoo_dynamics=zooplankton_default,
    palatability_matrix=nothing,
    assimilation_efficiency_matrix=nothing,
    sinking_tracers=nothing,
    grid=BoxModelGrid(),
    open_bottom::Bool=true,
)
    phyto_spec = diameter_specification(phyto_diameters)
    zoo_spec = diameter_specification(zoo_diameters)

    parameters = create_darwin_parameters(
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
        DIC_dynamics=DIC_dynamics,
        DIN_dynamics=DIN_dynamics,
        PO4_dynamics=PO4_dynamics,
        DOC_dynamics=DOC_dynamics,
        POC_dynamics=POC_dynamics,
        DON_dynamics=DON_dynamics,
        PON_dynamics=PON_dynamics,
        DOP_dynamics=DOP_dynamics,
        POP_dynamics=POP_dynamics,
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
    instantiate(
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
        open_bottom=true,
    )

A function to instantiate an object of `bgc_type` returned by `DARWIN.construct()`.

The type specifies the number of phytoplankton and zooplankton in the model and includes
default parameter values. The instantiate method is used to override the default values
of any of the model parameters or plankton diameters.

!!! tip

    Changing the parameter values of an existing DARWIN model type using `instantiate()` is useful in
    dynamic programming contexts such as `for` loops.


# Arguments
- `bgc_type`: subtype of Oceananigans.Biogeochemistry returned by `DARWIN.construct()`
   with a specified number of phytoplankton and zooplankton

# Keywords
- `FT`: floating point type used when (re-)expanding parameters
- `phyto_diameters`: diameter specification (e.g., `DiameterRangeSpecification`) from which phytoplankton diameters can be
    computed or a list of values to use (as many as the model expects)
- `zoo_diameters`: diameter specification (e.g., `DiameterRangeSpecification`) from which zooplankton diameters can be
    computed or a list of values to use (as many as the model expects)
- `phyto_pft_parameters`: phytoplankton PFT parameters (`DarwinPhytoPFTParameters`), for default values see
    `DARWIN.default_darwin_phyto_parameters(FT)`
- `zoo_pft_parameters`: zooplankton PFT parameters (`DarwinZooPFTParameters`), for default values see
    `DARWIN.default_darwin_zoo_parameters(FT)`
- `bgc_specification`: elemental cycling specification (`DarwinBiogeochemistrySpecification`), for default values see
    `DARWIN.default_darwin_bgc_specification(FT)`
- `palatability_matrix`: optional palatability matrix passed as an Array
- `assimilation_efficiency_matrix`: optional assimilation efficiency matrix passed as an Array
- `sinking_tracers`: optional NamedTuple of sinking speeds (passed as positive values) of
   the form (<tracer name expressed as symbol> = <speed>, ...)
- `grid`: optional Oceananigans grid object defining the geometry to build the model on, must
   be passed if `sinking_tracers` is defined, defaults to BoxModelGrid
- `open_bottom`: indicates whether the sinking velocity should be smoothly brought to zero
   at the bottom to prevent the tracers leaving the domain, defaults to `true`, which means
   the bottom is open and the tracers leave (i.e., no slowing of velocity to 0 is applied)

# Example
```julia
using Agate.Models: DARWIN
using Oceananigans.Units: day

darwin_2p_2z = DARWIN.construct()

# change one parameter value without touching the rest
phyto_pft = DARWIN.default_darwin_phyto_parameters(Float64)
phyto_pft = with(phyto_pft; maximum_growth_rate_a = 2/day)

darwin_2p_2z_model_obj = DARWIN.instantiate(darwin_2p_2z; phyto_pft_parameters = phyto_pft)

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
        default_darwin_phyto_parameters(FT)
    else
        phyto_pft_parameters
    end

    zoo_pft = if isnothing(zoo_pft_parameters)
        default_darwin_zoo_parameters(FT)
    else
        zoo_pft_parameters
    end

    bgc_spec = if isnothing(bgc_specification)
        DarwinBiogeochemistrySpecification{FT}(
            FT(p0.POC_remineralization),
            FT(p0.DOC_remineralization),
            FT(p0.PON_remineralization),
            FT(p0.DON_remineralization),
            FT(p0.POP_remineralization),
            FT(p0.DOP_remineralization),
            FT(p0.DOM_POM_fractionation),
            FT(p0.nitrogen_to_carbon),
            FT(p0.phosphorus_to_carbon),
        )
    else
        bgc_specification
    end

    parameters = create_darwin_parameters(
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
