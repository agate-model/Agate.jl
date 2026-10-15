using ...Construction
import ...Integrations

function Integrations.npd_configuration(::FrankenLOBSTERFamily, _runtime)
    return (;
        nutrient_tracers=(:NO₃, :NH₄), consumed_detritus=(:DOM,), dependencies=(:T,),
    )
end

function _construction_inputs(;
    size_structure=DEFAULT_SIZE_STRUCTURE,
    parameters::NamedTuple=(;),
    settings::NamedTuple=(;),
    grid=nothing,
    arch=nothing,
    scalar_type=nothing,
    sinking_tracers=nothing,
    open_bottom::Bool=true,
)
    family = FrankenLOBSTERFamily()
    realization = (;
        plankton_pfts=Construction.plankton_realization(family, size_structure),
        parameter_overrides=parameters,
        setting_overrides=settings,
        sinking_tracers,
        open_bottom,
    )
    return (; family, realization, execution=(; grid, arch, scalar_type))
end

"""
    construct(; size_structure=DEFAULT_SIZE_STRUCTURE, parameters=(;), settings=(;),
                grid=nothing, arch=nothing, scalar_type=nothing,
                sinking_tracers=nothing, open_bottom=true) -> NPDPlankton

Construct a size-structured phytoplankton, zooplankton, and heterotrophic bacterioplankton
community for composition with OceanBioME's `LOBSTER` nutrient and detritus model.
The returned `Agate.Integrations.NPDPlankton` wraps the compiled Agate runtime, while
OceanBioME manages nutrient (`NO₃`, `NH₄`), and dissolved organic matter (`DOM`), 
and small and large particulate organic matter (`sPOM`, `bPOM`).

The light-attenuation model supplies `PAR`.

The default `size_structure` has PFTs `P` (two sizes: 0.6 and 1.2 μm), `Z` (two
sizes: 6 and 12 μm), and `H` (one size: 0.6 μm).

Keywords
--------
- `size_structure`: a NamedTuple with `phytoplankton`, `zooplankton`, and
  `bacterioplankton` PFT mappings. Each PFT accepts an explicit diameter vector or
  a range `(n, min_esd, max_esd, spacing)`.
- `parameters=(;)`: biological parameter overrides (e.g. the partial
  `(maximum_growth_rate=(P_1=1e-5,),)` override).
- `settings=(;)`: model-level conversion-factor overrides. The defaults are
  `chlorophyll_ratio=1.31` and `carbon_ratio=6.56`.
- `grid=nothing`: Oceananigans grid; required when configuring `sinking_tracers`.
- `arch=nothing`, `scalar_type=nothing`: optional runtime architecture and scalar type.
- `sinking_tracers=nothing`, `open_bottom=true`: sinking configuration.

Pass the result as `plankton` to `LOBSTER(grid; plankton, ...)`.
"""
function construct(; kwargs...)
    inputs = _construction_inputs(; kwargs...)
    return Integrations.construct_npd_plankton(
        inputs.family; inputs.realization..., inputs.execution...
    )
end

"""
    construct_plus_recipe(; kwargs...) -> plankton, recipe

Construct FrankenLOBSTER with the same keywords as `construct` and return its
versioned `ModelRecipe`.
The recipe records canonical PFT realization, biological parameter and model-setting
overrides, and sinking configuration. Replay supplies the grid, architecture, and scalar
type separately.
"""
function construct_plus_recipe(; kwargs...)
    inputs = _construction_inputs(; kwargs...)
    return Integrations.construct_npd_plankton_plus_recipe(
        inputs.family; inputs.realization..., inputs.execution...
    )
end

"""
    construct(recipe::Construction.ModelRecipe; grid=nothing, arch=nothing,
              scalar_type=nothing) -> NPDPlankton

Replay a FrankenLOBSTER recipe into an OceanBioME-compatible plankton component.
The recipe family and scientific definition version are verified before construction.
Provide a `grid` when the recipe configures sinking tracers.
"""
function construct(
    recipe::Construction.ModelRecipe; grid=nothing, arch=nothing, scalar_type=nothing
)
    recipe.family == :FrankenLOBSTER || throw(ArgumentError(
        "FrankenLOBSTER.construct requires a FrankenLOBSTER recipe; got $(recipe.family)"
    ))
    return Integrations.construct_npd_plankton(recipe; grid, arch, scalar_type)
end
