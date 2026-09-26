using ...Construction
import ...Integrations

function Integrations.npd_configuration(::FrankenLOBSTERFamily, _runtime)
    return (; nutrient_tracers=(:NO₃, :NH₄), consumed_detritus=(:DOM,))
end

function _construction_inputs(;
    size_structure=DEFAULT_SIZE_STRUCTURE,
    parameters::NamedTuple=(;),
    settings::NamedTuple=(;),
    grid=nothing,
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
    return (; family, realization, execution=(; grid))
end

"""Construct the Agate plankton component for composition with OceanBioME `LOBSTER`."""
function construct(; kwargs...)
    inputs = _construction_inputs(; kwargs...)
    return Integrations.construct_npd_plankton(
        inputs.family; inputs.realization..., inputs.execution...
    )
end

"""Construct FrankenLOBSTER plankton and capture its versioned Agate recipe."""
function construct_plus_recipe(; kwargs...)
    inputs = _construction_inputs(; kwargs...)
    return Integrations.construct_npd_plankton_plus_recipe(
        inputs.family; inputs.realization..., inputs.execution...
    )
end

"""Replay a FrankenLOBSTER recipe into the Agate plankton component."""
function construct(recipe::Construction.ModelRecipe; grid=nothing)
    recipe.family == :FrankenLOBSTER || throw(ArgumentError(
        "FrankenLOBSTER.construct requires a FrankenLOBSTER recipe; got $(recipe.family)"
    ))
    return Integrations.construct_npd_plankton(recipe; grid)
end
