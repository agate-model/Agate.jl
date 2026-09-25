using ...Construction
import ...Integrations

const _CALCITE_DIAGNOSTIC_PROCESSES = (
    :nitrate_growth_P, :ammonia_growth_P, :grazing_Z_on_living, :mortality_P,
)

Integrations.npd_diagnostic_processes(::FrankenLOBSTERFamily) = _CALCITE_DIAGNOSTIC_PROCESSES

function Integrations.npd_configuration(
    ::FrankenLOBSTERFamily, runtime, settings::NamedTuple
)
    return (;
        owned_components=(:P, :Z, :H),
        phytoplankton_components=(:P,),
        nutrient_tracers=(:NO₃, :NH₄),
        exchange_tracers=(
            solid=:solid_waste, dissolved=:dissolved_waste, inorganic=:inorganic_waste,
        ),
        consumed_detritus=(:DOM,),
        dependencies=(:NO₃, :NH₄, :DOM, :T),
        traits=settings,
        coupling=FrankenLOBSTERCoupling(runtime.metadata.process_diagnostics),
    )
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
