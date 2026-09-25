using Adapt: adapt
import Adapt: adapt_structure

using ...Integrations: NPDPlankton, phytoplankton_tracers, process_tendency

using OceanBioME.Models.NutrientsPlanktonDetritusModels: NutrientsPlanktonDetritus
import OceanBioME.Models.NutrientsPlanktonDetritusModels.InorganicCarbonModels:
    biological_calcium_carbonate_dissolution,
    biological_calcium_carbonate_precipitation,
    particulate_calcium_carbonate_production

"""FrankenLOBSTER-specific diagnostics retained by the generic NPD integration wrapper."""
struct FrankenLOBSTERCoupling{Diagnostics}
    process_diagnostics::Diagnostics
end

@inline adapt_structure(to, coupling::FrankenLOBSTERCoupling) =
    FrankenLOBSTERCoupling(adapt(to, coupling.process_diagnostics))

const FrankenLOBSTERPlankton = NPDPlankton{<:FrankenLOBSTERCoupling}

@inline function _phytoplankton_process_tendency(
    plankton::FrankenLOBSTERPlankton,
    process::Val,
    i, j, k, grid, fields, auxiliary_fields,
)
    return process_tendency(
        plankton,
        plankton.coupling.process_diagnostics,
        phytoplankton_tracers(plankton),
        process,
        i, j, k, grid, fields, auxiliary_fields,
    )
end

@inline _calcite_scale(plankton) =
    plankton.traits.calcium_carbonate_rain_ratio * plankton.traits.carbon_ratio

@inline biological_calcium_carbonate_precipitation(
    i, j, k, grid, plankton::FrankenLOBSTERPlankton,
    ::NutrientsPlanktonDetritus, fields, auxiliary_fields,
) = _calcite_scale(plankton) * (
    _phytoplankton_process_tendency(
        plankton, Val(:nitrate_growth_P), i, j, k, grid, fields, auxiliary_fields
    ) + _phytoplankton_process_tendency(
        plankton, Val(:ammonia_growth_P), i, j, k, grid, fields, auxiliary_fields
    )
)

@inline function particulate_calcium_carbonate_production(
    i, j, k, grid, plankton::FrankenLOBSTERPlankton,
    ::NutrientsPlanktonDetritus, fields, auxiliary_fields,
)
    grazing = -_phytoplankton_process_tendency(
        plankton, Val(:grazing_Z_on_living), i, j, k, grid, fields, auxiliary_fields
    )
    mortality = -_phytoplankton_process_tendency(
        plankton, Val(:mortality_P), i, j, k, grid, fields, auxiliary_fields
    )
    dissolved = plankton.traits.zooplankton_calcium_carbonate_dissolution
    return _calcite_scale(plankton) * ((one(dissolved) - dissolved) * grazing + mortality)
end

@inline biological_calcium_carbonate_dissolution(
    i, j, k, grid, plankton::FrankenLOBSTERPlankton,
    ::NutrientsPlanktonDetritus, fields, auxiliary_fields,
) = _calcite_scale(plankton) * plankton.traits.zooplankton_calcium_carbonate_dissolution *
    -_phytoplankton_process_tendency(
        plankton, Val(:grazing_Z_on_living), i, j, k, grid, fields, auxiliary_fields
    )
