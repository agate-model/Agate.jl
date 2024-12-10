module Tracers

using NamedArrays

export typical_detritus,
    typical_nutrients, simplified_phytoplankton_growth, simplified_zooplankton_growth

function typical_nutrients(plankton_array)
    return :(
        net_linear_loss(
            [$(plankton_array...)], linear_mortality, mortality_export_fraction
        ) +
        net_quadratic_loss(
            [$(plankton_array...)], quadratic_mortality, mortality_export_fraction
        ) +
        remineralization(D, detritus_remineralization) - net_photosynthetic_growth(
            N,
            [$(plankton_array...)],
            PAR,
            maximum_growth_rate,
            nitrogen_half_saturation,
            alpha,
        )
    )
end

function typical_detritus(plankton_array)
    return :(
        net_linear_loss(
            [$(plankton_array...)], linear_mortality, 1 - mortality_export_fraction
        ) +
        net_predation_assimilation_loss(
            [$(plankton_array...)],
            holling_half_saturation,
            maximum_predation_rate,
            assimilation_efficiency,
            palatability,
        ) +
        net_quadratic_loss(
            [$(plankton_array...)], quadratic_mortality, 1 - mortality_export_fraction
        ) - remineralization(D, detritus_remineralization)
    )
end

function simplified_phytoplankton_growth(plankton_array, name)
    return :(phytoplankton_dt(
        $name,
        N,
        NamedArray(
            [$(plankton_array...)],
            String.($(plankton_array)),
            PAR,
            linear_mortality,
            quadratic_mortality,
            maximum_growth_rate,
            holling_half_saturation,
            nitrogen_half_saturation,
            alpha,
            maximum_predation_rate,
            palatability,
        ),
    ))
end

function simplified_zooplankton_growth(plankton_array, name)
    return :(zooplankton_dt(
        $name,
        NamedArray(
            [$(plankton_array...)],
            String.($(plankton_array)),
            linear_mortality,
            quadratic_mortality,
            holling_half_saturation,
            maximum_predation_rate,
            assimilation_efficiency,
            palatability,
        ),
    ))
end

end # module
