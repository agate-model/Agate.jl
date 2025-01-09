using Agate
using Oceananigans.Units

parameters = (
    maximum_growth_rate=[9.05e-06, 1.87e-05, 0, 0] / second,
    nitrogen_half_saturation=[6.73e-3, 0.12, 0, 0],
    detritus_remineralization=0.1213 / day,
    holling_half_saturation=[0, 0, 5.0, 5.0],
    linear_mortality=[8e-7, 8e-7, 8e-7, 8e-7] / second,
    quadratic_mortality=[0, 0, 1e-6, 1e-6] / second,
    maximum_predation_rate=[0, 0, 8.86e-5, 4.88e-5] / second,
    alpha=[0.1953, 0.1953, 1e-99, 1e-99] / day,
    feeding_export_poc_doc_fraction=0.5,
    mortality_export_fraction=0.5,
    palatability=[
        0 0 0 0 #P1
        0 0 0 0 #P2
        1 0.3 0 0 #Z1
        0.3 1 0 0 #Z2
    ],
    assimilation_efficiency=[
        0 0 0 0
        0 0 0 0
        0.32 0.32 0 0
        0.32 0.32 0 0
    ],
)
tracers = Dict(
    "N" => :(
        custom_net_linear_loss(
            [P1, P2, Z1, Z2], linear_mortality, mortality_export_fraction
        ) +
        custom_net_quadratic_loss(
            [P1, P2, Z1, Z2], quadratic_mortality, mortality_export_fraction
        ) +
        remineralization(D, detritus_remineralization) - net_photosynthetic_growth(
            N,
            [P1, P2, Z1, Z2],
            PAR,
            maximum_growth_rate,
            nitrogen_half_saturation,
            alpha,
        )
    ),
    "D" => :(
        custom_net_linear_loss(
            [P1, P2, Z1, Z2], linear_mortality, 1 - mortality_export_fraction
        ) +
        net_predation_assimilation_loss(
            [P1, P2, Z1, Z2],
            holling_half_saturation,
            maximum_predation_rate,
            assimilation_efficiency,
            palatability,
        ) +
        custom_net_quadratic_loss(
            [P1, P2, Z1, Z2], quadratic_mortality, 1 - mortality_export_fraction
        ) - remineralization(D, detritus_remineralization)
    ),
    "P1" => :(plankton_dt(
        1,
        N,
        [P1, P2, Z1, Z2],
        PAR,
        linear_mortality,
        quadratic_mortality,
        maximum_growth_rate,
        holling_half_saturation,
        nitrogen_half_saturation,
        alpha,
        maximum_predation_rate,
        assimilation_efficiency,
        palatability,
    )),
    "P2" => :(plankton_dt(
        2,
        N,
        [P1, P2, Z1, Z2],
        PAR,
        linear_mortality,
        quadratic_mortality,
        maximum_growth_rate,
        holling_half_saturation,
        nitrogen_half_saturation,
        alpha,
        maximum_predation_rate,
        assimilation_efficiency,
        palatability,
    )),
    "Z1" => :(plankton_dt(
        3,
        N,
        [P1, P2, Z1, Z2],
        PAR,
        linear_mortality,
        quadratic_mortality,
        maximum_growth_rate,
        holling_half_saturation,
        nitrogen_half_saturation,
        alpha,
        maximum_predation_rate,
        assimilation_efficiency,
        palatability,
    )),
    "Z2" => :(plankton_dt(
        4,
        N,
        [P1, P2, Z1, Z2],
        PAR,
        linear_mortality,
        quadratic_mortality,
        maximum_growth_rate,
        holling_half_saturation,
        nitrogen_half_saturation,
        alpha,
        maximum_predation_rate,
        assimilation_efficiency,
        palatability,
    )),
)

# `linear_loss`, `quadratic_loss` etc are functions defined in Agate.Library
# import remaining "custom" functions from file
N2P2ZD = define_tracer_functions(parameters, tracers; helper_functions="functions.jl")
