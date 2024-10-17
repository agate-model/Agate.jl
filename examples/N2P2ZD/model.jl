using Agate
using Oceananigans.Units

parameters = (
    maximum_growth_rate=[7.190e-6, 2.216e-5, 0, 0] / second,
    nitrogen_half_saturation=[6.73e-3, 0.12, 0, 0],
    detritus_remineralization=0.0102 / day,
    holling_half_saturation=[0, 0, 5.0, 5.0],
    linear_mortality=[8e-7, 8e-7, 8e-7, 8e-7] / second,
    quadratic_mortality=[0, 0, 1e-6, 1e-6] / second,
    maximum_predation_rate=[0, 0, 8.86e-5, 4.88e-5] / second,
    alpha=[0.1953, 0.1953, 0, 0],
    palatability=[
        0 0 0 0 #P1 
        0 0 0 0 #P2
        1 1 0 0 #Z1
        0 1 0 0
    ], #Z2 
    assimilation_efficiency=[
        0 0 0 0
        0 0 0 0
        0.68 0.68 0 0
        0 0.68 0 0
    ],
)
tracers = Dict(
    "N" => :(
        net_linear_loss([P1, P2, Z1, Z2], linear_mortality) +
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
        net_linear_loss([P1, P2, Z1, Z2], linear_mortality) +
        net_predation_assimilation_loss(
            [P1, P2, Z1, Z2],
            holling_half_saturation,
            maximum_predation_rate,
            assimilation_efficiency,
            palatability,
        ) +
        net_quadratic_loss([P1, P2, Z1, Z2], quadratic_mortality) -
        remineralization(D, detritus_remineralization)
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
aux_field_vars = [:PAR]
N2P2ZD = create_bgc_struct(:N2P2ZD, parameters)
add_bgc_methods(
    N2P2ZD, tracers; auxiliary_fields=aux_field_vars, helper_functions="functions.jl"
)
