# Test assertions for scalar outputs
include("functions.jl")
# Test assertions for scalar outputs

# Sample input values for testing (Replace these with actual values or function calls)
N = 1.0
P1, P2, Z1, Z2 = 1.0, 1.0, 1.0, 1.0
PAR = 100.0
maximum_growth_rate = [7.190e-6, 2.216e-5, 0, 0]
nitrogen_half_saturation = [6.73e-3, 0.12, 0, 0]
detritus_remineralization = [0.0102]
holling_half_saturation = [0, 0, 5.0, 5.0]
linear_mortality = [8e-7, 8e-7, 8e-7, 8e-7]
quadratic_mortality = [0, 0, 1e-6, 1e-6]
maximum_predation_rate = [0, 0, 8.86e-5, 4.88e-5]
palatability = [
    0 0 0 0
    0 0 0 0
    1 0.3 0 0
    0.3 1 0 0
]
assimilation_efficiency = [
    0 0 0 0
    0 0 0 0
    0.3 0.3 0 0
    0 0.3 0 0
]
alpha = [0.1, 0.1, 0.1, 0.1]
plankton_index = 1

# 1. Test photosynthetic_growth
@assert isa(
    photosynthetic_growth(
        N, P1, PAR, maximum_growth_rate[1], nitrogen_half_saturation[1], alpha[1]
    ),
    Number,
) "photosynthetic_growth is returning a non-scalar!"

# 2. Test net_photosynthetic_growth
@assert isa(
    net_photosynthetic_growth(
        N, [P1, P2], PAR, maximum_growth_rate, nitrogen_half_saturation, alpha
    ),
    Number,
) "net_photosynthetic_growth is returning a non-scalar!"

# 3. Test predation_loss
@assert isa(
    predation_loss(
        P1, Z1, maximum_predation_rate[3], holling_half_saturation[3], palatability[3, 1]
    ),
    Number,
) "predation_loss is returning a non-scalar!"

# 4. Test predation_gain
@assert isa(
    predation_gain(
        P1,
        Z1,
        assimilation_efficiency[3, 1],
        maximum_predation_rate[3],
        holling_half_saturation[3],
        palatability[3, 1],
    ),
    Number,
) "predation_gain is returning a non-scalar!"

# 5. Test summed_predation_loss
@assert isa(
    summed_predation_loss(
        1, [P1, P2, Z1, Z2], maximum_predation_rate, holling_half_saturation, palatability
    ),
    Number,
) "summed_predation_loss is returning a non-scalar!"

# 6. Test summed_predation_gain
@assert isa(
    summed_predation_gain(
        3,
        [P1, P2, Z1, Z2],
        assimilation_efficiency,
        maximum_predation_rate,
        holling_half_saturation,
        palatability,
    ),
    Number,
) "summed_predation_gain is returning a non-scalar!"

# 7. Test net_linear_loss
@assert isa(net_linear_loss([P1, P2, Z1, Z2], linear_mortality), Number) "net_linear_loss is returning a non-scalar!"

# 8. Test net_quadratic_loss
@assert isa(net_quadratic_loss([P1, P2, Z1, Z2], quadratic_mortality), Number) "net_quadratic_loss is returning a non-scalar!"

# 9. Test predation_assimilation_loss
@assert isa(
    predation_assimilation_loss(
        P1,
        Z1,
        assimilation_efficiency[3, 1],
        maximum_predation_rate[3],
        holling_half_saturation[3],
        palatability[3, 1],
    ),
    Number,
) "predation_assimilation_loss is returning a non-scalar!"

# 10. Test net_predation_assimilation_loss
@assert isa(
    net_predation_assimilation_loss(
        [P1, P2, Z1, Z2],
        holling_half_saturation,
        maximum_predation_rate,
        assimilation_efficiency,
        palatability,
    ),
    Number,
) "net_predation_assimilation_loss is returning a non-scalar!"

# 11. Test plankton_dt
@assert isa(
    plankton_dt(
        plankton_index,
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
    ),
    Number,
) "plankton_dt is returning a non-scalar!"

# Test menden_limitation
@assert isa(menden_limitation(N, nitrogen_half_saturation[1]), Number) "menden_limitation is returning a non-scalar!"

# Test smith_light_limitation
@assert isa(smith_light_limitation(PAR, alpha[1], maximum_growth_rate[1]), Number) "smith_light_limitation is returning a non-scalar!"

# Test holling_type_2
@assert isa(holling_type_2(P1, holling_half_saturation[1]), Number) "holling_type_2 is returning a non-scalar!"

# Test linear_loss
@assert isa(linear_loss(P1, linear_mortality[1]), Number) "linear_loss is returning a non-scalar!"

# Test quadratic_loss
@assert isa(quadratic_loss(P1, quadratic_mortality[1]), Number) "quadratic_loss is returning a non-scalar!"

# Test remineralization
@assert isa(remineralization(P1, detritus_remineralization[1]), Number) "remineralization is returning a non-scalar!"

# Test summed_predation_assimilation_loss
@assert isa(
    summed_predation_assimilation_loss(
        3,
        [P1, P2, Z1, Z2],
        assimilation_efficiency,
        maximum_predation_rate,
        holling_half_saturation,
        palatability,
    ),
    Number,
) "summed_predation_assimilation_loss is returning a non-scalar!"

# Test net_predation_assimilation_loss
@assert isa(
    net_predation_assimilation_loss(
        [P1, P2, Z1, Z2],
        holling_half_saturation,
        maximum_predation_rate,
        assimilation_efficiency,
        palatability,
    ),
    Number,
) "net_predation_assimilation_loss is returning a non-scalar!"

print("finished")
