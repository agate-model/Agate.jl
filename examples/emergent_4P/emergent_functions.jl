function dummy_emergent_growth(growth_a::Real, growth_b::Real, volume::Real)
    rate = 0.0
    if growth_a == 0
        return rate
    else
        if volume == 1
            rate = 7.190e-6
        elseif volume == 10
            rate = 2.216e-5
        end
    end

    return rate
end

function dummy_emergent_palat(
    prey_data,
    predator_data;
    volume_key::String="volume",
    optimum_ratio_key::String="optimum_predator_prey_ratio",
    protection_key::String="protection",
)
    prey_volume = prey_data[volume_key]
    predator_volume = predator_data[volume_key]
    optimum_predator_prey_ratio = predator_data[optimum_ratio_key]
    protection = prey_data[protection_key]

    ratio = predator_volume / prey_volume

    if optimum_predator_prey_ratio == 0
        palat = 0.0
    elseif ratio == optimum_predator_prey_ratio
        palat = 1 * (1 - protection)
    else
        palat = 0.3 * (1 - protection)
    end
    return palat
end

function dummy_emergent_predation_rate(
    predation_rate_a::Real, predation_rate_b::Real, volume::Real
)
    rate = 0
    if predation_rate_a == 0
        return rate = 0  # Early return if volume_a is zero
    else
        # Set rate based on the value of volume
        if volume == 10
            rate = 8.86e-5
        elseif volume == 100
            rate = 4.88e-5
        end
    end

    return rate
end

function dummy_emergent_nitrogen_half_saturation(
    nitrogen_half_saturation_a::Real, nitrogen_half_saturation_b::Real, volume::Real
)
    if nitrogen_half_saturation_a == 0
        return rate = 0  # Early return if volume_a is zero
    else
        # Set rate based on the value of volume
        if volume == 1
            rate = 6.73e-3
        elseif volume == 10
            rate = 0.12
        end
    end

    return rate
end

"""
Default fall-back function if no emergent function is defined
"""
function default_emergent(my_parameter::Real)
    return my_parameter
end

function dummy_emergent_assimilation_efficiency(
    prey_data, predator_data; assimilation_efficiency_key::String="assimilation_efficiency"
)
    assimilation_efficiency = 0
    # predators don't eat other predators
    if prey_data[assimilation_efficiency_key] == 0
        assimilation_efficiency = predator_data[assimilation_efficiency_key]
    end
    return assimilation_efficiency
end
