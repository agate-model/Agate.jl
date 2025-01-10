
function dummy_emergent_growth(growth_a::Real, growth_b::Real, diameter::Real)
    rate = 0.0
    if growth_a == 0
        return rate
    else
        if diameter == 2
            rate = 9.05e-06
        elseif diameter == 10
            rate = 1.87e-05
        end
    end
    return rate
end

function dummy_emergent_predation_rate(
    predation_rate_a::Real, predation_rate_b::Real, diameter::Real
)
    rate = 0
    if predation_rate_a == 0
        return rate = 0  # Early return if diameter_a is zero
    else
        # Set rate based on the value of diameter
        if diameter == 20
            rate = 8.86e-5
        elseif diameter == 100
            rate = 4.88e-5
        end
    end

    return rate
end

function dummy_emergent_nitrogen_half_saturation(
    nitrogen_half_saturation_a::Real, nitrogen_half_saturation_b::Real, diameter::Real
)
    if nitrogen_half_saturation_a == 0
        return rate = 0  # Early return if diameter_a is zero
    else
        # Set rate based on the value of diameter
        if diameter == 2
            rate = 6.73e-3
        elseif diameter == 10
            rate = 0.12
        end
    end

    return rate
end

function dummy_emergent_palat(prey_data, predator_data)
    prey_diameter = prey_data["diameters"]
    predator_diameter = predator_data["diameters"]
    optimum_predator_prey_ratio = predator_data["optimum_predator_prey_ratio"]
    protection = prey_data["protection"]

    ratio = predator_diameter / prey_diameter

    if optimum_predator_prey_ratio == 0
        palat = 0.0
    elseif ratio == optimum_predator_prey_ratio
        palat = 1 * (1 - protection)
    else
        palat = 0.3 * (1 - protection)
    end
    return palat
end

function dummy_emergent_assimilation_efficiency(prey_data, predator_data)
    assimilation_efficiency = 0
    # predators don't eat other predators
    if prey_data["assimilation_efficiency"] == 0
        assimilation_efficiency = predator_data["assimilation_efficiency"]
    end
    return assimilation_efficiency
end
