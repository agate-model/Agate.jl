function dummy_emergent_predation_rate(volume_a, volume_b, volume)
    if volume_a == 0
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

println(dummy_emergent_predation_rate(0, 0, 1)) #P1
println(dummy_emergent_predation_rate(0, 0, 10)) #P2
println(dummy_emergent_predation_rate(1, 1, 10)) #Z1
println(dummy_emergent_predation_rate(1, 1, 100)) #Z2
