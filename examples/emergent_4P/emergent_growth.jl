function dummy_emergent_growth(volume_a, volume_b, volume)
    if volume_a == 0
        return rate = 0  # Early return if volume_a is zero
    else
        # Set rate based on the value of volume
        if volume == 1
            rate = 7.190e-6
        elseif volume == 10
            rate = 2.216e-5
        end
    end

    return rate
end

println(dummy_emergent_growth(1, 1, 1)) #P1
println(dummy_emergent_growth(1, 1, 10)) #P2
println(dummy_emergent_growth(0, 0, 10)) #Z1
println(dummy_emergent_growth(0, 0, 100)) #Z2
