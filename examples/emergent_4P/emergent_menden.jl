function dummy_emergent_nitrogen_half_saturation(volume_a, volume_b, volume)
    
    if volume_a == 0
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

println(dummy_emergent_nitrogen_half_saturation(1, 1, 1)) #P1
println(dummy_emergent_nitrogen_half_saturation(1, 1, 10)) #P2
println(dummy_emergent_nitrogen_half_saturation(0, 0, 10)) #Z1
println(dummy_emergent_nitrogen_half_saturation(0, 0, 100)) #Z2