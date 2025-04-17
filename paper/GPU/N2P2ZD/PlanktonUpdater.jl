using CUDA

function instantiate_updater(arch="gpu")
    if arch=="gpu"
        array = CUDA.zeros(Float64, 4)
        function update_gpu!(P1, P2, Z1, Z2)            
            array[[1]] = P1
            array[[2]] = P2
            array[[3]] = Z1
            array[[4]] = Z2
            return array
        end
        return update_gpu!
    elseif arch=="cpu"
        array = zeros(Float64, 4)
        function update_cpu!(P1, P2, Z1, Z2)
            array[1] = P1
            array[2] = P2
            array[3] = Z1
            array[4] = Z2
            return array
        end
        return update_cpu!
    end
end

# use
# arch = "gpu" 
# update_plankton! = instantiate_updater(arch)
# update_plankton!(P1, P2, Z1, Z2)
