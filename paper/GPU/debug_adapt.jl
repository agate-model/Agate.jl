using CUDA, Adapt

#this does not work:

struct N2P2ZD
    maximum_growth_rate::Vector{Float64}
    nutrient_half_saturation::Vector{Float64}
end

Adapt.@adapt_structure N2P2ZD

# Create an instance on CPU
model = N2P2ZD([1.0, 2.0, 3.0], [0.1, 0.2, 0.3])

# Move to GPU
model_gpu = adapt(CuArray, model)

# Check the adapted fields
println(typeof(model_gpu.maximum_growth_rate))  # Vector{Float64}
println(typeof(model_gpu.nutrient_half_saturation))  # Vector{Float64}


#this is fine:

# Define a parametric struct
struct parametric_N2P2ZD{T}
    maximum_growth_rate::T
    nutrient_half_saturation::T
end

Adapt.@adapt_structure parametric_N2P2ZD

# Create an instance on CPU
parametric_model = parametric_N2P2ZD([1.0, 2.0, 3.0], [0.1, 0.2, 0.3])

# Move to GPU
parametric_model_gpu = adapt(CuArray, parametric_model)

# Check the adapted fields
println(typeof(parametric_model_gpu.maximum_growth_rate))  # CuArray{Float64, 1}
println(typeof(parametric_model_gpu.nutrient_half_saturation))  # CuArray{Float64, 1}


