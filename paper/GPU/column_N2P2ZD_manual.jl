using Agate
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields: FunctionField, ConstantField

const year = years = 365day

include(joinpath("N2P2ZD", "tracers_vectorized.jl"))

using Adapt, CUDA
CUDA.allowscalar(false)

Adapt.@adapt_structure N2P2ZD


adapted_instance = Adapt.adapt(CuArray, N2P2ZD())

biogeochemistry = Biogeochemistry(
    adapted_instance;
    #N2P2ZD();
    light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid(arch=GPU())),
)

#biogeochemistry = Adapt.adapt(CuArray, biogeochemistry)

const year = years = 365days
nothing #hide

@inline PAR⁰(x, y, t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

@inline H(t, t₀, t₁) = ifelse(t₀ < t < t₁, 1.0, 0.0)

@inline fmld1(t) = H(t, 50days, year) * (1 / (1 + exp(-(t - 100days) / 5days))) * (1 / (1 + exp((t - 330days) / 25days)))

@inline MLD(t) = - (10 + 340 * (1 - fmld1(year - eps(year)) * exp(-mod(t, year) / 25days) - fmld1(mod(t, year))))

@inline κₜ(x, y, z, t) = 1e-2 * (1 + tanh((z - MLD(t)) / 10)) / 2 + 1e-4

@inline temp(x, y, z, t) = 2.4 * cos(t * 2π / year + 50days) + 10


grid = RectilinearGrid(GPU(), size = (1, 1, 50), extent = (20meters, 20meters, 200meters))


clock = Clock(; time = 0.0)
T = FunctionField{Center, Center, Center}(temp, grid; clock)
S = ConstantField(35)

model = NonhydrostaticModel(; grid,
                              clock,
                              closure = ScalarDiffusivity(ν = κₜ, κ = κₜ),
                              biogeochemistry,
                              auxiliary_fields = (; T, S))

set!(model, P1=0.01, P2=0.01, Z1=0.05, Z2=0.05, N=7.0, D=1)

simulation = Simulation(model, Δt = 3minutes, stop_time = 100days)

filename = "column"
simulation.output_writers[:profiles] = JLD2OutputWriter(model, model.tracers,
                                                        filename = "$filename.jld2",
                                                        schedule = TimeInterval(1day),
                                                        overwrite_existing = true)


run!(simulation)
