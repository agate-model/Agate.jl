using OceanBioME, Oceananigans, Printf
using Oceananigans.Fields: FunctionField, ConstantField
using Oceananigans.Units

const year = years = 365days
nothing #hide

@inline PAR⁰(x, y, t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

@inline H(t, t₀, t₁) = ifelse(t₀ < t < t₁, 1.0, 0.0)

@inline fmld1(t) = H(t, 50days, year) * (1 / (1 + exp(-(t - 100days) / 5days))) * (1 / (1 + exp((t - 330days) / 25days)))

@inline MLD(t) = - (10 + 340 * (1 - fmld1(year - eps(year)) * exp(-mod(t, year) / 25days) - fmld1(mod(t, year))))

@inline κₜ(x, y, z, t) = 1e-2 * (1 + tanh((z - MLD(t)) / 10)) / 2 + 1e-4

@inline temp(x, y, z, t) = 2.4 * cos(t * 2π / year + 50days) + 10


grid = RectilinearGrid(GPU(), size = (1, 1, 50), extent = (20meters, 20meters, 200meters))

biogeochemistry = LOBSTER(; grid,
                            surface_photosynthetically_active_radiation = PAR⁰,
                            carbonates = true,
                            scale_negatives = true)


clock = Clock(; time = 0.0)
T = FunctionField{Center, Center, Center}(temp, grid; clock)
S = ConstantField(35)

model = NonhydrostaticModel(; grid,
                              clock,
                              closure = ScalarDiffusivity(ν = κₜ, κ = κₜ),
                              biogeochemistry,
                              auxiliary_fields = (; T, S))

set!(model, P = 0.03, Z = 0.03, NO₃ = 4.0, NH₄ = 0.05, DIC = 2239.8, Alk = 2409.0)

simulation = Simulation(model, Δt = 3minutes, stop_time = 100days)

filename = "column"
simulation.output_writers[:profiles] = JLD2OutputWriter(model, model.tracers,
                                                        filename = "$filename.jld2",
                                                        schedule = TimeInterval(1day),
                                                        overwrite_existing = true)


run!(simulation)
