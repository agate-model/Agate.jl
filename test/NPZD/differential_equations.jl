# This example shows how to integrate a Biogeochemistry (BGC) box model with DifferentialEquations.
using Agate.Library.Light

using DifferentialEquations
using Plots

using Oceananigans.Units
using Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers, required_biogeochemical_auxiliary_fields

const year = years = 365day

# ==================================================
# Define BGC model (NPZD)
# ==================================================

include("tracers.jl")
model = NPZD()

# ==================================================
# DifferentialEquations
# ==================================================

function model_ODEs(du, u, p, t)
    model = NPZD(p...)

    PAR = cyclical_PAR(-10, t)

    for (i, tracer) in enumerate(tracers)
        du[i] = model(Val(tracer), 0, 0, 0, t, u..., PAR)
    end

    return nothing
end

# make sure initial values are passed in right order (Z,P,N,D)
init_conditions = (N=7.0, P=0.01, Z=0.05, D=0.0)
tracers = required_biogeochemical_tracers(model)
u0 = [eval(:(init_conditions.$t)) for t in tracers]

# get model parameters
p = [getfield(model, f) for f in fieldnames(typeof(model))]

tspan = (0.0, 3years)

prob = ODEProblem(model_ODEs, u0, tspan, p)

sol = solve(prob, Tsit5())

# ==================================================
# Plotting
# ==================================================

# plot(sol)

# tracer order is Z,P,N,D
p = plot(sol.t, [sol.u[i][2] for i in range(1, length(sol))]; label="P")
plot!(p, sol.t, [sol.u[i][1] for i in range(1, length(sol))]; label="Z")
plot!(p, sol.t, [sol.u[i][4] for i in range(1, length(sol))]; label="D")
savefig(p, "NPZD_box_differential_equations.png")
using Agate
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units
using Plots
using Agate.Constructors: NiPiZD
using Agate.Models.Tracers
using Agate.Library.Photosynthesis

const year = years = 365day

# ==================================================
# Define BGC models
# ==================================================

# Default photosynthesis model
N2P2ZD_default_photosynthesis = NiPiZD.construct()
bgc_model_default_photosynthesis = Biogeochemistry(
    N2P2ZD_default_photosynthesis();
    light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid()),
)
full_model_default_photosynthesis = BoxModel(;
    biogeochemistry=bgc_model_default_photosynthesis
)
set!(full_model_default_photosynthesis; P1=0.01, P2=0.01, Z1=0.05, Z2=0.05, N=7.0, D=1)

# Geider photosynthesis model
N2P2ZD_geider = NiPiZD.construct(;
    phyto_args=NiPiZD.DEFAULT_PHYTO_GEIDER_ARGS,
    nutrient_dynamics=nutrients_geider_light,
    phyto_dynamics=phytoplankton_growth_single_nutrient_geider_light,
)
bgc_model_geider_photosynthesis = Biogeochemistry(
    N2P2ZD_geider(); light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid())
)
full_model_geider_photosynthesis = BoxModel(;
    biogeochemistry=bgc_model_geider_photosynthesis
)
set!(full_model_geider_photosynthesis; P1=0.01, P2=0.01, Z1=0.05, Z2=0.05, N=7.0, D=1)

# ==================================================
# Simulate models
# ==================================================

function run_simulation!(model, filename)
    simulation = Simulation(model; Δt=10minutes, stop_time=3years)
    simulation.output_writers[:fields] = JLD2OutputWriter(
        model,
        model.fields;
        filename=filename,
        schedule=TimeInterval(1day),
        overwrite_existing=true,
    )
    run!(simulation)
    return NamedTuple{keys(model.fields)}(
        FieldTimeSeries(filename, "$field")[1, 1, 1, :] for field in keys(model.fields)
    )
end

filename_default_photosynthesis = "box_n2p2zd_default_photosynthesis.jld2"
filename_geider_photosynthesis = "box_n2p2zd_geider_photosynthesis.jld2"

timeseries_default_photosynthesis = run_simulation!(
    full_model_default_photosynthesis, filename_default_photosynthesis
)
timeseries_geider_photosynthesis = run_simulation!(
    full_model_geider_photosynthesis, filename_geider_photosynthesis
)

# ==================================================
# Plotting
# ==================================================

fields = [:P1, :P2, :Z1, :Z2, :D, :N]
titles = [
    "Phytoplankton 1",
    "Phytoplankton 2",
    "Zooplankton 1",
    "Zooplankton 2",
    "Detritus",
    "Nutrient",
]

plots = [
    plot(; title=titles[i], xlabel="Time (days)", ylabel=string(fields[i])) for
    i in 1:length(fields)
]

for (i, field) in enumerate(fields)
    plot!(plots[i], timeseries_default_photosynthesis[field]; label="Default", color=:blue)
    plot!(plots[i], timeseries_geider_photosynthesis[field]; label="Geider", color=:green)
end

p = plot(plots...; layout=(2, 3), size=(900, 600))
savefig(p, "N2P2ZD_contrast_photosynthesis.png")
