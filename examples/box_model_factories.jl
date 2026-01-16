# # Box model factories
#
# This example demonstrates Agate's **factory-based API** in the simplest setting: a 0D box model.
#
# You will see how to:
#
# 1. Construct a default model from a factory.
# 2. Change community structure (number of size classes, diameters).
# 3. Override parameters (registry keywords, including allometry).
# 4. Swap components (dynamics functions).
# 5. Provide explicit interaction matrices.
#
# The same pattern applies to all factories (e.g., [`DarwinFactory`](@ref)).

using Agate
using Agate.Utils: parse_community
using Agate.Library.Light
using Agate.Library.Allometry: AllometricParam, PowerLaw
using Agate.Parameters: parameter_registry, update_registry

using OceanBioME
using OceanBioME: Biogeochemistry

using Oceananigans
using Oceananigans.Units

using CairoMakie

const year = years = 365days
nothing #hide

# ## 1. Start from a factory

factory = NiPiZDFactory()

# Inspect the factory's default parameter registry.
reg = parameter_registry(factory)
println(reg)

# `construct` returns the biogeochemistry instance directly.

bgc_default = construct(factory)

# ## 2. Pull defaults and override community structure

# Factories expose defaults for:
#
# - `plankton_dynamics`: how each plankton group is translated into tracer tendencies
# - `community`: how many size classes, what diameters, and optional per-PFT overrides
# - `biogeochem_dynamics`: non-plankton tracer tendencies (e.g., nutrients and detritus)
#
# Parameter defaults are **not** stored in these containers; they live in the model's
# parameter registry (single source of truth).

plankton_dynamics = Agate.Models.default_plankton_dynamics(factory)
community = Agate.Models.default_community(factory)
biogeochem_dynamics = Agate.Models.default_biogeochem_dynamics(factory)

# Reduce zooplankton to 1 size class at a fixed diameter.
# Increase phytoplankton to 3 size classes with log-spaced diameters.

community_custom = update_community(community, :Z; n=1, diameters=[60.0])
community_custom = update_community(community_custom, :P; n=3, diameters=(1.5, 20.0, :log_splitting))

# ## 3. Override parameter values (registry keywords)

# Parameter overrides are applied by updating the parameter registry.
# Vector parameters can be overridden with a per-group mapping (`(P=..., Z=...)`).

parameter_overrides = (
    detritus_remineralization = 0.18 / day,
    maximum_growth_rate = (P = AllometricParam(PowerLaw(); prefactor=3.0 / day, exponent=-0.15), Z = 0.0),
)

registry_custom = update_registry(parameter_registry(factory); parameter_overrides...)

# Inspect the updated registry after applying overrides.
println(registry_custom)

# ## 4. Swap components (dynamics)

# NiPiZD includes an optional Geider-style growth closure. We can swap in those dynamics
# functions, while keeping everything else unchanged.

using Agate.Models.NiPiZD.Tracers: nutrient_geider_light, phytoplankton_geider_light

plankton_dynamics_geider = update_dynamics(plankton_dynamics; P=phytoplankton_geider_light)
biogeochem_dynamics_geider = update_dynamics(biogeochem_dynamics; N=nutrient_geider_light)

# ## 5. Provide explicit interaction matrices

community_ctx = parse_community(
    Float64,
    community_custom;
    plankton_dynamics=plankton_dynamics_geider,
    biogeochem_dynamics=biogeochem_dynamics_geider,
)

n = community_ctx.n_total
pal = zeros(Float64, n, n)
assim = zeros(Float64, n, n)

# For illustration, make the single zooplankton class graze all phytoplankton equally.
z_idx = findall(==(:Z), community_ctx.group_symbols)
p_idx = findall(==(:P), community_ctx.group_symbols)

for i in z_idx, j in p_idx
    pal[i, j] = 1.0
    assim[i, j] = 0.3
end

interactions_custom = (
    palatability_matrix=pal,
    assimilation_matrix=assim,
)

# ## 6. Construct the customised model

bgc_custom = construct(
    factory;
    plankton_dynamics=plankton_dynamics_geider,
    biogeochem_dynamics=biogeochem_dynamics_geider,
    community=community_custom,
    registry=registry_custom,
    interactions=interactions_custom,
    sinking_tracers=(D=2.0 / day,),
)

# ## 7. Run a short box model simulation

light = FunctionFieldPAR(; grid=BoxModelGrid())
bgc_model = Biogeochemistry(bgc_custom; light_attenuation=light)

box = BoxModel(; biogeochemistry=bgc_model)

# With custom `community`, the tracer names (and count) can change.
println(tracer_names(bgc_custom))

set!(box; N=7.0, D=0.05, Z1=0.02, P1=0.01, P2=0.01, P3=0.01)

filename = "box_factories.jld2"

sim = Simulation(box; Δt=30minutes, stop_time=90days)
sim.output_writers[:fields] = JLD2Writer(
    box,
    box.fields;
    filename=filename,
    schedule=TimeInterval(1day),
    overwrite_existing=true,
)

run!(sim)

# ## 8. Plot a few tracers

times = FieldTimeSeries(filename, "N").times ./ days

P1 = FieldTimeSeries(filename, "P1")[1, 1, 1, :]
P2 = FieldTimeSeries(filename, "P2")[1, 1, 1, :]
P3 = FieldTimeSeries(filename, "P3")[1, 1, 1, :]
Z1 = FieldTimeSeries(filename, "Z1")[1, 1, 1, :]
N  = FieldTimeSeries(filename, "N")[1, 1, 1, :]

fig = Figure(; size=(900, 600), fontsize=18)

ax1 = Axis(fig[1, 1]; title="Plankton", xlabel="time (days)", ylabel="mmol N / m³")
lines!(ax1, times, P1; label="P1")
lines!(ax1, times, P2; label="P2")
lines!(ax1, times, P3; label="P3")
lines!(ax1, times, Z1; label="Z1")
axislegend(ax1; position=:rb)

ax2 = Axis(fig[2, 1]; title="Nutrient", xlabel="time (days)", ylabel="mmol N / m³")
lines!(ax2, times, N; label="N")

fig
