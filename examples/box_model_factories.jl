# # Box model factories
#
# This example demonstrates Agate's **factory-based API** in the simplest setting: a 0D box model.
#
# You will see how to:
#
# 1. Construct a default model from a factory.
# 2. Change community structure (number of size classes, diameters).
# 3. Override parameters (PFT and biogeochemical specification fields).
# 4. Swap components (dynamics functions).
# 5. Provide explicit interaction matrices.
#
# The same pattern applies to all factories (e.g., [`DarwinFactory`](@ref)).

using Agate
using Agate.Constructor: construct
using Agate.Models: NiPiZDFactory
using Agate.Constructor: PFTSpecification, BiogeochemistrySpecification, update_plankton_args, update_biogeochem_args, update_dynamics
using Agate.Utils: parse_community
using Agate.Library.Light

using OceanBioME
using OceanBioME: Biogeochemistry

using Oceananigans
using Oceananigans.Units

using CairoMakie

const year = years = 365days
nothing #hide

# ## 1. Start from a factory

factory = NiPiZDFactory()

# `construct` returns a **concrete** biogeochemistry type (a subtype of
# `AbstractContinuousFormBiogeochemistry`). You then instantiate it as usual.

bgc_type_default = construct(factory)
bgc_default = bgc_type_default()


# ## 2. Pull defaults and override community structure

# Factories expose defaults for:
#
# - `plankton_dynamics`: how each plankton group is translated into tracer tendencies
# - `plankton_args`: how many size classes, what diameters, and which PFT parameter sets
# - `biogeochem_dynamics`: non-plankton tracer tendencies (e.g., nutrients and detritus)
# - `biogeochem_args`: non-plankton model parameters (specification)

plankton_dynamics = Agate.Models.default_plankton_dynamics(factory)
plankton_args = Agate.Models.default_plankton_args(factory)
biogeochem_dynamics = Agate.Models.default_biogeochem_dynamics(factory)
biogeochem_args = Agate.Models.default_biogeochem_args(factory)

# `plankton_args` is a `NamedTuple` keyed by group prefix (e.g., `:Z`, `:P`).
# Each group specification must include a diameter specification (`diameters`) and, for
# non-explicit diameter specifications, a positive integer `n`.
#
# Here we:
# - reduce zooplankton to 1 size class at a fixed diameter
# - increase phytoplankton to 3 size classes with log-spaced diameters

plankton_args_custom = update_plankton_args(plankton_args, :Z; n=1, diameters=[60.0])
plankton_args_custom = update_plankton_args(plankton_args_custom, :P; n=3, diameters=(1.5, 20.0, :log_splitting))

# ## 3. Override parameter values

# PFT parameters are held in `PFTSpecification`, which stores a `NamedTuple` internally.
# You can create an overridden copy using keyword splatting.

# Apply a one-step PFT override via update_plankton_args
plankton_args_custom = update_plankton_args(plankton_args_custom, :P; maximum_growth_rate_a=3.0 / day)

# The non-plankton specification is similarly stored in `BiogeochemistrySpecification`.
biogeochem_args_fast = update_biogeochem_args(biogeochem_args; detritus_remineralization=0.18 / day)

# ## 4. Swap components (dynamics)

# NiPiZD includes an optional Geider-style growth closure. We can swap in those dynamics
# functions, while keeping everything else unchanged.

using Agate.Models.NiPiZD.Tracers: nutrient_geider_light, phytoplankton_geider_light

plankton_dynamics_geider = update_dynamics(plankton_dynamics; P=phytoplankton_geider_light)
biogeochem_dynamics_geider = update_dynamics(biogeochem_dynamics; N=nutrient_geider_light)

# ## 5. Provide explicit interaction matrices

# You can override interactions either by supplying:
#
# - a `NamedTuple` of matrices, or
# - a callable that returns such a `NamedTuple`.
#
# The matrices must be `n_total × n_total`, where `n_total` is the total number of plankton tracers
# (all groups concatenated in the order of `plankton_args`).

community = parse_community(
    Float64,
    plankton_args_custom;
    plankton_dynamics=plankton_dynamics_geider,
    biogeochem_dynamics=biogeochem_dynamics_geider,
)

n = community.n_total
pal = zeros(Float64, n, n)
assim = zeros(Float64, n, n)

# For illustration, make the single zooplankton class graze all phytoplankton equally.
z_idx = findall(==(:Z), community.group_symbols)
p_idx = findall(==(:P), community.group_symbols)

for i in z_idx, j in p_idx
    pal[i, j] = 1.0
    assim[i, j] = 0.3
end

interactions_custom = (
    palatability_matrix=pal,
    assimilation_efficiency_matrix=assim,
)

# ## 6. Construct the customised model type

bgc_type_custom = construct(
    factory;
    plankton_dynamics=plankton_dynamics_geider,
    plankton_args=plankton_args_custom,
    biogeochem_dynamics=biogeochem_dynamics_geider,
    biogeochem_args=biogeochem_args_fast,
    interactions=interactions_custom,
    sinking_tracers=(D=2.0 / day,),
)

bgc_custom = bgc_type_custom()

# ## 7. Run a short box model simulation

light = FunctionFieldPAR(; grid=BoxModelGrid())
bgc_model = Biogeochemistry(bgc_custom; light_attenuation=light)

box = BoxModel(; biogeochemistry=bgc_model)

# With custom `plankton_args`, the tracer names (and count) can change. For this example we know
# the configuration is: `N`, `D`, `Z1`, `P1`, `P2`, `P3`.
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
