# # [Comparing phytoplankton light strategies] (@id named_plankton_groups_example)

# Plankton are highly diverse, and representing some of this diversity can be
# important in ecosystem models. One important axis of variation is adaptation to
# the light environment, which can differ both between and within species.
# High- and low-light-adapted *Prochlorococcus* ecotypes provide a well-known
# example of this ecological differentiation.
#
# Here, we use this idea to construct two idealized phytoplankton functional types
# with the same cell size but different light-response traits. This allows us to
# isolate the effect of light adaptation and explore how the two functional types
# respond under low- and high-light conditions.
#
# The `low_light` type has a lower maximum growth rate but a steeper initial light-response
# slope than `high_light`. We place the same community in two otherwise identical box models,
# differing only in photosynthetically active radiation (PAR), and compare which type is more
# competitive in each light environment.

# ## Loading dependencies
# The example uses Agate.jl, Oceananigans.jl, and OceanBioME.jl for the ecosystem simulations.
# CairoMakie.jl is used for plotting.

using Agate
using Agate.Introspection: tracer_names
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Fields: ConstantField
using Oceananigans.Units
using CairoMakie

nothing #hide

# ## Ecosystem model

# Named groups are defined separately for the phytoplankton and zooplankton roles. Both
# phytoplankton types use the same diameter so their different behavior comes from the trait
# overrides rather than from size-dependent defaults.

size_structure = (;
    phytoplankton=(high_light=[0.7], low_light=[0.7]),
    zooplankton=(grazer=[30.0],),
)
nothing #hide

# The two phytoplankton types are given a simple growth trade-off. `high_light` has the larger
# maximum growth rate, while `low_light` has the larger `alpha`, the initial slope of the Smith
# photosynthesis-irradiance curve. The parameter values are illustrative and chosen so that
# `low_light` grows faster at low PAR and `high_light` grows faster at high PAR.

parameters = (;
    maximum_growth_rate=(high_light_1=2.0 / day, low_light_1=1.4 / day),
    alpha=(high_light_1=0.15 / day, low_light_1=0.30 / day),
)

bgc = Agate.Models.NiPiZD.construct(; size_structure, parameters)
nothing #hide

# The named groups are realized as deterministic tracer names.

println(tracer_names(bgc))

# ## Light treatments

# We compare a low-light treatment at PAR = 5 with a high-light treatment at PAR = 100.
# Every other part of the ecosystem setup is shared between the two simulations.

low_PAR = 5.0
high_PAR = 100.0

# ## Run the box models

function run_light_treatment(bgc, PAR, filename)
    grid = BoxModelGrid()
    light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(PAR))
    bgc_model = Biogeochemistry(bgc; light_attenuation)
    full_model = BoxModel(; grid, biogeochemistry=bgc_model)

    set!(
        full_model;
        N=20.0,
        D=0.0,
        grazer_1=0.0,
        high_light_1=0.01,
        low_light_1=0.01,
    )

    simulation = Simulation(full_model; Δt=60minutes, stop_time=3days)

    simulation.output_writers[:fields] = JLD2Writer(
        full_model,
        full_model.fields;
        filename,
        schedule=TimeInterval(1hour),
        overwrite_existing=true,
    )

    run!(simulation)

    return filename
end

low_filename = run_light_treatment(bgc, low_PAR, "named_groups_low_light.jld2")
high_filename = run_light_treatment(bgc, high_PAR, "named_groups_high_light.jld2")

nothing #hide

# ## Load the output time series

function phytoplankton_timeseries(filename)
    high_light = FieldTimeSeries(filename, "high_light_1")[1, 1, 1, :]
    low_light = FieldTimeSeries(filename, "low_light_1")[1, 1, 1, :]
    times = FieldTimeSeries(filename, "high_light_1").times
    return (; times, high_light, low_light)
end

low = phytoplankton_timeseries(low_filename)
high = phytoplankton_timeseries(high_filename)
nothing #hide

# ## Compare the light treatments

# Both simulations start with equal phytoplankton biomass. Under low PAR, the steeper initial
# light-response slope favors `low_light`. Under high PAR, the larger maximum growth rate favors
# `high_light`.

fig = Figure(; size=(900, 420), fontsize=14)

for (column, treatment, PAR) in ((1, low, low_PAR), (2, high, high_PAR))
    ax = Axis(
        fig[1, column];
        xlabel="Time (days)",
        ylabel="Concentration (mmol N / m³)",
        title="PAR = $(PAR)",
    )

    lines!(
        ax,
        treatment.times / day,
        treatment.high_light;
        label="high_light",
        linewidth=2,
    )
    lines!(
        ax,
        treatment.times / day,
        treatment.low_light;
        label="low_light",
        linewidth=2,
    )
    axislegend(ax; position=:lt)
end

Label(fig[0, 1:2], "Contrasting responses to low and high light"; fontsize=20)

save("named_plankton_groups.png", fig; px_per_unit=1)

fig
