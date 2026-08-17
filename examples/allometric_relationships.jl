# # [Allometric parameters] (@id allometric_relationships_example)

# Many plankton traits depend on organism size (allometry). Agate uses this to calculate the parameter 
# values of relevant parameters for each plankton size in the model. 
#
# We will first look at the  default [Agate.jl-NiPiZD](@ref NiPiZD) relationship for maximum growth rate. We change
# only its exponent, compare the growth rates that Agate assigns to each phytoplankton
# size, and then run the two versions in the same simple box model. We then replace the
# built-in relationship with one based on [Ward et al. (2017)](https://doi.org/10.1086/689992)
# and repeat the same steps.

# !!! info
#     The allometric relationships are defined with Agate. The short simulations later in
#     the example use [Oceananigans.jl](https://clima.github.io/OceananigansDocumentation/stable/)
#     and [OceanBioME.jl](https://oceanbiome.github.io/OceanBioME.jl/stable/) to run them in
#     a box model.

# ## Loading dependencies

using Agate.Models: NiPiZD
using Agate.Introspection: plankton_groups, plankton_tracers
using Agate.Library.Allometry: AllometricParam, PowerLaw
using Agate.Library.Light: FunctionFieldPAR
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units
using CairoMakie

# ## Choose the plankton sizes

# We use five phytoplankton sizes so that the effect of size on growth is easy to see.
# The numbers below are cell diameters in μm. Here `P` is the phytoplankton group and `Z`
# is the zooplankton group. The chosen sizes also span the rise and fall in growth rate
# described by Ward et al. later in the example. In the box-model experiments, zooplankton
# starts at zero so that the first few days mainly show phytoplankton growth.

size_structure = (;
    phytoplankton=(P=[0.8, 2.0, 6.0, 20.0, 60.0],),
    zooplankton=(Z=[100.0],),
)
nothing #hide


# NiPiZD already uses Agate's built-in `PowerLaw` for `maximum_growth_rate`. A power law
# has the form `prefactor × V^exponent`, where `V` is spherical cell volume calculated
# from diameter. The prefactor sets the overall scale, while the exponent controls how
# strongly the value changes with size.
# NiPiZD's defaults are `2 / day` and `-0.15`. We write them explicitly here so that the
# change in the next section is easy to follow. `AllometricParam` combines the relationship
# with its coefficients so that it can be supplied as a model parameter.

default_growth = AllometricParam(
    PowerLaw();
    prefactor=2 / day,
    exponent=-0.15,
)
nothing #hide

# Pass the relationship through `parameters` when constructing NiPiZD. Agate then
# calculates the corresponding `maximum_growth_rate` for each plankton size.

default_bgc = NiPiZD.construct(;
    size_structure,
    parameters=(; maximum_growth_rate=default_growth),
)
nothing #hide


# The model contains one maximum-growth value for each plankton class. Select the
# phytoplankton values so that we can plot them against the phytoplankton sizes.

phytoplankton_tracers = plankton_groups(default_bgc).P
phytoplankton_diameters = size_structure.phytoplankton.P
phytoplankton_indices = findall(in(phytoplankton_tracers), plankton_tracers(default_bgc))
default_growth_rates = default_bgc.parameters.maximum_growth_rate[phytoplankton_indices]
nothing #hide

# ## Change the built-in size relationship

# An exponent closer to zero makes the decrease in growth rate with size weaker. Keep
# the default prefactor and change only the exponent from `-0.15` to `-0.05`. When the
# model is constructed, Agate evaluates this relationship for every phytoplankton size.

shallower_growth = AllometricParam(
    PowerLaw();
    prefactor=2 / day,
    exponent=-0.05,
)

shallower_bgc = NiPiZD.construct(;
    size_structure,
    parameters=(; maximum_growth_rate=shallower_growth),
)
nothing #hide
shallower_growth_rates = shallower_bgc.parameters.maximum_growth_rate[phytoplankton_indices]
nothing #hide

# ### Compare the resulting growth rates

# Before running a simulation, compare the values that the two relationships give each
# phytoplankton size:
#
# `relationship + coefficients → construct → one parameter value for each size class`.

slope_fig = Figure(; size=(700, 440), fontsize=14)
slope_ax = Axis(
    slope_fig[1, 1];
    xlabel="Phytoplankton diameter (μm)",
    ylabel="Maximum growth rate (day⁻¹)",
    xscale=log10,
    title="Changing the PowerLaw exponent",
)

scatterlines!(
    slope_ax,
    phytoplankton_diameters,
    default_growth_rates .* day;
    label="default exponent = -0.15",
    linewidth=2,
)
scatterlines!(
    slope_ax,
    phytoplankton_diameters,
    shallower_growth_rates .* day;
    label="exponent = -0.05",
    linewidth=2,
)
axislegend(slope_ax; position=:rt)

slope_fig

# ### See the effect in a box model

# Now run the default and modified relationships in the same well-mixed box model. The
# two simulations use the same light, nutrients, starting phytoplankton biomass, and
# initial zooplankton. This keeps the comparison focused on the change in maximum growth
# rate with size.

# #### Run the simulations

# First, set up and run a short box model for each constructed Agate model. The helper
# below saves the model fields every hour so we can inspect them afterwards.

function run_growth_experiment(bgc, phytoplankton_tracers, filename)
    grid = BoxModelGrid()
    constant_PAR = t -> 100.0
    light_attenuation = FunctionFieldPAR(; grid, PAR_f=constant_PAR)
    bgc_model = Biogeochemistry(bgc; light_attenuation)
    model = BoxModel(; biogeochemistry=bgc_model)

    phytoplankton_initial_conditions = (;
        (tracer => 0.01 for tracer in phytoplankton_tracers)...
    )
    initial_conditions = merge(
        (; N=20.0, D=0.0, Z_1=0.0),
        phytoplankton_initial_conditions,
    )
    set!(model; initial_conditions...)

    simulation = Simulation(model; Δt=60minutes, stop_time=3days)
    simulation.output_writers[:fields] = JLD2Writer(
        model,
        model.fields;
        filename,
        schedule=TimeInterval(1hour),
        overwrite_existing=true,
    )

    run!(simulation)
    return filename
end

default_filename = run_growth_experiment(
    default_bgc,
    phytoplankton_tracers,
    "allometry_default_growth.jld2",
)
shallower_filename = run_growth_experiment(
    shallower_bgc,
    phytoplankton_tracers,
    "allometry_shallower_growth.jld2",
)

nothing #hide

# #### Load the results

# Next, load the saved phytoplankton biomass through time for each simulation.

function load_growth_experiment(filename, phytoplankton_tracers)
    series = (;
        (tracer => FieldTimeSeries(filename, string(tracer))[1, 1, 1, :] for tracer in phytoplankton_tracers)...
    )
    times = FieldTimeSeries(filename, string(first(phytoplankton_tracers))).times
    return times, series
end

times, default_series = load_growth_experiment(default_filename, phytoplankton_tracers)
_, shallower_series = load_growth_experiment(shallower_filename, phytoplankton_tracers)

nothing #hide

# #### Plot the comparison

# Finally, plot each phytoplankton size class so the effect of the changed allometric
# slope can be compared through time.

function plot_growth_comparison(
    times,
    reference_series,
    treatment_series,
    phytoplankton_tracers,
    phytoplankton_diameters;
    reference_label,
    treatment_label,
    title,
)
    fig = Figure(; size=(850, 760), fontsize=13)
    Label(fig[0, 1:2], title; fontsize=20)

    for (idx, (tracer, diameter)) in enumerate(zip(phytoplankton_tracers, phytoplankton_diameters))
        row = 1 + div(idx - 1, 2)
        col = 1 + mod(idx - 1, 2)
        ax = Axis(
            fig[row, col];
            xlabel="Time (days)",
            ylabel="Biomass (mmol N / m³)",
            title="$(tracer), $(diameter) μm",
        )

        lines!(
            ax,
            times / day,
            getproperty(reference_series, tracer);
            label=reference_label,
            linewidth=2,
        )
        lines!(
            ax,
            times / day,
            getproperty(treatment_series, tracer);
            label=treatment_label,
            linewidth=2,
            linestyle=:dash,
        )

        idx == 1 && axislegend(ax; position=:lt)
    end

    return fig
end

slope_box_fig = plot_growth_comparison(
    times,
    default_series,
    shallower_series,
    phytoplankton_tracers,
    phytoplankton_diameters;
    reference_label="default",
    treatment_label="shallower exponent",
    title="Box-model response to the PowerLaw slope",
)

slope_box_fig

# ## Use a custom relationship from Ward et al. (2017)

# A single `PowerLaw` is useful when a parameter changes steadily with size. Published
# relationships can also have other shapes. With `AllometricParam`, we can provide our
# own rule for calculating a parameter from diameter and then use it in the same way.
#
# Ward et al. (2017) describe maximum phytoplankton growth using three size-dependent
# quantities: maximum nitrogen uptake, ``\rho_{max}``; the minimum cellular nitrogen
# requirement, ``Q_{min}``; and the maximum metabolic growth rate, ``\mu_\infty``.
# Their equation (3) combines them as
#
# ```math
# \mu_{max} = \frac{\mu_\infty \rho_{max}}
#                   {\mu_\infty Q_{min} + \rho_{max}}.
# ```
#
# Ward et al. describe each of these quantities as a power law of cell volume, ``V``.
# Here we use the adjusted coefficients from their table 1:
#
# ```math
# \rho_{max} = 0.024 V^{1.10}, \qquad
# Q_{min} = 0.032 V^{0.76}, \qquad
# \mu_\infty = 4.7 V^{-0.26}.
# ```
#
# Here ``V`` is cell volume in μm³. ``\rho_{max}`` is in pg N cell⁻¹ day⁻¹,
# ``Q_{min}`` is in pg N cell⁻¹, and ``\mu_\infty`` is in day⁻¹. The nitrogen units
# cancel in equation (3), so the result is a growth rate. In the code below, rates given
# per day are divided by `day` so they use the same time units as the rest of the model.
#
# A custom relationship is a function of `(coeffs, diameter)`, in that order. Agate
# supplies each plankton diameter, while the keyword arguments passed to `AllometricParam`
# below are available inside the function through `coeffs`.

function ward_maximum_growth(coeffs, diameter)
    volume = (4 / 3) * pi * (diameter / 2)^3

    rho_max = coeffs.rho_prefactor * volume^coeffs.rho_exponent
    q_min = coeffs.q_min_prefactor * volume^coeffs.q_min_exponent
    mu_infinity = coeffs.mu_infinity_prefactor * volume^coeffs.mu_infinity_exponent

    return mu_infinity * rho_max / (mu_infinity * q_min + rho_max)
end

# Pair the Ward relationship with the coefficients from the paper.

ward_growth = AllometricParam(
    ward_maximum_growth;
    rho_prefactor=0.024 / day,
    rho_exponent=1.10,
    q_min_prefactor=0.032,
    q_min_exponent=0.76,
    mu_infinity_prefactor=4.7 / day,
    mu_infinity_exponent=-0.26,
)

ward_bgc = NiPiZD.construct(;
    size_structure,
    parameters=(; maximum_growth_rate=ward_growth),
)

ward_growth_rates = ward_bgc.parameters.maximum_growth_rate[phytoplankton_indices]
nothing #hide

# ### Compare the Ward relationship with the default

# Agate's default relationship decreases steadily as phytoplankton get larger. The Ward
# relationship instead increases from the smallest cells to an intermediate-size maximum,
# then decreases again for larger cells.

ward_fig = Figure(; size=(700, 440), fontsize=14)
ward_ax = Axis(
    ward_fig[1, 1];
    xlabel="Phytoplankton diameter (μm)",
    ylabel="Maximum growth rate (day⁻¹)",
    xscale=log10,
    title="Default and Ward et al. maximum growth",
)

scatterlines!(
    ward_ax,
    phytoplankton_diameters,
    default_growth_rates .* day;
    label="Agate default PowerLaw",
    linewidth=2,
)
scatterlines!(
    ward_ax,
    phytoplankton_diameters,
    ward_growth_rates .* day;
    label="Ward et al. (2017)",
    linewidth=2,
)
axislegend(ward_ax; position=:rt)

ward_fig

# ### See the effect in a box model

# Finally, put the Ward relationship through the same run, load, and plot workflow. The
# environmental conditions and starting values are unchanged, so differences come from
# the relationship used for `maximum_growth_rate`.

# #### Run the simulation

ward_filename = run_growth_experiment(
    ward_bgc,
    phytoplankton_tracers,
    "allometry_ward_growth.jld2",
)

nothing #hide

# #### Load the results

_, ward_series = load_growth_experiment(ward_filename, phytoplankton_tracers)

nothing #hide

# #### Plot the comparison

ward_box_fig = plot_growth_comparison(
    times,
    default_series,
    ward_series,
    phytoplankton_tracers,
    phytoplankton_diameters;
    reference_label="Agate default",
    treatment_label="Ward et al. (2017)",
    title="Box-model response to the custom growth relationship",
)

ward_box_fig
