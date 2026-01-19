# # Box model customisation
#
# This example demonstrates customising an Agate model in the simplest setting: a 0D box model.
#
# You will see how to:
#
# 1. Construct a default model.
# 2. Change community structure (number of size classes, diameters).
# 3. Override parameters.
# 4. Swap components (dynamics functions).
# 5. Provide explicit interaction matrices.
#
# The same pattern applies to other Agate model modules (e.g. `DARWIN.construct`).

using Agate
using Agate.Utils: parse_community
using Agate.Library.Light
using Agate.Library.Allometry: AllometricParam, PowerLaw

using OceanBioME
using OceanBioME: Biogeochemistry

using Oceananigans
using Oceananigans.Units

using CairoMakie

const year = years = 365days
nothing #hide

# ## 1. Start from a default model

bgc_default = NiPiZD.construct()

# ## 2. Override community structure

n_phyto = 3
n_zoo = 1

phyto_diameters = (1.5, 20.0, :log_splitting)
zoo_diameters = [60.0]

# ## 3. Override parameters

parameter_overrides = (
    detritus_remineralization = 0.18 / day,
    maximum_growth_rate = (P = AllometricParam(PowerLaw(); prefactor=3.0 / day, exponent=-0.15), Z = 0.0),
)

# ## 4. Swap components (dynamics)

using Agate.Models.NiPiZD.Tracers: nutrient_geider_light, phytoplankton_geider_light

# ## 5. Provide explicit interaction matrices

# We compute matrices in the *constructed community order*.
# `parse_community` is a lightweight helper that matches the constructor's own parsing.

community_ctx = parse_community(
    Float64,
    (;
        Z=(; n=n_zoo, diameters=zoo_diameters),
        P=(; n=n_phyto, diameters=phyto_diameters),
    );
    plankton_dynamics=(Z=Agate.Models.NiPiZD.Tracers.zooplankton_default,
                       P=phytoplankton_geider_light),
    biogeochem_dynamics=(N=nutrient_geider_light,
                         D=Agate.Models.NiPiZD.Tracers.detritus_default),
)

# Canonical interaction matrices are rectangular: (consumer, prey).
nc = length(community_ctx.consumer_indices)
np = length(community_ctx.prey_indices)

pal = zeros(Float64, nc, np)
assim = zeros(Float64, nc, np)

# For illustration, make the single zooplankton class graze all phytoplankton equally.
consumer_local = Dict(gidx => i for (i, gidx) in pairs(community_ctx.consumer_indices))
prey_local     = Dict(gidx => i for (i, gidx) in pairs(community_ctx.prey_indices))

z_global = findall(==(:Z), community_ctx.group_symbols)
p_global = findall(==(:P), community_ctx.group_symbols)

for i in z_global, j in p_global
    haskey(consumer_local, i) && haskey(prey_local, j) || continue
    ic = consumer_local[i]
    jp = prey_local[j]
    pal[ic, jp] = 1.0
    assim[ic, jp] = 0.3
end

# ## 6. Construct the customised model

bgc_custom = NiPiZD.construct(
    ;
    n_phyto,
    n_zoo,
    phyto_diameters,
    zoo_diameters,
    nutrient_dynamics=nutrient_geider_light,
    phyto_dynamics=phytoplankton_geider_light,
    parameters=parameter_overrides,
    palatability_matrix=pal,
    assimilation_matrix=assim,
    sinking_tracers=(D=2.0 / day,),
)

# ## 7. Run a short box model simulation

light = FunctionFieldPAR(; grid=BoxModelGrid())
bgc_model = Biogeochemistry(bgc_custom; light_attenuation=light)

box = BoxModel(; biogeochemistry=bgc_model)

println(tracer_names(bgc_custom))

set!(box; N=7.0, D=0.05, Z1=0.02, P1=0.01, P2=0.01, P3=0.01)

filename = "box_customisation.jld2"

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
