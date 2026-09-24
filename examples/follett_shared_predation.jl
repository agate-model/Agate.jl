# # [Shared predation and phytoplankton exclusion] (@id follett_shared_predation_example)
#
# Follett et al. (2022, PNAS, doi:10.1073/pnas.2110993118) proposed that
# *Prochlorococcus*-like phytoplankton can be excluded as nutrient supply increases because
# similarly sized heterotrophic bacteria share the same predators. Increasing nutrient supply
# supports larger bacterial size classes; shared predation then turns that bacterial increase
# into an indirect top-down pressure on similarly sized phytoplankton.
#
# Here we reproduce an **idealized version of that size-dependent mechanism** with
# FrankenLOBSTER rather than the paper's exact equations. Following the paper's zero-dimensional
# experiment, we use 10 size classes with the same logarithmic spacing as the first 10 classes of
# its 15-class 0.6--104 um global spectrum. P and H occupy matching prey sizes, while each Z class
# is centred on FrankenLOBSTER's default 10:1 predator:prey size optimum.
#
# The parameter sweep is expressed as a SciML `EnsembleProblem`. Each nutrient-supply value and
# treatment is one trajectory, and `EnsembleThreads` distributes those trajectories across Julia
# threads. This pattern extends directly to additional parameter axes.

# ## Loading dependencies

using Agate
using CairoMakie
using OrdinaryDiffEq: Tsit5, solve
using SciMLBase: EnsembleProblem, EnsembleThreads, ODEProblem, remake
using Statistics: mean

using OceanBioME: BoxModelGrid, PrescribedPhotosyntheticallyActiveRadiation
using OceanBioME.Models.NutrientsPlanktonDetritusModels: LOBSTER
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers
using Oceananigans.Fields: ConstantField
using Oceananigans.Units: day

const FrankenLOBSTER = Agate.Models.FrankenLOBSTER
nothing #hide

# ## A size spectrum that crosses the growth-rate optimum
#
# Follett et al. use 10 size classes in their zero-dimensional experiment, spaced as in the
# global model. The global spectrum has 15 logarithmically spaced classes between 0.6 and 104 um,
# so we use its first 10 classes here. This deliberately spans the ~3 um ESD breakpoint where the
# maximum-growth allometry changes from increasing to decreasing with size.

const FOLLETT_GRID = BoxModelGrid()
const FOLLETT_GLOBAL_SIZES = collect(10.0 .^ range(log10(0.6), log10(104.0); length=15))
const FOLLETT_PREY_SIZES = FOLLETT_GLOBAL_SIZES[1:10]
const FOLLETT_GRAZER_SIZES = 10 .* FOLLETT_PREY_SIZES
const FOLLETT_SIZE_BREAK = 3.0
const FOLLETT_SIZE_STRUCTURE = (
    phytoplankton=(phyto=FOLLETT_PREY_SIZES,),
    zooplankton=(grazer=FOLLETT_GRAZER_SIZES,),
    bacterioplankton=(bacteria=FOLLETT_PREY_SIZES,),
)

const N_SIZE_CLASSES = length(FOLLETT_PREY_SIZES)
const P_TRACERS = ntuple(i -> Symbol("phyto_$i"), N_SIZE_CLASSES)
const Z_TRACERS = ntuple(i -> Symbol("grazer_$i"), N_SIZE_CLASSES)
const H_TRACERS = ntuple(i -> Symbol("bacteria_$i"), N_SIZE_CLASSES)

# FrankenLOBSTER's canonical P/H defaults are fitted to its small (<3 um) classes and therefore
# use the positive V^0.28 branch throughout. Follett's size spectrum crosses 3 um, where maximum
# growth changes to a negative size dependence. For this example only, preserve the existing
# small-cell coefficients and make the relationship continuous at 3 um before switching to a
# Darwin-style large-cell V^-0.15 branch. Other FrankenLOBSTER allometries retain their normal defaults.

spherical_volume(diameter) = pi / 6 * diameter^3

function unimodal_rate(
    diameter,
    small_prefactor;
    breakpoint=FOLLETT_SIZE_BREAK,
    small_exponent=0.28,
    large_exponent=-0.15,
)
    volume = spherical_volume(diameter)
    break_volume = spherical_volume(breakpoint)
    if diameter <= breakpoint
        return small_prefactor * volume^small_exponent
    end
    rate_at_break = small_prefactor * break_volume^small_exponent
    return rate_at_break * (volume / break_volume)^large_exponent
end

named_values(names::Tuple, values) = NamedTuple{names}(Tuple(values))

const FOLLETT_PARAMETERS = (
    maximum_growth_rate=named_values(
        P_TRACERS,
        (unimodal_rate(diameter, 1.2066 / day) for diameter in FOLLETT_PREY_SIZES),
    ),
    bacterial_maximum_uptake_rate=named_values(
        H_TRACERS,
        (unimodal_rate(diameter, 1.836 / day) for diameter in FOLLETT_PREY_SIZES),
    ),
)

# The mechanistic contrast should come from community structure rather than a special detritus
# configuration. We therefore use OceanBioME's standard LOBSTER detritus and remineralization
# unchanged. FrankenLOBSTER's default 5% P exudation supplies one route to DOM, alongside the
# ordinary coupled detrital pathways.

const FOLLETT_LIGHT = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(100.0))

const FOLLETT_PLANKTON = FrankenLOBSTER.construct(
    ;
    size_structure=FOLLETT_SIZE_STRUCTURE,
    parameters=FOLLETT_PARAMETERS,
)
const FOLLETT_COUPLED = LOBSTER(
    FOLLETT_GRID;
    plankton=FOLLETT_PLANKTON,
    light_attenuation=FOLLETT_LIGHT,
    open_bottom=false,
)
const FOLLETT_BGC = FOLLETT_COUPLED.underlying_biogeochemistry
const FOLLETT_TRACERS = required_biogeochemical_tracers(FOLLETT_COUPLED)

nothing #hide

# ## A small SciML adapter for the OceanBioME box tendencies
#
# OceanBioME's NPD model evaluates tendencies from 1x1x1 tracer fields. `BoxCell` provides that
# interface without allocating a full field for every tracer at every ODE evaluation.

struct BoxCell{T} <: AbstractArray{T,3}
    value::T
end

Base.size(::BoxCell) = (1, 1, 1)
Base.IndexStyle(::Type{<:BoxCell}) = IndexCartesian()
@inline Base.getindex(cell::BoxCell, ::Int, ::Int, ::Int) = cell.value

@inline function box_fields(u)
    values = ntuple(i -> BoxCell(u[i]), length(FOLLETT_TRACERS))
    return NamedTuple{FOLLETT_TRACERS}(values)
end

const FOLLETT_AUXILIARY = (PAR=BoxCell(100.0),)

function tracer_index(name)
    index = findfirst(==(name), FOLLETT_TRACERS)
    isnothing(index) && error("Expected tracer :$name; got $(FOLLETT_TRACERS)")
    return index
end

const NO3_INDEX = tracer_index(:NO₃)
const P_INDICES = Tuple(tracer_index(name) for name in P_TRACERS)
const Z_INDICES = Tuple(tracer_index(name) for name in Z_TRACERS)
const H_INDICES = Tuple(tracer_index(name) for name in H_TRACERS)

const NITROGEN_TRACERS = (:NO₃, :NH₄, :DOM, :sPOM, :bPOM, P_TRACERS..., Z_TRACERS..., H_TRACERS...)
const NITROGEN_INDICES = Tuple(tracer_index(name) for name in NITROGEN_TRACERS)

expected_tracers = Set((NITROGEN_TRACERS..., :T))
Set(FOLLETT_TRACERS) == expected_tracers || error(
    "Unexpected FrankenLOBSTER tracer set for this example: $(FOLLETT_TRACERS)"
)

# The paper varies a constant inorganic-resource input. A continuously supplied box also needs an
# export term to keep total nitrogen bounded, so every N-bearing pool experiences the same parcel
# exchange rate:
#
# ```math
# \frac{dN_i}{dt} = F_i(\mathbf{N}) - D N_i, \qquad
# \frac{dNO_3}{dt} = F_{NO_3}(\mathbf{N}) + S_N - D NO_3.
# ```
#
# The SciML parameter object stores the trajectory-specific supply rate and common exchange rate.


function follett_rhs!(du, u, p, t)
    fields = box_fields(u)
    clock = (; time=t)

    for (i, tracer) in enumerate(FOLLETT_TRACERS)
        du[i] = FOLLETT_BGC(
            1, 1, 1, FOLLETT_GRID, Val(tracer), clock, fields, FOLLETT_AUXILIARY
        )
    end

    du[NO3_INDEX] += p.supply_rate

    for i in NITROGEN_INDICES
        du[i] -= p.dilution_rate * u[i]
    end

    return nothing
end
nothing #hide

# ## Nutrient-supply experiment
#
# Supply rates span the order of magnitude highlighted in the Follett et al. zero-dimensional
# experiments (~1e-7 mmol N m^-3 s^-1) and extend into both more oligotrophic and more productive
# conditions. Each supply rate is run twice: once with no bacterial seed and once with bacteria
# present. Since H growth is proportional to H biomass, a zero seed remains a bacteria-free
# control without altering the model equations.

const SUPPLY_RATES = 10.0 .^ range(-9.0, -5.5; length=30)
const DILUTION_RATE = 0.02 / day
const TREATMENTS = (:without_bacteria, :shared_predation)
const EXPERIMENTS = vec([
    (; supply_index, supply_rate=SUPPLY_RATES[supply_index], treatment)
    for supply_index in eachindex(SUPPLY_RATES), treatment in TREATMENTS
])

function initial_state(; bacteria_seed)
    u0 = zeros(length(FOLLETT_TRACERS))
    u0[NO3_INDEX] = 0.02
    u0[tracer_index(:NH₄)] = 0.0
    u0[tracer_index(:T)] = 20.0
    u0[tracer_index(:DOM)] = 0.01

    # Equal N biomass per size class keeps the initialization neutral with respect to size.
    for i in P_INDICES
        u0[i] = 0.002
    end
    for i in Z_INDICES
        u0[i] = 0.0005
    end
    for i in H_INDICES
        u0[i] = bacteria_seed
    end
    return u0
end

const STOP_TIME = 6 * 365day
const AVERAGING_WINDOW = 2 * 365day
const SAVE_INTERVAL = 10day
const SAVE_TIMES = (STOP_TIME - AVERAGING_WINDOW):SAVE_INTERVAL:STOP_TIME

base_problem = ODEProblem(
    follett_rhs!,
    initial_state(; bacteria_seed=0.0),
    (0.0, STOP_TIME),
    (; supply_rate=first(SUPPLY_RATES), dilution_rate=DILUTION_RATE),
)

function follett_prob_func(prob, context)
    experiment = EXPERIMENTS[context.sim_id]
    bacteria_seed = experiment.treatment === :shared_predation ? 0.0005 : 0.0
    return remake(
        prob;
        u0=initial_state(; bacteria_seed),
        p=(; supply_rate=experiment.supply_rate, dilution_rate=DILUTION_RATE),
    )
end

trailing_mean(sol, index) = mean(u[index] for u in sol.u)
trailing_spectrum(sol, indices) = [trailing_mean(sol, index) for index in indices]

function follett_output_func(sol, context)
    experiment = EXPERIMENTS[context.sim_id]
    result = (;
        experiment...,
        P=trailing_spectrum(sol, P_INDICES),
        H=trailing_spectrum(sol, H_INDICES),
    )
    return result, false
end

ensemble_problem = EnsembleProblem(
    base_problem;
    prob_func=follett_prob_func,
    output_func=follett_output_func,
    safetycopy=false, # prob_func remakes rather than mutates the shared template problem
)

ensemble = solve(
    ensemble_problem,
    Tsit5(),
    EnsembleThreads();
    trajectories=length(EXPERIMENTS),
    saveat=SAVE_TIMES,
    save_start=false,
    reltol=1e-7,
    abstol=1e-9,
)

nothing #hide

# ## Size-dependent exclusion along the supply gradient

function treatment_matrix(results, treatment, field)
    values = fill(NaN, N_SIZE_CLASSES, length(SUPPLY_RATES))
    for result in results
        result.treatment === treatment || continue
        values[:, result.supply_index] .= getproperty(result, field)
    end
    return values
end

P_without_bacteria = treatment_matrix(ensemble.u, :without_bacteria, :P)
P_shared = treatment_matrix(ensemble.u, :shared_predation, :P)
H_shared = treatment_matrix(ensemble.u, :shared_predation, :H)
relative_P = P_shared ./ max.(P_without_bacteria, eps(Float64))

# Plot supply in per-day units and biomass on a logarithmic color scale. The dashed horizontal
# line marks the 3 um maximum-growth breakpoint. A moving band of low `P shared / P control`
# indicates the Follett shared-predation exclusion mechanism extending to progressively larger
# prey as resource supply increases.

supply_per_day = SUPPLY_RATES .* day
biomass_floor = 1e-10
log_P_without = log10.(max.(P_without_bacteria, biomass_floor))
log_P_shared = log10.(max.(P_shared, biomass_floor))
log_H_shared = log10.(max.(H_shared, biomass_floor))

fig = Figure(; size=(1050, 760), fontsize=14)

function spectrum_axis(position, title)
    axis = Axis(
        position;
        xlabel="inorganic N supply (mmol N m^-3 day^-1)",
        ylabel="prey ESD (um)",
        xscale=log10,
        yscale=log10,
        title,
    )
    hlines!(axis, [FOLLETT_SIZE_BREAK]; linestyle=:dash)
    return axis
end

ax1 = spectrum_axis(fig[1, 1], "P biomass, bacteria absent")
hm1 = heatmap!(ax1, supply_per_day, FOLLETT_PREY_SIZES, log_P_without')
Colorbar(fig[1, 2], hm1; label="log10 biomass (mmol N m^-3)")

ax2 = spectrum_axis(fig[1, 3], "P biomass, shared predation")
hm2 = heatmap!(ax2, supply_per_day, FOLLETT_PREY_SIZES, log_P_shared')
Colorbar(fig[1, 4], hm2; label="log10 biomass (mmol N m^-3)")

ax3 = spectrum_axis(fig[2, 1], "H biomass, shared predation")
hm3 = heatmap!(ax3, supply_per_day, FOLLETT_PREY_SIZES, log_H_shared')
Colorbar(fig[2, 2], hm3; label="log10 biomass (mmol N m^-3)")

ax4 = spectrum_axis(fig[2, 3], "P with bacteria / P without bacteria")
hm4 = heatmap!(ax4, supply_per_day, FOLLETT_PREY_SIZES, relative_P')
Colorbar(fig[2, 4], hm4; label="relative P biomass")

output_path = joinpath(@__DIR__, "follett_shared_predation.png")
save(output_path, fig; px_per_unit=1)

fig
