# # Box model: add a detritivorous heterotroph group (H) to NiPiZD
#
# This example demonstrates how to extend an existing factory model (NiPiZD) with:
#
# - a **new plankton group** (`H`),
# - a **new plankton dynamics builder** for that group using the equation-based API,
# - a local **detritus (D) dynamics builder** that includes detrital uptake by heterotrophs,
# - and a simple mass balance test.
#
# Requirement: `H` should be eaten by `Z` but cannot eat.

using Agate
using Agate.Utils.Specifications: PFTSpecification
using Agate.Library.Allometry: AllometricParam, PowerLaw
using Agate.Library.Equations: Equation, sum_over
using Agate.Library.Mortality: linear_loss, linear_loss_sum, quadratic_loss_sum
using Agate.Library.Predation: grazing_loss, grazing_assimilation_loss
using Agate.Library.Light: FunctionFieldPAR
using Agate.Parameters: ParamSpec, parameter_registry, extend_registry

using OceanBioME
using OceanBioME: Biogeochemistry

using Oceananigans
using Oceananigans.Units

using Test

nothing #hide

# ## 1. Define a new plankton dynamics builder
#
# Dynamics builders return an `Equation`, assembled from symbolic building blocks.
# Here, `H` grows from detritus `D` with a Monod-style limitation,
# is grazed by predators (`Z`), and experiences linear mortality.
#
#+#+#+#+
# Parameter placeholders live under the canonical namespace `Agate.ParamVars`.
# When we extend the parameter registry below, `construct` will declare the
# required placeholders in `Agate.ParamVars` automatically.

const PV = Agate.ParamVars

@inline monod(D, k) = D / (D + k)

function heterotroph_growth(plankton_syms::AbstractVector{Symbol}, plankton_sym::Symbol, plankton_idx::Int)
    uptake = PV.maximum_detritus_uptake_rate[plankton_idx] *
             monod(:D, PV.detritus_half_saturation[plankton_idx]) *
             plankton_sym

    grazing = grazing_loss(plankton_sym, plankton_idx, plankton_syms)
    mort = linear_loss(plankton_sym, plankton_idx)

    return Equation(uptake - grazing - mort)
end

# ## 2. Define a detritus dynamics that includes heterotroph detrital uptake
#
# NiPiZD’s default detritus tendency does not include detrital uptake.
# Here we define a local tendency that adds a sink term:
#
# ```julia
# sum_over(plankton_syms) do sym, i
#     maximum_detritus_uptake_rate[i] * monod(D, detritus_half_saturation[i]) * sym
# end
# ```
#
# For groups that don’t consume detritus, these parameters are missing.
# We will declare them in the parameter registry with `scope=:zero_silent` so missing entries become zeros.

function detritus_with_heterotrophs(plankton_syms::AbstractVector{Symbol})
    linear_sum = linear_loss_sum(plankton_syms)
    quadratic_sum = quadratic_loss_sum(plankton_syms)
    assimilation_loss_sum = grazing_assimilation_loss(plankton_syms)

    export_frac = PV.mortality_export_fraction
    remin = PV.detritus_remineralization * :D

    uptake_sum = sum_over(plankton_syms) do sym, i
        PV.maximum_detritus_uptake_rate[i] *
        monod(:D, PV.detritus_half_saturation[i]) *
        sym
    end

    return Equation(
        (1 - export_frac) * linear_sum +
        assimilation_loss_sum +
        (1 - export_frac) * quadratic_sum -
        remin -
        uptake_sum,
    )
end

# ## 3. Extend the NiPiZD factory configuration with a new group `H`

factory = NiPiZDFactory()

# Optional: inspect the base registry before extending.
println(parameter_registry(factory))

plankton_dynamics   = Agate.Models.default_plankton_dynamics(factory)
community           = Agate.Models.default_community(factory)
biogeochem_dynamics = Agate.Models.default_biogeochem_dynamics(factory)

# Override detritus (`D`) dynamics to include heterotroph detrital uptake.
biogeochem_dynamics_H = update_dynamics(biogeochem_dynamics; D=detritus_with_heterotrophs)

# Add a new group entry to the community.
#
# Requirement: `H` cannot eat but can be eaten.
# Agate builds default palatability/assimilation matrices from these traits.

heterotroph_pft = PFTSpecification(
    # Interaction traits used by default allometric interactions.
    can_eat         = false,
    can_be_eaten    = true,
    optimum_predator_prey_ratio = 0.0,
    protection      = 0.0,
    specificity     = 0.0,
    assimilation_efficiency = 0.0,
)

community_H = extend_community(community; H=(; n=1, diameters=[6.0], pft=heterotroph_pft))

# Add matching dynamics for the new group.
plankton_dynamics_H = extend_dynamics(plankton_dynamics; H=heterotroph_growth)

# ## 4. Extend the parameter registry with the new parameters
#
# The resolver requires every parameter referenced in equations to have a `ParamSpec`.
# Here we add two new vector parameters used by heterotroph detrital uptake.

extra_specs = [
    ParamSpec(
        :maximum_detritus_uptake_rate,
        1 / day,
        "Maximum detritus uptake rate for detritivorous heterotrophs.",
        (H = AllometricParam(PowerLaw(); prefactor=1.5 / day, exponent=-0.15),);
        scope=:zero_silent,
    ),
    ParamSpec(
        :detritus_half_saturation,
        1,
        "Half-saturation constant for heterotroph detritus uptake.",
        (H = 0.04,);
        scope=:zero_silent,
    ),
]

base_reg = parameter_registry(factory)
extended_reg = extend_registry(base_reg, extra_specs...)

# Inspect the extended registry including the new heterotroph parameters.
println(extended_reg)

# Optional ergonomic override: tweak detritus half-saturation for H via the community PFT.
community_H = update_community(community_H, :H; detritus_half_saturation=0.05)

# ## 5. Construct a new concrete NiPiZDH model type

bgc = construct(
    factory;
    plankton_dynamics=plankton_dynamics_H,
    biogeochem_dynamics=biogeochem_dynamics_H,
    community=community_H,
    registry=extended_reg,
)

# The new runtime parameter vectors exist because they were required by the equations.
(
    H_maximum_detritus_uptake_rate = bgc.parameters.maximum_detritus_uptake_rate[end],
    H_detritus_half_saturation     = bgc.parameters.detritus_half_saturation[end],
)

# ## 6. Run a short box model simulation

light = FunctionFieldPAR(; grid=BoxModelGrid())
bgc_model = Biogeochemistry(bgc; light_attenuation=light)
box = BoxModel(; biogeochemistry=bgc_model)

# With the default NiPiZD community (`Z1`, `Z2`, `P1`, `P2`) plus our new group (`H1`),
# the plankton tracers are: `Z1`, `Z2`, `P1`, `P2`, `H1`.

set!(box; N=7.0, D=0.05, Z1=0.02, Z2=0.02, P1=0.01, P2=0.01, H1=0.01)

# ### Mass balance test
# For this closed box, total mass (N + D + total plankton) should be conserved
# (up to time-integration error).

function total_mass(box)
    tracers = (:N, :D, :Z1, :Z2, :P1, :P2, :H1)
    s = 0.0
    for name in tracers
        s += getproperty(box.fields, name)[1, 1, 1]
    end
    return s
end

m₀ = total_mass(box)

sim = Simulation(box; Δt=30minutes, stop_time=30days)
run!(sim)

m₁ = total_mass(box)
@test isapprox(m₁, m₀; rtol=1e-4, atol=1e-10)

# Inspect final values (at grid point 1,1,1)
(
    N  = box.fields.N[1, 1, 1],
    D  = box.fields.D[1, 1, 1],
    Z1 = box.fields.Z1[1, 1, 1],
    Z2 = box.fields.Z2[1, 1, 1],
    P1 = box.fields.P1[1, 1, 1],
    P2 = box.fields.P2[1, 1, 1],
    H1 = box.fields.H1[1, 1, 1],
)
