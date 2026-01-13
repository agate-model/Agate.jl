# # Box model: add a detritivorous heterotroph group (H) to NiPiZD
#
# This example demonstrates how to extend an existing factory model (NiPiZD) with:
#
# - a **new plankton group** (`H`),
# - a **new plankton dynamics builder** for that group using the *equation-based* API,
# - a local **detritus (D) dynamics builder** that includes detrital uptake by heterotrophs,
# - and a simple **mass balance test**.
#
# Requirement: `H` should be eaten by `Z` but cannot eat.

using Agate
using Agate.Models: NiPiZDFactory
using Agate.Constructor: construct, PFTSpecification, update_plankton_args
using Agate.Library.Allometry: AllometricParam, PowerLaw
using Agate.Library.Equations: Equation, Σ, group_param, community_param, bgc_param
using Agate.Library.Mortality: linear_loss, linear_loss_sum, quadratic_loss_sum
using Agate.Library.Predation: grazing_loss, grazing_assimilation_loss
using Agate.Library.Light: FunctionFieldPAR

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
# We intentionally use a **new parameter container** (`maximum_detritus_uptake_rate`) rather than
# reusing `maximum_growth_rate`, because NiPiZD’s default nutrient tendency subtracts only
# photosynthetic growth.

@inline monod(D, k) = D / (D + k)

function heterotroph_growth(plankton_syms::AbstractVector{Symbol}, plankton_sym::Symbol, plankton_idx::Int)
    uptake = group_param(:maximum_detritus_uptake_rate)[plankton_idx] *
             monod(:D, group_param(:detritus_half_saturation)[plankton_idx]) *
             plankton_sym

    grazing = grazing_loss(plankton_sym, plankton_idx, plankton_syms)
    mort = linear_loss(plankton_sym, plankton_idx)

    return Equation(uptake - grazing - mort)
end

# ## 2. Define a detritus dynamics that includes heterotroph detrital uptake
#
# NiPiZD’s default detritus tendency does not include detrital uptake.
# Here we define a local `detritus_with_heterotrophs` tendency that adds a sink term:
#
# ```
# ∑ᵢ maximum_detritus_uptake_rate[i] * monod(D, detritus_half_saturation[i]) * plankton[i]
# ```
#
# For groups that don’t consume detritus, these parameters may be missing or `nothing`.
# The constructor fills community-optional parameters with `0`, making those terms inactive.

function detritus_with_heterotrophs(plankton_syms::AbstractVector{Symbol})
    linear_sum = linear_loss_sum(plankton_syms)
    quadratic_sum = quadratic_loss_sum(plankton_syms)
    assimilation_loss_sum = grazing_assimilation_loss(plankton_syms)

    export_frac = bgc_param(:mortality_export_fraction)
    remin = bgc_param(:detritus_remineralization) * :D

    uptake_sum = Σ(plankton_syms) do sym, i
        community_param(:maximum_detritus_uptake_rate)[i] *
        monod(:D, community_param(:detritus_half_saturation)[i]) *
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

plankton_dynamics   = Agate.Models.default_plankton_dynamics(factory)
plankton_args       = Agate.Models.default_plankton_args(factory)
biogeochem_dynamics = Agate.Models.default_biogeochem_dynamics(factory)
biogeochem_args     = Agate.Models.default_biogeochem_args(factory)

# Override detritus (`D`) dynamics to include heterotroph detrital uptake.
biogeochem_dynamics_H = merge(biogeochem_dynamics, (; D = detritus_with_heterotrophs))

# Add a new group entry to `plankton_args`.
#
# Requirement: `H` **cannot eat** but **can be eaten**.
# Agate builds default palatability/assimilation matrices from these traits,
# so `Z` can graze `H` but `H` will not graze anything.

heterotroph_pft = PFTSpecification(
    # New parameters used by `heterotroph_growth`.
    maximum_detritus_uptake_rate = AllometricParam(PowerLaw(); prefactor=1.5 / day, exponent=-0.15),
    detritus_half_saturation     = 0.04,

    # Losses
    linear_mortality = 8e-7 / second,

    # Explicitly mark predator parameters as not applicable.
    maximum_predation_rate  = nothing,
    holling_half_saturation = nothing,

    # Interaction traits used by default allometric interactions.
    can_eat         = false,
    can_be_eaten    = true,
    optimum_predator_prey_ratio = 0.0,
    protection      = 0.0,
    specificity     = 0.0,
    assimilation_efficiency = 0.0,
)

plankton_args_H = merge(plankton_args, (; H=(; n=1, diameters=[6.0], pft=heterotroph_pft)))

# Add matching dynamics for the new group.
plankton_dynamics_H = merge(plankton_dynamics, (; H = heterotroph_growth))

# Optional ergonomic update (same API as other PFT keys):
plankton_args_H = update_plankton_args(plankton_args_H, :H; detritus_half_saturation=0.05)

# ## 4. Construct a new concrete NiPiZDH model type

bgc_type = construct(
    factory;
    plankton_dynamics=plankton_dynamics_H,
    plankton_args=plankton_args_H,
    biogeochem_dynamics=biogeochem_dynamics_H,
    biogeochem_args=biogeochem_args,
)

bgc = bgc_type()

# The new runtime parameter containers exist because they were required by the equations.
(
    H_maximum_detritus_uptake_rate = bgc.parameters.maximum_detritus_uptake_rate[end],
    H_detritus_half_saturation     = bgc.parameters.detritus_half_saturation[end],
)

# ## 5. Run a short box model simulation

light = FunctionFieldPAR(; grid=BoxModelGrid())
bgc_model = Biogeochemistry(bgc; light_attenuation=light)
box = BoxModel(; biogeochemistry=bgc_model)

# With the default NiPiZD community (`Z1`, `Z2`, `P1`, `P2`) plus our new group (`H1`),
# the plankton tracers are: `Z1`, `Z2`, `P1`, `P2`, `H1`.

set!(box; N=7.0, D=0.05, Z1=0.02, Z2=0.02, P1=0.01, P2=0.01, H1=0.01)

# ### Mass balance test
# For this closed box, total mass `N + D + Σ(plankton)` should be conserved
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
