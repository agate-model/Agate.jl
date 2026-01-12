# # Box model: add a detritivorous heterotroph group (H) to NiPiZD
#
# This example demonstrates how to extend an existing factory model (NiPiZD) with:
#
# - a **new plankton group** (`H`),
# - a **new plankton dynamics function** for that group,
# - **registering** the new dynamics so the constructor can allocate required runtime parameter containers,
# - adding **new PFT arguments** for the group,
# - updating the **detritus (D) dynamics** to include the detritus consumption term from heterotroph growth,
# - and a simple **mass balance test**.
#
# Requirement: `H` should be eaten by `Z` but cannot eat.

using Agate
using Agate.Models: NiPiZDFactory
using Agate.Constructor: construct, PFTSpecification, update_plankton_args
using Agate.Utils: @register_dynamics, sum_expr
using Agate.Library.Predation: predation_loss_sum, predation_assimilation_loss_sum
using Agate.Library.Nutrients: monod_limitation
using Agate.Library.Mortality: linear_loss_sum, quadratic_loss_sum
using Agate.Library.Light: FunctionFieldPAR

using OceanBioME
using OceanBioME: Biogeochemistry

using Oceananigans
using Oceananigans.Units

using Test

nothing #hide

# ## 1. Define a new plankton dynamics builder
#
# Dynamics builders return an `Expr` for the tracer tendency.
# Here, `H` grows from detritus `D` with a Monod limitation, is grazed by predators
# (primarily `Z`), and experiences linear mortality.
#
# We intentionally use a **new parameter container** (`maximum_detritus_uptake_rate`) rather than
# reusing `maximum_growth_rate`, because NiPiZD’s default **nutrient** tendency subtracts the
# photosynthetic growth summed over all plankton tracers.

function heterotroph_growth(plankton_syms::AbstractVector{Symbol}, plankton_sym::Symbol, plankton_idx::Int)
    growth = :(maximum_detritus_uptake_rate[$plankton_idx] *
               monod_limitation(D, detritus_half_saturation[$plankton_idx]) *
               $plankton_sym)

    grazing_loss = predation_loss_sum(plankton_syms, plankton_sym, plankton_idx)
    linear_mortality = :(linear_loss($plankton_sym, linear_mortality[$plankton_idx]))

    return :($growth - ($grazing_loss) - $linear_mortality)
end

# ## 2. Register the new dynamics (so the constructor can build parameter containers)
#
# Registering tells the constructor which runtime parameter containers must exist.
# Any non-scalar, non-matrix symbol becomes a **vector container** of length `n_total` in `bgc.parameters`.

@register_dynamics heterotroph_growth (
    :maximum_detritus_uptake_rate,
    :detritus_half_saturation,
    :maximum_predation_rate,
    :holling_half_saturation,
    :palatability_matrix,
    :linear_mortality,
)


# ## 3. Define a detritus dynamics that includes heterotroph detrital uptake
#
# NiPiZD’s default detritus tendency does not include detrital uptake.
# Here we define a local `detritus_with_heterotrophs` tendency that adds a sink term
# corresponding to the detritus-based growth of any detritivorous plankton (including `H`).
#
# The sink term is:
#
# ```
# ∑ᵢ maximum_detritus_uptake_rate[i] * monod_limitation(D, detritus_half_saturation[i]) * plankton[i]
# ```
#
# Groups that don’t consume detritus should set `maximum_detritus_uptake_rate = 0` (the defaults do this).

function detritus_with_heterotrophs(plankton_syms::AbstractVector{Symbol})
    linear_sum = linear_loss_sum(plankton_syms)
    quadratic_sum = quadratic_loss_sum(plankton_syms)
    assimilation_loss_sum = predation_assimilation_loss_sum(plankton_syms)

    # sum of heterotroph detrital uptake
    uptake_terms = Expr[]
    for (i, sym) in enumerate(plankton_syms)
        push!(
            uptake_terms,
            :(maximum_detritus_uptake_rate[$i] *
              monod_limitation(D, detritus_half_saturation[$i]) *
              $sym),
        )
    end
    heterotroph_uptake_sum = sum_expr(uptake_terms)

    return :(
        (1 - mortality_export_fraction) * ($linear_sum) +
        ($assimilation_loss_sum) +
        (1 - mortality_export_fraction) * ($quadratic_sum) -
        remineralization_idealized(D, detritus_remineralization) -
        ($heterotroph_uptake_sum)
    )
end

@register_dynamics detritus_with_heterotrophs (
    :linear_mortality,
    :quadratic_mortality,
    :maximum_predation_rate,
    :holling_half_saturation,
    :palatability_matrix,
    :assimilation_efficiency_matrix,
    :detritus_remineralization,
    :mortality_export_fraction,
    :maximum_detritus_uptake_rate,
    :detritus_half_saturation,
)


# ## 4\. Extend the NiPiZD factory configuration with a new group `H`

factory = NiPiZDFactory()

plankton_dynamics = Agate.Models.default_plankton_dynamics(factory)
plankton_args     = Agate.Models.default_plankton_args(factory)
biogeochem_dynamics = Agate.Models.default_biogeochem_dynamics(factory)
biogeochem_args     = Agate.Models.default_biogeochem_args(factory)

# Override detritus (`D`) dynamics to include heterotroph detrital uptake.
biogeochem_dynamics_H = merge(biogeochem_dynamics, (; D = detritus_with_heterotrophs))

# Add a new group entry to `plankton_args`.
#
# Requirement: `H` **cannot eat** but **can be eaten**.
# NiPiZD’s default interaction builder uses these flags when building palatability/assimilation matrices,
# so `Z` can graze `H` but `H` will not graze anything.

heterotroph_pft = PFTSpecification(
    # New parameters used by `heterotroph_growth`.
    maximum_detritus_uptake_rate_a = 1.5 / day,
    maximum_detritus_uptake_rate_b = -0.15,
    detritus_half_saturation = 0.04,

    # Losses
    linear_mortality = 8e-7 / second,

    # Interaction flags
    can_eat = false,
    can_be_eaten = true,

    # Optional interaction knobs used by NiPiZD's default allometric interactions
    optimum_predator_prey_ratio = 0.0,
    protection = 0.0,
    specificity = 0.0,
    assimilation_efficiency = 0.0,
)

plankton_args_H = merge(plankton_args, (; H=(; n=1, diameters=[6.0], pft=heterotroph_pft)))

# Add matching dynamics for the new group.
plankton_dynamics_H = merge(plankton_dynamics, (; H=heterotroph_growth))

# Optional ergonomic update (same API as other PFT keys):
plankton_args_H = update_plankton_args(plankton_args_H, :H; detritus_half_saturation=0.05)

# ## 5. Construct a new concrete NiPiZDH model type

bgc_type = construct(
    factory;
    plankton_dynamics=plankton_dynamics_H,
    plankton_args=plankton_args_H,
    biogeochem_dynamics=biogeochem_dynamics_H,
    biogeochem_args=biogeochem_args,
)

bgc = bgc_type()

# The new runtime parameter containers exist because we registered them.
bgc.parameters.detritus_half_saturation
bgc.parameters.maximum_detritus_uptake_rate

# ## 6. Run a short box model simulation

light = FunctionFieldPAR(; grid=BoxModelGrid())
bgc_model = Biogeochemistry(bgc; light_attenuation=light)
box = BoxModel(; biogeochemistry=bgc_model)

# With the default NiPiZD community (`Z1`, `Z2`, `P1`, `P2`) plus our new group (`H1`),
# the plankton tracers are: `Z1`, `Z2`, `P1`, `P2`, `H1`.

set!(box; N=7.0, D=0.05, Z1=0.02, Z2=0.02, P1=0.01, P2=0.01, H1=0.01)

# ### Mass balance test
# For this closed box, total mass `N + D + Σ(plankton)` should be conserved (up to time-integration error).

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
