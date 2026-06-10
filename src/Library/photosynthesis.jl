"""Photosynthesis, light limitation, and phytoplankton growth kernels."""

module Photosynthesis

using ..Nutrients: MonodLimitation, LiebigMinimum, SmoothLiebigMinimum

export smith_light_limitation,
    geider_light_limitation,
    liebig_nutrient_limitation,
    smooth_liebig_nutrient_limitation,
    nutrient_limitation,
    smith_growth,
    geider_growth

"""
    SmithLightLimitation(alpha, maximum_growth_0C)

Callable Smith (1936) light-limitation factor as a function of PAR.

!!! formulation
    ```math
    L_S(I) = \\frac{\\alpha I}{\\sqrt{\\mu_0^2 + (\\alpha I)^2}}
    ```

    where ``I`` is photosynthetically active radiation, ``\\alpha`` is the
    initial photosynthetic slope, and ``\\mu_0`` is the maximum growth rate at
    0 °C. The returned factor is dimensionless and is multiplied by ``\\mu_0``
    in `smith_growth`.
"""
struct SmithLightLimitation{T1,T2}
    alpha::T1
    maximum_growth_0C::T2
end

@inline function (f::SmithLightLimitation)(PAR)
    α = f.alpha
    μ₀ = f.maximum_growth_0C
    if α == zero(α) || μ₀ == zero(μ₀)
        return zero(α)
    end
    return α * PAR / sqrt(μ₀ * μ₀ + α * α * PAR * PAR)
end

"""
    GeiderLightLimitation(alpha, maximum_growth_rate, chlorophyll_to_carbon_ratio)

Callable Geider-style light-dependent carbon fixation rate.

!!! formulation
    ```math
    L_G(I) = P^C_{max} \\left[1 - \\exp\\left(-\\frac{\\alpha^{chl}\\theta^C I}{P^C_{max}}\\right)\\right]
    ```

    where ``I`` is photosynthetically active radiation, ``P^C_{max}`` is the
    maximum carbon-specific growth rate, ``\\alpha^{chl}`` is the chlorophyll-
    specific initial slope, and ``\\theta^C`` is the chlorophyll-to-carbon ratio.
    This follows the Geider et al. light response used by DARWIN-style growth.
"""
struct GeiderLightLimitation{T1,T2,T3}
    alpha::T1
    maximum_growth_rate::T2
    chlorophyll_to_carbon_ratio::T3
end

@inline function (f::GeiderLightLimitation)(PAR)
    α = f.alpha
    Pᶜₘₐₓ = f.maximum_growth_rate
    θᶜ = f.chlorophyll_to_carbon_ratio
    if Pᶜₘₐₓ == zero(Pᶜₘₐₓ)
        return zero(Pᶜₘₐₓ)
    end
    return Pᶜₘₐₓ * (one(Pᶜₘₐₓ) - exp((-α * θᶜ * PAR) / Pᶜₘₐₓ))
end

# -----------------------------------------------------------------------------
# Explicit function aliases (preferred developer UX).
# -----------------------------------------------------------------------------

"""
    smith_light_limitation(PAR, alpha, maximum_growth_0C)

Evaluate the Smith (1936) light-limitation factor.

!!! formulation
    ```math
    L_S(I) = \\frac{\\alpha I}{\\sqrt{\\mu_0^2 + (\\alpha I)^2}}
    ```

# Arguments
- `PAR`: photosynthetically active radiation ``I``.
- `alpha`: initial photosynthetic slope ``\\alpha``.
- `maximum_growth_0C`: maximum growth rate ``\\mu_0`` at 0 °C.
"""
@inline smith_light_limitation(PAR, alpha, maximum_growth_0C) = SmithLightLimitation(
    alpha, maximum_growth_0C
)(
    PAR
)

"""
    geider_light_limitation(PAR, alpha, maximum_growth_rate, chlorophyll_to_carbon_ratio)

Evaluate the Geider-style light-dependent growth rate.

!!! formulation
    ```math
    L_G(I) = P^C_{max} \\left[1 - \\exp\\left(-\\frac{\\alpha^{chl}\\theta^C I}{P^C_{max}}\\right)\\right]
    ```

# Arguments
- `PAR`: photosynthetically active radiation ``I``.
- `alpha`: chlorophyll-specific initial slope ``\\alpha^{chl}``.
- `maximum_growth_rate`: maximum carbon-specific growth rate ``P^C_{max}``.
- `chlorophyll_to_carbon_ratio`: chlorophyll-to-carbon ratio ``\\theta^C``.
"""
@inline geider_light_limitation(PAR, alpha, maximum_growth_rate, chlorophyll_to_carbon_ratio) = GeiderLightLimitation(
    alpha, maximum_growth_rate, chlorophyll_to_carbon_ratio
)(
    PAR
)

"""
    liebig_nutrient_limitation(resources, half_saturations, reference)

Compute a Liebig minimum over Monod limitation factors for any number of nutrient
resources.

!!! formulation
    ```math
    \\gamma = \\min_i \\frac{R_i}{K_i + R_i}
    ```

    where ``R_i`` is the concentration of nutrient resource `i` and ``K_i`` is
    its half-saturation parameter. `reference` supplies the scalar type and the
    initial value `one(reference)`.

# Arguments
- `resources`: tuple of nutrient concentrations ``R_i``.
- `half_saturations`: tuple of half-saturation constants ``K_i``.
- `reference`: numeric value used to initialize the limitation factor with the
  appropriate scalar type.
"""
@inline function monod_limitations(resources::Tuple, half_saturations::Tuple)
    return ntuple(length(resources)) do i
        MonodLimitation(half_saturations[i])(resources[i])
    end
end

@inline nutrient_limitation(::Val{:liebig}, limitations::Tuple, reference) = LiebigMinimum()((one(reference), limitations...))

@inline nutrient_limitation(::Val{:smooth_liebig}, limitations::Tuple, reference) = SmoothLiebigMinimum()((one(reference), limitations...))

@inline function liebig_nutrient_limitation(resources::Tuple, half_saturations::Tuple, reference)
    return nutrient_limitation(Val(:liebig), monod_limitations(resources, half_saturations), reference)
end

@inline function smooth_liebig_nutrient_limitation(resources::Tuple, half_saturations::Tuple, reference)
    return nutrient_limitation(Val(:smooth_liebig), monod_limitations(resources, half_saturations), reference)
end

@inline function nutrient_limitation(
    limitation::Val,
    resources::Tuple,
    half_saturations::Tuple,
    reference,
)
    return nutrient_limitation(limitation, monod_limitations(resources, half_saturations), reference)
end

"""
    smith_growth(resources, P, PAR, maximum_growth_0C, half_saturations, alpha)

Compute Smith-style phytoplankton biomass growth with Liebig nutrient limitation.

!!! formulation
    ```math
    G_S = \\mu_0\\,\\gamma\\,L_S(I)\\,P,
    \\qquad
    \\gamma = \\min_i \\frac{R_i}{K_i + R_i}
    ```

    where ``P`` is phytoplankton biomass, ``I`` is PAR, ``L_S`` is the Smith
    light-limitation factor, and ``\\gamma`` is the minimum Monod nutrient
    limitation across the supplied resources. A single nutrient is represented
    by a tuple of length one.

# Arguments
- `resources`: tuple of nutrient concentrations ``R_i``.
- `P`: phytoplankton biomass.
- `PAR`: photosynthetically active radiation ``I``.
- `maximum_growth_0C`: maximum growth rate ``\\mu_0`` at 0 °C.
- `half_saturations`: tuple of nutrient half-saturation constants ``K_i``.
- `alpha`: initial photosynthetic slope ``\\alpha``.
"""
@inline function smith_growth(
    resources::Tuple,
    P,
    PAR,
    maximum_growth_0C,
    half_saturations::Tuple,
    alpha,
)
    return smith_growth(Val(:liebig), resources, P, PAR, maximum_growth_0C, half_saturations, alpha)
end

@inline function smith_growth(
    limitation::Val,
    resources::Tuple,
    P,
    PAR,
    maximum_growth_0C,
    half_saturations::Tuple,
    alpha,
)
    γ = nutrient_limitation(limitation, resources, half_saturations, maximum_growth_0C)
    return maximum_growth_0C * γ * SmithLightLimitation(alpha, maximum_growth_0C)(PAR) * P
end

"""
    geider_growth(resources, P, PAR, maximum_growth_rate, half_saturations, alpha,
                  chlorophyll_to_carbon_ratio)

Compute Geider-style phytoplankton biomass growth with Liebig nutrient limitation.

!!! formulation
    ```math
    G_G = \\gamma\\,L_G(I)\\,P,
    \\qquad
    \\gamma = \\min_i \\frac{R_i}{K_i + R_i}
    ```

    where ``P`` is phytoplankton biomass, ``I`` is PAR, ``L_G`` is the Geider
    light-dependent growth rate, and ``\\gamma`` is the minimum Monod nutrient
    limitation across the supplied resources. A single nutrient is represented
    by a tuple of length one.

# Arguments
- `resources`: tuple of nutrient concentrations ``R_i``.
- `P`: phytoplankton biomass.
- `PAR`: photosynthetically active radiation ``I``.
- `maximum_growth_rate`: maximum carbon-specific growth rate ``P^C_{max}``.
- `half_saturations`: tuple of nutrient half-saturation constants ``K_i``.
- `alpha`: chlorophyll-specific initial slope ``\\alpha^{chl}``.
- `chlorophyll_to_carbon_ratio`: chlorophyll-to-carbon ratio ``\\theta^C``.
"""
@inline function geider_growth(
    resources::Tuple,
    P,
    PAR,
    maximum_growth_rate,
    half_saturations::Tuple,
    alpha,
    chlorophyll_to_carbon_ratio,
)
    return geider_growth(
        Val(:liebig),
        resources,
        P,
        PAR,
        maximum_growth_rate,
        half_saturations,
        alpha,
        chlorophyll_to_carbon_ratio,
    )
end

@inline function geider_growth(
    limitation::Val,
    resources::Tuple,
    P,
    PAR,
    maximum_growth_rate,
    half_saturations::Tuple,
    alpha,
    chlorophyll_to_carbon_ratio,
)
    γ = nutrient_limitation(limitation, resources, half_saturations, maximum_growth_rate)
    return γ *
           GeiderLightLimitation(alpha, maximum_growth_rate, chlorophyll_to_carbon_ratio)(PAR) *
           P
end

end # module
