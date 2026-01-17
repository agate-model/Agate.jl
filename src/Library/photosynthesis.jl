"""Functions related to phytoplankton light uptake."""

module Photosynthesis

using ..Nutrients: monod_limitation, liebig_minimum
using ...Equations: AExpr, req, merge_requirements, _to_aexpr
export light_limitation_smith,
    light_limitation_geider,
    photosynthetic_growth_single_nutrient,
    photosynthetic_growth_single_nutrient_geider_light,
    photosynthetic_growth_two_nutrients_geider_light

# Symbolic (construction-time) equation blocks
export growth_single_nutrient
export growth_single_nutrient_comm
export growth_single_nutrient_geider
export growth_single_nutrient_geider_comm
export growth_two_nutrients_geider
export growth_two_nutrients_geider_comm

@inline _vec(PV, key::Symbol, idx::Int) = getproperty(PV, key)[idx]

"""Build an `AExpr` for a runtime function call, merging argument requirements."""
function _call(fsym::Symbol, args...)
    merged = req()
    nodes = Any[]
    for a in args
        ae = _to_aexpr(a)
        merged = merge_requirements(merged, ae.req)
        push!(nodes, ae.node)
    end
    return AExpr(Expr(:call, fsym, nodes...), merged)
end

"""
    light_limitation_smith(PAR, initial_slope, maximum_growth_0C)
    
Smith 1936 formulation of light limitation (also see Evans and Parslow, 1985).

!!! formulation
    α * PAR / √(μ₀² + α² * PAR²)

    where:
    - PAR = photosynthetic active radiation
    - α = initial photosynthetic slope
    - μ₀ = maximum growth rate at 0 °C    

# Arguments
- `PAR`: photosynthetic active radiation
- `initial_slope`: initial photosynthetic slope
- `maximum_growth_0C`: maximum growth rate at T = 0 °C
"""
@inline function light_limitation_smith(PAR, initial_slope, maximum_growth_0C)
    if initial_slope == zero(initial_slope)
        return zero(initial_slope)
    end

    return initial_slope * PAR / sqrt(maximum_growth_0C^2 + initial_slope^2 * PAR^2)
end

"""
    light_limitation_geider(PAR, photosynthetic_slope, maximum_growth_rate, 
        chlorophyll_to_carbon_ratio)

A light limitation function which depends on the cellular ratio of chlorophyll to carbon.
This formulation is based on equation (4) from Geider et al., 1998.

!!! formulation
    Pᶜₘₐₓ[1-exp((-αᶜʰˡθᶜPAR)/Pᶜₘₐₓ)]
    
    where:
    - Pᶜₘₐₓ = maximum_growth rate
    - αᶜʰˡ = photosynthetic slope
    - PAR = photosynthetic active radiation
    - θᶜ = chlorophyll to carbon ratio

# Arguments
- `PAR`: photosynthetic active radiation
- `maximum_growth_rate`: maximum growth rate before nutrient limitation
- `photosynthetic_slope`: initial photosynthetic slope
- `chlorophyll_to_carbon_ratio`: ratio between cellular chlorophyll and carbon
"""
@inline function light_limitation_geider(
    PAR, photosynthetic_slope, maximum_growth_rate, chlorophyll_to_carbon_ratio
)
    if maximum_growth_rate == zero(maximum_growth_rate)
        return zero(maximum_growth_rate)
    end

    return maximum_growth_rate * (
        one(maximum_growth_rate) - exp(
            (-photosynthetic_slope * chlorophyll_to_carbon_ratio * PAR) /
            maximum_growth_rate,
        )
    )
end

"""
    photosynthetic_growth_single_nutrient(N, P, PAR, maximum_growth_0C, 
        nutrient_half_saturation, initial_slope)

Single nutrient monod smith photosynthetic growth (used, for example, in Kuhn 2015).

!!! formulation
    μ₀ * γᴿ * γᴾᴬᴿ * ``P``
    
    where:
    - μ₀ = maximum growth rate at 0 °C
    - γᴿ = `monod_limitation(R, kᵣ)`
    - R = nutrient concentration
    - kᵣ = nutrient half saturation
    - γᴾᴬᴿ = `light_limitation_smith(PAR, α, μ₀)`
    - PAR = photosynthetic active radiation
    - α = initial photosynthetic slope
    - P = plankton concentration

# Arguments
- `R`: nutrient concentration
- `P`: plankton concentration
- `PAR`: photosynthetic active radiation
- `maximum_growth_0C`: maximum growth rate at T = 0 °C
- `nutrient_half_saturation`: nutrient half saturation
- `initial_slope`: initial photosynthetic slope
"""
@inline function photosynthetic_growth_single_nutrient(
    R, P, PAR, maximum_growth_0C, nutrient_half_saturation, initial_slope
)
    return maximum_growth_0C *
           monod_limitation(R, nutrient_half_saturation) *
           light_limitation_smith(PAR, initial_slope, maximum_growth_0C) *
           P
end

"""
    photosynthetic_growth_single_nutrient_geider_light(R, P, PAR, maximum_growth_rate, 
        nutrient_half_saturation, photosynthetic_slope, chlorophyll_to_carbon_ratio)

Single nutrient geider photosynthetic growth.

!!! formulation
    γᴿ * γᴾᴬᴿ * P

    where:
    - γᴿ = `monod_limitation(R, kᵣ)`
    - R = nutrient concentration
    - kᵣ = nutrient half saturation
    - γᴾᴬᴿ = `light_limitation_geider(PAR, α, Pᶜₘₐₓ, θᶜ)`
    - PAR = photosynthetic active radiation
    - Pᶜₘₐₓ = maximum growth rate
    - θᶜ = chlorophyll to carbon ratio
    - α = initial photosynthetic slope
    - P = plankton concentration

# Arguments
- `N`: nutrient concentration
- `P`: phytoplankton concentration
- `PAR`: photosynthetic active radiation
- `maximum_growth_rate`: maximum growth rate before nutrient limitation
- `nutrient_half_saturation`: nutrient half saturation
- `photosynthetic_slope`: initial photosynthetic slope
- `chlorophyll_to_carbon_ratio`: ratio between cellular chlorophyll and carbon
"""
@inline function photosynthetic_growth_single_nutrient_geider_light(
    R,
    P,
    PAR,
    maximum_growth_rate,
    nutrient_half_saturation,
    photosynthetic_slope,
    chlorophyll_to_carbon_ratio,
)
    return monod_limitation(R, nutrient_half_saturation) *
           light_limitation_geider(
               PAR, photosynthetic_slope, maximum_growth_rate, chlorophyll_to_carbon_ratio
           ) *
           P
end

"""
    photosynthetic_growth_two_nutrients_geider_light(
        DIN,
        PO4,
        P,
        PAR,
        maximum_growth_rate,
        half_saturation_DIN,
        half_saturation_PO4,
        photosynthetic_slope,
        chlorophyll_to_carbon_ratio,
    )

Two nutrient geider photosynthetic growth.


!!! formulation
    Pᶜₘ *  γᴾᴬᴿ * P

    where:
    - Pᶜₘ = `liebig_minimum([γᴺ, γᴾ])` * Pᶜₘₐₓ
    - γᴾᴬᴿ = `light_limitation_geider(PAR, α, Pᶜₘ, θᶜ)`
    - γᴺ = `monod_limitation(DIN, Kₙ)`
    - γᴾ = `monod_limitation(PO₄, Kₚ)`
    - DIN = dissolved inorganic nitrogen concentration
    - PO₄ = phosphate concentration
    - Kₙ = nitrogen half saturation constant
    - Kₚ = phosphate half saturation constant
    - Pᶜₘₐₓ = maximum growth rate
    - θᶜ = chlorophyll to carbon ratio
    - α = initial photosynthetic slope
    - P = plankton concentration

!!! info
    Unlike the MITgcm-DARWIN formulation this function does not include light inhibition

# Arguments
- `DIN`: dissolved inorganic nitrogen concentration
- `PO4`: phosphate concentration
- `P`: phytoplankton concentration
- `PAR`: photosynthetic active radiation
- `maximum_growth_rate`: maximum growth rate before nutrient limitation
- `half_saturation_DIN`: nitrogen half saturation
- `half_saturation_PO4`: phosphate half saturation
- `photosynthetic_slope`: initial photosynthetic slope
- `chlorophyll_to_carbon_ratio`: ratio between cellular chlorophyll and carbon
"""
@inline function photosynthetic_growth_two_nutrients_geider_light(
    R1,
    R2,
    P,
    PAR,
    maximum_growth_rate,
    nutrient_half_saturation_1,
    nutrient_half_saturation_2,
    photosynthetic_slope,
    chlorophyll_to_carbon_ratio,
)
    nutrient_factor = liebig_minimum(
        monod_limitation(R1, nutrient_half_saturation_1),
        monod_limitation(R2, nutrient_half_saturation_2),
    )

    return nutrient_factor *
           light_limitation_geider(
               PAR, photosynthetic_slope, maximum_growth_rate, chlorophyll_to_carbon_ratio
           ) *
           P
end

# -----------------------------------------------------------------------------
# Symbolic (construction-time) equation blocks
# -----------------------------------------------------------------------------

"""\
    growth_single_nutrient(PV, nutrient_sym, plankton_sym, light_sym, idx; ...) -> AExpr

Symbolic single-nutrient photosynthetic growth for a *group's own* dynamics.
Uses group-owned parameters (missing in that group => error; explicit `nothing` => 0).
"""
function growth_single_nutrient(
    PV,
    nutrient_sym::Symbol,
    plankton_sym::Symbol,
    light_sym::Symbol,
    idx::Int;
    maximum_growth_rate::Symbol=:maximum_growth_rate,
    nutrient_half_saturation::Symbol=:nutrient_half_saturation,
    alpha::Symbol=:alpha,
)
    return _call(
        :photosynthetic_growth_single_nutrient,
        nutrient_sym,
        plankton_sym,
        light_sym,
        _vec(PV, maximum_growth_rate, idx),
        _vec(PV, nutrient_half_saturation, idx),
        _vec(PV, alpha, idx),
    )
end

"""\
    growth_single_nutrient_comm(PV, nutrient_sym, plankton_sym, light_sym, idx; ...) -> AExpr

Community-optional variant used inside community sums.
Missing or `nothing` parameters are treated as inactive (filled with 0).
"""
function growth_single_nutrient_comm(
    PV,
    nutrient_sym::Symbol,
    plankton_sym::Symbol,
    light_sym::Symbol,
    idx::Int;
    maximum_growth_rate::Symbol=:maximum_growth_rate,
    nutrient_half_saturation::Symbol=:nutrient_half_saturation,
    alpha::Symbol=:alpha,
)
    return _call(
        :photosynthetic_growth_single_nutrient,
        nutrient_sym,
        plankton_sym,
        light_sym,
        _vec(PV, maximum_growth_rate, idx),
        _vec(PV, nutrient_half_saturation, idx),
        _vec(PV, alpha, idx),
    )
end

"""Symbolic Geider-style single-nutrient growth for group-owned dynamics."""
function growth_single_nutrient_geider(
    PV,
    nutrient_sym::Symbol,
    plankton_sym::Symbol,
    light_sym::Symbol,
    idx::Int;
    maximum_growth_rate::Symbol=:maximum_growth_rate,
    nutrient_half_saturation::Symbol=:nutrient_half_saturation,
    photosynthetic_slope::Symbol=:photosynthetic_slope,
    chlorophyll_to_carbon_ratio::Symbol=:chlorophyll_to_carbon_ratio,
)
    return _call(
        :photosynthetic_growth_single_nutrient_geider_light,
        nutrient_sym,
        plankton_sym,
        light_sym,
        _vec(PV, maximum_growth_rate, idx),
        _vec(PV, nutrient_half_saturation, idx),
        _vec(PV, photosynthetic_slope, idx),
        _vec(PV, chlorophyll_to_carbon_ratio, idx),
    )
end

"""Community-optional Geider-style single-nutrient growth."""
function growth_single_nutrient_geider_comm(
    PV,
    nutrient_sym::Symbol,
    plankton_sym::Symbol,
    light_sym::Symbol,
    idx::Int;
    maximum_growth_rate::Symbol=:maximum_growth_rate,
    nutrient_half_saturation::Symbol=:nutrient_half_saturation,
    photosynthetic_slope::Symbol=:photosynthetic_slope,
    chlorophyll_to_carbon_ratio::Symbol=:chlorophyll_to_carbon_ratio,
)
    return _call(
        :photosynthetic_growth_single_nutrient_geider_light,
        nutrient_sym,
        plankton_sym,
        light_sym,
        _vec(PV, maximum_growth_rate, idx),
        _vec(PV, nutrient_half_saturation, idx),
        _vec(PV, photosynthetic_slope, idx),
        _vec(PV, chlorophyll_to_carbon_ratio, idx),
    )
end

"""Symbolic Geider-style two-nutrient growth for group-owned dynamics."""
function growth_two_nutrients_geider(
    PV,
    nutrient1_sym::Symbol,
    nutrient2_sym::Symbol,
    plankton_sym::Symbol,
    light_sym::Symbol,
    idx::Int;
    maximum_growth_rate::Symbol=:maximum_growth_rate,
    half_saturation_1::Symbol=:half_saturation_DIN,
    half_saturation_2::Symbol=:half_saturation_PO4,
    photosynthetic_slope::Symbol=:photosynthetic_slope,
    chlorophyll_to_carbon_ratio::Symbol=:chlorophyll_to_carbon_ratio,
)
    return _call(
        :photosynthetic_growth_two_nutrients_geider_light,
        nutrient1_sym,
        nutrient2_sym,
        plankton_sym,
        light_sym,
        _vec(PV, maximum_growth_rate, idx),
        _vec(PV, half_saturation_1, idx),
        _vec(PV, half_saturation_2, idx),
        _vec(PV, photosynthetic_slope, idx),
        _vec(PV, chlorophyll_to_carbon_ratio, idx),
    )
end

"""Community-optional Geider-style two-nutrient growth."""
function growth_two_nutrients_geider_comm(
    PV,
    nutrient1_sym::Symbol,
    nutrient2_sym::Symbol,
    plankton_sym::Symbol,
    light_sym::Symbol,
    idx::Int;
    maximum_growth_rate::Symbol=:maximum_growth_rate,
    half_saturation_1::Symbol=:half_saturation_DIN,
    half_saturation_2::Symbol=:half_saturation_PO4,
    photosynthetic_slope::Symbol=:photosynthetic_slope,
    chlorophyll_to_carbon_ratio::Symbol=:chlorophyll_to_carbon_ratio,
)
    return _call(
        :photosynthetic_growth_two_nutrients_geider_light,
        nutrient1_sym,
        nutrient2_sym,
        plankton_sym,
        light_sym,
        _vec(PV, maximum_growth_rate, idx),
        _vec(PV, half_saturation_1, idx),
        _vec(PV, half_saturation_2, idx),
        _vec(PV, photosynthetic_slope, idx),
        _vec(PV, chlorophyll_to_carbon_ratio, idx),
    )
end

end # module
