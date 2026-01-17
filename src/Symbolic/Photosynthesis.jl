module Photosynthesis

using ...Equations: AExpr, req, _to_aexpr, merge_requirements

export growth_single_nutrient,
    growth_single_nutrient_comm,
    growth_single_nutrient_geider,
    growth_single_nutrient_geider_comm,
    growth_two_nutrients_geider

@inline _vec(PV, key::Symbol, idx::Int) = getproperty(PV, key)[idx]

"""Build an `AExpr` for a function call, merging argument requirements."""
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
    growth_single_nutrient(PV, nutrient_sym, plankton_sym, light_sym, idx; ...)

Symbolic single-nutrient photosynthetic growth for a group.
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
    μ₀ = _vec(PV, maximum_growth_rate, idx)
    K = _vec(PV, nutrient_half_saturation, idx)
    α = _vec(PV, alpha, idx)

    nutrient_lim = nutrient_sym / (K + nutrient_sym)

    # NOTE: `:==` is not a valid symbol literal in Julia; use `Symbol("==")`.
    eqsym = Symbol("==")
    zeroα = _call(:zero, α)
    condμ = _call(eqsym, μ₀, _call(:zero, μ₀))
    condα = _call(eqsym, α, zeroα)

    denom = _call(:sqrt, μ₀ * μ₀ + α * α * light_sym * light_sym)
    inner = _call(:ifelse, condα, zeroα, (α * light_sym) / denom)
    light_lim = _call(:ifelse, condμ, zeroα, inner)

    return μ₀ * nutrient_lim * light_lim * plankton_sym
end

"""Community-optional variant used in community sums."""
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
    return growth_single_nutrient(
        PV,
        nutrient_sym,
        plankton_sym,
        light_sym,
        idx;
        maximum_growth_rate=maximum_growth_rate,
        nutrient_half_saturation=nutrient_half_saturation,
        alpha=alpha,
    )
end

"""Symbolic Geider-style single-nutrient photosynthetic growth."""
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
    Pᶜₘₐₓ = _vec(PV, maximum_growth_rate, idx)
    K = _vec(PV, nutrient_half_saturation, idx)
    α = _vec(PV, photosynthetic_slope, idx)
    θᶜ = _vec(PV, chlorophyll_to_carbon_ratio, idx)

    nutrient_lim = nutrient_sym / (K + nutrient_sym)

    eqsym = Symbol("==")
    cond = _call(eqsym, Pᶜₘₐₓ, _call(:zero, Pᶜₘₐₓ))
    exp_arg = (-α * θᶜ * light_sym) / Pᶜₘₐₓ
    light_lim = _call(
        :ifelse,
        cond,
        _call(:zero, Pᶜₘₐₓ),
        Pᶜₘₐₓ * (_call(:one, Pᶜₘₐₓ) - _call(:exp, exp_arg)),
    )

    return nutrient_lim * light_lim * plankton_sym
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
    return growth_single_nutrient_geider(
        PV,
        nutrient_sym,
        plankton_sym,
        light_sym,
        idx;
        maximum_growth_rate=maximum_growth_rate,
        nutrient_half_saturation=nutrient_half_saturation,
        photosynthetic_slope=photosynthetic_slope,
        chlorophyll_to_carbon_ratio=chlorophyll_to_carbon_ratio,
    )
end

"""Symbolic Geider-style two-nutrient photosynthetic growth."""
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
    Pᶜₘₐₓ = _vec(PV, maximum_growth_rate, idx)
    K1 = _vec(PV, half_saturation_1, idx)
    K2 = _vec(PV, half_saturation_2, idx)
    α = _vec(PV, photosynthetic_slope, idx)
    θᶜ = _vec(PV, chlorophyll_to_carbon_ratio, idx)

    lim1 = nutrient1_sym / (K1 + nutrient1_sym)
    lim2 = nutrient2_sym / (K2 + nutrient2_sym)
    nutrient_lim = _call(:min, lim1, lim2)

    eqsym = Symbol("==")
    cond = _call(eqsym, Pᶜₘₐₓ, _call(:zero, Pᶜₘₐₓ))
    exp_arg = (-α * θᶜ * light_sym) / Pᶜₘₐₓ
    light_lim = _call(
        :ifelse,
        cond,
        _call(:zero, Pᶜₘₐₓ),
        Pᶜₘₐₓ * (_call(:one, Pᶜₘₐₓ) - _call(:exp, exp_arg)),
    )

    return nutrient_lim * light_lim * plankton_sym
end

end # module
