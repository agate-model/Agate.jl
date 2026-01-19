"""Tracer tendency functors for the simplified DARWIN model.

Predation terms use the canonical rectangular interaction matrices stored in
`bgc.parameters.interactions` (consumer-by-prey). The legacy square
`palatability_matrix[predator_idx, prey_idx]` access pattern remains available
via a zero-padded view, but the DARWIN tracer kernels no longer rely on it.
"""

module Tracers

using ....Functors: CompiledEquation, req

using ....Library.Mortality: LinearLoss, QuadraticLoss
using ....Library.Photosynthesis: TwoNutrientGrowthGeider
using ....Library.Predation: PreferentialPredationLoss, PreferentialPredationGain, PreferentialPredationAssimilationLoss
using ....Library.Remineralization: LinearRemineralization

using ....Utils: sum_over

export DIC_geider_light,
    DIN_geider_light,
    PO4_geider_light,
    POC_default,
    DOC_default,
    PON_default,
    DON_default,
    POP_default,
    DOP_default,
    phytoplankton_growth_two_nutrients_geider_light,
    zooplankton_default

@inline function _state(bgc, args)
    K = keys(bgc.tracer_functions)
    N = length(K)
    vals = ntuple(i -> args[i], N)
    return NamedTuple{K}(vals), N
end

@inline function _aux_value(bgc, args, nstate::Int, name::Symbol)
    aux = bgc.auxiliary_fields
    @inbounds for i in eachindex(aux)
        if aux[i] === name
            return args[nstate + i]
        end
    end
    return nothing
end

@inline function _require_aux(bgc, args, nstate::Int, name::Symbol)
    v = _aux_value(bgc, args, nstate, name)
    v === nothing && throw(ArgumentError("missing required auxiliary field: $(name)"))
    return v
end

@inline function _uptake_sum_geider(p, state, plankton_syms, npl::Int, DIN, PO4, PAR)
    sum_over(npl, zero(DIN)) do i
        P = getproperty(state, plankton_syms[i])
        TwoNutrientGrowthGeider(
            p.maximum_growth_rate[i],
            p.half_saturation_DIN[i],
            p.half_saturation_PO4[i],
            p.photosynthetic_slope[i],
            p.chlorophyll_to_carbon_ratio[i],
        )(DIN, PO4, P, PAR)
    end
end

@inline function _mortality_loss_sum(p, state, plankton_syms, npl::Int, init)
    sum_over(npl, init) do i
        P = getproperty(state, plankton_syms[i])
        LinearLoss(p.linear_mortality[i])(P) + QuadraticLoss(p.quadratic_mortality[i])(P)
    end
end

@inline function _grazing_assimilation_loss_sum(p, state, plankton_syms, init)
    ints = p.interactions
    pal = ints.palatability
    assim = ints.assimilation
    consumer_global = ints.consumer_global
    prey_global = ints.prey_global

    sum_over(eachindex(consumer_global), init) do ic
        predator_idx = consumer_global[ic]
        predator = getproperty(state, plankton_syms[predator_idx])
        gmax = p.maximum_predation_rate[predator_idx]
        K = p.holling_half_saturation[predator_idx]

        sum_over(eachindex(prey_global), zero(init)) do ip
            prey_idx = prey_global[ip]
            prey = getproperty(state, plankton_syms[prey_idx])
            beta = assim[ic, ip]
            phi = pal[ic, ip]
            PreferentialPredationAssimilationLoss(beta, gmax, K, phi)(prey, predator)
        end
    end
end

@inline function _grazing_loss_sum(p, state, plankton_syms, prey, prey_idx::Int, init)
    ints = p.interactions
    ip = @inbounds ints.global_to_prey[prey_idx]
    ip == 0 && return init

    pal = ints.palatability
    consumer_global = ints.consumer_global

    sum_over(eachindex(consumer_global), init) do ic
        predator_idx = consumer_global[ic]
        predator = getproperty(state, plankton_syms[predator_idx])
        gmax = p.maximum_predation_rate[predator_idx]
        K = p.holling_half_saturation[predator_idx]
        phi = pal[ic, ip]
        PreferentialPredationLoss(gmax, K, phi)(prey, predator)
    end
end

@inline function _grazing_gain_sum(p, state, plankton_syms, predator, predator_idx::Int, init)
    ints = p.interactions
    ic = @inbounds ints.global_to_consumer[predator_idx]
    ic == 0 && return init

    pal = ints.palatability
    assim = ints.assimilation
    prey_global = ints.prey_global

    gmax = p.maximum_predation_rate[predator_idx]
    K = p.holling_half_saturation[predator_idx]

    sum_over(eachindex(prey_global), init) do ip
        prey_idx = prey_global[ip]
        prey = getproperty(state, plankton_syms[prey_idx])
        beta = assim[ic, ip]
        phi = pal[ic, ip]
        PreferentialPredationGain(beta, gmax, K, phi)(prey, predator)
    end
end

"""DIC tendency with Geider-style growth (carbon units)."""
function DIC_geider_light(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(
        scalars=(:DOC_remineralization, :POC_remineralization),
        vectors=(
            :maximum_growth_rate,
            :half_saturation_DIN,
            :half_saturation_PO4,
            :photosynthetic_slope,
            :chlorophyll_to_carbon_ratio,
        ),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        state, nstate = _state(bgc, args)
        DIN = state.DIN
        PO4 = state.PO4
        DOC = state.DOC
        POC = state.POC
        PAR = _require_aux(bgc, args, nstate, :PAR)

        uptake = _uptake_sum_geider(p, state, plankton_syms, npl, DIN, PO4, PAR)

        dic_remin = LinearRemineralization(p.DOC_remineralization)(DOC) +
                    LinearRemineralization(p.POC_remineralization)(POC)

        return dic_remin - uptake
    end

    return CompiledEquation(f, requirements)
end

"""DIN tendency assuming fixed stoichiometry (N:C)."""
function DIN_geider_light(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(
        scalars=(:DON_remineralization, :PON_remineralization, :nitrogen_to_carbon),
        vectors=(
            :maximum_growth_rate,
            :half_saturation_DIN,
            :half_saturation_PO4,
            :photosynthetic_slope,
            :chlorophyll_to_carbon_ratio,
        ),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        state, nstate = _state(bgc, args)
        DIN = state.DIN
        PO4 = state.PO4
        DON = state.DON
        PON = state.PON
        PAR = _require_aux(bgc, args, nstate, :PAR)

        uptake = _uptake_sum_geider(p, state, plankton_syms, npl, DIN, PO4, PAR)

        din_remin = LinearRemineralization(p.DON_remineralization)(DON) +
                    LinearRemineralization(p.PON_remineralization)(PON)

        return din_remin - p.nitrogen_to_carbon * uptake
    end

    return CompiledEquation(f, requirements)
end

"""PO4 tendency assuming fixed stoichiometry (P:C)."""
function PO4_geider_light(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(
        scalars=(:DOP_remineralization, :POP_remineralization, :phosphorus_to_carbon),
        vectors=(
            :maximum_growth_rate,
            :half_saturation_DIN,
            :half_saturation_PO4,
            :photosynthetic_slope,
            :chlorophyll_to_carbon_ratio,
        ),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        state, nstate = _state(bgc, args)
        DIN = state.DIN
        PO4 = state.PO4
        DOP = state.DOP
        POP = state.POP
        PAR = _require_aux(bgc, args, nstate, :PAR)

        uptake = _uptake_sum_geider(p, state, plankton_syms, npl, DIN, PO4, PAR)

        po4_remin = LinearRemineralization(p.DOP_remineralization)(DOP) +
                    LinearRemineralization(p.POP_remineralization)(POP)

        return po4_remin - p.phosphorus_to_carbon * uptake
    end

    return CompiledEquation(f, requirements)
end

# --- Organic matter ---------------------------------------------------------

"""DOC tendency from plankton losses and remineralization."""
function DOC_default(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(
        scalars=(:DOM_POM_fractionation, :DOC_remineralization),
        vectors=(:linear_mortality, :quadratic_mortality, :maximum_predation_rate, :holling_half_saturation),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        state, _ = _state(bgc, args)
        DOC = state.DOC

        M = _mortality_loss_sum(p, state, plankton_syms, npl, zero(DOC))
        g = _grazing_assimilation_loss_sum(p, state, plankton_syms, zero(DOC))
        R = LinearRemineralization(p.DOC_remineralization)(DOC)

        frac = one(p.DOM_POM_fractionation) - p.DOM_POM_fractionation
        return frac * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

"""POC tendency from plankton losses and remineralization."""
function POC_default(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(
        scalars=(:DOM_POM_fractionation, :POC_remineralization),
        vectors=(:linear_mortality, :quadratic_mortality, :maximum_predation_rate, :holling_half_saturation),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        state, _ = _state(bgc, args)
        POC = state.POC

        M = _mortality_loss_sum(p, state, plankton_syms, npl, zero(POC))
        g = _grazing_assimilation_loss_sum(p, state, plankton_syms, zero(POC))
        R = LinearRemineralization(p.POC_remineralization)(POC)

        return p.DOM_POM_fractionation * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

"""DON tendency assuming fixed stoichiometry (N:C)."""
function DON_default(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(
        scalars=(:DOM_POM_fractionation, :DON_remineralization, :nitrogen_to_carbon),
        vectors=(:linear_mortality, :quadratic_mortality, :maximum_predation_rate, :holling_half_saturation),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        state, _ = _state(bgc, args)
        DON = state.DON

        M = _mortality_loss_sum(p, state, plankton_syms, npl, zero(DON))
        g = _grazing_assimilation_loss_sum(p, state, plankton_syms, zero(DON))
        R = LinearRemineralization(p.DON_remineralization)(DON)

        frac = one(p.DOM_POM_fractionation) - p.DOM_POM_fractionation
        return frac * p.nitrogen_to_carbon * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

"""PON tendency assuming fixed stoichiometry (N:C)."""
function PON_default(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(
        scalars=(:DOM_POM_fractionation, :PON_remineralization, :nitrogen_to_carbon),
        vectors=(:linear_mortality, :quadratic_mortality, :maximum_predation_rate, :holling_half_saturation),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        state, _ = _state(bgc, args)
        PON = state.PON

        M = _mortality_loss_sum(p, state, plankton_syms, npl, zero(PON))
        g = _grazing_assimilation_loss_sum(p, state, plankton_syms, zero(PON))
        R = LinearRemineralization(p.PON_remineralization)(PON)

        return p.DOM_POM_fractionation * p.nitrogen_to_carbon * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

"""DOP tendency assuming fixed stoichiometry (P:C)."""
function DOP_default(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(
        scalars=(:DOM_POM_fractionation, :DOP_remineralization, :phosphorus_to_carbon),
        vectors=(:linear_mortality, :quadratic_mortality, :maximum_predation_rate, :holling_half_saturation),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        state, _ = _state(bgc, args)
        DOP = state.DOP

        M = _mortality_loss_sum(p, state, plankton_syms, npl, zero(DOP))
        g = _grazing_assimilation_loss_sum(p, state, plankton_syms, zero(DOP))
        R = LinearRemineralization(p.DOP_remineralization)(DOP)

        frac = one(p.DOM_POM_fractionation) - p.DOM_POM_fractionation
        return frac * p.phosphorus_to_carbon * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

"""POP tendency assuming fixed stoichiometry (P:C)."""
function POP_default(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(
        scalars=(:DOM_POM_fractionation, :POP_remineralization, :phosphorus_to_carbon),
        vectors=(:linear_mortality, :quadratic_mortality, :maximum_predation_rate, :holling_half_saturation),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        state, _ = _state(bgc, args)
        POP = state.POP

        M = _mortality_loss_sum(p, state, plankton_syms, npl, zero(POP))
        g = _grazing_assimilation_loss_sum(p, state, plankton_syms, zero(POP))
        R = LinearRemineralization(p.POP_remineralization)(POP)

        return p.DOM_POM_fractionation * p.phosphorus_to_carbon * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

# --- Plankton ---------------------------------------------------------------

"""Phytoplankton tendency with Geider-style, two-nutrient growth."""
function phytoplankton_growth_two_nutrients_geider_light(plankton_syms, plankton_sym::Symbol, plankton_idx::Int)
    requirements = req(
        vectors=(
            :maximum_growth_rate,
            :half_saturation_DIN,
            :half_saturation_PO4,
            :photosynthetic_slope,
            :chlorophyll_to_carbon_ratio,
            :linear_mortality,
            :maximum_predation_rate,
            :holling_half_saturation,
        ),
        matrices=(:palatability_matrix,),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        state, nstate = _state(bgc, args)
        DIN = state.DIN
        PO4 = state.PO4
        P = getproperty(state, plankton_sym)
        PAR = _require_aux(bgc, args, nstate, :PAR)

        growth = TwoNutrientGrowthGeider(
            p.maximum_growth_rate[plankton_idx],
            p.half_saturation_DIN[plankton_idx],
            p.half_saturation_PO4[plankton_idx],
            p.photosynthetic_slope[plankton_idx],
            p.chlorophyll_to_carbon_ratio[plankton_idx],
        )(DIN, PO4, P, PAR)

        grazing = _grazing_loss_sum(p, state, plankton_syms, P, plankton_idx, zero(P))

        mort = LinearLoss(p.linear_mortality[plankton_idx])(P)

        return growth - grazing - mort
    end

    return CompiledEquation(f, requirements)
end

"""Zooplankton tendency with preferential grazing gain."""
function zooplankton_default(plankton_syms, plankton_sym::Symbol, plankton_idx::Int)
    requirements = req(
        vectors=(:linear_mortality, :quadratic_mortality, :maximum_predation_rate, :holling_half_saturation),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        state, _ = _state(bgc, args)
        Z = getproperty(state, plankton_sym)

        gain = _grazing_gain_sum(p, state, plankton_syms, Z, plankton_idx, zero(Z))

        lin = LinearLoss(p.linear_mortality[plankton_idx])(Z)
        quad = QuadraticLoss(p.quadratic_mortality[plankton_idx])(Z)

        return gain - lin - quad
    end

    return CompiledEquation(f, requirements)
end

end # module
