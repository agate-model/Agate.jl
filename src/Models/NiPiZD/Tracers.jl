"""Tracer tendency functors for the NiPiZD model.

Predation terms use the canonical rectangular interaction matrices stored in
`bgc.parameters.interactions` (consumer-by-prey). The legacy square
`palatability_matrix[predator_idx, prey_idx]` access pattern remains available
via a zero-padded view, but the NiPiZD tracer kernels no longer rely on it.
"""

module Tracers

using ....Functors: CompiledEquation, req

using ....Library.Mortality: LinearLoss, QuadraticLoss
using ....Library.Photosynthesis: SingleNutrientGrowthSmith
using ....Library.Remineralization: LinearRemineralization

using ....Utils: sum_over

using ...PredationSums: _grazing_assimilation_loss_sum, _grazing_loss_sum, _grazing_gain_sum

export nutrient_default, detritus_default, phytoplankton_default, zooplankton_default

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

@inline function _mortality_loss_sum(p, state, plankton_syms, npl::Int, init)
    sum_over(npl, init) do i
        P = getproperty(state, plankton_syms[i])
        LinearLoss(p.linear_mortality[i])(P) + QuadraticLoss(p.quadratic_mortality[i])(P)
    end
end

@inline function _uptake_sum_smith(p, state, plankton_syms, npl::Int, N, PAR)
    sum_over(npl, zero(N)) do i
        P = getproperty(state, plankton_syms[i])
        SingleNutrientGrowthSmith(
            p.maximum_growth_rate[i], p.nutrient_half_saturation[i], p.alpha[i]
        )(
            N, P, PAR
        )
    end
end

"""Nutrient tendency with Smith growth and mortality/remineralization."""
function nutrient_default(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(;
        scalars=(:detritus_remineralization, :mortality_export_fraction),
        vectors=(
            :linear_mortality,
            :quadratic_mortality,
            :maximum_growth_rate,
            :nutrient_half_saturation,
            :alpha,
        ),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        state, nstate = _state(bgc, args)
        N = state.N
        D = state.D
        PAR = _require_aux(bgc, args, nstate, :PAR)

        lin_sum = sum_over(npl, zero(N)) do i
            P = getproperty(state, plankton_syms[i])
            LinearLoss(p.linear_mortality[i])(P)
        end

        quad_sum = sum_over(npl, zero(N)) do i
            P = getproperty(state, plankton_syms[i])
            QuadraticLoss(p.quadratic_mortality[i])(P)
        end

        uptake = _uptake_sum_smith(p, state, plankton_syms, npl, N, PAR)

        export_frac = p.mortality_export_fraction
        remin = LinearRemineralization(p.detritus_remineralization)(D)

        return export_frac * (lin_sum + quad_sum) + remin - uptake
    end

    return CompiledEquation(f, requirements)
end

"""Detritus tendency from mortality, sloppy feeding, and remineralization."""
function detritus_default(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(;
        scalars=(:detritus_remineralization, :mortality_export_fraction),
        vectors=(
            :linear_mortality,
            :quadratic_mortality,
            :maximum_predation_rate,
            :holling_half_saturation,
        ),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        state, _ = _state(bgc, args)
        D = state.D

        lin_sum = sum_over(npl, zero(D)) do i
            P = getproperty(state, plankton_syms[i])
            LinearLoss(p.linear_mortality[i])(P)
        end

        quad_sum = sum_over(npl, zero(D)) do i
            P = getproperty(state, plankton_syms[i])
            QuadraticLoss(p.quadratic_mortality[i])(P)
        end

        assim_loss = _grazing_assimilation_loss_sum(p, state, plankton_syms, zero(D))

        export_frac = p.mortality_export_fraction
        remin = LinearRemineralization(p.detritus_remineralization)(D)

        return (one(export_frac) - export_frac) * (lin_sum + quad_sum) + assim_loss - remin
    end

    return CompiledEquation(f, requirements)
end

"""Phytoplankton tendency with Smith growth, grazing loss, and linear mortality."""
function phytoplankton_default(plankton_syms, plankton_sym::Symbol, plankton_idx::Int)
    requirements = req(;
        vectors=(
            :maximum_growth_rate,
            :nutrient_half_saturation,
            :alpha,
            :linear_mortality,
            :maximum_predation_rate,
            :holling_half_saturation,
        ),
        matrices=(:palatability_matrix,),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        state, nstate = _state(bgc, args)
        N = state.N
        P = getproperty(state, plankton_sym)
        PAR = _require_aux(bgc, args, nstate, :PAR)

        growth = SingleNutrientGrowthSmith(
            p.maximum_growth_rate[plankton_idx],
            p.nutrient_half_saturation[plankton_idx],
            p.alpha[plankton_idx],
        )(
            N, P, PAR
        )

        grazing = _grazing_loss_sum(p, state, plankton_syms, P, plankton_idx, zero(P))

        mort = LinearLoss(p.linear_mortality[plankton_idx])(P)

        return growth - grazing - mort
    end

    return CompiledEquation(f, requirements)
end

"""Zooplankton tendency with preferential grazing gain and mortality losses."""
function zooplankton_default(plankton_syms, plankton_sym::Symbol, plankton_idx::Int)
    requirements = req(;
        vectors=(
            :linear_mortality,
            :quadratic_mortality,
            :maximum_predation_rate,
            :holling_half_saturation,
        ),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        state, _ = _state(bgc, args)
        Z = getproperty(state, plankton_sym)

        gmax = p.maximum_predation_rate[plankton_idx]
        K = p.holling_half_saturation[plankton_idx]

        gain = _grazing_gain_sum(p, state, plankton_syms, Z, plankton_idx, zero(Z))

        lin = LinearLoss(p.linear_mortality[plankton_idx])(Z)
        quad = QuadraticLoss(p.quadratic_mortality[plankton_idx])(Z)

        return gain - lin - quad
    end

    return CompiledEquation(f, requirements)
end

end # module
