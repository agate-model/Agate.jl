"""Construction specifications and runtime parameter containers.

Agate separates data containers into:

- *Specifications* (construction-only, may contain CPU-only metadata)
- *Parameters* (runtime, kernel-safe, Adapt-compatible)

For NiPiZD, user-facing flexibility is represented with `*Specification` types and
expanded into a single `NiPiZDParameters` runtime container.
"""

module Parameters

using Adapt

using Agate.Library.Allometry:
    PalatabilityPreyParameters,
    PalatabilityPredatorParameters,
    allometric_scaling_power,
    allometric_palatability_unimodal_protection

using Agate.Library.Predation:
    AssimilationPreyParameters,
    AssimilationPredatorParameters,
    assimilation_efficiency_emergent_binary

export AbstractDiameterSpecification,
    DiameterListSpecification,
    DiameterRangeSpecification,
    NiPiZDBiogeochemistrySpecification,
    PhytoPFTParameters,
    PhytoSpecification,
    ZooPFTParameters,
    ZooSpecification,
    NiPiZDParameters,
    create_nipizd_parameters,
    compute_nipizd_parameters

"""Abstract supertype for diameter specifications."""
abstract type AbstractDiameterSpecification end

"""A diameter specification defined by an explicit list of diameters."""
struct DiameterListSpecification{T, VT<:AbstractVector{T}} <: AbstractDiameterSpecification
    diameters::VT
end

"""A diameter specification defined by a range and a splitting method."""
struct DiameterRangeSpecification{T} <: AbstractDiameterSpecification
    min_diameter::T
    max_diameter::T
    splitting::Symbol
end

"""PFT-level constants for phytoplankton."""
struct PhytoPFTParameters{FT<:AbstractFloat}
    maximum_growth_rate_a::FT
    maximum_growth_rate_b::FT
    nutrient_half_saturation_a::FT
    nutrient_half_saturation_b::FT
    linear_mortality::FT
    alpha::FT
    photosynthetic_slope::FT
    chlorophyll_to_carbon_ratio::FT
    can_eat::Bool
    can_be_eaten::Bool
    optimum_predator_prey_ratio::FT
    protection::FT
    specificity::FT
    assimilation_efficiency::FT
end

"""PFT-level constants for zooplankton."""
struct ZooPFTParameters{FT<:AbstractFloat}
    maximum_predation_rate_a::FT
    maximum_predation_rate_b::FT
    linear_mortality::FT
    holling_half_saturation::FT
    quadratic_mortality::FT
    can_eat::Bool
    can_be_eaten::Bool
    optimum_predator_prey_ratio::FT
    protection::FT
    specificity::FT
    assimilation_efficiency::FT
end

"""Construction-time specification of phytoplankton size classes."""
struct PhytoSpecification{FT<:AbstractFloat, DS<:AbstractDiameterSpecification}
    n::Int
    diameters::DS
    pft::PhytoPFTParameters{FT}
end

"""Construction-time specification of zooplankton size classes."""
struct ZooSpecification{FT<:AbstractFloat, DS<:AbstractDiameterSpecification}
    n::Int
    diameters::DS
    pft::ZooPFTParameters{FT}
end

"""Construction-time constants for nutrient and detritus cycling."""
struct NiPiZDBiogeochemistrySpecification{FT<:AbstractFloat}
    detritus_remineralization::FT
    mortality_export_fraction::FT
end

"""Runtime NiPiZD parameters stored in the biogeochemistry object."""
struct NiPiZDParameters{FT<:AbstractFloat, VT<:AbstractVector{FT}, MT<:AbstractMatrix{FT}}
    n_P::Int
    n_Z::Int

    diameters::VT

    maximum_growth_rate::VT
    nutrient_half_saturation::VT
    maximum_predation_rate::VT

    linear_mortality::VT
    holling_half_saturation::VT
    quadratic_mortality::VT

    alpha::VT
    photosynthetic_slope::VT
    chlorophyll_to_carbon_ratio::VT

    palatability_matrix::MT
    assimilation_efficiency_matrix::MT

    detritus_remineralization::FT
    mortality_export_fraction::FT
end

Adapt.@adapt_structure NiPiZDParameters

@inline function _check_length(name::Symbol, expected::Int, got::Int)
    if expected != got
        throw(ArgumentError("$(name) must have length $(expected) but has length $(got)"))
    end
    return nothing
end

function _compute_diameters(::Type{FT}, n::Int, spec::DiameterRangeSpecification) where {FT<:AbstractFloat}
    min_d = FT(spec.min_diameter)
    max_d = FT(spec.max_diameter)

    if n == 1
        return FT[min_d]
    end

    diameters = Vector{FT}(undef, n)

    if spec.splitting === :log_splitting
        log_min = log(min_d)
        log_max = log(max_d)
        step = (log_max - log_min) / FT(n - 1)
        @inbounds for i in 1:n
            diameters[i] = exp(log_min + FT(i - 1) * step)
        end
    elseif spec.splitting === :linear_splitting
        step = (max_d - min_d) / FT(n - 1)
        @inbounds for i in 1:n
            diameters[i] = min_d + FT(i - 1) * step
        end
    else
        throw(ArgumentError("Unsupported splitting method: $(spec.splitting)"))
    end

    return diameters
end

function _compute_diameters(::Type{FT}, n::Int, spec::DiameterListSpecification) where {FT<:AbstractFloat}
    _check_length(:diameters, n, length(spec.diameters))
    diameters = Vector{FT}(undef, n)
    @inbounds for i in 1:n
        diameters[i] = FT(spec.diameters[i])
    end
    return diameters
end

@inline function _zoo_pft(zoo::ZooSpecification)
    return zoo.pft
end

@inline function _phyto_pft(phyto::PhytoSpecification)
    return phyto.pft
end

function compute_nipizd_parameters(
    ::Type{FT},
    phyto::PhytoSpecification{FT},
    zoo::ZooSpecification{FT},
    bgc::NiPiZDBiogeochemistrySpecification{FT};
    palatability_matrix=nothing,
    assimilation_efficiency_matrix=nothing,
) where {FT<:AbstractFloat}
    n_P = phyto.n
    n_Z = zoo.n
    n_plankton = n_Z + n_P

    zoo_diameters = _compute_diameters(FT, n_Z, zoo.diameters)
    phyto_diameters = _compute_diameters(FT, n_P, phyto.diameters)

    diameters = Vector{FT}(undef, n_plankton)
    @inbounds begin
        for i in 1:n_Z
            diameters[i] = zoo_diameters[i]
        end
        for i in 1:n_P
            diameters[n_Z + i] = phyto_diameters[i]
        end
    end

    maximum_growth_rate = zeros(FT, n_plankton)
    nutrient_half_saturation = zeros(FT, n_plankton)
    maximum_predation_rate = zeros(FT, n_plankton)

    linear_mortality = zeros(FT, n_plankton)
    holling_half_saturation = zeros(FT, n_plankton)
    quadratic_mortality = zeros(FT, n_plankton)

    alpha = zeros(FT, n_plankton)
    photosynthetic_slope = zeros(FT, n_plankton)
    chlorophyll_to_carbon_ratio = zeros(FT, n_plankton)

    @inbounds for i in 1:n_Z
        zpft = zoo.pft
        linear_mortality[i] = zpft.linear_mortality
        holling_half_saturation[i] = zpft.holling_half_saturation
        quadratic_mortality[i] = zpft.quadratic_mortality

        maximum_predation_rate[i] = allometric_scaling_power(
            zpft.maximum_predation_rate_a,
            zpft.maximum_predation_rate_b,
            diameters[i],
        )
    end

    @inbounds for i in 1:n_P
        idx = n_Z + i
        ppft = phyto.pft

        linear_mortality[idx] = ppft.linear_mortality
        alpha[idx] = ppft.alpha
        photosynthetic_slope[idx] = ppft.photosynthetic_slope
        chlorophyll_to_carbon_ratio[idx] = ppft.chlorophyll_to_carbon_ratio

        maximum_growth_rate[idx] = allometric_scaling_power(
            ppft.maximum_growth_rate_a,
            ppft.maximum_growth_rate_b,
            diameters[idx],
        )

        nutrient_half_saturation[idx] = allometric_scaling_power(
            ppft.nutrient_half_saturation_a,
            ppft.nutrient_half_saturation_b,
            diameters[idx],
        )
    end

    palatability = if isnothing(palatability_matrix)
        M = zeros(FT, n_plankton, n_plankton)
        @inbounds for pred in 1:n_plankton
            pred_is_zoo = pred <= n_Z
            pred_pft = pred_is_zoo ? zoo.pft : phyto.pft

            predator = PalatabilityPredatorParameters{FT}(
                pred_pft.can_eat,
                diameters[pred],
                pred_pft.optimum_predator_prey_ratio,
                pred_pft.specificity,
            )

            for prey in 1:n_plankton
                prey_is_zoo = prey <= n_Z
                prey_pft = prey_is_zoo ? zoo.pft : phyto.pft

                prey_params = PalatabilityPreyParameters{FT}(diameters[prey], prey_pft.protection)
                M[pred, prey] = allometric_palatability_unimodal_protection(prey_params, predator)
            end
        end
        M
    else
        Array{FT}(palatability_matrix)
    end

    assimilation = if isnothing(assimilation_efficiency_matrix)
        M = zeros(FT, n_plankton, n_plankton)
        @inbounds for pred in 1:n_plankton
            pred_is_zoo = pred <= n_Z
            pred_pft = pred_is_zoo ? zoo.pft : phyto.pft

            predator = AssimilationPredatorParameters{FT}(pred_pft.can_eat, pred_pft.assimilation_efficiency)

            for prey in 1:n_plankton
                prey_is_zoo = prey <= n_Z
                prey_pft = prey_is_zoo ? zoo.pft : phyto.pft

                prey_params = AssimilationPreyParameters(prey_pft.can_be_eaten)
                M[pred, prey] = assimilation_efficiency_emergent_binary(prey_params, predator)
            end
        end
        M
    else
        Array{FT}(assimilation_efficiency_matrix)
    end

    return NiPiZDParameters{FT, typeof(diameters), typeof(palatability)}(
        n_P,
        n_Z,
        diameters,
        maximum_growth_rate,
        nutrient_half_saturation,
        maximum_predation_rate,
        linear_mortality,
        holling_half_saturation,
        quadratic_mortality,
        alpha,
        photosynthetic_slope,
        chlorophyll_to_carbon_ratio,
        palatability,
        assimilation,
        bgc.detritus_remineralization,
        bgc.mortality_export_fraction,
    )
end


@inline function _cast_phyto_pft(::Type{FT}, p::PhytoPFTParameters) where {FT<:AbstractFloat}
    return PhytoPFTParameters{FT}(
        FT(p.maximum_growth_rate_a),
        FT(p.maximum_growth_rate_b),
        FT(p.nutrient_half_saturation_a),
        FT(p.nutrient_half_saturation_b),
        FT(p.linear_mortality),
        FT(p.alpha),
        FT(p.photosynthetic_slope),
        FT(p.chlorophyll_to_carbon_ratio),
        Bool(p.can_eat),
        Bool(p.can_be_eaten),
        FT(p.optimum_predator_prey_ratio),
        FT(p.protection),
        FT(p.specificity),
        FT(p.assimilation_efficiency),
    )
end

@inline function _cast_zoo_pft(::Type{FT}, p::ZooPFTParameters) where {FT<:AbstractFloat}
    return ZooPFTParameters{FT}(
        FT(p.maximum_predation_rate_a),
        FT(p.maximum_predation_rate_b),
        FT(p.linear_mortality),
        FT(p.holling_half_saturation),
        FT(p.quadratic_mortality),
        Bool(p.can_eat),
        Bool(p.can_be_eaten),
        FT(p.optimum_predator_prey_ratio),
        FT(p.protection),
        FT(p.specificity),
        FT(p.assimilation_efficiency),
    )
end

@inline function _cast_bgc_spec(::Type{FT}, s::NiPiZDBiogeochemistrySpecification) where {FT<:AbstractFloat}
    return NiPiZDBiogeochemistrySpecification{FT}(
        FT(s.detritus_remineralization),
        FT(s.mortality_export_fraction),
    )
end

"""
    create_nipizd_parameters(::Type{FT}; kwargs...) -> NiPiZDParameters

Convenience constructor that builds specifications and expands them to runtime parameters.
"""
function create_nipizd_parameters(
    ::Type{FT};
    n_phyto::Int,
    n_zoo::Int,
    phyto_diameters::AbstractDiameterSpecification,
    zoo_diameters::AbstractDiameterSpecification,
    phyto_pft_parameters::PhytoPFTParameters,
    zoo_pft_parameters::ZooPFTParameters,
    bgc_specification::NiPiZDBiogeochemistrySpecification;
    palatability_matrix=nothing,
    assimilation_efficiency_matrix=nothing,
) where {FT<:AbstractFloat}
    phyto_pft = _cast_phyto_pft(FT, phyto_pft_parameters)
    zoo_pft = _cast_zoo_pft(FT, zoo_pft_parameters)
    bgc_spec = _cast_bgc_spec(FT, bgc_specification)

    phyto = PhytoSpecification{FT, typeof(phyto_diameters)}(n_phyto, phyto_diameters, phyto_pft)
    zoo = ZooSpecification{FT, typeof(zoo_diameters)}(n_zoo, zoo_diameters, zoo_pft)

    return compute_nipizd_parameters(
        FT,
        phyto,
        zoo,
        bgc_spec;
        palatability_matrix=palatability_matrix,
        assimilation_efficiency_matrix=assimilation_efficiency_matrix,
    )
end

end # module
