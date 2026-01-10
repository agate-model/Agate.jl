module Parameters

using Adapt
using Oceananigans.Units

using Agate.Library.Allometry:
    PalatabilityPreyParameters,
    PalatabilityPredatorParameters,
    allometric_scaling_power,
    allometric_palatability_unimodal_protection

using Agate.Library.Predation:
    AssimilationPreyParameters,
    AssimilationPredatorParameters,
    assimilation_efficiency_emergent_binary

using Agate.Utils: param_check_square_matrix, param_cast_matrix, param_compute_diameters

export DarwinBiogeochemistrySpecification,
    DarwinPhytoPFTParameters,
    DarwinPhytoSpecification,
    DarwinZooPFTParameters,
    DarwinZooSpecification,
    DarwinParameterValues,
    create_darwin_parameters,
    compute_darwin_parameters,
    default_darwin_bgc_specification,
    default_darwin_phyto_parameters,
    default_darwin_zoo_parameters

using Agate.Utils:
    AbstractDiameterSpecification, DiameterListSpecification, DiameterRangeSpecification

"""PFT-level constants for phytoplankton."""
struct DarwinPhytoPFTParameters{FT<:AbstractFloat}
    maximum_growth_rate_a::FT
    maximum_growth_rate_b::FT
    half_saturation_DIN_a::FT
    half_saturation_DIN_b::FT
    half_saturation_PO4_a::FT
    half_saturation_PO4_b::FT
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
struct DarwinZooPFTParameters{FT<:AbstractFloat}
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
struct DarwinPhytoSpecification{FT<:AbstractFloat,DS<:AbstractDiameterSpecification}
    n::Int
    diameters::DS
    pft::DarwinPhytoPFTParameters{FT}
end

"""Construction-time specification of zooplankton size classes."""
struct DarwinZooSpecification{FT<:AbstractFloat,DS<:AbstractDiameterSpecification}
    n::Int
    diameters::DS
    pft::DarwinZooPFTParameters{FT}
end

"""Construction-time constants for the simplified DARWIN elemental cycling."""
struct DarwinBiogeochemistrySpecification{FT<:AbstractFloat}
    POC_remineralization::FT
    DOC_remineralization::FT
    PON_remineralization::FT
    DON_remineralization::FT
    POP_remineralization::FT
    DOP_remineralization::FT
    DOM_POM_fractionation::FT
    nitrogen_to_carbon::FT
    phosphorus_to_carbon::FT
end

"""
    DarwinBiogeochemistrySpecification{FT}(; kwargs...) where {FT<:AbstractFloat}

Keyword constructor for `DarwinBiogeochemistrySpecification{FT}`.

This constructor allows users to pass keyword arguments when generating DarwinBiogeochemistrySpecification.

# Keywords
- `POC_remineralization`: Remineralization rate for particulate organic carbon (POC).
- `DOC_remineralization`: Remineralization rate for dissolved organic carbon (DOC).
- `PON_remineralization`: Remineralization rate for particulate organic nitrogen (PON).
- `DON_remineralization`: Remineralization rate for dissolved organic nitrogen (DON).
- `POP_remineralization`: Remineralization rate for particulate organic phosphorus (POP).
- `DOP_remineralization`: Remineralization rate for dissolved organic phosphorus (DOP).
- `DOM_POM_fractionation=0.45`: Fractionation parameter controlling DOM/POM partitioning.
- `nitrogen_to_carbon=0.15`: N:C stoichiometric ratio.
- `phosphorus_to_carbon=0.009`: P:C stoichiometric ratio.

All keyword values are converted to `FT` internally via `FT(x)`.
"""
function DarwinBiogeochemistrySpecification{FT}(;
    POC_remineralization,
    DOC_remineralization,
    PON_remineralization,
    DON_remineralization,
    POP_remineralization,
    DOP_remineralization,
    DOM_POM_fractionation=0.45,
    nitrogen_to_carbon=0.15,
    phosphorus_to_carbon=0.009,
) where {FT<:AbstractFloat}
    return DarwinBiogeochemistrySpecification{FT}(
        FT(POC_remineralization),
        FT(DOC_remineralization),
        FT(PON_remineralization),
        FT(DON_remineralization),
        FT(POP_remineralization),
        FT(DOP_remineralization),
        FT(DOM_POM_fractionation),
        FT(nitrogen_to_carbon),
        FT(phosphorus_to_carbon),
    )
end

"""Default phytoplankton PFT parameter set for Geider-style growth (values chosen to match the original Agate baseline)."""
function default_darwin_phyto_parameters(::Type{FT}) where {FT<:AbstractFloat}
    return DarwinPhytoPFTParameters{FT}(
        FT(2 / day),
        FT(-0.15),
        FT(0.17),
        FT(0.27),
        FT(0.17),
        FT(0.27),        
        FT(8e-7 / second),
        FT(3 / day),
        FT(30),
        FT(0.04),
        false,
        true,
        zero(FT),
        zero(FT),
        zero(FT),
        zero(FT),
    )
end

"""Default zooplankton PFT parameter set (values chosen to match the original Agate baseline)."""
function default_darwin_zoo_parameters(::Type{FT}) where {FT<:AbstractFloat}
    return DarwinZooPFTParameters{FT}(
        FT(0.35 / day),
        FT(0.0),
        FT(0.05 / day),
        FT(0.5),
        FT(0.1 / day),
        true,
        false,
        FT(10),
        one(FT),
        FT(0.3),
        FT(0.32),
    )
end

"""Runtime parameters for the simplified DARWIN model."""
struct DarwinParameterValues{
    FT<:AbstractFloat,VT<:AbstractVector{FT},MT<:AbstractMatrix{FT}
}
    n_P::Int
    n_Z::Int

    diameters::VT

    maximum_growth_rate::VT
    half_saturation_DIN::VT
    half_saturation_PO4::VT
    maximum_predation_rate::VT

    linear_mortality::VT
    holling_half_saturation::VT
    quadratic_mortality::VT

    photosynthetic_slope::VT
    chlorophyll_to_carbon_ratio::VT

    palatability_matrix::MT
    assimilation_efficiency_matrix::MT

    POC_remineralization::FT
    DOC_remineralization::FT
    PON_remineralization::FT
    DON_remineralization::FT
    POP_remineralization::FT
    DOP_remineralization::FT

    DOM_POM_fractionation::FT
    nitrogen_to_carbon::FT
    phosphorus_to_carbon::FT
end

Adapt.@adapt_structure DarwinParameterValues

@inline function _darwin_zoo_pft(zoo::DarwinZooSpecification)
    return zoo.pft
end

@inline function _darwin_phyto_pft(phyto::DarwinPhytoSpecification)
    return phyto.pft
end

@inline function _cast_darwin_phyto_pft(
    ::Type{FT}, p::DarwinPhytoPFTParameters
) where {FT<:AbstractFloat}
    return DarwinPhytoPFTParameters{FT}(
        FT(p.maximum_growth_rate_a),
        FT(p.maximum_growth_rate_b),
        FT(p.half_saturation_DIN_a),
        FT(p.half_saturation_DIN_b),
        FT(p.half_saturation_PO4_a),
        FT(p.half_saturation_PO4_a),        
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

@inline function _cast_darwin_zoo_pft(
    ::Type{FT}, p::DarwinZooPFTParameters
) where {FT<:AbstractFloat}
    return DarwinZooPFTParameters{FT}(
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

function default_darwin_bgc_specification(::Type{FT}) where {FT<:AbstractFloat}
    return DarwinBiogeochemistrySpecification{FT}(;
        POC_remineralization=FT(0.1213 / day),
        DOC_remineralization=FT(0.1213 / day),
        PON_remineralization=FT(0.1213 / day),
        DON_remineralization=FT(0.1213 / day),
        POP_remineralization=FT(0.1213 / day),
        DOP_remineralization=FT(0.1213 / day),
        DOM_POM_fractionation=FT(0.45),
        nitrogen_to_carbon=FT(0.15),
        phosphorus_to_carbon=FT(0.009),
    )
end

@inline function _cast_darwin_bgc_spec(
    ::Type{FT}, s::DarwinBiogeochemistrySpecification
) where {FT<:AbstractFloat}
    return DarwinBiogeochemistrySpecification{FT}(
        FT(s.POC_remineralization),
        FT(s.DOC_remineralization),
        FT(s.PON_remineralization),
        FT(s.DON_remineralization),
        FT(s.POP_remineralization),
        FT(s.DOP_remineralization),
        FT(s.DOM_POM_fractionation),
        FT(s.nitrogen_to_carbon),
        FT(s.phosphorus_to_carbon),
    )
end

"""Expand DARWIN specifications into a runtime `DarwinParameterValues` container."""
function compute_darwin_parameters(
    ::Type{FT},
    phyto::DarwinPhytoSpecification{FT},
    zoo::DarwinZooSpecification{FT},
    bgc::DarwinBiogeochemistrySpecification{FT};
    palatability_matrix=nothing,
    assimilation_efficiency_matrix=nothing,
) where {FT<:AbstractFloat}
    n_P = phyto.n
    n_Z = zoo.n
    n_plankton = n_Z + n_P

    zoo_diameters = param_compute_diameters(FT, n_Z, zoo.diameters)
    phyto_diameters = param_compute_diameters(FT, n_P, phyto.diameters)

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
    half_saturation_DIN = zeros(FT, n_plankton)
    half_saturation_PO4 = zeros(FT, n_plankton)
    maximum_predation_rate = zeros(FT, n_plankton)

    linear_mortality = zeros(FT, n_plankton)
    holling_half_saturation = zeros(FT, n_plankton)
    quadratic_mortality = zeros(FT, n_plankton)

    photosynthetic_slope = zeros(FT, n_plankton)
    chlorophyll_to_carbon_ratio = zeros(FT, n_plankton)

    @inbounds for i in 1:n_Z
        zpft = zoo.pft
        linear_mortality[i] = zpft.linear_mortality
        holling_half_saturation[i] = zpft.holling_half_saturation
        quadratic_mortality[i] = zpft.quadratic_mortality

        maximum_predation_rate[i] = allometric_scaling_power(
            zpft.maximum_predation_rate_a, zpft.maximum_predation_rate_b, diameters[i]
        )
    end

    @inbounds for i in 1:n_P
        idx = n_Z + i
        ppft = phyto.pft

        linear_mortality[idx] = ppft.linear_mortality
        photosynthetic_slope[idx] = ppft.photosynthetic_slope
        chlorophyll_to_carbon_ratio[idx] = ppft.chlorophyll_to_carbon_ratio

        maximum_growth_rate[idx] = allometric_scaling_power(
            ppft.maximum_growth_rate_a, ppft.maximum_growth_rate_b, diameters[idx]
        )

        hs_DIN = allometric_scaling_power(
            ppft.half_saturation_DIN_a, ppft.half_saturation_DIN_b, diameters[idx]
        )

        hs_PO4 = allometric_scaling_power(
            ppft.half_saturation_PO4_a, ppft.half_saturation_PO4_b, diameters[idx]
        )

        half_saturation_DIN[idx] = hs_DIN
        half_saturation_PO4[idx] = hs_PO4
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

                prey_params = PalatabilityPreyParameters{FT}(
                    diameters[prey], prey_pft.protection
                )
                M[pred, prey] = allometric_palatability_unimodal_protection(
                    prey_params, predator
                )
            end
        end
        M
    else
        param_check_square_matrix(:palatability_matrix, n_plankton, palatability_matrix)
        param_cast_matrix(FT, palatability_matrix)
    end

    assimilation = if isnothing(assimilation_efficiency_matrix)
        M = zeros(FT, n_plankton, n_plankton)
        @inbounds for pred in 1:n_plankton
            pred_is_zoo = pred <= n_Z
            pred_pft = pred_is_zoo ? zoo.pft : phyto.pft

            predator = AssimilationPredatorParameters{FT}(
                pred_pft.can_eat, pred_pft.assimilation_efficiency
            )

            for prey in 1:n_plankton
                prey_is_zoo = prey <= n_Z
                prey_pft = prey_is_zoo ? zoo.pft : phyto.pft

                prey_params = AssimilationPreyParameters(prey_pft.can_be_eaten)
                M[pred, prey] = assimilation_efficiency_emergent_binary(
                    prey_params, predator
                )
            end
        end
        M
    else
        param_check_square_matrix(
            :assimilation_efficiency_matrix, n_plankton, assimilation_efficiency_matrix
        )
        param_cast_matrix(FT, assimilation_efficiency_matrix)
    end

    return DarwinParameterValues{FT,typeof(diameters),typeof(palatability)}(
        n_P,
        n_Z,
        diameters,
        maximum_growth_rate,
        half_saturation_DIN,
        half_saturation_PO4,
        maximum_predation_rate,
        linear_mortality,
        holling_half_saturation,
        quadratic_mortality,
        photosynthetic_slope,
        chlorophyll_to_carbon_ratio,
        palatability,
        assimilation,
        bgc.POC_remineralization,
        bgc.DOC_remineralization,
        bgc.PON_remineralization,
        bgc.DON_remineralization,
        bgc.POP_remineralization,
        bgc.DOP_remineralization,
        bgc.DOM_POM_fractionation,
        bgc.nitrogen_to_carbon,
        bgc.phosphorus_to_carbon,
    )
end

"""
    create_darwin_parameters(::Type{FT}; kwargs...) -> DarwinParameterValues

Convenience constructor that builds specifications and expands them to runtime parameters.
"""
function create_darwin_parameters(
    ::Type{FT};
    n_phyto::Int=2,
    n_zoo::Int=2,
    phyto_diameters::AbstractDiameterSpecification=DiameterRangeSpecification(
        2, 10, :log_splitting
    ),
    zoo_diameters::AbstractDiameterSpecification=DiameterRangeSpecification(
        20, 100, :linear_splitting
    ),
    phyto_pft_parameters::DarwinPhytoPFTParameters=default_darwin_phyto_parameters(FT),
    zoo_pft_parameters::DarwinZooPFTParameters=default_darwin_zoo_parameters(FT),
    bgc_specification::DarwinBiogeochemistrySpecification=default_darwin_bgc_specification(
        FT
    ),
    palatability_matrix=nothing,
    assimilation_efficiency_matrix=nothing,
) where {FT<:AbstractFloat}
    phyto_pft = _cast_darwin_phyto_pft(FT, phyto_pft_parameters)
    zoo_pft = _cast_darwin_zoo_pft(FT, zoo_pft_parameters)
    bgc_spec = _cast_darwin_bgc_spec(FT, bgc_specification)

    phyto = DarwinPhytoSpecification{FT,typeof(phyto_diameters)}(
        n_phyto, phyto_diameters, phyto_pft
    )
    zoo = DarwinZooSpecification{FT,typeof(zoo_diameters)}(n_zoo, zoo_diameters, zoo_pft)

    return compute_darwin_parameters(
        FT,
        phyto,
        zoo,
        bgc_spec;
        palatability_matrix=palatability_matrix,
        assimilation_efficiency_matrix=assimilation_efficiency_matrix,
    )
end

end # module
