using Oceananigans.Biogeochemistry:
    required_biogeochemical_auxiliary_fields, required_biogeochemical_tracers
using Test

struct LegacyNiPiZDRegressionFactory <: Agate.Factories.AbstractBGCFactory end

const LEGACY_NIPIZD_TENDENCIES = Agate.Tendencies.TendencyConfig(;
    growth=:smith,
    organic_cycling=:simple_detritus,
    zooplankton=:preferential_grazing,
    nutrient_limitation=:liebig,
    nutrients=(
        Agate.Tendencies.nutrient_coupling(
            :N,
            :nutrient_half_saturation;
            remineralization=((:D, :detritus_remineralization),),
        ),
    ),
)

Agate.Factories.parameter_definitions(::LegacyNiPiZDRegressionFactory) =
    Agate.Factories.parameter_definitions(Agate.Models.NiPiZD.NiPiZDFactory())
Agate.Factories.default_community(::LegacyNiPiZDRegressionFactory) =
    Agate.Factories.default_community(Agate.Models.NiPiZD.NiPiZDFactory())
Agate.Factories.default_plankton_dynamics(::LegacyNiPiZDRegressionFactory) = (
    Z=idx -> Agate.Tendencies.zooplankton_tendency(LEGACY_NIPIZD_TENDENCIES; plankton_idx=idx),
    P=idx -> Agate.Tendencies.phytoplankton_tendency(LEGACY_NIPIZD_TENDENCIES; plankton_idx=idx),
)
Agate.Factories.default_biogeochem_dynamics(::LegacyNiPiZDRegressionFactory) = (
    N=() -> Agate.Tendencies.inorganic_tendency(LEGACY_NIPIZD_TENDENCIES; target=:N),
    D=() -> Agate.Tendencies.detritus_tendency(LEGACY_NIPIZD_TENDENCIES),
)

function legacy_nipizd_regression_bgc()
    return Agate.Construction.construct_factory(
        LegacyNiPiZDRegressionFactory();
        interaction_roles=(consumers=(:Z,), prey=(:P,)),
        parameter_roles=(producers=(:P,), consumers=(:Z,)),
        auxiliary_fields=(:PAR,),
        grid=dummy_grid(Float64),
    )
end

@testset "NiPiZD canonical process construction" begin
    bgc, recipe = Agate.Models.NiPiZD.construct_plus_recipe(; grid=dummy_grid(Float64))

    @test recipe isa Agate.Construction.ProcessModelRecipe
    @test required_biogeochemical_tracers(bgc) == (:N, :D, :Z_1, :Z_2, :P_1, :P_2)
    @test required_biogeochemical_auxiliary_fields(bgc) == (:PAR,)
    @test all(
        equation -> equation isa Agate.Compilation.StaticContributionEquation,
        values(bgc.tracer_functions),
    )
    for name in (:NIPIZD_TENDENCIES, :phytoplankton_nipizd, :zooplankton_nipizd)
        @test !isdefined(Agate.Models.NiPiZD, name)
    end

    _, manifest = Agate.Construction.construct_factory_plus_manifest(
        recipe; grid=dummy_grid(Float64)
    )
    @test isempty(manifest.parameter_role_indices)
    @test manifest.ecological_roles == (phytoplankton=(:P,), zooplankton=(:Z,))

    @test count(x -> !iszero(x), bgc.parameters.maximum_growth_rate) == 2
    @test count(x -> !iszero(x), bgc.parameters.maximum_predation_rate) == 2
    @test bgc.parameters.protection == [1.0, 1.0, 0.0, 0.0]
end

@testset "NiPiZD canonical path matches v0.11 construction" begin
    canonical = Agate.Models.NiPiZD.construct(; grid=dummy_grid(Float64))
    reference = legacy_nipizd_regression_bgc()

    @test canonical.parameters == reference.parameters
    for tracer in NIPIZD_PROCESS_TRACER_ORDER
        @test process_compiler_isapprox(
            canonical(Val(tracer), NIPIZD_PROCESS_ARGS...),
            reference(Val(tracer), NIPIZD_PROCESS_ARGS...),
        )
    end
end
