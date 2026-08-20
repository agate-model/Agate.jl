using Oceananigans.Biogeochemistry:
    required_biogeochemical_auxiliary_fields, required_biogeochemical_tracers
using Test

const LEGACY_NIPIZD_TENDENCIES = LegacyTendencies.TendencyConfig(;
    growth=:smith,
    organic_cycling=:simple_detritus,
    zooplankton=:preferential_grazing,
    nutrient_limitation=:liebig,
    nutrients=(
        LegacyTendencies.nutrient_coupling(
            :N,
            :nutrient_half_saturation;
            remineralization=((:D, :detritus_remineralization),),
        ),
    ),
)

function legacy_nipizd_regression_bgc(parameters)
    context = Agate.Configuration.parse_community(
        Float64,
        default_nipizd_community();
        biogeochem_tracers=(:N, :D),
        interaction_roles=(consumers=(:Z,), prey=(:P,)),
    )
    tracer_order = (:N, :D, Tuple(context.plankton_symbols)...)
    tracers = (
        N=LegacyTendencies.inorganic_tendency(LEGACY_NIPIZD_TENDENCIES; target=:N),
        D=LegacyTendencies.detritus_tendency(LEGACY_NIPIZD_TENDENCIES),
        Z_1=LegacyTendencies.zooplankton_tendency(
            LEGACY_NIPIZD_TENDENCIES; plankton_idx=1
        ),
        Z_2=LegacyTendencies.zooplankton_tendency(
            LEGACY_NIPIZD_TENDENCIES; plankton_idx=2
        ),
        P_1=LegacyTendencies.phytoplankton_tendency(
            LEGACY_NIPIZD_TENDENCIES; plankton_idx=3
        ),
        P_2=LegacyTendencies.phytoplankton_tendency(
            LEGACY_NIPIZD_TENDENCIES; plankton_idx=4
        ),
    )
    tracer_index = Agate.Runtime.build_tracer_index(
        context, tracer_order, (:PAR,); n_biogeochem_tracers=2
    )
    factory = Agate.Construction.define_tracer_functions(
        parameters, tracers; auxiliary_fields=(:PAR,), tracer_index
    )
    return factory(parameters; plankton_diameters=Tuple(context.diameters))
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
    @test manifest.ecological_roles == (phytoplankton=(:P,), zooplankton=(:Z,))

    @test count(x -> !iszero(x), bgc.parameters.maximum_growth_rate) == 2
    @test count(x -> !iszero(x), bgc.parameters.maximum_predation_rate) == 2
    @test bgc.parameters.protection == [1.0, 1.0, 0.0, 0.0]
end

@testset "NiPiZD process compiler matches legacy equations" begin
    canonical = Agate.Models.NiPiZD.construct(; grid=dummy_grid(Float64))
    reference = legacy_nipizd_regression_bgc(canonical.parameters)

    for tracer in NIPIZD_PROCESS_TRACER_ORDER
        @test process_compiler_isapprox(
            canonical(Val(tracer), NIPIZD_PROCESS_ARGS...),
            reference(Val(tracer), NIPIZD_PROCESS_ARGS...),
        )
    end
end
