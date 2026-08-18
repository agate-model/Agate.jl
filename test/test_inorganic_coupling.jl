using Agate
using Test
using Logging

@testset "Inorganic coupling" begin
    config = MULTI_NUTRIENT_LIEBIG
    sources = ((:DOC, :organic_remineralization), (:POC, :organic_remineralization))

    @test_logs min_level=Logging.Warn Agate.Tendencies.inorganic_tendency(
        config; target=:DIC, remineralization=sources
    )

    @test Agate.Tendencies.inorganic_target_coupling(
        config, :DIC, sources, nothing
    ) == (sources, :one)

    @test Agate.Tendencies.inorganic_target_coupling(
        config, :DIN, nothing, nothing
    ) == (
        ((:DON, :organic_remineralization), (:PON, :organic_remineralization)),
        :nitrogen_to_carbon,
    )

    @test_throws ArgumentError Agate.Tendencies.inorganic_tendency(config; target=:DIC)
end
