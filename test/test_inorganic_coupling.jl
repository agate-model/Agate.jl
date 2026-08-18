using Agate
using Test

@testset "Inorganic coupling" begin
    config = MULTI_NUTRIENT_LIEBIG
    sources = ((:DOC, :organic_remineralization), (:POC, :organic_remineralization))

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
