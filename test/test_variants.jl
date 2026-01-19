@testset "Model variants" begin
    id = Agate.Models.ModelId(:DARWIN, :citation2026, :A)

    # Registry should contain the example variant (loaded when DARWIN is included).
    @test id in Agate.Models.list_variants(family = :DARWIN)

    spec = Agate.Models.variant(id; n_phyto = 3, n_zoo = 2)

    @test spec.id == id
    @test spec.community.P.n == 3
    @test spec.community.Z.n == 2
end
