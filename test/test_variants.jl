@testset "Model variants" begin
    id = Agate.Models.ModelId(:DARWIN, :citation2026, :A)

    # Registry should contain the example variants (loaded when DARWIN is included).
    @test id in Agate.Models.list_variants(; family=:DARWIN)

    spec = Agate.Models.variant(id; n_phyto=3, n_zoo=2)

    @test spec.id == id
    @test spec.community.P.n == 3
    @test spec.community.Z.n == 2

    # Variants should construct concrete interaction matrices without requiring user
    # interaction overrides.
    bgc = Agate.Models.construct(spec; grid=dummy_grid(Float32))
    @test size(bgc.parameters.palatability_matrix) == (2, 3)
    @test size(bgc.parameters.assimilation_matrix) == (2, 3)

    # Trait vectors used for derivations must have eltype FT (the grid float type).
    n_total = 3 + 2
    bad_traits = fill(10.0, n_total) # Float64 literals

    try
        Agate.Models.construct(
            spec;
            grid=dummy_grid(Float32),
            parameters=(; optimum_predator_prey_ratio=bad_traits),
        )
        @test false
    catch e
        @test e isa ArgumentError
        msg = sprint(showerror, e)
        @test occursin("expected trait vector", msg)
        @test occursin("Vector{Float32}", msg)
    end

    # Explicit interaction matrices are never overwritten by derived recomputation.
    explicit_assim = fill(0.123f0, 2, 3)
    good_traits = fill(0.5f0, n_total)

    bgc2 = Agate.Models.construct(
        spec;
        grid=dummy_grid(Float32),
        parameters=(; assimilation_efficiency=good_traits),
        interaction_overrides=(; assimilation_matrix=explicit_assim),
    )
    @test bgc2.parameters.assimilation_matrix == explicit_assim
end
