using Test

@testset "Legacy execution architecture removed" begin
    @test !isdefined(Agate, :Tendencies)
    @test all(
        name -> !isdefined(Agate.Factories, name),
        (:default_plankton_dynamics, :default_biogeochem_dynamics, :default_community),
    )
    @test all(
        name -> !isdefined(Agate.Construction, name),
        (:ModelRecipe, :capture_model_recipe),
    )
end
