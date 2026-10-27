# Run from an environment with Agate installed from a remotely accessible Git commit.
using Agate
using Agate.Construction: export_model_release, encode_recipe
using Agate.Models: NiPiZD
using Test

@testset "Portable release and fresh-process replay" begin
    mktempdir() do workspace
        _, recipe = NiPiZD.construct_plus_recipe(; grid=nothing)
        release = joinpath(workspace, "release")
        @test export_model_release(release, recipe) == release
        @test Set(readdir(release)) == Set(("recipe.json", "Project.toml", "Manifest.toml"))

        # The receiving Julia process installs the copied environment and reconstructs the model.
        replay = """
        using Pkg
        Pkg.instantiate()
        using Agate
        using Agate.Models: NiPiZD
        using Agate.Construction: import_recipe, encode_recipe
        recipe = import_recipe(joinpath(ARGS[1], "recipe.json"))
        @assert encode_recipe(recipe)["content_hash"] == ARGS[2]
        @assert NiPiZD.construct(recipe; grid=nothing).parameters ==
                NiPiZD.construct(; grid=nothing).parameters
        """
        hash = encode_recipe(recipe)["content_hash"]
        run(`$(Base.julia_cmd()) --startup-file=no --project=$release -e $replay $release $hash`)
    end
end
