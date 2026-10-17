# Ensure dependency-update PRs exercise the newest declared supported series.
# Each compat entry is a comma-separated list of Julia minor-version floors.
using Pkg, TOML, Test

compat = TOML.parsefile(joinpath(@__DIR__, "..", "..", "Project.toml"))["compat"]
installed = Dict(info.name => info.version for info in values(Pkg.dependencies()) if !isnothing(info.version))

@testset "GPU dependency versions" begin
    for name in ("OceanBioME", "Oceananigans")
        newest = maximum(VersionNumber(strip(part)) for part in split(compat[name], ','))
        actual = installed[name]
        @info "GPU dependency" name actual newest_declared=newest
        @test (actual.major, actual.minor) == (newest.major, newest.minor)
        @test actual >= newest
    end
end
