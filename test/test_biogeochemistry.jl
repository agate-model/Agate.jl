using Agate
using Test

using Agate.Construction: define_tracer_functions
using Agate.Equations: CompiledEquation
using Oceananigans.Biogeochemistry:
    required_biogeochemical_auxiliary_fields, required_biogeochemical_tracers

struct LVParameters{T<:AbstractFloat}
    alpha::T
    beta::T
    delta::T
    gamma::T
end

@testset "Biogeochemistry" begin
    parameters = LVParameters{Float64}(2 / 3, 4 / 3, 1.0, 1.0)
    tracers = (
        R=CompiledEquation((bgc, x, y, z, t, R, F, PAR) -> begin
            p = bgc.parameters
            p.alpha * R - p.beta * R * F
        end),
        F=CompiledEquation((bgc, x, y, z, t, R, F, PAR) -> begin
            p = bgc.parameters
            -p.gamma * F + p.delta * R * F
        end),
    )
    LV = define_tracer_functions(parameters, tracers; auxiliary_fields=(:PAR,))
    model = LV(parameters)
    modified = LV(LVParameters{Float64}(1.0, 1.0, 2.0, 2.0))

    @test required_biogeochemical_tracers(model) == (:R, :F)
    @test required_biogeochemical_tracers(typeof(model)) == (:R, :F)
    @test required_biogeochemical_auxiliary_fields(model) == (:PAR,)
    @test required_biogeochemical_auxiliary_fields(typeof(model)) == (:PAR,)
    @test model(Val(:R), 0, 0, 0, 0, 10, 2, 0) == 2 / 3 * 10 - 4 / 3 * 10 * 2
    @test model(Val(:F), 0, 0, 0, 0, 10, 2, 0) == -1 * 2 + 1 * 10 * 2
    @test modified(Val(:R), 0, 0, 0, 0, 10, 2, 0) == 1 * 10 - 1 * 10 * 2
    @test modified(Val(:F), 0, 0, 0, 0, 10, 2, 0) == -2 * 2 + 2 * 10 * 2
end
