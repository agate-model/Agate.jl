# Standalone test used by Buildkite. Also included by test/runtests.jl when
# AGATE_TEST_CUDA=1 to keep the full-suite opt-in available.
using Agate
using CUDA
using Test
using OceanBioME: Biogeochemistry, PrescribedPhotosyntheticallyActiveRadiation
using OceanBioME.Models.NutrientsPlanktonDetritusModels: LOBSTER
using Oceananigans: RectilinearGrid, NonhydrostaticModel, Clock, Center, set!, time_step!
using Oceananigans.Architectures: GPU, array_type
using Oceananigans.Biogeochemistry: biogeochemical_drift_velocity
using Oceananigans.Fields: FunctionField, ZeroField, interior
using Oceananigans.Grids: Periodic, Bounded

@inline gpu_smoke_PAR(x, y, z, t) = 100f0 * exp(0.2f0 * z)

@testset "CUDA smoke: direct Agate model and FrankenLOBSTER" begin
    functional = CUDA.functional()
    @test functional
    if functional
        grid = RectilinearGrid(
            GPU(), Float32;
            topology=(Periodic, Periodic, Bounded),
            size=(2, 2, 4),
            x=(0f0, 2f0), y=(0f0, 2f0), z=(-4f0, 0f0),
        )

        # A small directly authored model exercises the generic GPU compilation path.
        components = (
            N=Agate.Components.Pool(:nitrogen),
            P=Agate.Components.Plankton(;
                states=(nitrogen=:nitrogen,), reference_state=:nitrogen,
                size_structure=[2.0],
            ),
        )
        processes = (growth=Agate.Processes.Growth(;
            plankton=:P, reference_resource=:N,
            bindings=(maximum_rate=:maximum_growth_rate,),
            factors=(light=Agate.Processes.Light(Agate.Processes.Smith(); driver=:PAR),),
        ),)
        parameters = (
            maximum_growth_rate=Agate.Parameters.Parameter(2e-5),
            alpha=Agate.Parameters.Parameter(2e-6),
        )
        definition = Agate.ModelDefinition(; components, processes, parameters)
        sinking_rate = 2.5f0 / 86400f0
        direct = Agate.Construction.construct(definition;
            grid, sinking_tracers=(P_1=sinking_rate,),
        )
        @test direct.parameters.maximum_growth_rate isa array_type(GPU())
        drift = biogeochemical_drift_velocity(direct, Val(:P_1)).w
        @test parent(drift) isa array_type(GPU())
        @test any(==(-sinking_rate), Array(parent(drift)))
        @test biogeochemical_drift_velocity(direct, Val(:N)).w == ZeroField()

        clock = Clock(; time=zero(grid))
        PAR = FunctionField{Center,Center,Center}(gpu_smoke_PAR, grid; clock)
        bgc = Biogeochemistry(
            direct; light_attenuation=PrescribedPhotosyntheticallyActiveRadiation(PAR)
        )
        model = NonhydrostaticModel(grid; clock, biogeochemistry=bgc)
        set!(model; N=7f0, P_1=0.01f0)
        time_step!(model, 60f0)
        @test model.clock.iteration == 1
        @test all(isfinite, Array(interior(model.tracers.P_1)))

        # FrankenLOBSTER: GPU sinking through the composed OceanBioME NPD model.
        franken_sinking_rate = 0.2f0 / 86400f0
        plankton = Agate.Models.FrankenLOBSTER.construct(;
            grid, scalar_type=Float32,
            sinking_tracers=(P_1=franken_sinking_rate,), open_bottom=true,
        )
        @test plankton.runtime.parameters.palatability_matrix isa array_type(GPU())
        @test plankton.runtime.parameters.bacterial_dom_half_saturation isa array_type(GPU())
        @test plankton.settings == plankton.runtime.metadata.model_settings

        franken_clock = Clock(; time=zero(grid))
        franken_PAR = FunctionField{Center,Center,Center}(
            gpu_smoke_PAR, grid; clock=franken_clock
        )
        franken_bgc = LOBSTER(
            grid; plankton,
            light_attenuation=PrescribedPhotosyntheticallyActiveRadiation(franken_PAR),
        )
        franken_npd = franken_bgc.underlying_biogeochemistry
        franken_drift = biogeochemical_drift_velocity(franken_npd, Val(:P_1)).w
        @test parent(franken_drift) isa array_type(GPU())
        @test any(v -> v < 0, Array(parent(franken_drift)))
        @test biogeochemical_drift_velocity(franken_npd, Val(:Z_1)).w == ZeroField()

        franken_model = NonhydrostaticModel(
            grid; clock=franken_clock, biogeochemistry=franken_bgc
        )
        set!(
            franken_model;
            NO₃=2f0, NH₄=0.2f0, DOM=0.1f0, sPOM=0.1f0, bPOM=0.1f0,
            P_1=0.02f0, P_2=0.01f0, Z_1=0.005f0, Z_2=0.005f0,
            H_1=0.01f0, T=20f0,
        )
        time_step!(franken_model, 60f0)
        @test franken_model.clock.iteration == 1
        for tracer in (:P_1, :Z_1, :H_1, :NO₃, :NH₄, :DOM, :sPOM, :bPOM)
            @test all(isfinite, Array(interior(getproperty(franken_model.tracers, tracer))))
        end
    end
end
