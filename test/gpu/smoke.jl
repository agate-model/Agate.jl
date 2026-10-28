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

@testset "CUDA smoke: direct Agate model and NPD coupling" begin
    functional = CUDA.functional()
    @test functional
    if functional
        grid = RectilinearGrid(
            GPU(), Float32;
            topology=(Periodic, Periodic, Bounded),
            size=(2, 2, 4),
            x=(0f0, 2f0), y=(0f0, 2f0), z=(-4f0, 0f0),
        )

        # A small directly authored model exercises generic GPU compilation and NPD coupling.
        components = (
            NO₃=Agate.Components.Pool(:nitrogen),
            P=Agate.Components.Plankton(;
                states=(nitrogen=:nitrogen,), reference_state=:nitrogen,
                size_structure=[2.0],
            ),
        )
        processes = (growth=Agate.Processes.Growth(;
            plankton=:P, reference_resource=:NO₃,
            bindings=(maximum_rate=:maximum_growth_rate,),
            factors=(light=Agate.Processes.Light(Agate.Processes.Smith(); driver=:PAR),),
        ),)
        parameters = (
            maximum_growth_rate=Agate.Parameters.Parameter(2e-5),
            alpha=Agate.Parameters.Parameter(2e-6),
        )
        # Include a grazing interaction with a derived matrix to cover GPU adaptation
        # and consumption/product-routing kernels without relying on an external model.
        grazing_components = merge(components, (
            P=Agate.Components.Plankton(;
                states=(nitrogen=:nitrogen,), reference_state=:nitrogen,
                size_structure=[2.0, 4.0],
            ),
            Z=Agate.Components.Plankton(;
                states=(nitrogen=:nitrogen,), reference_state=:nitrogen,
                size_structure=[20.0],
            ),
        ))
        grazing_processes = merge(processes, (
            mortality=Agate.Processes.Mortality(Agate.Processes.LinearMortality();
                plankton=:P, bindings=(rate=:mortality_rate,),
                products=Agate.Processes.Products(:NO₃),
            ),
            grazing=Agate.Processes.Consumption(Agate.Processes.PreferentialGrazing();
                consumers=:Z, resources=:P,
                bindings=(
                    maximum_rate=:maximum_predation_rate,
                    half_saturation=:holling_half_saturation,
                    palatability=:palatability_matrix,
                    assimilation=:assimilation_matrix,
                ),
                unassimilated_products=:NO₃,
            ),
        ))
        CP, D = Agate.Parameters.ConstructionParameter, Agate.Parameters.DerivedDefault
        grazing_parameters = merge(parameters, (
            mortality_rate=Agate.Parameters.Parameter(1e-6),
            maximum_predation_rate=Agate.Parameters.Parameter(5e-5),
            holling_half_saturation=Agate.Parameters.Parameter(0.1),
            optimum_predator_prey_ratio=CP(10.0; axes=:plankton),
            specificity=CP(0.3; axes=:plankton),
            protection=CP(0.0; axes=:plankton),
            assimilation_efficiency=CP(0.7; axes=:plankton),
            palatability_matrix=Agate.Parameters.Parameter(D(
                Agate.Parameters.AllometricPalatability();
                deps=(:optimum_predator_prey_ratio, :specificity, :protection),
            )),
            assimilation_matrix=Agate.Parameters.Parameter(D(
                Agate.Parameters.BroadcastConsumerAcrossResources();
                deps=(:assimilation_efficiency,),
            )),
        ))
        definition = Agate.ModelDefinition(;
            components=grazing_components, processes=grazing_processes,
            parameters=grazing_parameters,
        )
        sinking_rate = 2.5f0 / 86400f0
        direct = Agate.Construction.construct(definition;
            grid, sinking_tracers=(P_1=sinking_rate,),
        )
        @test direct.parameters.maximum_growth_rate isa array_type(GPU())
        @test direct.parameters.palatability_matrix isa array_type(GPU())
        @test direct.parameters.assimilation_matrix isa array_type(GPU())
        @test size(direct.parameters.palatability_matrix) == (1, 2)
        drift = biogeochemical_drift_velocity(direct, Val(:P_1)).w
        @test parent(drift) isa array_type(GPU())
        @test any(==(-sinking_rate), Array(parent(drift)))
        @test biogeochemical_drift_velocity(direct, Val(:NO₃)).w == ZeroField()

        clock = Clock(; time=zero(grid))
        PAR = FunctionField{Center,Center,Center}(gpu_smoke_PAR, grid; clock)
        light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(PAR)
        bgc = Biogeochemistry(direct; light_attenuation)
        model = NonhydrostaticModel(grid; clock, biogeochemistry=bgc)
        set!(model; NO₃=7f0, P_1=0.01f0, P_2=0.01f0, Z_1=0.05f0)
        time_step!(model, 60f0)
        @test model.clock.iteration == 1
        @test all(isfinite, Array(interior(model.tracers.P_1)))
        @test all(isfinite, Array(interior(model.tracers.P_2)))
        @test all(Array(interior(model.tracers.Z_1)) .> 0.05f0)

        # A temperature-dependent variant exercises the composed OceanBioME NPD adapter.
        npd_processes = (growth=Agate.Processes.Growth(;
            plankton=:P, reference_resource=:NO₃,
            bindings=(maximum_rate=:maximum_growth_rate,),
            factors=(
                light=Agate.Processes.Light(Agate.Processes.Smith(); driver=:PAR),
                temperature=Agate.Processes.Temperature(
                    Agate.Processes.Q10(:plankton); driver=:T,
                    bindings=(q10=:temperature_q10, reference_temperature=:reference_temperature),
                ),
            ),
        ),)
        npd_parameters = merge(parameters, (
            temperature_q10=Agate.Parameters.Parameter(2f0),
            reference_temperature=Agate.Parameters.Parameter(20f0),
        ))
        npd_definition = Agate.ModelDefinition(;
            components, processes=npd_processes, parameters=npd_parameters,
        )
        npd_runtime = Agate.Construction.construct(npd_definition;
            grid, sinking_tracers=(P_1=sinking_rate,),
        )
        plankton = Agate.Integrations.NPDPlankton(
            npd_runtime;
            owned_components=(:P,), phytoplankton_components=(:P,),
            nutrient_tracers=(:NO₃,), dependencies=(:T,),
            settings=(carbon_ratio=6.56f0, chlorophyll_ratio=1.31f0),
        )
        npd_clock = Clock(; time=zero(grid))
        npd_PAR = FunctionField{Center,Center,Center}(
            gpu_smoke_PAR, grid; clock=npd_clock
        )
        npd_bgc = LOBSTER(
            grid; plankton,
            light_attenuation=PrescribedPhotosyntheticallyActiveRadiation(npd_PAR),
        )
        npd = npd_bgc.underlying_biogeochemistry
        npd_drift = biogeochemical_drift_velocity(npd, Val(:P_1)).w
        @test parent(npd_drift) isa array_type(GPU())
        @test any(v -> v < 0, Array(parent(npd_drift)))

        npd_model = NonhydrostaticModel(
            grid; clock=npd_clock, biogeochemistry=npd_bgc
        )
        initial = (
            NO₃=2f0, NH₄=0.2f0, DOM=0.1f0, sPOM=0.1f0, bPOM=0.1f0,
            P_1=0.02f0,
        )
        set!(npd_model; initial..., T=20f0)
        time_step!(npd_model, 60f0)
        @test npd_model.clock.iteration == 1
        for tracer in (:P_1, :NO₃, :NH₄, :DOM, :sPOM, :bPOM, :T)
            @test all(isfinite, Array(interior(getproperty(npd_model.tracers, tracer))))
        end

        P20 = Array(interior(npd_model.tracers.P_1))
        set!(npd_model; initial..., T=30f0)
        time_step!(npd_model, 60f0)
        P30 = Array(interior(npd_model.tracers.P_1))
        @test all(P30 .> P20)
    end
end
