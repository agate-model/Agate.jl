using Agate

const NiPiZD = Agate.Models.NiPiZD

using Test

using Oceananigans.Units
using Agate.Introspection: interaction_matrix, plankton_diameters
using Agate.Library.Allometry: AllometricParam, ConstantParam, PowerLaw
using Oceananigans.Fields: ZeroField
using Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers, required_biogeochemical_auxiliary_fields,
    biogeochemical_drift_velocity

@inline gpu_smoke_PAR(x, y, z, t) = 100f0 * exp(0.2f0 * z)

@testset "Public model constructors" begin

    @testset "NiPiZD defaults" begin
        bgc = NiPiZD.construct(; grid=dummy_grid(Float32))
        tracers = (:N, :D, :P_1, :P_2, :Z_1, :Z_2)
        args = nipizd_runtime_args(Float32)

        @test required_biogeochemical_tracers(bgc) == tracers
        @test required_biogeochemical_tracers(typeof(bgc)) == tracers
        @test required_biogeochemical_auxiliary_fields(bgc) == (:PAR,)
        @test required_biogeochemical_auxiliary_fields(typeof(bgc)) == (:PAR,)
        @test isfinite(@inferred(bgc(Val(:P_1), args...)))
        @test all(tracer -> isfinite(bgc(Val(tracer), args...)), (:N, :D, :P_1, :Z_1))
    end

    @testset "NiPiZD size structure" begin
        phyto_diameters = [2.0, 10.0]
        zoo_diameters = [20.0, 100.0]

        default_size_structure = NiPiZD.construct(;
            size_structure=(;
                phytoplankton=(P=phyto_diameters,),
                zooplankton=(Z=zoo_diameters,),
            ),
            grid=dummy_grid(Float32),
        )
        @test required_biogeochemical_tracers(default_size_structure) ==
            (:N, :D, :P_1, :P_2, :Z_1, :Z_2)

        unsized = (;
            phytoplankton=(P=(n=0,),),
            zooplankton=(Z=(n=0,),),
        )
        allometric_message = argument_error_message(() -> NiPiZD.construct(;
            size_structure=unsized, grid=dummy_grid(Float32)
        ))
        @test occursin("maximum_growth_rate", allometric_message)
        @test occursin("no diameter metadata", allometric_message)

        palatability_message = argument_error_message(() -> NiPiZD.construct(;
            size_structure=unsized,
            grid=dummy_grid(Float32),
            parameters=(;
                maximum_growth_rate=(P=2 / day,),
                nutrient_half_saturation=(P=0.2,),
                maximum_predation_rate=(Z=1 / day,),
            ),
        ))
        @test occursin("AllometricPalatability", palatability_message)
        @test occursin("diameter metadata for SizeClass", palatability_message)

        unsized_bgc = NiPiZD.construct(;
            size_structure=unsized,
            grid=dummy_grid(Float32),
            parameters=(;
                maximum_growth_rate=(P=2 / day,),
                nutrient_half_saturation=(P=0.2,),
                maximum_predation_rate=(Z=1 / day,),
            ),
            palatability_matrix=reshape(Float32[0.5], 1, 1),
        )
        @test required_biogeochemical_tracers(unsized_bgc) == (:N, :D, :P, :Z)
        @test plankton_diameters(unsized_bgc) == [nothing, nothing]

        mixed_bgc = NiPiZD.construct(;
            size_structure=(;
                phytoplankton=(plain=(n=0,), sized=(n=1, min_esd=5.0, max_esd=5.0, spacing=:log)),
                zooplankton=(Z=[20.0],),
            ),
            grid=dummy_grid(Float32),
            parameters=(;
                maximum_growth_rate=(plain=2 / day,),
                nutrient_half_saturation=(plain=0.2,),
            ),
            palatability_matrix=reshape(Float32[0.5, 0.5], 1, 2),
        )
        @test required_biogeochemical_tracers(mixed_bgc) ==
              (:N, :D, :plain, :sized_1, :Z_1)
        @test mixed_bgc.parameters.maximum_growth_rate[1] ≈ Float32(2 / day)
        @test isfinite(mixed_bgc.parameters.maximum_growth_rate[2])
        @test plankton_diameters(mixed_bgc) == [nothing, Float32(5), Float32(20)]

        named = NiPiZD.construct(;
            size_structure=nipizd_named_size_structure(), grid=dummy_grid(Float32)
        )
        tracers = required_biogeochemical_tracers(named)
        @test tracers == (
            :N,
            :D,
            :diat_1,
            :diat_2,
            :diat_3,
            :dino_1,
            :dino_2,
            :mesozoo_1,
            :microzoo_1,
            :microzoo_2,
        )
        @test size(named.parameters.palatability_matrix) == (3, 5)
        @test size(named.parameters.assimilation_matrix) == (3, 5)

        invalid_size_structures = (
            1,
            (; phytoplankton=(diat=[2.0],)),
            (;
                phytoplankton=(diat=[2.0],),
                zooplankton=(microzoo=[30.0],),
                detritus=(small=[1.0],),
            ),
            (; phytoplankton=[2.0], zooplankton=(microzoo=[30.0],)),
            (; phytoplankton=(;), zooplankton=(microzoo=[30.0],)),
            (; phytoplankton=(shared=[2.0],), zooplankton=(shared=[30.0],)),
            (; phytoplankton=(diat=Float64[],), zooplankton=(microzoo=[30.0],)),
            (; phytoplankton=(diat=(n=0, min_esd=1.0, max_esd=2.0, spacing=:log),), zooplankton=(microzoo=[30.0],)),
        )
        for size_structure in invalid_size_structures
            @test_throws ArgumentError NiPiZD.construct(; size_structure)
        end
    end

    @testset "NiPiZD PFT parameter semantics" begin
        phyto_diameters = [2.0, 5.0, 10.0]
        zoo_diameters = [30.0, 60.0, 100.0]
        named_size_structure = (;
            phytoplankton=(diat=phyto_diameters[1:2], dino=phyto_diameters[3:3]),
            zooplankton=(microzoo=zoo_diameters[1:2], mesozoo=zoo_diameters[3:3]),
        )

        named = NiPiZD.construct(;
            size_structure=named_size_structure, grid=dummy_grid(Float32)
        )
        growth = copy(named.parameters.maximum_growth_rate)
        predation = copy(named.parameters.maximum_predation_rate)
        specificity = fill(Float32(0.3), length(plankton_diameters(named)))
        growth[1] = Float32(1.2 / day)
        predation[3] = Float32(0.7 / day)
        # Canonical PFT/SizeClass order is identity-based: diat, dino, mesozoo, microzoo.
        specificity[5] = 2.0f0

        named_overrides = NiPiZD.construct(;
            size_structure=named_size_structure,
            grid=dummy_grid(Float32),
            parameters=(;
                maximum_growth_rate=(diat_1=1.2 / day,),
                maximum_predation_rate=(microzoo_2=0.7 / day,),
                specificity=(microzoo_1=2.0,),
            ),
        )
        positional_overrides = NiPiZD.construct(;
            size_structure=named_size_structure,
            grid=dummy_grid(Float32),
            parameters=(;
                maximum_growth_rate=growth,
                maximum_predation_rate=predation,
                specificity=specificity,
            ),
        )
        @test named_overrides.parameters == positional_overrides.parameters
        @test !hasproperty(named_overrides.parameters, :specificity)

        palatability = reshape(Float32.(1:9), 3, 3)
        assimilation = reshape(Float32.(11:19) ./ 20, 3, 3)
        explicit_interactions = NiPiZD.construct(;
            size_structure=named_size_structure,
            grid=dummy_grid(Float32),
            palatability_matrix=palatability,
            assimilation_matrix=assimilation,
        )
        @test explicit_interactions.parameters.palatability_matrix == palatability
        @test explicit_interactions.parameters.assimilation_matrix == assimilation
        @test_throws ArgumentError NiPiZD.construct(;
            size_structure=named_size_structure, parameters=(palatability_matrix=palatability,),
            palatability_matrix=palatability,
        )

        @test_throws ArgumentError NiPiZD.construct(;
            size_structure=named_size_structure,
            grid=dummy_grid(Float32),
            palatability_matrix=zeros(Float32, 4, 4),
        )
    end

    @testset "NiPiZD interaction parameter overrides" begin
        bgc = NiPiZD.construct(; grid=dummy_grid(Float32))
        axes = interaction_matrix(bgc, :palatability_matrix)
        n_total = length(plankton_diameters(bgc))
        n_cons = length(axes.rows)
        n_prey = length(axes.columns)

        # Rectangular consumer-by-prey matrices are stored as-is.
        rect = reshape(Float32.(1:(n_cons * n_prey)), n_cons, n_prey)
        assimilation = rect ./ (length(rect) + 1)
        bgc_rect = NiPiZD.construct(;
            grid=dummy_grid(Float32), palatability_matrix=rect, assimilation_matrix=assimilation
        )
        @test bgc_rect.parameters.palatability_matrix == rect
        @test bgc_rect.parameters.assimilation_matrix == assimilation

        # Full-square matrices are not accepted for role-aware interaction parameters.
        full = zeros(Float32, n_total, n_total)
        @test_throws ArgumentError NiPiZD.construct(;
            grid=dummy_grid(Float32),
            palatability_matrix=full,
            assimilation_matrix=full,
        )

        # Provider/callable values are not parameter values.
        rect_provider(_) = fill(Float32(9), 1, 1)
        message = argument_error_message(() -> NiPiZD.construct(;
            grid=dummy_grid(Float32),
            palatability_matrix=rect_provider,
            assimilation_matrix=rect_provider,
        ))
        @test occursin("must be a matrix", message)
    end

    @testset "Derived interaction matrices" begin
        # If a model exposes interaction traits, overriding one of those traits
        # should regenerate the derived matrices (unless the matrix itself is
        # explicitly overridden).

        bgc0 = NiPiZD.construct(; grid=dummy_grid(Float32))
        pal0 = bgc0.parameters.palatability_matrix
        specificity = fill(Float32(3), length(plankton_diameters(bgc0)))

        bgc1 = NiPiZD.construct(;
            grid=dummy_grid(Float32), parameters=(; specificity=specificity)
        )
        pal1 = bgc1.parameters.palatability_matrix
        @test any(pal1 .!= pal0)

        rect = fill(Float32(11), size(pal0))
        bgc2 = NiPiZD.construct(;
            grid=dummy_grid(Float32),
            parameters=(; specificity=specificity),
            palatability_matrix=rect,
        )
        @test all(bgc2.parameters.palatability_matrix .== rect)
    end


    @testset "Parameter-law overrides" begin
        phyto_diameters = [2.0, 8.0]
        zoo_diameters = [20.0, 100.0]
        size_structure = (;
            phytoplankton=(P=phyto_diameters,), zooplankton=(Z=zoo_diameters,)
        )

        powerlaw_value(T, prefactor, exponent, diameter) = begin
            d = T(diameter)
            volume = (T(4) / T(3)) * T(π) * (d / T(2))^T(3)
            T(prefactor) * volume^T(exponent)
        end

        growth_prefactor = 2 / day
        growth_exponent = -0.15
        predation_prefactor = 30.84 / day
        predation_exponent = -0.16

        bgc = NiPiZD.construct(;
            size_structure,
            grid=dummy_grid(Float32),
            parameters=(;
                maximum_growth_rate=AllometricParam(
                    PowerLaw(); prefactor=growth_prefactor, exponent=growth_exponent
                ),
                maximum_predation_rate=AllometricParam(
                    PowerLaw(); prefactor=predation_prefactor, exponent=predation_exponent
                ),
            ),
        )

        expected_growth = Float32[
            powerlaw_value(Float32, growth_prefactor, growth_exponent, phyto_diameters[1]),
            powerlaw_value(Float32, growth_prefactor, growth_exponent, phyto_diameters[2]),
        ]
        expected_predation = Float32[
            powerlaw_value(Float32, predation_prefactor, predation_exponent, zoo_diameters[1]),
            powerlaw_value(Float32, predation_prefactor, predation_exponent, zoo_diameters[2]),
        ]

        @test bgc.parameters.maximum_growth_rate ≈ expected_growth
        @test bgc.parameters.maximum_predation_rate ≈ expected_predation
        @test eltype(bgc.parameters.maximum_growth_rate) === Float32
        @test eltype(bgc.parameters.maximum_predation_rate) === Float32

        mortality_prefactor = 0.5 / day
        mortality_exponent = -0.1
        bgc_laws = NiPiZD.construct(;
            size_structure,
            grid=dummy_grid(Float32),
            parameters=(;
                linear_mortality=AllometricParam(
                    PowerLaw();
                    prefactor=mortality_prefactor,
                    exponent=mortality_exponent,
                ),
                maximum_growth_rate=ConstantParam(1.5 / day),
            ),
        )

        all_diameters = [phyto_diameters; zoo_diameters]
        expected_mortality = Float32[
            powerlaw_value(Float32, mortality_prefactor, mortality_exponent, diameter) for
            diameter in all_diameters
        ]
        expected_constant_growth = Float32[1.5 / day, 1.5 / day]

        @test bgc_laws.parameters.linear_mortality ≈ expected_mortality
        @test bgc_laws.parameters.maximum_growth_rate == expected_constant_growth

        bad_relationship = (coeffs, diameter, extra) -> 0.0

        @test_throws ArgumentError(
            "AllometricParam relationship must be callable as relationship(coeffs, diameter)"
        ) NiPiZD.construct(;
            size_structure,
            grid=dummy_grid(Float32),
            parameters=(;
                maximum_growth_rate=AllometricParam(
                    bad_relationship; prefactor=growth_prefactor
                ),
            ),
        )

        message = argument_error_message(() -> NiPiZD.construct(;
            grid=dummy_grid(Float32),
            parameters=(; detritus_remineralization=AllometricParam(
                PowerLaw(); prefactor=1.0, exponent=0.0
            )),
        ))
        @test occursin("detritus_remineralization", message)
        @test occursin("diameter-indexed vector", message)
    end

    @testset "NiPiZD sinking" begin
        sinking_tracers = (P_1=0.1 / day, P_2=0.5 / day, D=2.7489 / day)
        bgc = NiPiZD.construct(; sinking_tracers)

        @test biogeochemical_drift_velocity(bgc, Val(:P_1)).w.data[1, 1, 1] == -0.1 / day
        @test biogeochemical_drift_velocity(bgc, Val(:P_2)).w.data[1, 1, 1] == -0.5 / day
        @test biogeochemical_drift_velocity(bgc, Val(:D)).w.data[1, 1, 1] == -2.7489 / day
        @test biogeochemical_drift_velocity(bgc, Val(:Z_1)).w == ZeroField()

        for invalid in ((not_a_tracer=1.0,), (D=-1.0,), (D=Inf,), (D=true,))
            @test_throws ArgumentError NiPiZD.construct(; sinking_tracers=invalid)
        end
    end

    # Loading CUDA can fail hard in misconfigured environments, so GPU execution is
    # an explicit opt-in test via AGATE_TEST_CUDA=1 rather than part of the default suite.
    if lowercase(get(ENV, "AGATE_TEST_CUDA", "0")) in ("1", "true", "yes")
        @testset "GPU smoke test" begin
            @eval using CUDA
            @eval using OceanBioME: Biogeochemistry, PrescribedPhotosyntheticallyActiveRadiation
            @eval using Oceananigans: RectilinearGrid, NonhydrostaticModel, Clock, Center
            @eval using Oceananigans: set!, time_step!
            @eval using Oceananigans.Fields: FunctionField
            @eval using Oceananigans.Architectures: GPU, array_type
            @eval using Oceananigans.Grids: Periodic, Bounded

            cuda_functional = CUDA.functional()
            @test cuda_functional
            if cuda_functional
                bgc_cpu = NiPiZD.construct(; grid=dummy_grid(Float32))
                bgc_gpu = NiPiZD.construct(; grid=dummy_grid(Float32; arch=GPU()))

                @test required_biogeochemical_tracers(bgc_gpu) ==
                    required_biogeochemical_tracers(bgc_cpu)
                @test bgc_gpu.parameters.palatability_matrix isa array_type(GPU())
                @test bgc_gpu.parameters.maximum_predation_rate isa array_type(GPU())

                grid = RectilinearGrid(
                    GPU(), Float32;
                    topology=(Periodic, Periodic, Bounded),
                    size=(2, 2, 4),
                    x=(0f0, 2f0),
                    y=(0f0, 2f0),
                    z=(-4f0, 0f0),
                )
                sinking_rate = 2.5f0 / 86400f0
                bgc_sinking = NiPiZD.construct(;
                    grid, sinking_tracers=(D=sinking_rate,)
                )
                drift = biogeochemical_drift_velocity(bgc_sinking, Val(:D)).w

                @test parent(drift) isa array_type(GPU())
                @test any(==(-sinking_rate), Array(parent(drift)))
                @test biogeochemical_drift_velocity(bgc_sinking, Val(:Z_1)).w == ZeroField()

                clock = Clock(; time=zero(grid))
                PAR = FunctionField{Center,Center,Center}(gpu_smoke_PAR, grid; clock)
                light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(PAR)
                bgc_model = Biogeochemistry(bgc_sinking; light_attenuation)
                model = NonhydrostaticModel(grid; clock, biogeochemistry=bgc_model)
                set!(
                    model;
                    N=7f0,
                    D=0.01f0,
                    P_1=0.01f0,
                    P_2=0.01f0,
                    Z_1=0.05f0,
                    Z_2=0.05f0,
                )
                time_step!(model, 60f0)
                @test model.clock.iteration == 1
            end
        end
    end

    @testset "Input validation" begin
        message = argument_error_message(() -> NiPiZD.construct(;
            grid=dummy_grid(Float32),
            parameters=(; optimum_predator_prey_ratio=(Z_3=5.0,)),
        ))
        @test occursin("Unknown key `Z_3`", message)
        @test occursin("P_1, P_2, Z_1, Z_2", message)

        for parameters in (
            (; detritus_remineralization=(Z_1=1.0,)),
            (; palatability_matrix=(Z_1=(P_1=1.0,),)),
        )
            @test_throws ArgumentError NiPiZD.construct(; grid=dummy_grid(Float32), parameters)
        end

        # Grid determines precision unless an explicit scalar type is supplied.
        @test NiPiZD.construct(; grid=dummy_grid(Float32)).parameters.detritus_remineralization isa Float32
        bgc_f32 = NiPiZD.construct(; grid=dummy_grid(Float64), scalar_type=Float32)
        @test bgc_f32.parameters.detritus_remineralization isa Float32
        @test eltype(bgc_f32.parameters.maximum_growth_rate) === Float32

        for scalar_type in (Real, ComplexF64, 1.0)
            @test_throws ArgumentError NiPiZD.construct(; scalar_type)
        end
    end
end
