using Agate

const NiPiZD = Agate.Models.NiPiZD
const DARWIN = Agate.Models.DARWIN

using Test

using OceanBioME: BoxModelGrid

using Oceananigans.Units
using Agate.Library.Allometry: AllometricParam, PowerLaw
using Oceananigans.Fields: ZeroField
using Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers, biogeochemical_drift_velocity

@testset "Public model constructors" begin
    @testset "NiPiZD defaults" begin
        bgc = NiPiZD.construct(; grid=dummy_grid(Float32))

        # Guardrail for GPU compilation: tracer callables must be concretely typed.
        @test !any(t -> t === Any, fieldtypes(typeof(bgc.tracer_functions)))

        @test required_biogeochemical_tracers(bgc) == (:N, :D, :Z1, :Z2, :P1, :P2)

        P1 = 0.01f0
        P2 = 0.01f0
        Z1 = 0.05f0
        Z2 = 0.05f0
        N = 7.0f0
        D = 1.0f0
        PAR = 100.0f0

        tracer_vals(sym) =
            if sym === :P1
                P1
            elseif sym === :P2
                P2
            elseif sym === :Z1
                Z1
            elseif sym === :Z2
                Z2
            elseif sym === :N
                N
            else
                D
            end

        ordered = [tracer_vals(s) for s in required_biogeochemical_tracers(bgc)]

        @test isfinite(bgc(Val(:N), 0, 0, 0, 0, ordered..., PAR))
        @test isfinite(bgc(Val(:D), 0, 0, 0, 0, ordered..., PAR))
        @test isfinite(bgc(Val(:P1), 0, 0, 0, 0, ordered..., PAR))
        @test isfinite(bgc(Val(:Z1), 0, 0, 0, 0, ordered..., PAR))
    end

    @testset "NiPiZD interaction overrides" begin
        bgc = NiPiZD.construct(; grid=dummy_grid(Float32))
        ints0 = bgc.parameters.interactions
        n_total = length(ints0.global_to_prey)
        n_cons = length(ints0.consumer_global)
        n_prey = length(ints0.prey_global)

        wrong = zeros(Float32, 3, 3)
        @test_throws ArgumentError NiPiZD.construct(;
            grid=dummy_grid(Float32), palatability_matrix=wrong, assimilation_matrix=wrong
        )

        # Rectangular consumer-by-prey matrices are stored as-is.
        rect = reshape(Float32.(1:(n_cons * n_prey)), n_cons, n_prey)
        bgc_rect = NiPiZD.construct(;
            grid=dummy_grid(Float32), palatability_matrix=rect, assimilation_matrix=rect
        )
        M = bgc_rect.parameters.palatability_matrix
        @test size(M) == (n_cons, n_prey)
        @test all(M .== rect)

        # Non-axis-sized matrices are not accepted; provide an explicit
        # axis-sized matrix (n_consumer, n_prey) instead.
        axis_block = reshape(Float32[7], 1, 1)
        @test_throws ArgumentError NiPiZD.construct(;
            grid=dummy_grid(Float32),
            palatability_matrix=axis_block,
            assimilation_matrix=axis_block,
        )

        # Full-square matrices are not accepted for role-aware interaction overrides.
        correct = zeros(Float32, n_total, n_total)
        @test_throws ArgumentError NiPiZD.construct(;
            grid=dummy_grid(Float32),
            palatability_matrix=correct,
            assimilation_matrix=correct,
        )

        # Provider/callable values are not supported for interaction overrides.
        rect_provider(ctx) = fill(
            Float32(9), length(ctx.consumer_indices), length(ctx.prey_indices)
        )
        err = try
            NiPiZD.construct(;
                grid=dummy_grid(Float32),
                palatability_matrix=rect_provider,
                assimilation_matrix=rect_provider,
            )
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("providers are not supported", sprint(showerror, err))
    end

    @testset "Derived interaction matrices" begin
        # If a model exposes interaction traits, overriding one of those traits
        # should regenerate the derived matrices (unless the matrix itself is
        # explicitly overridden).

        bgc0 = NiPiZD.construct(; grid=dummy_grid(Float32))
        pal0 = bgc0.parameters.interactions.palatability
        n_total = length(bgc0.parameters.interactions.global_to_prey)

        specificity = zeros(Float32, n_total)
        specificity[bgc0.parameters.interactions.consumer_global] .= 3.0f0

        bgc1 = NiPiZD.construct(;
            grid=dummy_grid(Float32), parameters=(; specificity=specificity)
        )
        pal1 = bgc1.parameters.interactions.palatability
        @test any(pal1 .!= pal0)

        rect = fill(Float32(11), size(pal0))
        bgc2 = NiPiZD.construct(;
            grid=dummy_grid(Float32),
            parameters=(; specificity=specificity),
            palatability_matrix=rect,
        )
        @test all(bgc2.parameters.interactions.palatability .== rect)

        dar0 = DARWIN.construct(; grid=dummy_grid(Float32))
        dar_pal0 = dar0.parameters.interactions.palatability
        dar_n_total = length(dar0.parameters.interactions.global_to_prey)
        dar_spec = zeros(Float32, dar_n_total)
        dar_spec[dar0.parameters.interactions.consumer_global] .= 2.0f0
        dar1 = DARWIN.construct(;
            grid=dummy_grid(Float32), parameters=(; specificity=dar_spec)
        )
        dar_pal1 = dar1.parameters.interactions.palatability
        @test any(dar_pal1 .!= dar_pal0)
    end


    @testset "Named parameter vector overrides" begin
        bgc_default = NiPiZD.construct(; grid=dummy_grid(Float32))
        vopt = copy(bgc_default.parameters.optimum_predator_prey_ratio)
        vopt[1] = 5.0f0
        vopt[2] = 5.0f0
        growth = copy(bgc_default.parameters.maximum_growth_rate)
        growth[3] = Float32(1.2 / day)

        bgc_named = NiPiZD.construct(;
            grid=dummy_grid(Float32),
            parameters=(;
                optimum_predator_prey_ratio=(Z1=5.0, Z2=5.0),
                maximum_growth_rate=(P1=1.2 / day,),
            ),
        )
        bgc_positional = NiPiZD.construct(;
            grid=dummy_grid(Float32),
            parameters=(; optimum_predator_prey_ratio=vopt, maximum_growth_rate=growth),
        )

        @test bgc_named.parameters.optimum_predator_prey_ratio == vopt
        @test bgc_named.parameters.maximum_growth_rate == growth
        @test bgc_named.parameters.optimum_predator_prey_ratio ==
            bgc_positional.parameters.optimum_predator_prey_ratio
        @test bgc_named.parameters.maximum_growth_rate ==
            bgc_positional.parameters.maximum_growth_rate
        @test bgc_named.parameters.interactions.palatability ==
            bgc_positional.parameters.interactions.palatability

        dar_default = DARWIN.construct(; grid=dummy_grid(Float32))
        dar_spec = copy(dar_default.parameters.specificity)
        dar_spec[1] = 2.0f0
        dar_spec[2] = 2.5f0
        dar_din = copy(dar_default.parameters.half_saturation_DIN)
        dar_din[4] = 0.3f0

        dar_named = DARWIN.construct(;
            grid=dummy_grid(Float32),
            parameters=(; specificity=(Z1=2.0, Z2=2.5), half_saturation_DIN=(P2=0.3,)),
        )
        dar_positional = DARWIN.construct(;
            grid=dummy_grid(Float32),
            parameters=(; specificity=dar_spec, half_saturation_DIN=dar_din),
        )

        @test dar_named.parameters.specificity == dar_spec
        @test dar_named.parameters.half_saturation_DIN == dar_din
        @test dar_named.parameters.interactions.palatability ==
            dar_positional.parameters.interactions.palatability

        err = try
            NiPiZD.construct(;
                grid=dummy_grid(Float32),
                parameters=(; optimum_predator_prey_ratio=(Z3=5.0,)),
            )
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("Unknown key `Z3`", sprint(showerror, err))
        @test occursin("Z1, Z2, P1, P2", sprint(showerror, err))

        @test_throws ArgumentError NiPiZD.construct(;
            grid=dummy_grid(Float32), parameters=(; detritus_remineralization=(Z1=1.0,))
        )
        @test_throws ArgumentError NiPiZD.construct(;
            grid=dummy_grid(Float32), parameters=(; palatability_matrix=(Z1=(P1=1.0,),))
        )
    end

    @testset "Allometric parameter overrides" begin
        phyto_diameters = [2.0, 8.0]
        zoo_diameters = [20.0, 100.0]

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
            phyto_size_structure=phyto_diameters,
            zoo_size_structure=zoo_diameters,
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
            0,
            0,
            powerlaw_value(Float32, growth_prefactor, growth_exponent, phyto_diameters[1]),
            powerlaw_value(Float32, growth_prefactor, growth_exponent, phyto_diameters[2]),
        ]
        expected_predation = Float32[
            powerlaw_value(Float32, predation_prefactor, predation_exponent, zoo_diameters[1]),
            powerlaw_value(Float32, predation_prefactor, predation_exponent, zoo_diameters[2]),
            0,
            0,
        ]

        @test bgc.parameters.maximum_growth_rate ≈ expected_growth
        @test bgc.parameters.maximum_predation_rate ≈ expected_predation
        @test eltype(bgc.parameters.maximum_growth_rate) === Float32
        @test eltype(bgc.parameters.maximum_predation_rate) === Float32

        bgc_full_vector = NiPiZD.construct(;
            phyto_size_structure=phyto_diameters,
            zoo_size_structure=zoo_diameters,
            grid=dummy_grid(Float32),
            parameters=(; maximum_growth_rate=expected_growth),
        )
        bgc_named = NiPiZD.construct(;
            phyto_size_structure=phyto_diameters,
            zoo_size_structure=zoo_diameters,
            grid=dummy_grid(Float32),
            parameters=(; maximum_growth_rate=(P1=1.2 / day,)),
        )

        @test bgc_full_vector.parameters.maximum_growth_rate == expected_growth
        @test bgc_named.parameters.maximum_growth_rate[3] == Float32(1.2 / day)

        err = try
            NiPiZD.construct(;
                grid=dummy_grid(Float32),
                parameters=(;
                    detritus_remineralization=AllometricParam(
                        PowerLaw(); prefactor=1.0, exponent=0.0
                    ),
                ),
            )
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("detritus_remineralization", sprint(showerror, err))
        @test occursin("diameter-indexed vector", sprint(showerror, err))
    end


    @testset "NiPiZD community structure overrides" begin
        bgc = NiPiZD.construct(; phyto_size_structure=[3.0], grid=dummy_grid(Float32))
        @test required_biogeochemical_tracers(bgc) == (:N, :D, :Z1, :Z2, :P1)
    end

    @testset "NiPiZD sinking" begin
        bgc = NiPiZD.construct(;
            sinking_tracers=(P1=0.2551 / day, P2=0.2551 / day, D=2.7489 / day)
        )

        @test biogeochemical_drift_velocity(bgc, Val(:P1)).w.data[1, 1, 1] == -0.2551 / day
        @test biogeochemical_drift_velocity(bgc, Val(:D)).w.data[1, 1, 1] == -2.7489 / day
        @test biogeochemical_drift_velocity(bgc, Val(:Z1)).w == ZeroField()
    end

    @testset "DARWIN defaults" begin
        bgc = DARWIN.construct(; grid=dummy_grid(Float32))

        @test !any(t -> t === Any, fieldtypes(typeof(bgc.tracer_functions)))

        @test required_biogeochemical_tracers(bgc)[1:9] ==
            (:DIC, :DIN, :PO4, :DOC, :POC, :DON, :PON, :DOP, :POP)

        bgc_explicit = DARWIN.construct(; scalar_type=Float32)
        @test bgc_explicit.parameters.DOC_remineralization isa Float32
        @test eltype(bgc_explicit.parameters.maximum_growth_rate) === Float32
    end

    @testset "GPU smoke test" begin
        # NOTE: Loading CUDA can crash Julia in misconfigured environments (e.g. mixed system/toolkit libs).
        # To keep the default test suite robust, this test only runs when explicitly enabled.
        if lowercase(get(ENV, "AGATE_TEST_CUDA", "0")) in ("1", "true", "yes")
            @eval using CUDA
            @eval using Oceananigans.Architectures: GPU, array_type

            if CUDA.functional()
                bgc_cpu = NiPiZD.construct(; grid=dummy_grid(Float32))
                bgc_gpu = NiPiZD.construct(; grid=dummy_grid(Float32; arch=GPU()))

                @test required_biogeochemical_tracers(bgc_gpu) ==
                    required_biogeochemical_tracers(bgc_cpu)
                @test bgc_gpu.parameters.interactions.palatability isa array_type(GPU())
                @test bgc_gpu.parameters.maximum_predation_rate isa array_type(GPU())
            else
                @test true
            end
        else
            @test true
        end
    end

    @testset "Input validation" begin
        @test_throws ArgumentError NiPiZD.construct(;
            phyto_size_structure=(n=0, min_esd=2, max_esd=10, splitting=:log_splitting)
        )
        @test_throws ArgumentError NiPiZD.construct(;
            zoo_size_structure=(n=0, min_esd=20, max_esd=100, splitting=:linear_splitting)
        )

        # Grid determines precision unless an explicit scalar type is supplied.
        bgc_f32 = NiPiZD.construct(; grid=dummy_grid(Float32))
        @test bgc_f32.parameters.detritus_remineralization isa Float32

        bgc_explicit_f32 = NiPiZD.construct(; grid=dummy_grid(Float64), scalar_type=Float32)
        @test bgc_explicit_f32.parameters.detritus_remineralization isa Float32
        @test eltype(bgc_explicit_f32.parameters.maximum_growth_rate) === Float32

        @test_throws ArgumentError NiPiZD.construct(; scalar_type=Real)
        @test_throws ArgumentError NiPiZD.construct(; scalar_type=ComplexF64)
        @test_throws ArgumentError NiPiZD.construct(; scalar_type=1.0)

        # Wrong interaction matrix sizes should error.
        bgc = NiPiZD.construct(; grid=dummy_grid(Float64))
        n_cons = size(bgc.parameters.palatability_matrix, 1)
        n_prey = size(bgc.parameters.palatability_matrix, 2)
        wrong = zeros(Float64, n_cons + 1, n_prey + 1)
        @test_throws ArgumentError NiPiZD.construct(;
            grid=dummy_grid(Float64), palatability_matrix=wrong, assimilation_matrix=wrong
        )
    end
end
