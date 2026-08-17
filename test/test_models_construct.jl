using Agate

const NiPiZD = Agate.Models.NiPiZD
const DARWIN = Agate.Models.DARWIN

using Test

using OceanBioME: BoxModelGrid

using Oceananigans.Units
using Agate.Library.Allometry: AllometricParam, ConstantParam, PowerLaw
using Oceananigans.Fields: ZeroField
using Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers, biogeochemical_drift_velocity

struct ThreeInteractionMatrixFactory <: Agate.Factories.AbstractBGCFactory end

function Agate.Factories.parameter_definitions(::ThreeInteractionMatrixFactory)
    return (
        Agate.Factories.ParameterDefinition(
            Agate.Factories.ParameterSpec(
                :encounter_matrix,
                :matrix;
                axes=(:consumer, :prey),
                doc="Test encounter interaction.",
            ),
            Agate.Factories.NoDefault(),
        ),
        Agate.Factories.ParameterDefinition(
            Agate.Factories.ParameterSpec(
                :capture_efficiency_matrix,
                :matrix;
                axes=(:consumer, :prey),
                doc="Test capture-efficiency interaction.",
            ),
            Agate.Factories.NoDefault(),
        ),
        Agate.Factories.ParameterDefinition(
            Agate.Factories.ParameterSpec(
                :handling_time_matrix,
                :matrix;
                axes=(:consumer, :prey),
                doc="Test handling-time interaction.",
            ),
            Agate.Factories.NoDefault(),
        ),
    )
end

@testset "Public model constructors" begin
    @testset "NiPiZD construction input classification" begin
        recipe_fields = Set(NiPiZD.RECIPE_INPUT_FIELDS)
        environment_fields = Set(NiPiZD.ENVIRONMENT_INPUT_FIELDS)

        @test isempty(intersect(recipe_fields, environment_fields))
        @test union(recipe_fields, environment_fields) ==
              Set(fieldnames(NiPiZD.NiPiZDConstructionOptions))
    end

    @testset "DARWIN construction input classification" begin
        recipe_fields = Set(DARWIN.RECIPE_INPUT_FIELDS)
        environment_fields = Set(DARWIN.ENVIRONMENT_INPUT_FIELDS)

        @test isempty(intersect(recipe_fields, environment_fields))
        @test union(recipe_fields, environment_fields) ==
              Set(fieldnames(DARWIN.DARWINConstructionOptions))
    end

    @testset "NiPiZD defaults" begin
        bgc = NiPiZD.construct(; grid=dummy_grid(Float32))

        # Guardrail for GPU compilation: tracer callables must be concretely typed.
        @test !any(t -> t === Any, fieldtypes(typeof(bgc.tracer_functions)))

        @test required_biogeochemical_tracers(bgc) == (:N, :D, :Z_1, :Z_2, :P_1, :P_2)

        P_1 = 0.01f0
        P_2 = 0.01f0
        Z_1 = 0.05f0
        Z_2 = 0.05f0
        N = 7.0f0
        D = 1.0f0
        PAR = 100.0f0

        tracer_vals(sym) =
            if sym === :P_1
                P_1
            elseif sym === :P_2
                P_2
            elseif sym === :Z_1
                Z_1
            elseif sym === :Z_2
                Z_2
            elseif sym === :N
                N
            else
                D
            end

        ordered = [tracer_vals(s) for s in required_biogeochemical_tracers(bgc)]

        @test isfinite(bgc(Val(:N), 0, 0, 0, 0, ordered..., PAR))
        @test isfinite(bgc(Val(:D), 0, 0, 0, 0, ordered..., PAR))
        @test isfinite(bgc(Val(:P_1), 0, 0, 0, 0, ordered..., PAR))
        @test isfinite(bgc(Val(:Z_1), 0, 0, 0, 0, ordered..., PAR))
    end

    @testset "NiPiZD default recipe" begin
        _, recipe = NiPiZD.construct_with_recipe()

        @test recipe.ecological_roles == (phytoplankton=(:P,), zooplankton=(:Z,))
        @test recipe.interaction_roles == (consumers=(:Z,), prey=(:P,))
        @test recipe.parameter_roles == (producers=(:P,), consumers=(:Z,))
        @test isempty(recipe.parameter_overrides)
        @test isnothing(recipe.interaction_overrides)
        @test keys(recipe.interaction_definitions) ==
              (:palatability_matrix, :assimilation_matrix)
        @test isnothing(recipe.sinking_tracers)
        @test recipe.open_bottom
        @test recipe.scalar_type === Float64
        @test !hasproperty(recipe, :grid)
        @test !hasproperty(recipe, :arch)
    end

    @testset "NiPiZD authored recipe inputs" begin
        phyto_diameters = [2.0, 8.0]
        size_structure = (;
            phytoplankton=(diat=phyto_diameters,),
            zooplankton=(;
                microzoo=(n=2, min_esd=30.0, max_esd=90.0, splitting=:log_splitting),
            ),
        )
        growth_override = (diat_2=1.25 / day,)
        mortality_law = AllometricParam(
            PowerLaw(); prefactor=0.05 / day, exponent=-0.1
        )
        alpha_definition = ConstantParam(0.2 / day)
        parameters = (;
            maximum_growth_rate=growth_override,
            linear_mortality=mortality_law,
            alpha=alpha_definition,
        )
        palatability = [0.8 0.2; 0.3 0.7]
        sinking_tracers = (D=2.5 / day,)

        _, recipe = NiPiZD.construct_with_recipe(;
            size_structure,
            scalar_type=Float32,
            parameters,
            palatability_matrix=palatability,
            sinking_tracers,
            open_bottom=false,
        )

        @test recipe.scalar_type === Float32
        @test recipe.sinking_tracers == sinking_tracers
        @test recipe.open_bottom === false
        @test recipe.community.diat.diameters isa Agate.Configuration.DiameterListSpecification
        @test recipe.community.microzoo.diameters isa
              Agate.Configuration.DiameterRangeSpecification
        @test recipe.community.microzoo.diameters.splitting === :log_splitting
        @test recipe.parameter_overrides == parameters
        @test keys(recipe.interaction_overrides) == (:palatability_matrix,)
        @test recipe.interaction_overrides.palatability_matrix == palatability

        phyto_diameters[1] = 999.0
        palatability[1, 1] = 999.0
        @test recipe.community.diat.diameters.diameters[1] == 2.0
        @test recipe.interaction_overrides.palatability_matrix[1, 1] == 0.8
    end

    @testset "NiPiZD exact in-memory replay" begin
        _, recipe, realization = NiPiZD.construct_with_realization()
        _, replayed_realization = NiPiZD.construct_with_realization(recipe)

        @test replayed_realization == realization
        @test (
            realization.group_tracers,
            realization.tracer_order,
            realization.interaction_matrix_sources,
        ) == (
            (Z=(:Z_1, :Z_2), P=(:P_1, :P_2)),
            (:N, :D, :Z_1, :Z_2, :P_1, :P_2),
            (palatability_matrix=:derived, assimilation_matrix=:derived),
        )

        size_structure = (;
            phytoplankton=(diat=[2.0, 8.0],),
            zooplankton=(;
                microzoo=(n=2, min_esd=30.0, max_esd=90.0, splitting=:log_splitting),
            ),
        )
        parameters = (;
            maximum_growth_rate=(diat_2=1.25 / day,),
            linear_mortality=AllometricParam(
                PowerLaw(); prefactor=0.05 / day, exponent=-0.1
            ),
        )

        _, authored_recipe, authored_realization = NiPiZD.construct_with_realization(;
            size_structure,
            scalar_type=Float32,
            parameters,
            palatability_matrix=Float32[0.8 0.2; 0.3 0.7],
            sinking_tracers=(D=2.5f0 / day,),
            open_bottom=false,
        )
        _, replayed_authored_realization = NiPiZD.construct_with_realization(authored_recipe)

        @test replayed_authored_realization == authored_realization
        @test authored_realization.interaction_matrix_sources == (
            palatability_matrix=:explicit, assimilation_matrix=:derived
        )
    end

    @testset "NiPiZD size structure" begin
        phyto_diameters = [2.0, 10.0]
        zoo_diameters = [20.0, 100.0]

        default_groups = NiPiZD.construct(;
            size_structure=(;
                phytoplankton=(P=phyto_diameters,),
                zooplankton=(Z=zoo_diameters,),
            ),
            grid=dummy_grid(Float32),
        )
        @test required_biogeochemical_tracers(default_groups) ==
            (:N, :D, :Z_1, :Z_2, :P_1, :P_2)

        named = NiPiZD.construct(;
            size_structure=(;
                phytoplankton=(diat=[2.0, 5.0, 10.0], dino=[8.0, 20.0]),
                zooplankton=(microzoo=[30.0, 60.0], mesozoo=[100.0]),
            ),
            grid=dummy_grid(Float32),
        )
        tracers = required_biogeochemical_tracers(named)
        @test tracers == (
            :N,
            :D,
            :microzoo_1,
            :microzoo_2,
            :mesozoo_1,
            :diat_1,
            :diat_2,
            :diat_3,
            :dino_1,
            :dino_2,
        )
        @test size(named.parameters.interactions.palatability) == (3, 5)
        @test size(named.parameters.interactions.assimilation) == (3, 5)
        @test count(x -> x != 0, named.parameters.maximum_growth_rate) == 5
        @test count(x -> x != 0, named.parameters.maximum_predation_rate) == 3

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
        )
        for size_structure in invalid_size_structures
            @test_throws ArgumentError NiPiZD.construct(; size_structure)
        end
    end

    @testset "NiPiZD grouped parameter semantics" begin
        phyto_diameters = [2.0, 5.0, 10.0]
        zoo_diameters = [30.0, 60.0, 100.0]
        named_size_structure = (;
            phytoplankton=(diat=phyto_diameters[1:2], dino=phyto_diameters[3:3]),
            zooplankton=(microzoo=zoo_diameters[1:2], mesozoo=zoo_diameters[3:3]),
        )

        function same_parameter_values(a, b)
            keys(a) == keys(b) || return false
            return all(key -> getproperty(a, key) == getproperty(b, key), keys(a))
        end

        flat = NiPiZD.construct(;
            size_structure=(;
                phytoplankton=(P=phyto_diameters,),
                zooplankton=(Z=zoo_diameters,),
            ),
            grid=dummy_grid(Float32),
        )
        named = NiPiZD.construct(;
            size_structure=named_size_structure, grid=dummy_grid(Float32)
        )
        @test same_parameter_values(named.parameters, flat.parameters)

        growth = copy(named.parameters.maximum_growth_rate)
        predation = copy(named.parameters.maximum_predation_rate)
        specificity = copy(named.parameters.specificity)
        growth[4] = Float32(1.2 / day)
        predation[2] = Float32(0.7 / day)
        specificity[1] = 2.0f0

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
        @test same_parameter_values(named_overrides.parameters, positional_overrides.parameters)

        palatability = reshape(Float32.(1:9), 3, 3)
        assimilation = reshape(Float32.(11:19), 3, 3)
        explicit_interactions = NiPiZD.construct(;
            size_structure=named_size_structure,
            grid=dummy_grid(Float32),
            palatability_matrix=palatability,
            assimilation_matrix=assimilation,
        )
        @test explicit_interactions.parameters.interactions.palatability == palatability
        @test explicit_interactions.parameters.interactions.assimilation == assimilation

        @test_throws ArgumentError NiPiZD.construct(;
            size_structure=named_size_structure,
            grid=dummy_grid(Float32),
            palatability_matrix=zeros(Float32, 4, 4),
        )
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

    @testset "InteractionMatrices structural equality" begin
        interactions = NiPiZD.construct(; grid=dummy_grid(Float32)).parameters.interactions

        matrix_names = keys(interactions.matrices)
        copied_matrices = NamedTuple{matrix_names}(
            Tuple(copy(getproperty(interactions.matrices, name)) for name in matrix_names)
        )
        copied = typeof(interactions)(
            copied_matrices,
            copy(interactions.consumer_global),
            copy(interactions.prey_global),
            copy(interactions.global_to_consumer),
            copy(interactions.global_to_prey),
        )
        @test copied == interactions

        function replace_field(interactions, field::Symbol, value)
            names = fieldnames(typeof(interactions))
            values = ntuple(
                i -> names[i] === field ? value : getfield(interactions, i), length(names)
            )
            return typeof(interactions)(values...)
        end

        for field in fieldnames(typeof(interactions))
            changed = if field === :matrices
                name = first(keys(interactions.matrices))
                matrix = copy(getproperty(interactions.matrices, name))
                matrix[1, 1] += one(eltype(matrix))
                merge(interactions.matrices, NamedTuple{(name,)}((matrix,)))
            else
                value = copy(getproperty(interactions, field))
                value[1] = value[1] == 0 ? 1 : 0
                value
            end
            @test replace_field(interactions, field, changed) != interactions
        end
    end

    @testset "Generic interaction matrix collection" begin
        context = Agate.Configuration.CommunityContext(
            Float32,
            3,
            Float32[10, 2, 5],
            [Agate.Configuration.PFTSpecification() for _ in 1:3],
            [:consumer_1, :prey_1, :prey_2],
            [:consumer, :prey, :prey],
            [1, 1, 2],
            Dict(:consumer => [1], :prey => [2, 3]),
            [1],
            [2, 3],
            (producers=[2, 3], consumers=[1]),
            (;),
            (;),
        )
        params = (;
            encounter_matrix=Float32[1 2],
            capture_efficiency_matrix=Float32[0.5 0.75],
            handling_time_matrix=Float32[3 4],
        )

        finalized = Agate.Configuration.finalize_interaction_parameters(
            ThreeInteractionMatrixFactory(), context, params
        )
        interactions = finalized.interactions

        @test keys(interactions.matrices) == (:encounter, :capture_efficiency, :handling_time)
        @test interactions.encounter === finalized.encounter_matrix
        @test interactions.capture_efficiency === finalized.capture_efficiency_matrix
        @test interactions.handling_time === finalized.handling_time_matrix
        @test interactions.consumer_global == [1]
        @test interactions.prey_global == [2, 3]
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
                optimum_predator_prey_ratio=(Z_1=5.0, Z_2=5.0),
                maximum_growth_rate=(P_1=1.2 / day,),
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
            parameters=(; specificity=(Z_1=2.0, Z_2=2.5), half_saturation_DIN=(P_2=0.3,)),
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
                parameters=(; optimum_predator_prey_ratio=(Z_3=5.0,)),
            )
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("Unknown key `Z_3`", sprint(showerror, err))
        @test occursin("Z_1, Z_2, P_1, P_2", sprint(showerror, err))

        @test_throws ArgumentError NiPiZD.construct(;
            grid=dummy_grid(Float32), parameters=(; detritus_remineralization=(Z_1=1.0,))
        )
        @test_throws ArgumentError NiPiZD.construct(;
            grid=dummy_grid(Float32), parameters=(; palatability_matrix=(Z_1=(P_1=1.0,),))
        )
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

        all_diameters = [zoo_diameters; phyto_diameters]
        expected_mortality = Float32[
            powerlaw_value(Float32, mortality_prefactor, mortality_exponent, diameter) for
            diameter in all_diameters
        ]
        expected_constant_growth = Float32[0, 0, 1.5 / day, 1.5 / day]

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

        bgc_full_vector = NiPiZD.construct(;
            size_structure,
            grid=dummy_grid(Float32),
            parameters=(; maximum_growth_rate=expected_growth),
        )
        bgc_named = NiPiZD.construct(;
            size_structure,
            grid=dummy_grid(Float32),
            parameters=(; maximum_growth_rate=(P_1=1.2 / day,)),
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
        bgc = NiPiZD.construct(;
            size_structure=(;
                phytoplankton=(P=[3.0],), zooplankton=(Z=[20.0, 100.0],)
            ),
            grid=dummy_grid(Float32),
        )
        @test required_biogeochemical_tracers(bgc) == (:N, :D, :Z_1, :Z_2, :P_1)
    end

    @testset "NiPiZD sinking" begin
        sinking_tracers = (P_1=0.2551 / day, P_2=0.2551 / day, D=2.7489 / day)
        bgc = NiPiZD.construct(; sinking_tracers)

        @test biogeochemical_drift_velocity(bgc, Val(:P_1)).w.data[1, 1, 1] == -0.2551 / day
        @test biogeochemical_drift_velocity(bgc, Val(:D)).w.data[1, 1, 1] == -2.7489 / day
        @test biogeochemical_drift_velocity(bgc, Val(:Z_1)).w == ZeroField()
    end

    @testset "DARWIN exact in-memory migration replay" begin
        for T in (Float64, Float32)
            grid = dummy_grid(T)
            legacy, recipe, realization = DARWIN.construct_with_realization(; grid, scalar_type=T)
            replayed, replayed_realization = DARWIN.construct_with_realization(recipe; grid)

            @test replayed_realization == realization
            @test replayed.parameters == legacy.parameters
            @test required_biogeochemical_tracers(replayed) ==
                  required_biogeochemical_tracers(legacy)
            @test Agate.Introspection.plankton_diameters(replayed) ==
                  Agate.Introspection.plankton_diameters(legacy)
        end

        grid = BoxModelGrid()
        parameters = (;
            maximum_growth_rate=AllometricParam(
                PowerLaw(); prefactor=1.75f0 / day, exponent=-0.12f0
            ),
            specificity=(Z_1=1.5f0,),
        )
        palatability = Float32[0.8 0.2 0.1; 0.3 0.7 0.4]
        sinking_tracers = (POC=2.5f0 / day,)

        legacy, recipe, realization = DARWIN.construct_with_realization(;
            grid,
            scalar_type=Float32,
            phyto_size_structure=Float32[1.5, 5, 20],
            zoo_size_structure=(;
                n=2, min_esd=20.0f0, max_esd=100.0f0, splitting=:log_splitting
            ),
            parameters,
            palatability_matrix=palatability,
            sinking_tracers,
            open_bottom=false,
        )
        replayed, replayed_realization = DARWIN.construct_with_realization(recipe; grid)

        @test replayed_realization == realization
        @test replayed.parameters == legacy.parameters
        @test realization.interaction_matrix_sources == (
            palatability_matrix=:explicit, assimilation_matrix=:derived
        )
        @test realization.tracer_order == (
            :DIC, :DIN, :PO4, :DOC, :POC, :DON, :PON, :DOP, :POP,
            :Z_1, :Z_2, :P_1, :P_2, :P_3,
        )
        @test biogeochemical_drift_velocity(replayed, Val(:POC)).w.data[1, 1, 1] ==
              biogeochemical_drift_velocity(legacy, Val(:POC)).w.data[1, 1, 1]
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
            size_structure=(;
                phytoplankton=(P=(n=0, min_esd=2, max_esd=10, splitting=:log_splitting),),
                zooplankton=(Z=[20.0, 100.0],),
            )
        )
        @test_throws ArgumentError NiPiZD.construct(;
            size_structure=(;
                phytoplankton=(P=[2.0, 10.0],),
                zooplankton=(Z=(n=0, min_esd=20, max_esd=100, splitting=:linear_splitting),),
            )
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
