using Agate
using Agate.Introspection:
    tracer_names,
    parameter_names,
    plankton_groups,
    plankton_tracers,
    plankton_diameters,
    nonplankton_tracers,
    tracer_groups,
    interaction_matrix
using Test

using Agate.Construction: define_tracer_functions
using Agate.Equations: CompiledEquation

@testset "Public introspection helpers" begin
    @testset "Model-constructed instance" begin
        bgc = Agate.Models.NiPiZD.construct(; grid=dummy_grid(Float32))

        @test tracer_names(bgc) == [:N, :D, :Z_1, :Z_2, :P_1, :P_2]

        groups = tracer_groups(bgc)
        @test groups.all == tracer_names(bgc)
        @test groups.by_group.Z == [:Z_1, :Z_2]
        @test groups.by_group.P == [:P_1, :P_2]
        @test groups.plankton == [:Z_1, :Z_2, :P_1, :P_2]
        @test groups.nonplankton == [:N, :D]
        @test plankton_groups(bgc) == groups.by_group
        @test plankton_tracers(bgc) == groups.plankton
        @test length(plankton_diameters(bgc)) == length(groups.plankton)
        @test plankton_diameters(bgc) == collect(bgc.plankton_diameters)
        @test eltype(plankton_diameters(bgc)) === Float32
        @test nonplankton_tracers(bgc) == groups.nonplankton

        pars = parameter_names(bgc)
        @test !isempty(pars)
        @test :data ∉ pars
        @test :palatability_matrix in pars
        @test :assimilation_matrix in pars
        @test :interactions ∉ pars

        pal = interaction_matrix(bgc, :palatability_matrix)
        assim = interaction_matrix(bgc, :assimilation_matrix)

        @test pal.parameter == :palatability_matrix
        @test assim.parameter == :assimilation_matrix
        @test pal.rows == [:Z_1, :Z_2]
        @test pal.columns == [:P_1, :P_2]
        @test assim.rows == pal.rows
        @test assim.columns == pal.columns
        @test pal.row_axis == :consumer
        @test pal.column_axis == :prey
        @test size(pal.matrix) == (length(pal.rows), length(pal.columns))
        @test size(assim.matrix) == (length(assim.rows), length(assim.columns))
        @test all(row in tracer_names(bgc) for row in pal.rows)
        @test all(col in tracer_names(bgc) for col in pal.columns)

        synthetic_bgc = (;
            parameters=(encounter_matrix=zeros(1, 1),),
            interaction_axes=(;
                parameters=(:encounter_matrix,), consumers=(:Z_1,), prey=(:P_1,)
            ),
        )
        encounter = interaction_matrix(synthetic_bgc, :encounter_matrix)
        @test encounter.parameter == :encounter_matrix
        @test encounter.matrix === synthetic_bgc.parameters.encounter_matrix

        @test_throws ArgumentError interaction_matrix(bgc, :consumer_global)
        @test_throws ArgumentError interaction_matrix(bgc, :unknown)
        try
            interaction_matrix(bgc, :unknown)
        catch err
            @test err isa ArgumentError
            @test occursin("palatability", sprint(showerror, err))
            @test occursin("assimilation", sprint(showerror, err))
        end


        phyto_diameters = [2.0, sqrt(20.0), 10.0]
        zoo_diameters = [20.0, 100.0]
        sized_bgc = Agate.Models.NiPiZD.construct(;
            size_structure=(;
                phytoplankton=(P=phyto_diameters,), zooplankton=(Z=zoo_diameters,)
            )
        )
        @test plankton_tracers(sized_bgc) == [:Z_1, :Z_2, :P_1, :P_2, :P_3]
        @test plankton_diameters(sized_bgc) ≈ [zoo_diameters; phyto_diameters]

        named_bgc = Agate.Models.NiPiZD.construct(;
            size_structure=(;
                phytoplankton=(diat=[2.0, 5.0, 10.0], dino=[8.0, 20.0]),
                zooplankton=(microzoo=[30.0, 60.0], mesozoo=[100.0]),
            ),
            grid=dummy_grid(Float32),
        )
        @test tracer_names(named_bgc) == [
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
        ]

        named_groups = plankton_groups(named_bgc)
        @test keys(named_groups) == (:microzoo, :mesozoo, :diat, :dino)
        @test named_groups.microzoo == [:microzoo_1, :microzoo_2]
        @test named_groups.mesozoo == [:mesozoo_1]
        @test named_groups.diat == [:diat_1, :diat_2, :diat_3]
        @test named_groups.dino == [:dino_1, :dino_2]
        @test plankton_tracers(named_bgc) == [
            :microzoo_1,
            :microzoo_2,
            :mesozoo_1,
            :diat_1,
            :diat_2,
            :diat_3,
            :dino_1,
            :dino_2,
        ]
        @test plankton_diameters(named_bgc) ≈
            [30.0, 60.0, 100.0, 2.0, 5.0, 10.0, 8.0, 20.0]

        named_pal = interaction_matrix(named_bgc, :palatability_matrix)
        @test named_pal.rows == [:microzoo_1, :microzoo_2, :mesozoo_1]
        @test named_pal.columns == [:diat_1, :diat_2, :diat_3, :dino_1, :dino_2]

        broken_bgc = (;
            parameters=merge(bgc.parameters, (palatability_matrix=ones(1, 1),)),
            interaction_axes=bgc.interaction_axes,
        )
        @test_throws ArgumentError interaction_matrix(broken_bgc, :palatability_matrix)
    end

    @testset "Generated model without groups" begin
        parameters = (rate=1.0,)
        tracers = (X=CompiledEquation((bgc, x, y, z, t, X) -> -bgc.parameters.rate * X),)
        model = define_tracer_functions(parameters, tracers)(parameters)
        groups = tracer_groups(model)

        @test tracer_names(model) == [:X]
        @test groups == (all=[:X], plankton=Symbol[], nonplankton=[:X], by_group=NamedTuple())
        @test plankton_groups(model) == NamedTuple()
        @test isempty(plankton_diameters(model))
        @test parameter_names(model) == [:rate]
        @test_throws ArgumentError interaction_matrix(model, :palatability_matrix)
    end
end
