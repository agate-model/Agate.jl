using Agate
using Agate.Introspection:
    tracer_names,
    parameter_names,
    parameter_domains,
    pfts,
    plankton_tracers,
    plankton_diameters,
    nonplankton_tracers,
    tracer_groups,
    interaction_matrix,
    model_summary
using Test


@testset "Public introspection helpers" begin
    @testset "Model-constructed instance" begin
        bgc = Agate.Models.NiPiZD.construct(; grid=dummy_grid(Float32))

        @test tracer_names(bgc) == [:N, :D, :P_1, :P_2, :Z_1, :Z_2]

        groups = tracer_groups(bgc)
        @test groups.all == tracer_names(bgc)
        @test groups.entities_by_pft.Z == [:Z_1, :Z_2]
        @test groups.entities_by_pft.P == [:P_1, :P_2]
        @test groups.plankton == [:P_1, :P_2, :Z_1, :Z_2]
        @test groups.nonplankton == [:N, :D]
        @test pfts(bgc) == groups.entities_by_pft
        @test plankton_tracers(bgc) == groups.plankton
        @test length(plankton_diameters(bgc)) == length(groups.plankton)
        @test eltype(plankton_diameters(bgc)) === Float32
        @test nonplankton_tracers(bgc) == groups.nonplankton
        @test !model_summary(bgc).has_sinking_velocities
        sinking_bgc = Agate.Models.NiPiZD.construct(; sinking_tracers=(D=1.0,))
        @test model_summary(sinking_bgc).has_sinking_velocities

        pars = parameter_names(bgc)
        @test !isempty(pars)
        @test :data ∉ pars
        @test :palatability_matrix in pars
        @test :assimilation_matrix in pars
        @test :interactions ∉ pars
        @test parameter_domains(bgc, :assimilation_matrix) == (:unit_interval,)

        pal = interaction_matrix(bgc, :palatability_matrix)
        assim = interaction_matrix(bgc, :assimilation_matrix)

        @test pal.parameter == :palatability_matrix
        @test assim.parameter == :assimilation_matrix
        @test pal.rows == [:Z_1, :Z_2]
        @test pal.columns == [:P_1, :P_2]
        @test assim.rows == pal.rows
        @test assim.columns == pal.columns
        @test pal.row_axis == :consumer
        @test pal.column_axis == :resource
        @test size(pal.matrix) == (length(pal.rows), length(pal.columns))
        @test size(assim.matrix) == (length(assim.rows), length(assim.columns))
        @test all(row in tracer_names(bgc) for row in pal.rows)
        @test all(col in tracer_names(bgc) for col in pal.columns)

        @test_throws ArgumentError interaction_matrix(bgc, :consumer_global)
        message = argument_error_message(() -> interaction_matrix(bgc, :unknown))
        @test occursin("palatability", message)
        @test occursin("assimilation", message)

        phyto_diameters = [2.0, sqrt(20.0), 10.0]
        zoo_diameters = [20.0, 100.0]
        sized_bgc = Agate.Models.NiPiZD.construct(;
            size_structure=(;
                phytoplankton=(P=phyto_diameters,), zooplankton=(Z=zoo_diameters,)
            )
        )
        @test plankton_tracers(sized_bgc) == [:P_1, :P_2, :P_3, :Z_1, :Z_2]
        @test plankton_diameters(sized_bgc) ≈ [phyto_diameters; zoo_diameters]

        named_bgc = Agate.Models.NiPiZD.construct(;
            size_structure=nipizd_named_size_structure(), grid=dummy_grid(Float32)
        )

        named_pfts = pfts(named_bgc)
        @test keys(named_pfts) == (:diat, :dino, :mesozoo, :microzoo)
        @test named_pfts.microzoo == [:microzoo_1, :microzoo_2]
        @test named_pfts.mesozoo == [:mesozoo_1]
        @test named_pfts.diat == [:diat_1, :diat_2, :diat_3]
        @test named_pfts.dino == [:dino_1, :dino_2]
        @test plankton_diameters(named_bgc) ≈
            [2.0, 5.0, 10.0, 8.0, 20.0, 100.0, 30.0, 60.0]

        named_pal = interaction_matrix(named_bgc, :palatability_matrix)
        @test named_pal.rows == [:mesozoo_1, :microzoo_1, :microzoo_2]
        @test named_pal.columns == [:diat_1, :diat_2, :diat_3, :dino_1, :dino_2]

    end

    @testset "Authored model without plankton PFTs" begin
        components = (X=Agate.Components.Pool(:carbon), Y=Agate.Components.Pool(:carbon))
        processes = (transfer=Agate.Processes.Remineralization(
            Agate.Processes.LinearRemineralization();
            sources=:X, destination=:Y, bindings=(rate=(X=:rate,),),
        ),)
        parameters = (rate=Agate.Parameters.Parameter(Agate.Parameters.ConstantDefault(1.0)),)
        model = Agate.Construction.construct(Agate.Processes.ModelDefinition(; components, processes, parameters))
        groups = tracer_groups(model)

        @test tracer_names(model) == [:X, :Y]
        @test (groups.all, groups.plankton, groups.nonplankton, groups.entities_by_pft) ==
              ([:X, :Y], Symbol[], [:X, :Y], NamedTuple())
        @test pfts(model) == NamedTuple()
        @test isempty(plankton_diameters(model))
        @test parameter_names(model) == [:rate]
        @test_throws ArgumentError interaction_matrix(model, :palatability_matrix)
    end
end
