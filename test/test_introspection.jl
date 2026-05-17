using Agate
using Agate.Introspection:
    tracer_names,
    parameter_names,
    plankton_groups,
    plankton_tracers,
    nonplankton_tracers,
    tracer_groups,
    interaction_matrices,
    interaction_matrix,
    interaction_axes,
    interaction_table
using Test

@testset "Public introspection helpers" begin
    @testset "Model-constructed instance" begin
        bgc = Agate.Models.NiPiZD.construct(; grid=dummy_grid(Float32))

        @test tracer_names(bgc) == [:N, :D, :Z1, :Z2, :P1, :P2]

        groups = tracer_groups(bgc)
        @test groups.all == tracer_names(bgc)
        @test groups.by_group.Z == [:Z1, :Z2]
        @test groups.by_group.P == [:P1, :P2]
        @test groups.plankton == [:Z1, :Z2, :P1, :P2]
        @test groups.nonplankton == [:N, :D]
        @test plankton_groups(bgc) == groups.by_group
        @test plankton_tracers(bgc) == groups.plankton
        @test nonplankton_tracers(bgc) == groups.nonplankton

        pars = parameter_names(bgc)
        @test !isempty(pars)
        @test :data ∉ pars

        matrices = interaction_matrices(bgc)
        @test propertynames(matrices) == (:palatability, :assimilation)

        synthetic_interactions = (;
            palatability=ones(1, 1),
            encounter=zeros(1, 1),
            consumer_global=[1],
            prey_global=[1],
            global_to_consumer=[1],
            global_to_prey=[1],
        )
        synthetic_bgc = (; parameters=(; interactions=synthetic_interactions))
        synthetic_matrices = interaction_matrices(synthetic_bgc)
        @test propertynames(synthetic_matrices) == (:palatability, :encounter)
        @test synthetic_matrices.encounter === synthetic_interactions.encounter
        @test interaction_matrix(synthetic_bgc, :encounter) === synthetic_interactions.encounter

        axes = interaction_axes(bgc)
        @test axes.rows == [:Z1, :Z2]
        @test axes.columns == [:P1, :P2]
        @test axes.row_axis == :consumer
        @test axes.column_axis == :prey

        pal = interaction_table(bgc, :palatability)
        assim = interaction_table(bgc, :assimilation)

        @test pal.kind == :palatability
        @test assim.kind == :assimilation
        @test pal.matrix === interaction_matrix(bgc, :palatability)
        @test assim.matrix === interaction_matrix(bgc, :assimilation)
        @test pal.rows == axes.rows
        @test pal.columns == axes.columns
        @test assim.rows == axes.rows
        @test assim.columns == axes.columns
        @test pal.row_axis == :consumer
        @test pal.column_axis == :prey
        @test size(pal.matrix) == (length(pal.rows), length(pal.columns))
        @test size(assim.matrix) == (length(assim.rows), length(assim.columns))
        @test all(row in tracer_names(bgc) for row in pal.rows)
        @test all(col in tracer_names(bgc) for col in pal.columns)
        @test_throws ArgumentError interaction_matrix(bgc, :unknown)
    end

    @testset "Generated model (define_tracer_functions)" begin
        include(joinpath(@__DIR__, "NPZD", "tracers.jl"))

        model = AgateNPZD(parameters)
        @test tracer_names(model) == [:N, :D, :P, :Z]
        @test plankton_groups(model) == NamedTuple()
        @test plankton_tracers(model) == Symbol[]
        @test nonplankton_tracers(model) == tracer_names(model)

        groups = tracer_groups(model)
        @test groups.all == tracer_names(model)
        @test groups.by_group == NamedTuple()
        @test groups.plankton == Symbol[]
        @test groups.nonplankton == tracer_names(model)
        @test !isempty(parameter_names(model))
        @test_throws ArgumentError interaction_matrices(model)
        @test_throws ArgumentError interaction_axes(model)
        @test_throws ArgumentError interaction_table(model, :palatability)
    end
end
