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
        @test length(plankton_diameters(bgc)) == length(groups.plankton)
        @test plankton_diameters(bgc) == collect(bgc.plankton_diameters)
        @test eltype(plankton_diameters(bgc)) === Float32
        @test nonplankton_tracers(bgc) == groups.nonplankton

        pars = parameter_names(bgc)
        @test !isempty(pars)
        @test :data ∉ pars

        pal = interaction_matrix(bgc, :palatability)
        assim = interaction_matrix(bgc, :assimilation)

        @test pal.kind == :palatability
        @test assim.kind == :assimilation
        @test pal.rows == [:Z1, :Z2]
        @test pal.columns == [:P1, :P2]
        @test assim.rows == pal.rows
        @test assim.columns == pal.columns
        @test pal.row_axis == :consumer
        @test pal.column_axis == :prey
        @test size(pal.matrix) == (length(pal.rows), length(pal.columns))
        @test size(assim.matrix) == (length(assim.rows), length(assim.columns))
        @test all(row in tracer_names(bgc) for row in pal.rows)
        @test all(col in tracer_names(bgc) for col in pal.columns)

        synthetic_bgc = (;
            parameters=(;
                interactions=(; encounter=zeros(1, 1), consumer_global=[1], prey_global=[1]),
            ),
            tracers=bgc.tracers,
        )
        encounter = interaction_matrix(synthetic_bgc, :encounter)
        @test encounter.kind == :encounter
        @test encounter.matrix === synthetic_bgc.parameters.interactions.encounter

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
            phyto_size_structure=phyto_diameters, zoo_size_structure=zoo_diameters
        )
        @test plankton_tracers(sized_bgc) == [:Z1, :Z2, :P1, :P2, :P3]
        @test plankton_diameters(sized_bgc) ≈ [zoo_diameters; phyto_diameters]

        darwin_bgc = Agate.Models.DARWIN.construct(; grid=dummy_grid(Float32))
        @test length(plankton_diameters(darwin_bgc)) == length(plankton_tracers(darwin_bgc))
        @test plankton_diameters(darwin_bgc) == collect(darwin_bgc.plankton_diameters)

        params = bgc.parameters
        interactions = params.interactions
        interaction_names = propertynames(interactions)
        broken_interactions = NamedTuple{interaction_names}(
            map(interaction_names) do name
                name === :palatability ? ones(1, 1) : getproperty(interactions, name)
            end,
        )
        broken_bgc = (; parameters=(; interactions=broken_interactions), tracers=bgc.tracers)
        @test_throws ArgumentError interaction_matrix(broken_bgc, :palatability)
    end

    @testset "Generated model (define_tracer_functions)" begin
        include(joinpath(@__DIR__, "NPZD", "tracers.jl"))

        model = AgateNPZD(parameters)
        @test tracer_names(model) == [:N, :D, :P, :Z]
        @test plankton_groups(model) == NamedTuple()
        @test plankton_tracers(model) == Symbol[]
        @test isempty(plankton_diameters(model))
        @test nonplankton_tracers(model) == tracer_names(model)

        groups = tracer_groups(model)
        @test groups.all == tracer_names(model)
        @test groups.by_group == NamedTuple()
        @test groups.plankton == Symbol[]
        @test groups.nonplankton == tracer_names(model)
        @test !isempty(parameter_names(model))
        @test_throws ArgumentError interaction_matrix(model, :palatability)
    end
end
