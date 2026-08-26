using Agate
using Test

using Agate.Configuration:
    DiameterRangeSpecification, Population, Pool, population_state, realize_components,
    realize_model_layout, normalize_diameters
using Agate.Runtime: class, resolve_class, resolve_state, class_count, build_tracer_index, Tracers

@testset "Independent interaction roles" begin
    components = (
        P=Population(:carbon; size_structure=[1.0]),
        B=Population(:carbon; size_structure=[1.0]),
        M=Population(:carbon; size_structure=[1.0]),
        Z=Population(:carbon; size_structure=[1.0]),
    )
    layout = realize_model_layout(
        components; interaction_roles=(consumers=(:Z, :M), prey=(:P, :B, :M))
    )

    @test layout.consumer_indices == (3, 4)
    @test layout.prey_indices == (1, 2, 3)
end

@testset "ClassRef + Tracers accessors" begin
    layout = default_nipizd_layout(Float64; auxiliary_fields=(:PAR,))

    @test class_count(layout, :Z) == 2
    @test class_count(layout, :P) == 2
    @test resolve_class(layout, class(:Z, 1)) == 3
    @test resolve_class(layout, class(:Z, 2)) == 4
    @test resolve_class(layout, class(:P, 1)) == 5
    @test resolve_class(layout, class(:P, 2)) == 6

    tracers = Tracers(build_tracer_index(layout))
    args = (10.0, 20.0, 1.0, 2.0, 3.0, 4.0, 42.0)

    @test tracers.N(args) == 10.0
    @test tracers.D(args) == 20.0
    @test Tuple(tracers.plankton(args, i) for i in 1:4) == (1.0, 2.0, 3.0, 4.0)
    @test tracers.Z(args, 2) == 2.0
    @test tracers.P(args, 1) == 3.0
    @test tracers.PAR(args) == 42.0
    @test @inferred(tracers.P(args, 2)) == 4.0
end

@testset "Generic component ClassRef" begin
    layout = realize_components((
        B=Population(:carbon),
        POM=Pool(:carbon; size_structure=[0.5, 5.0, 50.0]),
    ))

    @test class_count(layout, :B) == 1
    @test class_count(layout, :POM) == 3
    @test resolve_class(layout, class(:B, 1)) == 4
    @test resolve_class(layout, class(:POM, 1)) == 1
    @test resolve_class(layout, class(:POM, 3)) == 3
    @test_throws ArgumentError resolve_class(layout, class(:POM, 4))

    multistate = realize_components((
        P=Population(; states=(carbon=:carbon, nitrogen=:nitrogen), size_structure=[1.0, 2.0]),
    ))
    @test class_count(multistate, :P) == 2
    @test_throws ArgumentError resolve_class(multistate, class(:P, 1))
    @test resolve_state(multistate, class(:P, 1), population_state(:P, :carbon)) == 1
    @test resolve_state(multistate, class(:P, 2), population_state(:P, :nitrogen)) == 4
end

@testset "Diameter input normalization" begin
    components = (
        Z=Population(:carbon),
        P=Population(:carbon),
    )
    layout = realize_model_layout(
        components;
        population_groups=(Z=(:Z,), P=(:P,)),
        group_diameters=(
            Z=DiameterRangeSpecification(2, 20.0, 100.0, :linear_splitting),
            P=(n=3, min_esd=2.0, max_esd=10.0, splitting=:log_splitting),
        ),
    )

    @test class_count(layout, :Z) == 2
    @test class_count(layout, :P) == 3
    @test layout.diameters[1:2] == (20.0, 100.0)
    @test collect(layout.diameters[3:5]) ≈ [2.0, sqrt(20.0), 10.0]

    bad_path = "population group :P diameters"
    for invalid in (
        Float64[], [1.0, 0.0], [1.0, Inf], [true],
        (n=0, min_esd=1.0, max_esd=2.0, splitting=:log_splitting),
        (n=true, min_esd=1.0, max_esd=2.0, splitting=:log_splitting),
        (n=2, min_esd=0.0, max_esd=2.0, splitting=:log_splitting),
        (n=2, min_esd=2.0, max_esd=1.0, splitting=:log_splitting),
        (n=2, min_esd=1.0, max_esd=2.0, splitting=:unsupported),
    )
        err = try
            realize_model_layout(
                (P=Population(:carbon),);
                population_groups=(P=(:P,),), group_diameters=(P=invalid,),
            )
            nothing
        catch caught
            caught
        end
        @test err isa ArgumentError
        @test occursin(bad_path, sprint(showerror, err))
    end

    @test normalize_diameters([1.0, 2.0]).n == 2
end
