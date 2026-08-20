using Agate
using Agate.Introspection: parameter_names
using Test

import Agate.Parameters:
    ConstDefault,
    DerivedDefault,
    FillDefault,
    ParameterDefinition,
    derive_default,
    parameter_definitions,
    parameter_directory

struct DerivedDefaultFixture end
struct AddOneDefault end
struct DoubleDefault end

derive_default(::AddOneDefault, ::DerivedDefaultFixture, ::Any, params::NamedTuple) =
    params.base + one(params.base)
derive_default(::DoubleDefault, ::DerivedDefaultFixture, ::Any, params::NamedTuple) =
    2 * params.middle

parameter_definitions(::DerivedDefaultFixture) = (
    ParameterDefinition(:base, ConstDefault(2.0)),
    ParameterDefinition(:top, DerivedDefault(DoubleDefault(); deps=(:middle,))),
    ParameterDefinition(:middle, DerivedDefault(AddOneDefault(); deps=(:base,))),
)

struct CyclicDerivedDefaultFixture end
parameter_definitions(::CyclicDerivedDefaultFixture) = (
    ParameterDefinition(:a, DerivedDefault(AddOneDefault(); deps=(:b,))),
    ParameterDefinition(:b, DerivedDefault(DoubleDefault(); deps=(:a,))),
)

@testset "Parameter directory" begin
    @testset "NiPiZD" begin
        family = Agate.Models.NiPiZD.NiPiZDFamily()
        dir = parameter_directory(family)
        @test !isempty(dir)

        bgc = Agate.Models.NiPiZD.construct(; grid=dummy_grid(Float32))
        dir_names = Set(spec.name for spec in dir)

        # All constructed parameters should be declared in the directory.
        for k in parameter_names(bgc)
            @test k in dir_names
        end

        specmap = Dict(spec.name => spec for spec in dir)
        definitions = Dict(def.spec.name => def for def in parameter_definitions(family))
        @test definitions[:linear_mortality].default isa FillDefault
        @test definitions[:palatability_matrix].default isa DerivedDefault
        @test definitions[:palatability_matrix].default.deps == (
            :optimum_predator_prey_ratio, :specificity, :protection
        )
        @test definitions[:assimilation_matrix].default isa DerivedDefault
        @test definitions[:assimilation_matrix].default.deps == (:assimilation_efficiency,)
        @test all(!isempty(spec.provides) for spec in values(specmap))
        @test length(specmap[:linear_mortality].provides) == 2
        @test length(specmap[:mortality_export_fraction].provides) == 3
        @test specmap[:detritus_remineralization].shape == :scalar
        @test specmap[:maximum_growth_rate].shape == :vector
        @test specmap[:maximum_growth_rate].axes == :plankton
        @test !isnothing(specmap[:maximum_growth_rate].materialization)
        @test specmap[:maximum_growth_rate].materialization.fill_value == 0
        @test !isnothing(specmap[:linear_mortality].materialization)
        @test specmap[:linear_mortality].materialization.fill_value == 0
        @test specmap[:palatability_matrix].shape == :matrix
        @test specmap[:palatability_matrix].axes == (:consumer, :prey)
        @test specmap[:assimilation_matrix].axes == (:consumer, :prey)
        @test (
            specmap[:maximum_growth_rate].runtime_path,
            specmap[:palatability_matrix].runtime_path,
            specmap[:assimilation_matrix].runtime_path,
        ) == (
            (:maximum_growth_rate,),
            (:interactions, :palatability),
            (:interactions, :assimilation),
        )
    end

    @testset "Derived default dependency resolution" begin
        source = DerivedDefaultFixture()
        context = (; scalar_type=Float64)
        @test Agate.Construction.validate_parameter_directory(source) == (:base, :top, :middle)

        defaults = Agate.Construction.build_process_parameter_defaults(
            source, nothing, nothing, context, Float64
        )
        @test defaults == (base=2.0,)

        resolve(overrides=(;)) = Agate.Construction.resolve_parameter_defaults(
            source, context, merge(defaults, overrides), Tuple(keys(overrides))
        )

        resolved = resolve()
        @test (resolved.base, resolved.middle, resolved.top) == (2.0, 3.0, 6.0)

        resolved_base = resolve((base=4.0,))
        @test (resolved_base.middle, resolved_base.top) == (5.0, 10.0)

        resolved_middle = resolve((middle=9.0,))
        @test (resolved_middle.middle, resolved_middle.top) == (9.0, 18.0)

        resolved_top = resolve((base=4.0, top=99.0))
        @test (resolved_top.middle, resolved_top.top) == (5.0, 99.0)

        @test_throws ArgumentError Agate.Construction.resolve_parameter_defaults(
            CyclicDerivedDefaultFixture(), context, (;), ()
        )
    end
end
