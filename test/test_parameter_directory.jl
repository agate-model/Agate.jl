using Agate
using Agate.Introspection: parameter_names
using Test

import Agate.Parameters:
    ConstantDefault,
    DerivedDefault,
    Parameter,
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
    base=Parameter(ConstantDefault(2.0)),
    top=Parameter(DerivedDefault(DoubleDefault(); deps=(:middle,))),
    middle=Parameter(DerivedDefault(AddOneDefault(); deps=(:base,))),
)

struct CyclicDerivedDefaultFixture end
parameter_definitions(::CyclicDerivedDefaultFixture) = (
    a=Parameter(DerivedDefault(AddOneDefault(); deps=(:b,))),
    b=Parameter(DerivedDefault(DoubleDefault(); deps=(:a,))),
)

@testset "Parameter directory" begin
    @testset "NiPiZD" begin
        family = Agate.Models.NiPiZD.NiPiZDFamily()
        dir = parameter_directory(family)
        @test !isempty(dir)

        bgc = Agate.Models.NiPiZD.construct(; grid=dummy_grid(Float32))
        dir_names = Set(keys(dir))

        # All constructed parameters should be declared in the directory.
        for k in parameter_names(bgc)
            @test k in dir_names
        end

        definitions = parameter_definitions(family)
        @test !hasproperty(dir.maximum_growth_rate, :name)
        @test definitions.linear_mortality.default isa Agate.Parameters.DiameterIndexedVectorDefault
        @test definitions.palatability_matrix.default isa DerivedDefault
        @test definitions.palatability_matrix.default.deps == (
            :optimum_predator_prey_ratio, :specificity, :protection
        )
        @test definitions.assimilation_matrix.default isa DerivedDefault
        @test definitions.assimilation_matrix.default.deps == (:assimilation_efficiency,)
        @test propertynames(dir.maximum_growth_rate) == (:axes,)
        @test dir.detritus_remineralization.axes === nothing
        @test dir.mortality_export_fraction.axes === nothing
        @test dir.maximum_growth_rate.axes == :plankton
        @test definitions.maximum_growth_rate.default isa Agate.Parameters.DiameterIndexedVectorDefault
        @test definitions.maximum_growth_rate.default.default == 0
        @test definitions.linear_mortality.default.default == 0
        @test definitions.mortality_export_fraction.default.value == 0.2
        @test dir.palatability_matrix.axes == (:consumer, :prey)
        @test dir.assimilation_matrix.axes == (:consumer, :prey)
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
