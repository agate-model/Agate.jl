using Agate
using Agate.Introspection: parameter_names
using Test

import Agate.Parameters:
    DerivedDefault,
    ConstructionParameter,
    Parameter,
    derive_default,
    parameter_definitions

struct AddOneDefault end

derive_default(
    ::AddOneDefault, ::Agate.Processes.ModelDefinition, ::Any, params::NamedTuple
) = params.base + one(params.base)

function derived_default_definition()
    components = (D=Agate.Configuration.Pool(:nitrogen), N=Agate.Configuration.Pool(:nitrogen))
    processes = (remineralization=Agate.Processes.Remineralization(
        Agate.Processes.LinearRemineralization();
        sources=:D,
        destination=:N,
        bindings=(rate=(D=:top,),),
    ),)
    parameters = (
        base=ConstructionParameter(2.0),
        top=Parameter(DerivedDefault(AddOneDefault(); deps=(:base,))),
    )
    return Agate.Processes.ModelDefinition(; components, processes, parameters)
end

@testset "Parameter definitions and derived defaults" begin
    @testset "NiPiZD" begin
        family = Agate.Models.NiPiZD.NiPiZDFamily()
        definitions = parameter_definitions(family)
        bgc = Agate.Models.NiPiZD.construct(; grid=dummy_grid(Float32))

        @test all(name -> hasproperty(definitions, name), parameter_names(bgc))
    end

    @testset "Derived default dependency resolution" begin
        definition = derived_default_definition()
        resolved_top(overrides=(;)) = Agate.Construction.construct(
            definition; parameter_overrides=overrides
        ).parameters.top

        @test resolved_top() == 3.0
        @test resolved_top((base=4.0,)) == 5.0
        @test resolved_top((top=99.0,)) == 99.0
        @test resolved_top((base=4.0, top=99.0)) == 99.0
        @test !hasproperty(
            Agate.Construction.construct(definition).parameters, :base
        )
    end
end
