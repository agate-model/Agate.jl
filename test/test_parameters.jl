using Agate
using Agate.Introspection: parameter_names
using Test

import Agate.Parameters:
    DerivedDefault,
    MetaParameter,
    Parameter,
    derive_default,
    parameter_definitions

struct DerivedDefaultFixture end
struct AddOneDefault end

derive_default(::AddOneDefault, ::DerivedDefaultFixture, ::Any, params::NamedTuple) =
    params.base + one(params.base)

parameter_definitions(::DerivedDefaultFixture) = (
    base=MetaParameter(2.0),
    top=Parameter(DerivedDefault(AddOneDefault(); deps=(:base,))),
)

@testset "Parameter definitions and planning" begin
    @testset "NiPiZD" begin
        family = Agate.Models.NiPiZD.NiPiZDFamily()
        definitions = parameter_definitions(family)
        bgc = Agate.Models.NiPiZD.construct(; grid=dummy_grid(Float32))

        @test all(name -> hasproperty(definitions, name), parameter_names(bgc))
        @test definitions.maximum_growth_rate isa Parameter
        @test definitions.linear_mortality.default isa Agate.Parameters.DiameterIndexedVectorDefault
        @test definitions.palatability_matrix.default isa DerivedDefault
        @test definitions.palatability_matrix.default.deps == (
            :optimum_predator_prey_ratio, :specificity, :protection
        )
        @test definitions.assimilation_matrix.default isa DerivedDefault
        @test definitions.assimilation_matrix.default.deps == (:assimilation_efficiency,)
        @test definitions.maximum_growth_rate.default.default == 0
        @test definitions.linear_mortality.default.default == 0
        @test definitions.mortality_export_fraction.default.value == 0.2
        @test definitions.specificity isa MetaParameter
        @test definitions.specificity.axes === :plankton
        @test definitions.assimilation_efficiency isa MetaParameter
    end

    @testset "Realized ParameterPlan" begin
        family = Agate.Models.NiPiZD.NiPiZDFamily()
        normalized = Agate.Processes.normalize_model(Agate.Processes.ModelDefinition(family))
        layout = default_nipizd_layout()
        plan = Agate.Processes.build_parameter_plan(normalized, layout)

        growth = plan.parameters.maximum_growth_rate
        @test (growth.rank, growth.axes, growth.storage_shape, growth.storage_labels) ==
            (1, (:population,), (2,), ((:P_1, :P_2),))
        @test growth.storage_diameters == layout.diameters[3:4]
        @test growth.runtime_bound

        palatability = plan.parameters.palatability_matrix
        @test (palatability.rank, palatability.storage_shape) == (2, (2, 2))
        @test palatability.storage_labels == ((:Z_1, :Z_2), (:P_1, :P_2))
        @test palatability.runtime_bound

        specificity = plan.parameters.specificity
        @test specificity.axes == (:plankton,)
        @test Agate.Construction.evaluate_process_default(
            Agate.Parameters.ConstantDefault(2), specificity, Float32
        ) == fill(2f0, 4)
        @test !specificity.runtime_bound
        @test :specificity ∉ plan.runtime_names
        @test :palatability_matrix ∈ plan.runtime_names

        remineralization = plan.parameters.detritus_remineralization
        @test (remineralization.rank, remineralization.storage_shape, remineralization.storage_labels) ==
            (0, (), ())

        mortality = plan.parameters.linear_mortality
        @test mortality.storage_labels == ((:Z_1, :Z_2, :P_1, :P_2),)
        mortality_slots = Tuple(
            Agate.Processes.planned_parameter_slot(plan, binding)
            for binding in normalized.parameter_bindings if binding.parameter === :linear_mortality
        )
        @test length(mortality_slots) == 2
        @test Set(mortality_slots) == Set((((1, 2),), ((3, 4),)))

        metadata = Agate.Processes.parameter_plan_metadata(plan)
        @test metadata.maximum_growth_rate.labels == growth.storage_labels
        @test metadata.palatability_matrix.labels == palatability.storage_labels
        @test metadata.specificity.derived_runtime_parameters == (:palatability_matrix,)

        growth_binding = only(filter(normalized.parameter_bindings) do binding
            binding.parameter === :maximum_growth_rate
        end)
        @test Agate.Processes.planned_parameter_slot(plan, growth_binding) == ((1, 2),)
    end

    @testset "Derived default dependency resolution" begin
        source = DerivedDefaultFixture()
        components = (D=Agate.Configuration.Pool(:nitrogen),)
        normalized = Agate.Processes.normalize_model(Agate.Processes.ModelDefinition(;
            components,
            processes=(remineralization=Agate.Processes.Remineralization(
                Agate.Processes.LinearRemineralization();
                sources=:D,
                destination=:D,
                bindings=(rate=(D=:top,),),
            ),),
            parameters=parameter_definitions(source),
        ))
        layout = Agate.Configuration.realize_model_layout(components; scalar_type=Float64)

        plan = Agate.Processes.build_parameter_plan(normalized, layout)
        defaults = Agate.Construction.build_process_parameter_defaults(plan, Float64)
        @test defaults == (base=2.0,)

        resolve(overrides=(;)) = Agate.Construction.resolve_parameter_defaults(
            plan, layout, merge(defaults, overrides); derivation_owner=source
        )

        @test resolve() == (base=2.0, top=3.0)
        @test resolve((base=4.0,)) == (base=4.0, top=5.0)
        @test resolve((top=99.0,)) == (base=2.0, top=99.0)
        @test resolve((base=4.0, top=99.0)) == (base=4.0, top=99.0)
    end
end
