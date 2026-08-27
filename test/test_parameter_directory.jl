using Agate
using Agate.Introspection: parameter_names
using Test

import Agate.Parameters:
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
    base=Parameter(2.0),
    top=Parameter(DerivedDefault(DoubleDefault(); deps=(:middle,))),
    middle=Parameter(DerivedDefault(AddOneDefault(); deps=(:base,))),
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

    @testset "Realized ParameterPlan" begin
        family = Agate.Models.NiPiZD.NiPiZDFamily()
        normalized = Agate.Processes.normalize_model(Agate.Processes.ModelDefinition(family))
        layout = default_nipizd_layout()
        plan = Agate.Processes.build_parameter_plan(normalized, layout)

        growth = plan.parameters.maximum_growth_rate
        @test (growth.rank, growth.storage_shape, growth.storage_labels) ==
            (1, (4,), ((:Z_1, :Z_2, :P_1, :P_2),))
        @test growth.applicable_indices == ((3, 4),)
        @test growth.storage_diameters == layout.diameters
        @test growth.runtime_bound

        palatability = plan.parameters.palatability_matrix
        @test (palatability.rank, palatability.storage_shape) == (2, (2, 2))
        @test palatability.storage_labels == ((:Z_1, :Z_2), (:P_1, :P_2))
        @test palatability.runtime_bound

        specificity = plan.parameters.specificity
        @test specificity.applicable_indices == ((1, 2, 3, 4),)
        @test Agate.Construction.evaluate_process_default(
            Agate.Parameters.ConstantDefault(2), specificity, Float32
        ) == fill(2f0, 4)
        @test !specificity.runtime_bound
        @test :specificity ∉ plan.runtime_names
        @test :palatability_matrix ∈ plan.runtime_names

        remineralization = plan.parameters.detritus_remineralization
        @test (remineralization.rank, remineralization.storage_shape, remineralization.storage_labels) ==
            (0, (), ())

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
        growth_slot = Agate.Processes.planned_parameter_slot(plan, growth_binding)
        @test growth_slot == ((3, 4),)
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
        @test plan.derived_order == (:middle, :top)
        defaults = Agate.Construction.build_process_parameter_defaults(plan, Float64)
        @test defaults == (base=2.0,)

        resolve(overrides=(;)) = Agate.Construction.resolve_parameter_defaults(
            plan, layout, merge(defaults, overrides), Tuple(keys(overrides));
            derivation_owner=source,
        )

        resolved = resolve()
        @test (resolved.base, resolved.middle, resolved.top) == (2.0, 3.0, 6.0)

        resolved_base = resolve((base=4.0,))
        @test (resolved_base.middle, resolved_base.top) == (5.0, 10.0)

        resolved_middle = resolve((middle=9.0,))
        @test (resolved_middle.middle, resolved_middle.top) == (9.0, 18.0)

        resolved_top = resolve((base=4.0, top=99.0))
        @test (resolved_top.middle, resolved_top.top) == (5.0, 99.0)
    end
end
