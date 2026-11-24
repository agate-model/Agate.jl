using Agate
using Enzyme
using Agate.Introspection: parameter_names
using Test

import Agate.Parameters:
    DerivedDefault,
    ConstructionParameter,
    Parameter,
    derive_default

struct AddOneDefault end

derive_default(
    ::AddOneDefault, ::Agate.Processes.ModelDefinition, ::Any, params::NamedTuple
) = params.base + one(params.base)

function derived_default_definition()
    components = (D=Agate.Components.Pool(:nitrogen), N=Agate.Components.Pool(:nitrogen))
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
    @testset "Direct model" begin
        definitions = runtime_npz_definition().parameters
        bgc = runtime_npz_model(; grid=dummy_grid(Float32))

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
        @test_throws ArgumentError resolved_top((top=nothing,))
        @test !hasproperty(
            Agate.Construction.construct(definition).parameters, :base
        )
    end
end

# Use the ordinary setup path to prepare the fixed numerical realization.
function prepared_parameter_fixture(definition; scalar_type=Float64, parameter_overrides=(;))
    canonical = Agate.Processes.canonicalize_model(definition)
    layout = Agate.Construction._realize_process_definition(
        canonical, scalar_type; auxiliary_fields=Agate.Processes.driver_identities(canonical)
    )
    plan = Agate.Processes.build_parameter_plan(canonical, layout)
    return Agate.Construction._prepare_parameter_realization(
        canonical, layout, plan, parameter_overrides; derivation_owner=definition
    )
end

@testset "Prepared scalar trait realization" begin
    definition = derived_default_definition()
    for (T, overrides) in ((Float32, (;)), (Float64, (;)),
                           (Float64, (base=4.0,)), (Float32, (base=4.0, top=99.0)))
        bgc = Agate.Construction.construct(definition; scalar_type=T, parameter_overrides=overrides)
        prepared = prepared_parameter_fixture(definition; scalar_type=T, parameter_overrides=overrides)
        active = Agate.Runtime.active_parameters(bgc; top=true)
        selected = Agate.Construction._prepare_trait_realization(prepared, (:base,), active)
        baseline = (; base=prepared.values.base)
        updated = (; base=T(6))
        expected = Agate.Construction.construct(
            definition; scalar_type=T, parameter_overrides=merge(overrides, updated)
        ).parameters
        realize(updates) = Agate.Construction._realize_trait_parameters(selected, updates)
        @test realize(baseline) == bgc.parameters
        @test realize(updated) == expected
        @test realize(baseline) == bgc.parameters
        p = Agate.Construction._realize_trait_vector(T, selected.numerics, updated, selected.gather)
        @test p == T[expected.top] && eltype(p) === T
        @test keys(selected.numerics.derivations) == (hasproperty(overrides, :top) ? () : (:top,))
    end
end

# A weighted two-output map gives a cheap analytic oracle for the Enzyme VJP,
# without compiling another ODE adjoint. The direct NPZ fixture is off-optimum.
function weighted_prepared_traits(tau, numerics, gather, baseline)
    updates = (;
        specificity=[baseline[1], tau[1]],
        assimilation_efficiency=[baseline[2], tau[2]],
    )
    p = Agate.Construction._realize_trait_vector(Float64, numerics, updates, gather)
    return 2p[1] + 3p[2]
end

@testset "Prepared named SizeClass traits" begin
    definition = direct_npz_definition()
    initial = (;
        optimum_predator_prey_ratio=(Z_1=8.0,),
        specificity=(Z_1=0.4,),
        assimilation_efficiency=(Z_1=0.8,),
    )
    for T in (Float32, Float64), fixed in ((;), (palatability_matrix=fill(0.25, 1, 1),))
        overrides = merge(initial, fixed)
        bgc = Agate.Construction.construct(definition; scalar_type=T, parameter_overrides=overrides)
        prepared = prepared_parameter_fixture(definition; scalar_type=T, parameter_overrides=overrides)
        # Deliberately choose a different order from the canonical parameter names.
        active = Agate.Runtime.active_parameters(bgc;
            palatability_matrix=((:Z_1, :P_1),), assimilation_matrix=((:Z_1, :P_1),)
        )
        selected = Agate.Construction._prepare_trait_realization(
            prepared, (:specificity, :assimilation_efficiency), active
        )
        baseline = (; specificity=initial.specificity, assimilation_efficiency=initial.assimilation_efficiency)
        updated = (; specificity=(Z_1=0.6,), assimilation_efficiency=(Z_1=0.5,))
        expected = Agate.Construction.construct(
            definition; scalar_type=T, parameter_overrides=merge(overrides, updated)
        ).parameters
        realize(updates) = Agate.Construction._realize_trait_parameters(selected, updates)
        @test realize(baseline) == bgc.parameters
        actual = realize(updated)
        @test actual == expected
        @test all(eltype(getproperty(actual, name)) === T for name in (:palatability_matrix, :assimilation_matrix))
        @test Agate.Construction._gather_runtime_vector(T, actual, selected.gather) ==
              T[only(expected.palatability_matrix), only(expected.assimilation_matrix)]
        @test keys(selected.numerics.derivations) ==
              (isempty(fixed) ? (:assimilation_matrix, :palatability_matrix) : (:assimilation_matrix,))
        # The checked result owns its arrays; neither trait updates nor callers
        # modifying returned storage can change the prepared baseline.
        actual.palatability_matrix .= 0
        @test realize(baseline) == bgc.parameters
        @test_throws ArgumentError realize(merge(updated, (assimilation_efficiency=(Z_1=1.1,),)))
        @test_throws ArgumentError realize(merge(updated, (specificity=[T(0.3)],)))

        if T === Float64 && isempty(fixed)
            values = prepared.values
            classes = only(prepared.plan.parameters.specificity.storage_labels)
            @test classes == (:P_1, :Z_1)
            # Check that unrelated derived outputs are retained as baseline values.
            single = Agate.Construction._prepare_trait_realization(prepared, (:specificity,), active)
            @test keys(single.numerics.derivations) == (:palatability_matrix,)
            tau = [values.specificity[2], values.assimilation_efficiency[2]]
            inactive = (values.specificity[1], values.assimilation_efficiency[1])
            gradient = zeros(2)
            Enzyme.autodiff(
                Enzyme.set_runtime_activity(Enzyme.Reverse),
                Enzyme.Const(weighted_prepared_traits), Enzyme.Active,
                Enzyme.Duplicated(copy(tau), gradient), Enzyme.Const(selected.numerics),
                Enzyme.Const(selected.gather), Enzyme.Const(inactive),
            )
            @test gradient ≈ [-2log(5) * 5^(-tau[1]), 3.0]
        end
    end
end
