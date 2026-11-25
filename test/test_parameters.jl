using Agate
using Enzyme
using Agate.Introspection: parameter_names
using Agate.Construction: active_traits, realize_trait_parameters
using Test
using Oceananigans.Architectures: GPU

import Agate.Parameters: DerivedDefault, ConstructionParameter, Parameter, derive_default

struct AddOneDefault end

function derive_default(
    ::AddOneDefault, ::Agate.Processes.ModelDefinition, ::Any, params::NamedTuple
)
    return params.base + one(params.base)
end

struct AddTraitsDefault end

function derive_default(
    ::AddTraitsDefault, ::Agate.Processes.ModelDefinition, ::Any, params::NamedTuple
)
    return params.base + params.shift
end

function derived_default_definition()
    components = (
        D=Agate.Components.Pool(:nitrogen),
        E=Agate.Components.Pool(:nitrogen),
        N=Agate.Components.Pool(:nitrogen),
    )
    processes = (
        remineralization=Agate.Processes.Remineralization(
            Agate.Processes.LinearRemineralization();
            sources=:D,
            destination=:N,
            bindings=(rate=(D=:top,),),
        ),
        other_remineralization=Agate.Processes.Remineralization(
            Agate.Processes.LinearRemineralization();
            sources=:E,
            destination=:N,
            bindings=(rate=(E=:other,),),
        ),
    )
    parameters = (
        base=ConstructionParameter(2.0),
        shift=ConstructionParameter(1.0),
        top=Parameter(DerivedDefault(AddOneDefault(); deps=(:base,))),
        other=Parameter(DerivedDefault(AddTraitsDefault(); deps=(:base, :shift))),
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
        resolved_top(overrides=(;)) =
            Agate.Construction.construct(definition; parameter_overrides=overrides).parameters.top

        @test resolved_top() == 3.0
        @test resolved_top((base=4.0,)) == 5.0
        @test resolved_top((top=99.0,)) == 99.0
        @test resolved_top((base=4.0, top=99.0)) == 99.0
        @test_throws ArgumentError resolved_top((top=nothing,))
        @test !hasproperty(Agate.Construction.construct(definition).parameters, :base)
    end
end

# Differentiate a weighted runtime vector without compiling an ODE adjoint.
function weighted_trait_parameters(tau, parameter_map, weights)
    p = parameter_map(tau)
    total = zero(eltype(tau))
    for i in eachindex(p, weights)
        total += p[i] * weights[i]
    end
    return total
end

function trait_map_gradient(tau, parameter_map, weights)
    gradient = zeros(eltype(tau), length(tau))
    Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.Reverse),
        Enzyme.Const(weighted_trait_parameters),
        Enzyme.Active,
        Enzyme.Duplicated(copy(tau), gradient),
        Enzyme.Const(parameter_map),
        Enzyme.Const(weights),
    )
    return gradient
end

# One trait feeds two outputs; both selected traits feed :other. Reuse this
# fixture to check input override semantics, output overrides and ordering.
@testset "Active scalar construction traits" begin
    definition = derived_default_definition()
    for (T, overrides) in
        ((Float32, (;)), (Float64, (base=4.0,)), (Float32, (base=4.0, top=99.0)))
        traits = active_traits(
            definition;
            scalar_type=T,
            parameter_overrides=overrides,
            trait_selection=(shift=true, base=true),
        )
        active = traits.active_parameters
        targets = hasproperty(overrides, :top) ? (:other,) : (:other, :top)
        @test traits.labels == ("shift", "base")
        @test active.labels == string.(targets)
        @test traits.explicit_override_keys == keys(overrides)
        tau0 = copy(traits.values)
        @test tau0 == T[1, get(overrides, :base, 2.0)]
        tau = T[2, 6]
        expected = Agate.Construction.construct(
            definition;
            scalar_type=T,
            parameter_overrides=merge(overrides, (base=tau[2], shift=tau[1])),
        )
        p = realize_trait_parameters(traits, tau)
        @test p == T[getproperty(expected.parameters, name) for name in targets]
        @test traits.parameter_map(tau) == p
        @test realize_trait_parameters(traits, tau0) == active.values
        @test traits.bgc.parameters == Agate.Construction.construct(
            definition; scalar_type=T, parameter_overrides=overrides
        ).parameters
        @test eltype(p) === T
        if T === Float64
            # tau = (shift, base); output order = (other, top).
            @test trait_map_gradient(tau, traits.parameter_map, [3.0, 2.0]) ≈ [3.0, 5.0]
            for invalid in (tau[1:1], [tau; 0.0])
                @test_throws DimensionMismatch traits.parameter_map(invalid)
            end
        end
    end

    # The architecture guard runs before allocating any device storage.
    arch = GPU(nothing)
    for options in ((; arch), (; grid=dummy_grid(Float64; arch)))
        @test_throws ArgumentError active_traits(
            definition; trait_selection=(base=true,), options...
        )
    end
    for selection in ((;), (missing=true,), (top=true,), (base=false,), (base=(:D,),))
        @test_throws ArgumentError active_traits(definition; trait_selection=selection)
    end
    @test_throws ArgumentError active_traits(
        definition;
        trait_selection=(base=true, shift=true),
        parameter_overrides=(other=1.0,),
    ) # :shift has no live output
    @test_throws ArgumentError active_traits(
        definition; trait_selection=(base=true,), parameter_overrides=(top=1.0, other=1.0)
    )
end

@testset "Active named SizeClass traits" begin
    definition = direct_npz_definition()
    initial = (;
        optimum_predator_prey_ratio=(Z_1=8.0,),
        specificity=(Z_1=0.4,),
        assimilation_efficiency=(Z_1=0.8,),
    )
    for T in (Float32, Float64), fixed in ((;), (palatability_matrix=fill(0.25, 1, 1),))
        overrides = merge(initial, fixed)
        selection = if isempty(fixed)
            (specificity=(:Z_1,), assimilation_efficiency=(:Z_1,))
        else
            (assimilation_efficiency=(:Z_1,),)
        end
        traits = active_traits(
            definition;
            scalar_type=T,
            parameter_overrides=overrides,
            trait_selection=selection,
        )
        active = traits.active_parameters
        baseline = deepcopy(traits.bgc.parameters)
        tau0 = copy(traits.values)
        @test tau0 == (isempty(fixed) ? T[0.4, 0.8] : T[0.8])
        targets = if isempty(fixed)
            (:assimilation_matrix, :palatability_matrix)
        else
            (:assimilation_matrix,)
        end
        @test active.labels == map(name -> "$(name)[Z_1, P_1]", targets)
        updates = if isempty(fixed)
            (specificity=(Z_1=0.6,), assimilation_efficiency=(Z_1=0.5,))
        else
            (assimilation_efficiency=(Z_1=0.5,),)
        end
        expected = Agate.Construction.construct(
            definition; scalar_type=T, parameter_overrides=merge(overrides, updates)
        )
        tau = isempty(fixed) ? T[0.6, 0.5] : T[0.5]
        p = realize_trait_parameters(traits, tau)
        @test p == T[only(getproperty(expected.parameters, name)) for name in targets]
        @test eltype(p) === T
        @test traits.parameter_map(tau) == p
        @test Agate.Runtime.parameterized(traits.bgc, p; active_parameters=active).parameters.palatability_matrix[
            1, 1
        ] == expected.parameters.palatability_matrix[1, 1]
        p .= 0
        @test realize_trait_parameters(traits, tau0) == active.values
        @test traits.bgc.parameters == baseline

        if T === Float64 && isempty(fixed)
            for invalid in (Float32.(tau), Float64[], [tau[1], 1.1], [tau[1], NaN])
                @test_throws ArgumentError realize_trait_parameters(traits, invalid)
            end
            for selection in (
                (specificity=true,),
                (specificity=(:unknown,),),
                (specificity=(:Z_1, :Z_1),),
                (specificity=(),),
            )
                @test_throws ArgumentError active_traits(
                    definition; trait_selection=selection
                )
            end
            # Plan order: assimilation, palatability.
            gradient = trait_map_gradient(tau0, traits.parameter_map, [3.0, 2.0])
            @test gradient ≈ [-2log(5) * 5^(-tau0[1]), 3.0]

            # Check the same runtime map against an independently rebuilt RHS.
            u0 = [1.0, 0.1, 0.05]
            p = realize_trait_parameters(traits, tau)
            prob = Agate.Runtime.ode_problem(
                traits.bgc,
                u0,
                (0.0, 1.0);
                p,
                active_parameters=active,
                auxiliary=(PAR=100.0,),
            )
            reference = Agate.Runtime.ode_problem(
                expected, u0, (0.0, 1.0); auxiliary=(PAR=100.0,)
            )
            du, du_ref = similar(u0), similar(u0)
            prob.f(du, u0, p, 0.0)
            reference.f(du_ref, u0, reference.p, 0.0)
            @test du ≈ du_ref
        end
    end
end

# Exercise mixed scalar/matrix outputs and enough runtime entries to expose
# boxing from a long tuple map; two selected traits share the matrix output.
@testset "Enzyme trait map with mixed parameter storage" begin
    npz, scalar = direct_npz_definition(), derived_default_definition()
    definition = Agate.Processes.ModelDefinition(;
        components=merge(npz.components, scalar.components),
        processes=merge(npz.processes, scalar.processes),
        parameters=merge(npz.parameters, scalar.parameters),
    )
    traits = active_traits(
        definition;
        plankton_pfts=(
            P=(P=[1.0, 1.5, 2.0, 3.0, 4.0, 6.0],),
            Z=(Z=[10.0, 15.0, 20.0, 30.0, 40.0, 60.0],),
        ),
        parameter_overrides=(optimum_predator_prey_ratio=(Z_1=8.0,),),
        trait_selection=(
            base=true, specificity=(:Z_1,), optimum_predator_prey_ratio=(:Z_1,)
        ),
    )
    tau = copy(traits.values)
    f = traits.parameter_map
    @test f(tau) == realize_trait_parameters(traits, tau) == traits.active_parameters.values
    @test length(traits.active_parameters) == 38 # two scalars and a 6×6 matrix
    weights = collect(range(1.0, 2.0; length=length(traits.active_parameters)))
    gradient = trait_map_gradient(tau, f, weights)
    fd = map(eachindex(tau)) do i
        plus, minus = copy(tau), copy(tau)
        h = 1e-5
        plus[i] += h
        minus[i] -= h
        (
            weighted_trait_parameters(plus, f, weights) -
            weighted_trait_parameters(minus, f, weights)
        ) / (2h)
    end
    @test all(abs.(fd) .> 1e-5)
    @test isapprox(gradient, fd; rtol=1e-6, atol=1e-8)
end
