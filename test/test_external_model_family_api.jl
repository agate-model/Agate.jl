using Test
using Agate
using Agate.Construction:
    capture_model_recipe,
    construct_factory,
    construct_factory_plus_manifest,
    decode_recipe,
    encode_recipe,
    replay_factory,
    resolve_construction_scalar_type
import Agate.Construction: recipe_family, recipe_factory
using Agate.Equations: CompiledEquation
using Agate.Factories:
    AbstractBGCFactory,
    ConstDefault,
    ParameterDefinition,
    ParameterSpec
import Agate.Factories:
    default_biogeochem_dynamics,
    default_community,
    default_plankton_dynamics,
    parameter_definitions
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers
using Oceananigans.Units: day

struct ExternalDecayFactory <: AbstractBGCFactory end

recipe_family(::ExternalDecayFactory) = :ExternalDecay
recipe_factory(::Val{:ExternalDecay}) = ExternalDecayFactory()

parameter_definitions(::ExternalDecayFactory) = (
    ParameterDefinition(ParameterSpec(:decay_rate, :scalar), ConstDefault(0.1 / day)),
)
default_community(::ExternalDecayFactory) = (;)
default_plankton_dynamics(::ExternalDecayFactory) = (;)
external_decay_tendency() = CompiledEquation(
    (bgc, x, y, z, t, X) -> -bgc.parameters.decay_rate * X
)
default_biogeochem_dynamics(::ExternalDecayFactory) = (X=external_decay_tendency,)

const EXTERNAL_MODEL_PUBLIC_IMPORTS = (
    :capture_model_recipe,
    :construct_factory,
    :construct_factory_plus_manifest,
    :decode_recipe,
    :encode_recipe,
    :replay_factory,
    :resolve_construction_scalar_type,
    :recipe_family,
    :recipe_factory,
    :CompiledEquation,
    :AbstractBGCFactory,
    :ConstDefault,
    :ParameterDefinition,
    :ParameterSpec,
    :default_biogeochem_dynamics,
    :default_community,
    :default_plankton_dynamics,
    :parameter_definitions,
)

const EXTERNAL_MODEL_PUBLIC_CONCEPT_GROUPS = (
    :factory_anchor,
    :parameter_definition,
    :community_default,
    :plankton_dynamics_default,
    :biogeochem_dynamics_default,
    :construction,
    :manifest_construction,
    :recipe_capture,
    :recipe_family_identity,
    :scalar_resolution,
    :serialization_and_replay,
)

const EXTERNAL_MODEL_EXTENSION_HOOKS = (
    :recipe_family,
    :recipe_factory,
    :parameter_definitions,
    :default_community,
    :default_plankton_dynamics,
    :default_biogeochem_dynamics,
)

@testset "External model-family API" begin
    factory = ExternalDecayFactory()
    grid = dummy_grid(Float32)
    parameter_overrides = (decay_rate=Float32(0.2 / day),)
    interaction_roles = (consumers=(), prey=())
    parameter_roles = (;)

    @test length(EXTERNAL_MODEL_PUBLIC_IMPORTS) == 18
    @test length(EXTERNAL_MODEL_PUBLIC_CONCEPT_GROUPS) == 11
    @test length(EXTERNAL_MODEL_EXTENSION_HOOKS) == 6
    @test resolve_construction_scalar_type(grid, nothing) === Float32
    @test recipe_family(factory) === :ExternalDecay
    @test replay_factory(
        capture_model_recipe(
            factory;
            community=(;),
            parameter_overrides,
            ecological_roles=(;),
            interaction_roles,
            parameter_roles,
            auxiliary_fields=(),
            scalar_type=Float32,
        ),
    ) isa ExternalDecayFactory

    recipe = capture_model_recipe(
        factory;
        community=(;),
        parameter_overrides,
        ecological_roles=(;),
        interaction_roles,
        parameter_roles,
        auxiliary_fields=(),
        scalar_type=Float32,
    )
    decoded = decode_recipe(encode_recipe(recipe))
    @test decoded == recipe
    @test recipe_family(replay_factory(decoded)) === :ExternalDecay

    construction_kwargs = (;
        community=(;),
        parameters=parameter_overrides,
        interaction_roles,
        parameter_roles,
        auxiliary_fields=(),
        grid,
        scalar_type=Float32,
    )
    bgc = construct_factory(factory; construction_kwargs...)
    bgc_with_manifest, manifest = construct_factory_plus_manifest(
        factory; construction_kwargs...
    )
    replayed = construct_factory(replay_factory(decoded); construction_kwargs...)

    @test required_biogeochemical_tracers(bgc) == (:X,)
    @test required_biogeochemical_tracers(replayed) == (:X,)
    @test manifest.scalar_type === Float32
    @test bgc.parameters == bgc_with_manifest.parameters == replayed.parameters
    @test bgc(Val(:X), 0, 0, 0, 0, 2f0) == -2f0 * parameter_overrides.decay_rate
end
