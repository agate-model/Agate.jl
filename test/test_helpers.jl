"""Small test-only helpers.

These exist to keep unit tests independent of specific Oceananigans grid
constructors. In Agate, **grid element type decides precision**, so tests can
use a minimal grid object that exposes `eltype(::grid)` and
`Oceananigans.Architectures.architecture(::grid)`.
"""

import Oceananigans.Architectures: architecture, CPU, GPU

"""A minimal grid stand-in for testing constructor precision/architecture inference.

`Oceananigans` determines the active architecture from the grid. For CPU architectures,
`architecture(grid)` is typically a singleton like `CPU()`, but GPU architectures
carry backend state (e.g. a CUDA backend) and are not nullary-constructible.

This test grid therefore stores an *architecture instance* and returns it directly.
"""
struct DummyGrid{T,Arch}
    arch::Arch
end

Base.eltype(::DummyGrid{T,Arch}) where {T,Arch} = T
architecture(g::DummyGrid) = g.arch

"""Construct a `DummyGrid` that behaves like an Oceananigans grid."""
dummy_grid(::Type{T}; arch=CPU()) where {T<:AbstractFloat} = DummyGrid{T,typeof(arch)}(arch)

struct MultiNutrientTestFactory <: Agate.Factories.AbstractBGCFactory end

const MULTI_NUTRIENT_LIEBIG = Agate.Tendencies.TendencyConfig(;
    growth=:smith,
    organic_cycling=:dom_pom,
    nutrient_limitation=:liebig,
    nutrients=(
        Agate.Tendencies.nutrient_coupling(
            :DIN,
            :half_saturation_DIN;
            stoichiometry=:nitrogen_to_carbon,
            remineralization=((:DON, :organic_remineralization), (:PON, :organic_remineralization)),
        ),
        Agate.Tendencies.nutrient_coupling(
            :PO4,
            :half_saturation_PO4;
            stoichiometry=:phosphorus_to_carbon,
            remineralization=((:DOP, :organic_remineralization), (:POP, :organic_remineralization)),
        ),
    ),
)

const MULTI_NUTRIENT_FRANK = Agate.Tendencies.TendencyConfig(;
    growth=:smith,
    organic_cycling=:dom_pom,
    nutrient_limitation=Agate.Library.Nutrients.FrankTNorm(50),
    nutrients=MULTI_NUTRIENT_LIEBIG.nutrients,
)

multi_nutrient_phyto(plankton_idx::Int) = Agate.Tendencies.phytoplankton_tendency(
    MULTI_NUTRIENT_LIEBIG; plankton_idx
)
multi_nutrient_zoo(plankton_idx::Int) = Agate.Tendencies.zooplankton_tendency(
    MULTI_NUTRIENT_LIEBIG; plankton_idx
)

function Agate.Factories.default_plankton_dynamics(::MultiNutrientTestFactory)
    return (Z=multi_nutrient_zoo, P=multi_nutrient_phyto)
end

function Agate.Factories.default_community(::MultiNutrientTestFactory)
    empty_pft = Agate.Configuration.PFTSpecification()
    return (
        Z=(; diameters=[20.0], pft=empty_pft),
        P=(; diameters=[2.0], pft=empty_pft),
    )
end

function Agate.Factories.default_biogeochem_dynamics(::MultiNutrientTestFactory)
    config = MULTI_NUTRIENT_LIEBIG
    return (
        DIC=() -> Agate.Tendencies.inorganic_tendency(
            config;
            target=:DIC,
            remineralization=((:DOC, :organic_remineralization), (:POC, :organic_remineralization)),
            stoichiometry=:one,
        ),
        DIN=() -> Agate.Tendencies.inorganic_tendency(config; target=:DIN),
        PO4=() -> Agate.Tendencies.inorganic_tendency(config; target=:PO4),
        DOC=() -> Agate.Tendencies.organic_matter_tendency(
            config; target=:DOC, remineralization=:organic_remineralization, fraction=:DOM
        ),
        POC=() -> Agate.Tendencies.organic_matter_tendency(
            config; target=:POC, remineralization=:organic_remineralization, fraction=:POM
        ),
        DON=() -> Agate.Tendencies.organic_matter_tendency(
            config;
            target=:DON,
            remineralization=:organic_remineralization,
            fraction=:DOM,
            stoichiometry=:nitrogen_to_carbon,
        ),
        PON=() -> Agate.Tendencies.organic_matter_tendency(
            config;
            target=:PON,
            remineralization=:organic_remineralization,
            fraction=:POM,
            stoichiometry=:nitrogen_to_carbon,
        ),
        DOP=() -> Agate.Tendencies.organic_matter_tendency(
            config;
            target=:DOP,
            remineralization=:organic_remineralization,
            fraction=:DOM,
            stoichiometry=:phosphorus_to_carbon,
        ),
        POP=() -> Agate.Tendencies.organic_matter_tendency(
            config;
            target=:POP,
            remineralization=:organic_remineralization,
            fraction=:POM,
            stoichiometry=:phosphorus_to_carbon,
        ),
    )
end

function Agate.Factories.parameter_definitions(::MultiNutrientTestFactory)
    F = Agate.Factories
    producer = F.DiameterIndexedMaterialization(:producers; fill_value=0)
    consumer = F.DiameterIndexedMaterialization(:consumers; fill_value=0)
    return (
        F.ParameterDefinition(
            F.ParameterSpec(:organic_remineralization, :scalar),
            F.ConstDefault(0.1213 / 86400),
        ),
        F.ParameterDefinition(
            F.ParameterSpec(:nitrogen_to_carbon, :scalar), F.ConstDefault(16 / 106)
        ),
        F.ParameterDefinition(
            F.ParameterSpec(:phosphorus_to_carbon, :scalar), F.ConstDefault(1 / 106)
        ),
        F.ParameterDefinition(
            F.ParameterSpec(:DOM_POM_fractionation, :scalar), F.ConstDefault(0.5)
        ),
        F.ParameterDefinition(
            F.ParameterSpec(:linear_mortality, :vector; axes=:plankton), F.FillDefault(8e-7)
        ),
        F.ParameterDefinition(
            F.ParameterSpec(
                :quadratic_mortality, :vector; axes=:plankton, materialization=consumer
            ),
            F.DiameterIndexedVectorDefault(1e-6, :consumers; default=0),
        ),
        F.ParameterDefinition(
            F.ParameterSpec(
                :maximum_growth_rate, :vector; axes=:plankton, materialization=producer
            ),
            F.DiameterIndexedVectorDefault(2 / 86400, :producers; default=0),
        ),
        F.ParameterDefinition(
            F.ParameterSpec(
                :half_saturation_DIN, :vector; axes=:plankton, materialization=producer
            ),
            F.DiameterIndexedVectorDefault(0.5, :producers; default=0),
        ),
        F.ParameterDefinition(
            F.ParameterSpec(
                :half_saturation_PO4, :vector; axes=:plankton, materialization=producer
            ),
            F.DiameterIndexedVectorDefault(0.5, :producers; default=0),
        ),
        F.ParameterDefinition(
            F.ParameterSpec(:alpha, :vector; axes=:plankton, materialization=producer),
            F.DiameterIndexedVectorDefault(0.1 / 86400, :producers; default=0),
        ),
        F.ParameterDefinition(
            F.ParameterSpec(
                :maximum_predation_rate, :vector; axes=:plankton, materialization=consumer
            ),
            F.DiameterIndexedVectorDefault(0.5 / 86400, :consumers; default=0),
        ),
        F.ParameterDefinition(
            F.ParameterSpec(
                :holling_half_saturation, :vector; axes=:plankton, materialization=consumer
            ),
            F.DiameterIndexedVectorDefault(1.0, :consumers; default=0),
        ),
        F.ParameterDefinition(
            F.ParameterSpec(:palatability_matrix, :matrix; axes=(:consumer, :prey)),
            F.FillDefault(1.0),
        ),
        F.ParameterDefinition(
            F.ParameterSpec(:assimilation_matrix, :matrix; axes=(:consumer, :prey)),
            F.FillDefault(0.32),
        ),
    )
end

function multi_nutrient_test_model(; grid=nothing)
    return Agate.Construction.construct_factory(
        MultiNutrientTestFactory();
        ecological_roles=(phytoplankton=(:P,), zooplankton=(:Z,)),
        interaction_roles=(consumers=(:Z,), prey=(:P,)),
        parameter_roles=(producers=(:P,), consumers=(:Z,)),
        grid,
    )
end
