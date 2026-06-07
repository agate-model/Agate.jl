const SUPPORTED_GROWTH_FORMULATIONS = (:smith, :geider)
const SUPPORTED_ORGANIC_CYCLING_FORMULATIONS = (:simple_detritus, :dom_pom)
const SUPPORTED_ZOOPLANKTON_FORMULATIONS = (:preferential_grazing,)
const SUPPORTED_NUTRIENT_LIMITATIONS = (:liebig,)

struct TendencyConfig{Growth,OrganicCycling,Zooplankton,NutrientLimitation,Nutrients}
    nutrients::Nutrients
end

"""
    TendencyConfig(; growth, organic_cycling, zooplankton=:preferential_grazing,
                     nutrient_limitation=:liebig, nutrients=())

Small configuration object used by reusable tendency builders.

`growth` selects the phytoplankton growth formulation. Supported values are
`:smith` and `:geider`. `organic_cycling` selects how plankton losses are routed
through organic matter. Supported values are `:simple_detritus` and `:dom_pom`.
`zooplankton` remains a separate selector because grazing formulations can vary
independently of growth and organic-matter cycling. `nutrient_limitation` selects
how nutrient limitation terms are aggregated by nutrient-coupled growth
formulations.
"""
function TendencyConfig(;
    growth::Symbol,
    organic_cycling::Symbol,
    zooplankton::Symbol=:preferential_grazing,
    nutrient_limitation::Symbol=:liebig,
    nutrients::Tuple=(),
)
    growth in SUPPORTED_GROWTH_FORMULATIONS ||
        error("Unsupported growth formulation: $growth")
    organic_cycling in SUPPORTED_ORGANIC_CYCLING_FORMULATIONS ||
        error("Unsupported organic cycling formulation: $organic_cycling")
    zooplankton in SUPPORTED_ZOOPLANKTON_FORMULATIONS ||
        error("Unsupported zooplankton formulation: $zooplankton")
    nutrient_limitation in SUPPORTED_NUTRIENT_LIMITATIONS ||
        error("Unsupported nutrient limitation rule: $nutrient_limitation")

    return TendencyConfig{
        growth,
        organic_cycling,
        zooplankton,
        nutrient_limitation,
        typeof(nutrients),
    }(nutrients)
end
