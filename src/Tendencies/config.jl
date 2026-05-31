const SUPPORTED_TENDENCY_FORMULATIONS = (:smith_detritus, :geider_dom_pom)
const SUPPORTED_ZOOPLANKTON_FORMULATIONS = (:preferential_grazing,)
const SUPPORTED_NUTRIENT_LIMITATIONS = (:liebig,)

struct TendencyConfig{Formulation,Zooplankton,NutrientLimitation,Nutrients}
    nutrients::Nutrients
end

"""
    TendencyConfig(; formulation, zooplankton=:preferential_grazing,
                     nutrient_limitation=:liebig, nutrients=())
    TendencyConfig(formulation; zooplankton=:preferential_grazing,
                   nutrient_limitation=:liebig, nutrients=())

Small configuration object used by reusable tendency builders.

`formulation` selects a coherent bundle of coupled phytoplankton, inorganic,
organic-matter, and detrital tendency assumptions. Supported formulations are
`:smith_detritus` and `:geider_dom_pom`. `zooplankton` remains a separate
selector because grazing formulations can vary independently of nutrient and
organic-matter cycling. `nutrient_limitation` selects how nutrient limitation
terms are aggregated by nutrient-coupled growth formulations.
"""
function TendencyConfig(;
    formulation::Symbol,
    zooplankton::Symbol=:preferential_grazing,
    nutrient_limitation::Symbol=:liebig,
    nutrients::Tuple=(),
)
    formulation in SUPPORTED_TENDENCY_FORMULATIONS ||
        error("Unsupported tendency formulation: $formulation")
    zooplankton in SUPPORTED_ZOOPLANKTON_FORMULATIONS ||
        error("Unsupported zooplankton formulation: $zooplankton")
    nutrient_limitation in SUPPORTED_NUTRIENT_LIMITATIONS ||
        error("Unsupported nutrient limitation rule: $nutrient_limitation")

    return TendencyConfig{
        formulation,
        zooplankton,
        nutrient_limitation,
        typeof(nutrients),
    }(nutrients)
end

TendencyConfig(
    formulation::Symbol;
    zooplankton::Symbol=:preferential_grazing,
    nutrient_limitation::Symbol=:liebig,
    nutrients::Tuple=(),
) = TendencyConfig(;
    formulation,
    zooplankton,
    nutrient_limitation,
    nutrients,
)
