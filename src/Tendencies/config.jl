const SUPPORTED_GROWTH_FORMULATIONS = (:smith, :geider)
const SUPPORTED_ORGANIC_CYCLING_FORMULATIONS = (:simple_detritus, :dom_pom)
const SUPPORTED_ZOOPLANKTON_FORMULATIONS = (:preferential_grazing,)

struct TendencyConfig{Growth,OrganicCycling,Zooplankton,NutrientLimitation,Nutrients}
    nutrient_limitation::NutrientLimitation
    nutrients::Nutrients
end

function nutrient_limitation_operator(limitation::Symbol)
    limitation === :liebig ||
        throw(ArgumentError("Unsupported nutrient limitation formulation: $limitation"))
    return LiebigMinimum()
end

nutrient_limitation_operator(limitation::Union{LiebigMinimum,FrankTNorm}) = limitation

function nutrient_limitation_operator(limitation)
    throw(
        ArgumentError(
            "Unsupported nutrient limitation operator: $(typeof(limitation))"
        )
    )
end

"""
    TendencyConfig(; growth, organic_cycling, zooplankton=:preferential_grazing,
                     nutrient_limitation=:liebig, nutrients=())

Small configuration object used by reusable tendency builders.

`growth` selects the phytoplankton growth formulation. Supported values are
`:smith` and `:geider`. `organic_cycling` selects how plankton losses are routed
through organic matter. Supported values are `:simple_detritus` and `:dom_pom`.
`zooplankton` remains a separate selector because grazing formulations can vary
independently of growth and organic-matter cycling. `nutrient_limitation` accepts
`:liebig` or a limitation operator such as `FrankTNorm(sharpness)`.
"""
function TendencyConfig(;
    growth::Symbol,
    organic_cycling::Symbol,
    zooplankton::Symbol=:preferential_grazing,
    nutrient_limitation=:liebig,
    nutrients::Tuple=(),
)
    growth in SUPPORTED_GROWTH_FORMULATIONS ||
        throw(ArgumentError("Unsupported growth formulation: $growth"))
    organic_cycling in SUPPORTED_ORGANIC_CYCLING_FORMULATIONS ||
        throw(ArgumentError("Unsupported organic cycling formulation: $organic_cycling"))
    zooplankton in SUPPORTED_ZOOPLANKTON_FORMULATIONS ||
        throw(ArgumentError("Unsupported zooplankton formulation: $zooplankton"))

    limitation = nutrient_limitation_operator(nutrient_limitation)

    return TendencyConfig{
        growth,
        organic_cycling,
        zooplankton,
        typeof(limitation),
        typeof(nutrients),
    }(limitation, nutrients)
end
