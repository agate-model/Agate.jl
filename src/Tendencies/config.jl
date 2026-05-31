const SUPPORTED_GROWTH = (:smith, :geider)
const SUPPORTED_LIMITATION = (:liebig,)
const SUPPORTED_CYCLING = (:simple_detritus, :dom_pom)

struct TendencyConfig{Growth,Limitation,Cycling,Nutrients}
    nutrients::Nutrients
end

"""
    TendencyConfig(; growth, nutrient_limitation=:liebig, cycling, nutrients=())

Small configuration object used by reusable tendency builders. It keeps public
builder names stable while allowing growth, nutrient-limitation, and detrital
cycling choices to vary internally.
"""
function TendencyConfig(;
    growth::Symbol,
    nutrient_limitation::Symbol=:liebig,
    cycling::Symbol,
    nutrients::Tuple=(),
)
    growth in SUPPORTED_GROWTH || error("Unsupported growth formulation: $growth")
    nutrient_limitation in SUPPORTED_LIMITATION ||
        error("Unsupported nutrient limitation rule: $nutrient_limitation")
    cycling in SUPPORTED_CYCLING || error("Unsupported detrital cycling: $cycling")

    return TendencyConfig{growth,nutrient_limitation,cycling,typeof(nutrients)}(nutrients)
end
