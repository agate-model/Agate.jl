"
    α * PAR / sqrt(μ₀ ^ 2 + α ^ 2 * PAR ^ 2)

Smith 1936 formulation of light limitation (also see Evans and Parslow, 1985).

Where: 
α = Initial photosynthetic slope
PAR = Photosynthetic Active Radiation
μ₀ = Maximum growth rate at T = 0 °C (this seems weird?, from Kuhn 2015)

"
function smith_light_limitation(PAR, α, μ₀)
    α * PAR / sqrt(μ₀ ^ 2 + α ^ 2 * PAR ^ 2)
end