"""Evaluate the canonical scientific rate for a mortality formulation."""
@inline process_rate(::LinearMortality, biomass, coefficient) = linear_loss(biomass, coefficient)
@inline process_rate(::QuadraticMortality, biomass, coefficient) = quadratic_loss(biomass, coefficient)
@inline process_rate(process::Mortality, biomass, coefficient) =
    process_rate(formulation(process), biomass, coefficient)
