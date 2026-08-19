include("production_oracle_helpers.jl")

reference_model, _ = build_nipizd_production_box(Float64)
set_nipizd_production_state!(reference_model, Float64)
reference = production_trajectory(reference_model, NIPIZD_PRODUCTION_TRACERS)

perturbed_model, _ = build_nipizd_production_box(Float64)
set_nipizd_production_state!(perturbed_model, Float64; ulp_perturb=true)
perturbed = production_trajectory(perturbed_model, NIPIZD_PRODUCTION_TRACERS)
sensitivity = max_trajectory_difference(reference, perturbed, NIPIZD_PRODUCTION_TRACERS)

reference_path = joinpath(@__DIR__, "reference", "nipizd_production_v0.10.1.csv")
write_trajectory_reference(
    reference_path,
    reference,
    NIPIZD_PRODUCTION_TRACERS;
    metadata=(;
        trajectory_rtol=sensitivity.max_rel,
        trajectory_atol=sensitivity.max_abs,
        ulp_sensitivity_max_rel=sensitivity.max_rel,
        ulp_sensitivity_max_abs=sensitivity.max_abs,
    ),
)
println(reference_path)
println("ulp_sensitivity_max_rel=", sensitivity.max_rel)
println("ulp_sensitivity_max_abs=", sensitivity.max_abs)
