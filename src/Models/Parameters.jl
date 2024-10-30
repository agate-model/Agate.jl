module Parameters

"""
A function to dynamically create a tracer dictionary for n plankton based on a template and merge with
existing biogeochemical_tracers dictionary.

# Example
```Julia

biogeochemical_tracers = Dict(
    "N" => :(
        net_linear_loss(plankton_list, linear_mortality, mortality_export_fraction) +
        net_quadratic_loss(
            plankton_list, quadratic_mortality, mortality_export_fraction
        ) +
        remineralization(D, detritus_remineralization) - net_photosynthetic_growth(
            N,
            plankton_list,
            PAR,
            maximum_growth_rate,
            nitrogen_half_saturation,
            alpha,
        )
    ),
    "D" => :(
        net_linear_loss(plankton_list, linear_mortality, 1 - mortality_export_fraction) +
        net_predation_assimilation_loss(
            plankton_list,
            holling_half_saturation,
            maximum_predation_rate,
            assimilation_efficiency,
            palatability,
        ) +
        net_quadratic_loss(
            plankton_list, quadratic_mortality, 1 - mortality_export_fraction
        ) - remineralization(D, detritus_remineralization)
    )
)

plankton_dt_template = :(plankton_dt(
    idx,
    N,
    plankton_list,
    PAR,
    linear_mortality,
    quadratic_mortality,
    maximum_growth_rate,
    holling_half_saturation,
    nitrogen_half_saturation,
    alpha,
    maximum_predation_rate,
    assimilation_efficiency,
    palatability,
))

tracers = create_tracers(Dict("P" => 2, "Z" => 2), biogeochemical_tracers, plankton_template)
```
"""
function add_plankton_tracers(
    plankton_counts::Dict{String,Int}, biogeochemical_tracers::Dict, plankton_template::Expr
)
    # Construct a list of plankton symbols (e.g., P1, P2, Z1, etc.)
    plankton_list = [
        Symbol("$type$i") for (type, count) in plankton_counts for i in 1:count
    ]

    # Initialize the tracers dictionary with base tracers, substituting `plankton_list`
    tracers = Dict(
        key => replace(value, :plankton_list => plankton_list) for
        (key, value) in biogeochemical_tracers
    )

    # Create dynamic entries for each plankton type
    idx = 1  # Keeps track of the global index for each plankton entry
    for (type, count) in plankton_counts
        for i in 1:count
            key = "$type$i"  # e.g., "P1", "Z3", "cocco5", etc.
            # Replace template variables with actual values
            expr = replace(
                plankton_template, Dict(:idx => idx, :plankton_list => plankton_list)
            )
            tracers[key] = expr
            idx += 1
        end
    end

    return tracers
end

end #module
