struct NutrientCoupling{Tracer,HalfSaturation,Stoichiometry,Remineralization} end

"""
    nutrient_coupling(tracer, half_saturation; stoichiometry=:one, remineralization=())

Describe how one dissolved nutrient tracer is coupled to growth limitation,
stoichiometric uptake, and remineralization sources.
"""
nutrient_coupling(
    tracer::Symbol,
    half_saturation::Symbol;
    stoichiometry::Symbol=:one,
    remineralization::Tuple=(),
) = NutrientCoupling{tracer,half_saturation,stoichiometry,remineralization}()

@inline function tracer_name(::NutrientCoupling{Tracer,HalfSaturation,Stoichiometry,Remineralization}) where {Tracer,HalfSaturation,Stoichiometry,Remineralization}
    return Tracer
end

@inline function stoichiometry_name(::NutrientCoupling{Tracer,HalfSaturation,Stoichiometry,Remineralization}) where {Tracer,HalfSaturation,Stoichiometry,Remineralization}
    return Stoichiometry
end

@inline function remineralization_sources(::NutrientCoupling{Tracer,HalfSaturation,Stoichiometry,Remineralization}) where {Tracer,HalfSaturation,Stoichiometry,Remineralization}
    return Remineralization
end

function target_coupling(nutrients::Tuple, target::Symbol)
    @inbounds for nutrient in nutrients
        tracer_name(nutrient) === target && return nutrient
    end
    error("No nutrient coupling found for target tracer: $target")
end
