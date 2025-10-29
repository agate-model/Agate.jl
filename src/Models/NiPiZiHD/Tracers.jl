module Tracers

export nutrients, phytoplankton, zooplankton, heterotrophs, detritus

# -------------------------------------------------------- #
#                     Base functions  - remove at some point
# -------------------------------------------------------- #

function photosynthetic_growth(
    Rp::Real,
    P::Real,
    PAR::Real,
    maximum_growth_rate::Real,  
    nutrient_half_saturation::Real,
    alpha::Real,
)
    
    return P * maximum_growth_rate * monod_limitation(Rp, nutrient_half_saturation) * light_limitation_smith(PAR, alpha, maximum_growth_rate)
         
end

function heterotrophic_growth(
    Rh::Real,
    H::Real,
    maximum_growth_rate::Real,
    detritus_half_saturation::Real,
    remin_frac_to_N::Real
)    
    return (1 - remin_frac_to_N) * H * maximum_growth_rate * monod_limitation(Rh, detritus_half_saturation)
end

function heterotrophic_remin(
    Rh::Real,
    H::Real,
    maximum_growth_rate::Real,
    detritus_half_saturation::Real,
    remin_frac_to_N::Real,
)
    return remin_frac_to_N * heterotrophic_growth(Rh, H, maximum_growth_rate, detritus_half_saturation, 0.)
end 

function grazing_loss(
    P::Real,
    Z::Real,
    maximum_grazing_rate::Real,
    holling_half_saturation::Real
)    
    return Z * maximum_grazing_rate * holling_type_2(P, holling_half_saturation)
end

function grazing_gain(
    P::Real,
    Z::Real,
    maximum_grazing_rate::Real,
    holling_half_saturation::Real,
	assimilation_efficiency::Real,
)    
    return Z * maximum_grazing_rate * holling_type_2(P, holling_half_saturation) * assimilation_efficiency
end

function messy_grazing(
	P::Real,
    Z::Real,
    maximum_grazing_rate::Real,
    holling_half_saturation::Real,
	assimilation_efficiency::Real
)
	
	return (1 - assimilation_efficiency) * grazing_loss(P, Z, maximum_grazing_rate, holling_half_saturation)
end

function linear_loss(x, rate)
    return x * rate
end

function quadratic_loss(x, rate)
    return x^2 * rate
end


# -------------------------------------------------------- #
#                     Dynamics
# -------------------------------------------------------- #

"""
"""
function nutrients(plankton_array, plankton_name, plankton_idx)
	
	plankton_symbol = Symbol(plankton_name)
	
	return :(
		sum(
			
			heterotrophic_remin.(
				D,
				$(plankton_symbol),
				maximum_growth_rate[$plankton_idx],
				detritus_half_saturation[$plankton_idx],
				remin_frac_to_N
			)
			
		) - 
		sum(
			
			photosynthetic_growth.(
				N,
				$(plankton_symbol),
				PAR,
				maximum_growth_rate[$plankton_idx],
				nutrient_half_saturation[$plankton_idx],
				alpha[$plankton_idx]
			)
			
		)
	)
end

"""
"""
function phytoplankton(plankton_array, plankton_name, plankton_idx)
	
	plankton_symbol = Symbol(plankton_name)
		
	return :(
		photosynthetic_growth.(
			
			N,
			$(plankton_symbol),
			PAR,
			maximum_growth_rate[$plankton_idx],
			nutrient_half_saturation[$plankton_idx],
			alpha[$plankton_idx]
			
		) - sum(
			
			grazing_loss.(
				$(plankton_symbol),
				[$(plankton_array...)],
				maximum_grazing_rate[$plankton_idx],
				holling_half_saturation,
				
			)
		) - sum(
			
			linear_loss.(
				$(plankton_symbol), linear_mortality[$plankton_idx]
			)
			
		)
	)
end

"""
"""
function zooplankton(plankton_array, plankton_name, plankton_idx)
	
	plankton_symbol = Symbol(plankton_name)
	
	return :(
		
		sum(
			
			grazing_gain.(
				[$(plankton_array...)]$(plankton_symbol),
				$(plankton_symbol),
				maximum_grazing_rate[$plankton_idx],
				holling_half_saturation,
				assimilation_efficiency
			)
			
		) - sum(
			
				linear_loss.(
					$(plankton_symbol), linear_mortality[$plankton_idx]
			)
			
		) - sum(
			
			quadratic_loss.(
				$(plankton_symbol), quadratic_mortality[$plankton_idx]
			)
			
		)
		
	)
end

"""
"""
function heterotrophs(plankton_array, plankton_name, plankton_idx)
	
	plankton_symbol = Symbol(plankton_name)
	
	return :(
		
		heterotrophic_growth.(
			D,
			$(plankton_symbol),
			maximum_growth_rate[$plankton_idx],
			detritus_half_saturation,
			remin_frac_to_N	
		) - sum(
			
			grazing_loss.(
				$(plankton_symbol),
				[$(plankton_array...)],
				maximum_grazing_rate[$plankton_idx],
				holling_half_saturation,		
			)
			
		) - sum(
			
			linear_loss.(
				[$(plankton_array...)], linear_mortality
			)
			
		)
		
	)
end

"""
"""
function detritus(plankton_array, plankton_name, plankton_idx)
	
	plankton_symbol = Symbol(plankton_name)
	
	return :(
		sum(
			
			linear_loss.([$(plankton_array...)], linear_mortality)
			
		) + 
		sum(
			
			quadratic_loss.([$(plankton_array...)], quadratic_mortality)
			
		) +
		sum(
			
			messy_grazing.(
				[$(plankton_array...)]$(plankton_symbol),
				$(plankton_symbol),
				maximum_grazing_rate[$plankton_idx],
				holling_half_saturation,
				assimilation_efficiency
			)
			
		) -
		sum(
			
			heterotrophic_growth.(
				D,
				$(plankton_symbol),
				maximum_growth_rate[$plankton_idx],
				detritus_half_saturation,
				0.
			)
			
		)
	)
end


end # module
