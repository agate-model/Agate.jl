# Photosynthetic growth

## Default 

```math
PC_j = PC_j^{max} \cdot \gamma^{nut}_j \cdot \gamma^{light}_j \cdot f^{temp}_j \cdot \gamma^{CO2}_j
```
   
```@docs
Agate.default_PCⱼ(PCⱼᵐᵃˣ, γⱼⁿᵘᵗ,  γⱼˡⁱᵍʰᵗ, fⱼᵗᵉᵐᵖ,  γⱼᶜᵒ²)
```

```math
\gamma^{light}_j = (1 - e^{(k_j^{sat} I)}) \cdot e^{k_{j}^{inh}} \cdot n_j^{light}
```

```@docs
Agate.γⱼˡⁱᵍʰᵗ(I, kⱼˢᵃᵗ, kⱼⁱⁿʰ, nⱼˡⁱᵍʰᵗ)
```


## Geider

The Geider growth formulation differs from the default DARWIN by accounting for the Chlorophyll content of the cell.


Four flavours:

non-spectral and fixed Chl quota,

spectral and fixed Chl quota,

non-spectral and variable Chl quota,

spectral and variable Chl quota.

geider_light_limitation * geider_light_inhibition * geider_light_saturated_growth

### Non-spectral

In this version of the formulation a mean spectral absorption coefficient and the total PAR is used.

```math

```

```@docs
Agate.Library.Light.non_spectral_carbon_yield()
```

### Spectral

In this version of the formulation wavelength-specific spectral absorption coefficient and the total irradiance at
each respective wavelength is used. These values are then summed to estimate maximum carbon yield. 

```math

```

```@docs
Agate.Library.Light.spectral_carbon_yield()
```


### Fixed Chlorophyll quota

In this version of the formulation, Chlorophyll is not included as a tracer in the model,
but a chlorophyll acclimation ratio is estimated based on irradiance and minimum and maximum 
chlorophyll:carbon ratios.

```math

```

```@docs
Agate.Library.Light.chlorophyll_carbon_ratio()
```


### Variable Chlorophyll quota

In this version of the formulation, Chlorophyll is included as a tracer in the model,
and the Chlorophyll:Carbon ratio is estimated based on the tracer values:


```math

```

```@docs
Agate.Library.Light.acclimated_chlorophyll_carbon_ratio()
```
