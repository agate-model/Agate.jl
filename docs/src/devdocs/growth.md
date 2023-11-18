# Growth

## Default

### Realised photosynthetic rate
```math
PC_j = PC_j^{max} \cdot \gamma^{nut}_j \cdot \gamma^{light}_j \cdot f^{temp}_j \cdot \gamma^{CO2}_j
```
   
```@docs
Agate.PCⱼ(PCⱼᵐᵃˣ, γⱼⁿᵘᵗ,  γⱼˡⁱᵍʰᵗ, fⱼᵗᵉᵐᵖ,  γⱼᶜᵒ²)
```

### Light inhibition (default)

```math
\gamma^{light}_j = (1 - e^{(k_j^{sat} I)}) \cdot e^{k_{j}^{inh}} \cdot n_j^{light}
```

```@docs
Agate.default_γⱼˡⁱᵍʰᵗ(I, kⱼˢᵃᵗ, kⱼⁱⁿʰ, nⱼˡⁱᵍʰᵗ)
```

## Geider

### Light inhibition (geider)
An alternative formulation of light limitation which considers effect of Chl:C ratios and Fe limitation on photosynthesis.


```math
\gamma^{light}_j = (1 - e^{(k_j^{sat} I)}) \cdot e^{k_{j}^{inh}} \cdot n_j^{light}
```

```@docs
Agate.geider_γⱼˡⁱᵍʰᵗ(γⱼⁱⁿʰ, Iₜₒₜ, Iₘᵢₙ, Chl¨Cⱼ, γⱼᶜᶠᵉ, fⱼᵃᴵ)
```