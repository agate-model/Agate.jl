# Growth

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
