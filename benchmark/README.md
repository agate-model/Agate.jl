# Agate performance benchmarks

Set up the benchmark environment from the repository root:

```bash
julia --project=benchmark -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
```

Run the CPU scaling benchmark:

```bash
julia --project=benchmark benchmark/grazing_scaling.jl
```

Run the same benchmark on an NVIDIA GPU:

```bash
julia --project=benchmark benchmark/grazing_scaling.jl --gpu
```

Use `--quick` for a short CI/development smoke benchmark (2P/2Z and 4P/4Z only). The default benchmark uses 2, 4, and 8 matched phytoplankton/zooplankton size classes. Use `--extended` to add 16P/16Z; this is intentionally opt-in because pathological compiler scaling can make that case memory-intensive.

The benchmark reports both model construction plus the first warm time step and the median warm time-step cost. It measures complete Oceananigans time steps rather than the scalar Agate rate helper, so it captures the generated-code and runtime scaling users actually experience.
