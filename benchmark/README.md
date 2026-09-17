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

Use `--quick` for a short CI/development smoke benchmark (2P/2Z and 4P/4Z only). The full benchmark uses 2, 4, 8, and 16 matched phytoplankton/zooplankton size classes.

The benchmark measures complete Oceananigans time steps rather than the scalar Agate rate helper. This intentionally includes the actual tracer-tendency execution strategy used by Oceananigans and therefore detects whether repeated shared-food expressions remain a practical scaling bottleneck after compiler optimization/fusion.
