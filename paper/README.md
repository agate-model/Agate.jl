# Paper and figure reproduction

This directory contains scripts and container recipes used to reproduce figures
and results associated with Agate.

## GPU reproduction

`paper/GPU/` contains scripts intended to run on an NVIDIA GPU (via CUDA.jl),
plus a Podman/Docker recipe (`cuda.Podmanfile`) that builds a CUDA-capable Julia
environment.

See `paper/GPU/README.md` for the exact run commands.

## Notes

Figure-reproduction scripts may require additional packages such as Oceananigans,
OceanBioME, and CairoMakie beyond Agate's library dependencies.
