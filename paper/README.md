# Paper and figure reproduction

This directory contains scripts and container recipes used to reproduce figures
and results associated with Agate.

## GPU reproduction

`paper/GPU/` contains scripts intended to run on an NVIDIA GPU (via CUDA.jl),
plus a Podman/Docker recipe (`cuda.Podmanfile`) that builds a CUDA-capable Julia
environment.

See `paper/GPU/README.md` for the exact run commands.

## Notes

These scripts are not part of the package API and may depend on additional
packages (e.g. Oceananigans, OceanBioME, CairoMakie) beyond the minimal runtime
needed to use Agate as a library.
