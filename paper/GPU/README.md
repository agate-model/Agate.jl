# GPU reproduction (Podman)

This directory contains scripts intended to run on an NVIDIA GPU via CUDA.jl.

## Build the container

Run these commands from the **repository root** (so the Podmanfile can copy the
Agate project into the image):

```bash
podman build -f paper/GPU/cuda.Podmanfile -t agate-gpu .
```

## Run a script

Mount the repository into the container so outputs are written to your working
directory, but use the image's project environment:

```bash
podman run --rm --device=nvidia.com/gpu=all \
  -v "$(pwd)":/scripts -w /scripts \
  agate-gpu julia --project=/opt/julia_env paper/GPU/column_N2P2ZD_agate.jl
```

## Cleanup

```bash
podman container prune
podman system prune -a -f --volumes
```
