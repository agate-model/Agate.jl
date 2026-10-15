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

## Run the GPU test suite

The default package tests do not load CUDA. From the repository root, the
lightweight GPU-only smoke test (also used by Buildkite) is:

```bash
JULIA_PKG_PRECOMPILE_AUTO=0 julia --project=test/gpu -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
julia --project=test/gpu test/gpu/smoke.jl
```

To run the full CPU/AD suite with the same GPU checks, set `AGATE_TEST_CUDA=1`:

```bash
AGATE_TEST_CUDA=1 julia --project=. -e 'using Pkg; Pkg.test()'
```

When using the Podman image, mount the repository and run its project directly:

```bash
podman run --rm --device=nvidia.com/gpu=all \
  -v "$(pwd)":/scripts -w /scripts \
  agate-gpu env AGATE_TEST_CUDA=1 \
  julia --project=/scripts -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

## Cleanup

```bash
podman container prune
podman system prune -a -f --volumes
```
