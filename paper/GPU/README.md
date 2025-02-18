# Reproducing the plots

To run the 1_degree_simulation inside the container:
`podman run --rm --device=nvidia.com/gpu=all -v $(pwd):/scripts -w /scripts localhost/my-gpu-julia-image julia 1_degree_simulation.jl`

To build the Podman container:
`podman build -f paper/GPU/gpu.Podmanfile -t my-container .`
