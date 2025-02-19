# Reproducing the plots

To run the 1_degree_simulation inside the container:
`podman run --rm --device=nvidia.com/gpu=all -v $(pwd):/scripts -w /scripts localhost/my-gpu-julia-image julia 1_degree_simulation.jl`

Note:
Podman should be ran from `./Agate.jl` not `./Agate.jl/paper/GPU`

To build the Podman container:
`podman build -f paper/GPU/gpu.Podmanfile -t my-container .`


To run the column:

`podman run --rm --device=nvidia.com/gpu=all -v $(pwd):/scripts -w /scripts localhost/my-gpu-julia-image julia paper/GPU/column_agate.jl`


To delete old containers:

`podman container prune`


To delete `tmp` cache:

`sudo rm -rf /tmp/*`