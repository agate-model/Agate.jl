# Reproducing the plots

Note:
Podman should be ran from `./Agate.jl` not `./Agate.jl/paper/GPU`

To build the Podman container:
`podman build -f paper/GPU/cuda.Podmanfile -t my-container .`

To run the column:

`podman run --rm --device=nvidia.com/gpu=all -v $(pwd):/scripts -w /scripts localhost/my-container julia paper/GPU/column_N2P2ZD_agate.jl`

To delete old containers:

`podman container prune`

Force remove everything:

`podman system prune -a -f --volumes`

To delete `tmp` cache:

`sudo rm -rf /tmp/*`
