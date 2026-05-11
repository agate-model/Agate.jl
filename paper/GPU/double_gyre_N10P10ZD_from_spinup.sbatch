#!/bin/bash
#SBATCH --time=0-4:00:00
#SBATCH --nodes=1
#SBATCH --mem=12000M
#SBATCH --cpus-per-task=4
#SBATCH --account=GEOG024542
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH --job-name=agate-double-gyre-N10P10ZD
#SBATCH --output=agate_double_gyre_N10P10ZD_%j.out
#SBATCH --error=agate_double_gyre_N10P10ZD_%j.err
module load apptainer/1.3.1

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

WORKDIR=/user/work/$(whoami)/Agate.jl
DEPOT=/user/work/$(whoami)/julia_depot_agate

mkdir -p "$DEPOT"

singularity exec \
  --nv \
  --cleanenv \
  --env JULIA_NUM_THREADS=$JULIA_NUM_THREADS \
  --env JULIA_PROJECT=/opt/agate \
  --env JULIA_DEPOT_PATH=$DEPOT:/opt/julia_depot \
  -B /user/work/$(whoami):/user/work/$(whoami) \
  "$WORKDIR/agate_cuda.sif" \
  julia --project=/opt/agate "$WORKDIR/paper/GPU/double_gyre_N10P10ZD_from_spinup.jl"
