#!/bin/bash -l 
#SBATCH --job-name=distributed_simulation
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=25
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --output=logs/sim_output_%j.log
#SBATCH --error=logs/sim_error_%j.log

# 1. Basic environment --------------------------------------------------------
module purge
# module load python/3.8
# module load gnu9/9.3.0
# module load openmpi4/4.0.5  # Add this if available on  cluster
module load gnu15/15.2.0
module load openmpi5/5.0.10

# Fix the conda activation command
# source ~/.local/miniconda3/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
conda activate simuscale_env  # Use the correct environment name
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

mkdir -p logs


echo "=== job $SLURM_JOB_ID started on $(date) at $(hostname) ==="
echo "nodes : $SLURM_NODELIST"
echo "ranks : $SLURM_NTASKS  ( $SLURM_NNODES × $SLURM_TASKS_PER_NODE )"
echo "python : $(which python3)"

# 2. Run the MPI program ------------------------------------------------------
# Use a supported MPI type (pmi2) or use mpirun instead
srun --mpi=pmi2 python3 -u notebook_server_mpi.py 

echo "=== job $SLURM_JOB_ID finished at $(date) ==="
