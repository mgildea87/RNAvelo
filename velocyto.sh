#!/bin/bash -l
#SBATCH -J velocyto
#SBATCH --mem=128000
#SBATCH --cpus-per-task=20
#SBATCH -p cpu_short
#SBATCH --export=ALL
#SBATCH --time=10:00:00
#SBATCH --array=1

module unload miniconda3/4.6.14
module load miniconda3/cpu/4.9.2
conda activate /gpfs/data/fisherlab/conda_envs/scVelo/

command=$(sed -n "$SLURM_ARRAY_TASK_ID"p velocyto_commands_CV7292.txt)
srun $command