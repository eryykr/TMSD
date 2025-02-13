#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=240:00:00
#SBATCH --partition=long
#SBATCH --clusters=arc
#SBATCH --job-name=ffs_fast
#SBATCH --error=ffs_%A_%a.err
#SBATCH --output=ffs_%A_%a.out


# Print this sub-job's task ID
# echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

./run_replica.sh

