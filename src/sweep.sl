#!/bin/bash -e
#SBATCH --job-name=0.001     # job name (shows up in the queue)
#SBATCH --time=400:00:00      # Walltime (HH:MM:SS)
#SBATCH --mem=2GB         # Memory in MB
#SBATCH --array=0-8
#SBATCH --cpus-per-task=24


module load Python/3.10.5-gimkl-2022a
python working.py ${SLURM_ARRAY_TASK_ID}
