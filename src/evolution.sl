#!/bin/bash -e
#SBATCH --job-name=milanc     # job name (shows up in the queue)
#SBATCH --time=23:59:00      # Walltime (HH:MM:SS)
#SBATCH --mem=2048MB          # Memory in MB
#SBATCH --cpus-per-task=2 


module load GCCcore/7.4.0
g++ -O3 fission.cpp -o out
./out
