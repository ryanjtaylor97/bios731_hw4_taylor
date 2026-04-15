#!/bin/bash
#SBATCH --array=1-3
#SBATCH --job-name=run_mix_job
#SBATCH --partition=wrobel
#SBATCH --output=run_mixture_sim.out
#SBATCH --error=run_mixture_sim.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript run_mixture_sim.R $JOBID


