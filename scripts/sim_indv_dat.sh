#!/bin/bash
#SBATCH --job-name=sim_indv_dat
#SBATCH --array=1-8
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=6
#SBATCH --output=../results/logs/%x_%A_%a.out

# Change directory and load R
cd ~/projects/def-mlegault/gf591137/BayesFuseMR/scripts/
module load r

# Pass the array index to R script
Rscript sim_indv_dat.R ${SLURM_ARRAY_TASK_ID}

# Run sequentially
# for i in {1..8}
# do
#   Rscript sim_indv_dat.R $i
# done