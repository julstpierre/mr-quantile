#!/bin/bash
#SBATCH --account=def-mlegault
#SBATCH --job-name=sim_indv_dat
#SBATCH --array=1-42
#SBATCH --time=05:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=../data/simulations/logs/%x_%A_%a.out

# -----------------------------
# Change directory and load R
# -----------------------------
cd ~/projects/def-mlegault/gf591137/mr-quantile/scripts/
module load r

# -----------------------------
# Scenario ID from array index
# -----------------------------
SCENARIO_ID=${SLURM_ARRAY_TASK_ID}

echo "Running scenario D_S${SCENARIO_ID}"

# -----------------------------
# Run simulation
# -----------------------------
Rscript sim_indv_dat.R ${SCENARIO_ID}
