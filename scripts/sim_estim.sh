#!/bin/bash
#SBATCH --job-name=mr_all
#SBATCH --array=0-79
#SBATCH --time=03:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=../results/logs/%x_%A_%a.out

cd ~/projects/def-mlegault/gf591137/BayesFuseMR/scripts/
module load r

SCENARIOS=("D_S1" "D_S2" "D_S3" "D_S4" "D_S5" "D_S6" "D_S7" "D_S8")

scenario_id=$((SLURM_ARRAY_TASK_ID / 10))    # 80 jobs / 8 scenarios = 10 chunks per scenario
chunk_index=$((SLURM_ARRAY_TASK_ID % 10))   # 0-based chunk index

SCENARIO=${SCENARIOS[$scenario_id]}
CHUNK_SIZE=500

START_CHUNK=$((chunk_index * CHUNK_SIZE + 1))
END_CHUNK=$((START_CHUNK + CHUNK_SIZE - 1))

echo "Running $SCENARIO, chunks $START_CHUNK to $END_CHUNK"
Rscript sim_estim.R $SCENARIO $START_CHUNK $END_CHUNK
