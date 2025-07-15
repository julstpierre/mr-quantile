#!/bin/bash
#SBATCH --job-name=mr_all
#SBATCH --array=0-29
#SBATCH --time=04:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=../results/logs/%x_%A_%a.out

cd ~/projects/def-mlegault/gf591137/BayesFuseMR/scripts/
module load r

SCENARIOS=("D_S1_null" "D_S1" "D_S2_null" "D_S2" "D_S3" "D_S4")

scenario_id=$((SLURM_ARRAY_TASK_ID / 5))    # 30 jobs / 6 scenarios = 5 chunks per scenario
chunk_index=$((SLURM_ARRAY_TASK_ID % 5))   # 0-based chunk index

SCENARIO=${SCENARIOS[$scenario_id]}
CHUNK_SIZE=200

START_CHUNK=$((chunk_index * CHUNK_SIZE + 1))
END_CHUNK=$((START_CHUNK + CHUNK_SIZE - 1))

echo "Running $SCENARIO, chunks $START_CHUNK to $END_CHUNK"
Rscript sim_estim.R $SCENARIO $START_CHUNK $END_CHUNK

# time Rscript run_chunk_all_methods.R D_S4 1 10
# /usr/bin/time -v Rscript run_chunk_all_methods.R D_S4 1 10
