#!/bin/bash
#SBATCH --account=def-mlegault
#SBATCH --job-name=mr_all
#SBATCH --array=0-279
#SBATCH --time=10:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=../results/logs/%x_%A_%a.out

cd ~/projects/def-mlegault/gf591137/MRFuse/scripts/
module load r

SCENARIOS=($(seq -f "D_S%g" 1 56))

scenario_id=$((SLURM_ARRAY_TASK_ID / 5))   # 280 jobs / 56 scenarios = 5 chunks per scenario
chunk_index=$((SLURM_ARRAY_TASK_ID % 5))   # 0-based chunk index

SCENARIO=${SCENARIOS[$scenario_id]}
CHUNK_SIZE=400

START_CHUNK=$((chunk_index * CHUNK_SIZE + 1))
END_CHUNK=$((START_CHUNK + CHUNK_SIZE - 1))

echo "Running $SCENARIO, chunks $START_CHUNK to $END_CHUNK"
Rscript sim_estim.R $SCENARIO $START_CHUNK $END_CHUNK
