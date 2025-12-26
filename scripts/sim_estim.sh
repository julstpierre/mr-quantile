#!/bin/bash
#SBATCH --account=def-mlegault
#SBATCH --job-name=mr_all
#SBATCH --array=0-1679
#SBATCH --time=10:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=../results/logs/%x_%A_%a.out

cd ~/projects/def-mlegault/gf591137/MRFuse/scripts/
module load r

SCENARIOS=($(seq -f "D_S%g" 1 168))

scenario_id=$((SLURM_ARRAY_TASK_ID / 10))   # 1680 jobs / 168 scenarios = 10 chunks per scenario
chunk_index=$((SLURM_ARRAY_TASK_ID % 10))   # 0-based chunk index

SCENARIO=${SCENARIOS[$scenario_id]}
CHUNK_SIZE=200

START_CHUNK=$((chunk_index * CHUNK_SIZE + 1))
END_CHUNK=$((START_CHUNK + CHUNK_SIZE - 1))

echo "Running $SCENARIO, chunks $START_CHUNK to $END_CHUNK"
Rscript sim_estim.R $SCENARIO $START_CHUNK $END_CHUNK
