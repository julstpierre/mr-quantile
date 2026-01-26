#!/bin/bash
#SBATCH --account=def-mlegault
#SBATCH --job-name=mr_main_fast
#SBATCH --array=0-489
#SBATCH --time=05:00:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH --output=../results/logs/%x_%A_%a.out

cd ~/projects/def-mlegault/gf591137/mr-quantile/scripts/
module load r

# -----------------------------
# Scenario definitions
# -----------------------------
SCENARIOS=($(seq -f "D_S%g" 1 42))
THETA0_SCENARIOS=(D_S3 D_S4 D_S9 D_S10 D_S15 D_S16 D_S21 D_S22 D_S27 D_S28 D_S33 D_S34 D_S39 D_S40)

CHUNK_SIZE=400
TASK_ID=$SLURM_ARRAY_TASK_ID
TASK_OFFSET=0

# -----------------------------
# Resolve scenario + chunk
# -----------------------------
for SCEN in "${SCENARIOS[@]}"; do

  if [[ " ${THETA0_SCENARIOS[@]} " =~ " ${SCEN} " ]]; then
    N_ITER=10000
  else
    N_ITER=2000
  fi

  N_CHUNKS=$((N_ITER / CHUNK_SIZE))

  if (( TASK_ID < TASK_OFFSET + N_CHUNKS )); then
    SCENARIO=$SCEN
    chunk_index=$((TASK_ID - TASK_OFFSET))
    break
  fi

  TASK_OFFSET=$((TASK_OFFSET + N_CHUNKS))
done

# -----------------------------
# Chunk bounds
# -----------------------------
START_CHUNK=$((chunk_index * CHUNK_SIZE + 1))
END_CHUNK=$((START_CHUNK + CHUNK_SIZE - 1))

# -----------------------------
# Run
# -----------------------------
echo "Running $SCENARIO, chunks $START_CHUNK to $END_CHUNK, fast methods"

Rscript sim_estim.R "main" $SCENARIO $START_CHUNK $END_CHUNK "fast"
