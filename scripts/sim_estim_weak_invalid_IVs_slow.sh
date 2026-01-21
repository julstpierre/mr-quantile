#!/bin/bash
#SBATCH --account=def-mlegault
#SBATCH --job-name=mr_weak_invalid_IVs_slow
#SBATCH --array=0-439
#SBATCH --time=10:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=../results/logs/%x_%A_%a.out

cd ~/projects/def-mlegault/gf591137/MRFuse/scripts/
module load r

# -----------------------------
# Scenario definitions
# -----------------------------
SCENARIOS=($(seq -f "D_S%g" 1 56))
THETA0_SCENARIOS=(D_S4 D_S11 D_S18 D_S25 D_S32 D_S39 D_S46 D_S53)

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
echo "Running $SCENARIO, chunks $START_CHUNK to $END_CHUNK, slow methods"

Rscript sim_estim.R "weak_invalid_IVs" $SCENARIO $START_CHUNK $END_CHUNK "slow"
