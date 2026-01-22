#!/bin/bash
#SBATCH --account=def-mlegault
#SBATCH --job-name=mr_all
#SBATCH --time=05:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=../results/logs/%x_%j.out

cd ~/projects/def-mlegault/gf591137/mr-quantile/scripts/
module load r

SCENARIO=$1
START_CHUNK=$2
END_CHUNK=$3

echo "Running $SCENARIO, chunks $START_CHUNK to $END_CHUNK"
Rscript sim_estim.R "weak_invalid_IVs" "$SCENARIO" "$START_CHUNK" "$END_CHUNK" "slow"
