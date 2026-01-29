#!/bin/bash
#SBATCH --account=def-mlegault
#SBATCH --job-name=mr_all
#SBATCH --time=10:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=../results/logs/%x_%j.out

cd ~/projects/def-mlegault/gf591137/mr-quantile/scripts/
module load r

SCENARIO_NAME=$1
SCENARIO=$2
START_CHUNK=$3
END_CHUNK=$4
METHODS=$5

echo "Running $SCENARIO, chunks $START_CHUNK to $END_CHUNK, $METHODS methods"
Rscript sim_estim.R "$SCENARIO_NAME" "$SCENARIO" "$START_CHUNK" "$END_CHUNK" "$METHODS"
