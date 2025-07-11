#!/bin/bash
#SBATCH --job-name=mr_D_S2_null_chunk2
#SBATCH --time=04:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/%x_%j.out

cd ~/projects/def-mlegault/gf591137/BayesFuseMR/scripts/
module load r

SCENARIO="D_S2_null"
START_CHUNK=201
END_CHUNK=400

echo "Running $SCENARIO, chunks $START_CHUNK to $END_CHUNK"
Rscript sim_estim.R $SCENARIO $START_CHUNK $END_CHUNK
