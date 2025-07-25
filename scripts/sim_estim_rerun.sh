#!/bin/bash
#SBATCH --job-name=mr_D_S5
#SBATCH --time=04:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=../results/logs/%x_%j.out

cd ~/projects/def-mlegault/gf591137/BayesFuseMR/scripts/
module load r

SCENARIO="D_S5"
START_CHUNK=1
END_CHUNK=1000

echo "Running $SCENARIO, chunks $START_CHUNK to $END_CHUNK"
Rscript sim_estim.R $SCENARIO $START_CHUNK $END_CHUNK
