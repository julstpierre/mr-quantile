#!/bin/bash

# List of scenarios
SCENARIOS=($(seq -f "D_S%g" 1 168))

CHUNK_SIZE=400

# Navigate to the scripts folder
cd ~/Dropbox/Postdoc/Projects/MRFuse/scripts

# Function to run a single chunk
run_chunk() {
    local scenario=$1
    local chunk_index=$2
    local start=$((chunk_index * CHUNK_SIZE + 1))
    local end=$((start + CHUNK_SIZE - 1))
    echo "Running $scenario, chunks $start to $end"
    Rscript sim_estim.R "$scenario" "$start" "$end"
}

# Export the function so parallel can see it
export -f run_chunk
export CHUNK_SIZE

# Generate a list of all scenario + chunk combinations and run in parallel
parallel -j 5 --line-buffer run_chunk {1} {2} ::: "${SCENARIOS[@]}" ::: 0 1 2 3 4