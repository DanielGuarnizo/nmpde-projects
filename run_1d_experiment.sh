#!/bin/bash

# This script replicates the 1D experiment from Weickenmeier et al. (2019)
# It is designed to be run from the main project directory, one level above 'build'.

echo "Starting 1D parameter sweep..."

# --- Define the path to the executable ---
EXECUTABLE="./build/neuro_disease_1D"

# --- Define the path to the mesh file we generated ---
MESH_FILE="meshes/mesh-1D-centered.msh"

# --- Define fixed simulation parameters ---
TOTAL_TIME=20.0
TIME_STEP=0.1
POLY_DEGREE=1
C_0_VALUE=0.4 # The initial peak concentration

# --- Define the parameters to test via sweeping ---
BASE_ALPHA=1.0
BASE_D_AXN=0.0001 # This is the 'd' from the paper's Fig 3
D_EXT=0.0

# --- Create a main directory for the results ---
OUTPUT_BASE_DIR="experiment_1d_results"
mkdir -p $OUTPUT_BASE_DIR

# Check if the executable exists before starting
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable not found at $EXECUTABLE"
    echo "Please compile the project first by running your build script."
    exit 1
fi

# --- Main experiment loop ---
for alpha_multiplier in 1 2 4
do
  for d_multiplier in 1 2 4
  do
    CURRENT_ALPHA=$(echo "$BASE_ALPHA * $alpha_multiplier" | bc)
    CURRENT_D_AXN=$(echo "$BASE_D_AXN * $d_multiplier" | bc)
    OUTPUT_DIR="${OUTPUT_BASE_DIR}/alpha_${alpha_multiplier}a_d_${d_multiplier}d/"
    FILENAME_PREFIX="solution"
    mkdir -p $OUTPUT_DIR

    echo "---------------------------------------------------------"
    echo "RUNNING: alpha = $CURRENT_ALPHA, d_axn = $CURRENT_D_AXN"
    echo "Output will be in: $OUTPUT_DIR"
    echo "---------------------------------------------------------"

    $EXECUTABLE \
      -m $MESH_FILE \
      -T $TOTAL_TIME \
      -t $TIME_STEP \
      -g $POLY_DEGREE \
      -a $CURRENT_ALPHA \
      -x $CURRENT_D_AXN \
      -e $D_EXT \
      -d $OUTPUT_DIR \
      -o $FILENAME_PREFIX
  done
done

echo "Experiment finished. Results are in the '${OUTPUT_BASE_DIR}' directory."