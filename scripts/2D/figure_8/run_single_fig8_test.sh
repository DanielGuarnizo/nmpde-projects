#!/bin/bash

# This script runs a single, specific test case from the Figure 8 sensitivity analysis:
# The "baseline" case.

echo "--- Running Single Test: Figure 8 (baseline) ---"

# --- Step 1: Load the necessary environment for the executable to run ---
echo "Loading environment modules..."
module load gcc-glibc dealii

# --- Step 2: Define all parameters for this specific test case ---

# Paths
EXECUTABLE="./build/neuro_disease_2D"
# MESH_FILE="meshes/slice_generated.msh"
# --- THE FIX: Point to the new, robustly created .ucd mesh file ---
MESH_FILE="meshes/slice_with_hole.msh"

# Simulation control parameters
TOTAL_TIME=48.0
TIME_STEP=0.24
POLY_DEGREE=1
C_0_VALUE=0.2
GRAY_MATTER_THRESH=4.0

# Seeding and Fiber types for the Figure 8 experiment
SEEDING_TYPE=2  # Corresponds to SeedingRegionType::Tau
FIBER_TYPE=2    # Corresponds to FiberFieldType::AxonBased

# --- THE CHANGE: Physical parameters for the "baseline" case ---
# Set all values to their baseline (1x) level.
ALPHA_VALUE=0.6
D_EXT_VALUE=1.5
D_AXN_VALUE=3.0

# --- Step 3: Define the output location ---
# --- THE CHANGE: Point to a new, specific output directory ---
OUTPUT_DIR="results/2D_tests/single_test_baseline/"
FILENAME_PREFIX="solution"

# Create the output directory, making sure it's clean first
rm -rf $OUTPUT_DIR
mkdir -p $OUTPUT_DIR

# Check if the executable exists before running
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable not found at $EXECUTABLE"
    echo "Please compile the project first by running ./build.sh"
    exit 1
fi

# --- Step 4: Run the simulation ---
echo "Executing simulation with the following parameters:"
echo "  alpha = $ALPHA_VALUE"
echo "  d_ext = $D_EXT_VALUE"
echo "  d_axn = $D_AXN_VALUE"
echo "Output will be in: $OUTPUT_DIR"
echo "---------------------------------------------------------"

# Get the number of available processors for MPI
n=$(nproc)

mpirun -n $n $EXECUTABLE \
  -m $MESH_FILE \
  -T $TOTAL_TIME \
  -t $TIME_STEP \
  -g $POLY_DEGREE \
  -a $ALPHA_VALUE \
  -x $D_AXN_VALUE \
  -e $D_EXT_VALUE \
  -c $C_0_VALUE \
  -H $GRAY_MATTER_THRESH \
  -s $SEEDING_TYPE \
  -f $FIBER_TYPE \
  -d $OUTPUT_DIR \
  -o $FILENAME_PREFIX

echo "--- Simulation Finished ---"
echo "Output files are in the '${OUTPUT_DIR}' directory."
echo "You can open the output files in ParaView to see the results."