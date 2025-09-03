#!/bin/bash

# This script runs the sensitivity analysis for extracellular diffusion (d_ext),
# axonal transport (d_axn), and growth rate (alpha), to replicate Figure 8.

echo "--- Starting Sensitivity Analysis for Figure 8 ---"
echo "Loading environment modules..."
module load gcc-glibc dealii

# --- Common Parameters for all simulations in this experiment ---
EXECUTABLE="./build/neuro_disease_2D"
MESH_FILE="meshes/new_slice_generated.msh"

# --- THE FIX: Add the dimension parameter ---
DIMENSION=2

TOTAL_TIME=48.0
TIME_STEP=0.24
POLY_DEGREE=1
C_0_VALUE=0.2
GRAY_MATTER_THRESH=4.0
SEEDING_TYPE=2
FIBER_TYPE=2

# --- Baseline parameters from the paper ---
BASELINE_DEXT=1.5
BASELINE_DAXN=3.0
BASELINE_ALPHA=0.6

# --- Base output directory ---
OUTPUT_BASE_DIR="results/2D_tests/figure_8"
mkdir -p $OUTPUT_BASE_DIR

# --- Reusable run function ---
run_simulation() {
    local dext=$1
    local daxn=$2
    local alpha=$3
    local output_sub_dir=$4
    local full_output_dir="${OUTPUT_BASE_DIR}/${output_sub_dir}/"
    mkdir -p $full_output_dir

    echo "---------------------------------------------------------"
    echo "RUNNING CASE: ${output_sub_dir}"
    echo "Parameters: d_ext=${dext}, d_axn=${daxn}, alpha=${alpha}"
    echo "Output will be in: ${full_output_dir}"
    echo "---------------------------------------------------------"

    n=$(nproc)
    # --- THE FIX: Add the '-D' flag for the dimension ---
    mpirun -n $n $EXECUTABLE \
        -m $MESH_FILE \
        -D $DIMENSION \
        -T $TOTAL_TIME \
        -t $TIME_STEP \
        -g $POLY_DEGREE \
        -e $dext \
        -x $daxn \
        -a $alpha \
        -c $C_0_VALUE \
        -H $GRAY_MATTER_THRESH \
        -s $SEEDING_TYPE \
        -f $FIBER_TYPE \
        -d $full_output_dir \
        -o "solution"
}

# --- Run the 4 simulations for Figure 8 ---
run_simulation $BASELINE_DEXT $BASELINE_DAXN $BASELINE_ALPHA "baseline"
run_simulation 6.0 $BASELINE_DAXN $BASELINE_ALPHA "4x_d_ext"
run_simulation $BASELINE_DEXT 24.0 $BASELINE_ALPHA "8x_d_axn"
run_simulation $BASELINE_DEXT $BASELINE_DAXN 1.2 "2x_alpha"

echo "--- Figure 8 Sensitivity Analysis Finished ---"