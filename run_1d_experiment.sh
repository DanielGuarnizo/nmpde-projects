# #!/bin/bash

# # This script replicates the 1D experiment from Weickenmeier et al. (2019)
# # using the corrected, more dynamic base parameters from the successful single test.

# echo "Starting 1D parameter sweep..."

# # --- Load the necessary environment ---
# echo "Loading environment modules..."
# module load gcc-glibc dealii

# # --- Define paths and parameters ---
# EXECUTABLE="./build/neuro_disease_1D"
# MESH_FILE="meshes/mesh-1D-centered.msh"
# TOTAL_TIME=20.0
# TIME_STEP=0.1
# POLY_DEGREE=1

# # --- THE FIX: Use the parameters from your successful 'run_single_test.sh' as the base ---
# BASE_ALPHA=1.8
# BASE_D_AXN=0.2
# C_0_VALUE=0.4

# # --- Create and clean the results directory ---
# OUTPUT_BASE_DIR="experiment_1d_results"
# rm -rf $OUTPUT_BASE_DIR
# mkdir -p $OUTPUT_BASE_DIR

# # Check if the executable exists
# if [ ! -f "$EXECUTABLE" ]; then
#     echo "Error: Executable not found at $EXECUTABLE"
#     exit 1
# fi

# # --- Main experiment loop ---
# for alpha_multiplier in 1 2 4
# do
#   for d_multiplier in 1 2 4
#   do
#     # Note: bc cannot handle multipliers like '1'. It must be '1.0' for float math.
#     # We will use awk for more robust floating point multiplication.
#     CURRENT_ALPHA=$(awk -v base="$BASE_ALPHA" -v mult="$alpha_multiplier" 'BEGIN{print base * mult}')
#     CURRENT_D_AXN=$(awk -v base="$BASE_D_AXN" -v mult="$d_multiplier" 'BEGIN{print base * mult}')

#     OUTPUT_DIR="${OUTPUT_BASE_DIR}/alpha_${alpha_multiplier}a_d_${d_multiplier}d/"
#     FILENAME_PREFIX="solution"
#     mkdir -p $OUTPUT_DIR

#     echo "---------------------------------------------------------"
#     echo "RUNNING: alpha = $CURRENT_ALPHA, d_axn = $CURRENT_D_AXN"
#     echo "Output will be in: $OUTPUT_DIR"
#     echo "---------------------------------------------------------"

#     $EXECUTABLE \
#       -m $MESH_FILE \
#       -T $TOTAL_TIME \
#       -t $TIME_STEP \
#       -g $POLY_DEGREE \
#       -a $CURRENT_ALPHA \
#       -x $CURRENT_D_AXN \
#       -e 0.0 \
#       -c $C_0_VALUE \
#       -d $OUTPUT_DIR \
#       -o $FILENAME_PREFIX
#   done
# done

# echo "Experiment finished. Results are in the '${OUTPUT_BASE_DIR}' directory."

#!/bin/bash

# This script replicates the 1D experiment from Weickenmeier et al. (2019)
# using parameters designed to show a traveling wave of activation.

echo "Starting 1D parameter sweep..."
echo "Loading environment modules..."
module load gcc-glibc dealii

EXECUTABLE="./build/neuro_disease_1D"
MESH_FILE="meshes/mesh-1D-centered.msh"
TOTAL_TIME=20.0
TIME_STEP=0.1
POLY_DEGREE=1

# --- THE FIX: Use the original paper's "slow and small" parameters ---
BASE_ALPHA=1.0
BASE_D_AXN=0.0001 
C_0_VALUE=0.1      # A small initial seed, BELOW the activation threshold

OUTPUT_BASE_DIR="experiment_1d_results"
rm -rf $OUTPUT_BASE_DIR
mkdir -p $OUTPUT_BASE_DIR

if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable not found at $EXECUTABLE"
    exit 1
fi

for alpha_multiplier in 1 2 4
do
  for d_multiplier in 1 2 4
  do
    CURRENT_ALPHA=$(awk -v base="$BASE_ALPHA" -v mult="$alpha_multiplier" 'BEGIN{print base * mult}')
    CURRENT_D_AXN=$(awk -v base="$BASE_D_AXN" -v mult="$d_multiplier" 'BEGIN{print base * mult}')
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
      -e 0.0 \
      -c $C_0_VALUE \
      -d $OUTPUT_DIR \
      -o $FILENAME_PREFIX
  done
done

echo "Experiment finished. Results are in the '${OUTPUT_BASE_DIR}' directory."