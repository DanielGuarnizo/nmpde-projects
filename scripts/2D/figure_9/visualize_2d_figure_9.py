import os
import subprocess
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as mcolors
import numpy as np
import sys

# This script orchestrates the visualization process and assembles the final plot.

if __name__ == "__main__":
    # --- Configuration ---
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..', '..'))

    RESULTS_BASE_DIR = os.path.join(PROJECT_ROOT, "results/2D_tests/figure_9")
    COMPUTE_SCRIPT_PATH = os.path.join(PROJECT_ROOT, "scripts/2D/compute_activation.py")
    PLOT_SCRIPT_PATH = os.path.join(PROJECT_ROOT, "scripts/2D/generate_single_plot.py")
    
    diseases = ["alpha_synuclein", "amyloid_beta", "tau", "tdp43"]
    fiber_models = ["isotropic", "circumferential", "radial", "axon_based"]
    
    # --- Stage 1 & 2: Pre-process data and generate a single plot for each case ---
    print("--- STAGE 1 & 2: Generating individual plots for all 16 cases ---")
    for disease_name in diseases:
        for fiber_model_name in fiber_models:
            case_name = f"{disease_name} with {fiber_model_name}"
            full_path = os.path.join(RESULTS_BASE_DIR, disease_name, fiber_model_name)
            
            input_vtu_for_plotting = os.path.join(full_path, "activation_time.vtu")
            output_png = os.path.join(full_path, "activation_plot.png")

            print(f"\n--- Checking case: {case_name} ---")
            
            if os.path.exists(output_png):
                print(f"Final plot '{os.path.basename(output_png)}' already exists. Skipping.")
                continue

            if not os.path.isdir(full_path):
                print(f"Warning: Directory not found. Skipping.")
                continue

            if os.path.exists(input_vtu_for_plotting):
                print(f"Found existing '{os.path.basename(input_vtu_for_plotting)}'. Skipping computation.")
            else:
                print(f"File not found. Running pre-processing...")
                compute_command = f'bash -c "source /opt/miniforge3/bin/activate pv-env; pvpython {COMPUTE_SCRIPT_PATH} ."'
                compute_result = subprocess.run(compute_command, shell=True, capture_output=True, text=True, cwd=full_path)
                if compute_result.returncode != 0:
                    print(f"--- PVPYTHON (compute) FAILED for {case_name} ---"); print("STDERR:", compute_result.stderr)
                    continue
                else:
                    print(compute_result.stdout)
            
            print(f"Generating final plot...")
            plot_command = f'bash -c "source /opt/miniforge3/bin/activate pv-env; pvpython {PLOT_SCRIPT_PATH} activation_time.vtu activation_plot.png"'
            plot_result = subprocess.run(plot_command, shell=True, capture_output=True, text=True, cwd=full_path)
            if plot_result.returncode != 0:
                print(f"--- PVPYTHON (plot) FAILED for {case_name} ---"); print("STDERR:", plot_result.stderr)
            else:
                print(plot_result.stdout)

    # --- Stage 3: Assemble the final 4x4 figure ---
    print("\n--- STAGE 3: Assembling final 4x4 figure ---")
    
    fig, axes = plt.subplots(4, 4, figsize=(12, 12), sharex=True, sharey=True) # Adjusted size
    fig.suptitle("2D Fiber & Seeding Sensitivity Analysis (Replication of Figure 9)", fontsize=16)
    
    for i, disease_name in enumerate(diseases):
        for j, fiber_model_name in enumerate(fiber_models):
            ax = axes[i, j]
            image_path = os.path.join(RESULTS_BASE_DIR, disease_name, fiber_model_name, "activation_plot.png")
            
            if os.path.exists(image_path):
                img = mpimg.imread(image_path)
                ax.imshow(img)
            else:
                ax.text(0.5, 0.5, 'Image\nnot found', ha='center', va='center', fontsize=9)
            
            ax.axis('off')
            
            # Set titles for the first row (the fiber models)
            if i == 0:
                ax.set_title(fiber_model_name.replace('_', ' ').capitalize(), fontsize=12)

            # --- THE FIX #1: Add disease names as row labels on the left ---
            if j == 0:
                # Use a text object for better placement control than ylabel
                fig.text(0.05, 0.77 - i*0.2, disease_name.replace('_', ' ').capitalize(), 
                         ha='center', va='center', rotation='vertical', fontsize=14)

    # --- THE FIX #2: Remove the colorbar ---
    # The block for creating the colorbar has been completely removed.
    
    # Adjust layout to prevent title overlap
    plt.tight_layout(rect=[0.05, 0, 1, 0.96])

    output_filename = os.path.join(RESULTS_BASE_DIR, "figure9_2D_replication.png")
    plt.savefig(output_filename, dpi=300)
    print(f"\nFinal composite plot saved as {output_filename}")
    
    # plt.show()