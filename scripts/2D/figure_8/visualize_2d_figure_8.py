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

    RESULTS_BASE_DIR = os.path.join(PROJECT_ROOT, "results/2D_tests/figure_8")
    COMPUTE_SCRIPT_PATH = os.path.join(PROJECT_ROOT, "scripts/2D/compute_activation.py")
    PLOT_SCRIPT_PATH = os.path.join(PROJECT_ROOT, "scripts/2D/generate_single_plot.py")
    
    cases = ["baseline", "4x_d_ext", "8x_d_axn", "2x_alpha"]
    
    # --- Stage 1: Generate a single plot for each case using pvpython ---
    print("--- STAGE 1: Generating individual plots with pvpython ---")
    for case_name in cases:
        full_path = os.path.join(RESULTS_BASE_DIR, case_name)
        input_vtu_for_plotting = os.path.join(full_path, "activation_time.vtu")
        output_png = os.path.join(full_path, "activation_plot.png")
        
        if not os.path.isdir(full_path):
            print(f"Warning: Directory not found for case '{case_name}'. Skipping.")
            continue

        # First, run the compute script to generate the activation_time.vtu
        compute_command = f'bash -c "source /opt/miniforge3/bin/activate pv-env; pvpython {COMPUTE_SCRIPT_PATH} {full_path}"'
        print(f"\nExecuting pre-processing for {case_name}...")
        compute_result = subprocess.run(compute_command, shell=True, capture_output=True, text=True)
        if compute_result.returncode != 0:
            print(f"--- PVPYTHON (compute) FAILED for {case_name} ---"); print("STDERR:", compute_result.stderr)
            continue
        else:
            print(compute_result.stdout)

        # Then, run the plotting script on the file that was just created.
        print(f"Executing plotting for {case_name}...")
        plot_command = f'bash -c "source /opt/miniforge3/bin/activate pv-env; pvpython {PLOT_SCRIPT_PATH} {input_vtu_for_plotting} {output_png}"'
        plot_result = subprocess.run(plot_command, shell=True, capture_output=True, text=True)
        if plot_result.returncode != 0:
            print(f"--- PVPYTHON (plot) FAILED for {case_name} ---"); print("STDERR:", plot_result.stderr)
        else:
            print(plot_result.stdout)

    # --- Stage 2: Assemble the final figure from the generated plots ---
    print("\n--- STAGE 2: Assembling final 2x2 figure ---")
    
    fig, axes = plt.subplots(2, 2, figsize=(10, 10)) # Adjusted size for no colorbar
    fig.suptitle("2D Sensitivity Analysis (Replication of Figure 8)", fontsize=16)
    
    case_map = { "baseline": (0, 0), "4x_d_ext": (0, 1), "8x_d_axn": (1, 0), "2x_alpha": (1, 1) }
    titles = { "baseline": "baseline simulation", "4x_d_ext": "4-fold diffusion d_ext", "8x_d_axn": "8-fold axonal transport d_axn", "2x_alpha": "2-fold growth rate α" }

    for case_name, (r, c) in case_map.items():
        ax = axes[r, c]
        image_path = os.path.join(RESULTS_BASE_DIR, case_name, "activation_plot.png")
        
        if os.path.exists(image_path):
            img = mpimg.imread(image_path)
            ax.imshow(img)
        else:
            print(f"Warning: Image not found for case '{case_name}' at {image_path}")
            ax.text(0.5, 0.5, 'Image not found', ha='center', va='center')
        
        ax.set_title(titles[case_name])
        ax.axis('off')

    # Use tight_layout to automatically adjust spacing
    fig.tight_layout(rect=[0, 0, 1, 0.96]) # Adjust for suptitle

    output_filename = os.path.join(PROJECT_ROOT, "results/2D_tests/figure_8/figure8_2D_replication.png")
    plt.savefig(output_filename, dpi=300)
    print(f"\nFinal composite plot saved as {output_filename}")
    
    # plt.show() # Comment out for non-GUI environment