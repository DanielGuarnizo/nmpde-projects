import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import meshio

# --- Configuration ---
RESULTS_BASE_DIR = "experiment_1d_results"
TOTAL_TIME = 20.0
TIME_STEP = 0.1
ACTIVATION_THRESHOLD = 0.5

# --- Main function to process one simulation run ---
def process_one_simulation(directory_path):
    """
    Reads all .pvtu files in a directory using the correct format hint,
    calculates activation times, and returns data for plotting.
    """
    # Look for the master .pvtu files for each time step
    pvtu_files = sorted(glob.glob(os.path.join(directory_path, "*.pvtu")))
    
    if not pvtu_files:
        print(f"Warning: No .pvtu files found in {directory_path}")
        return None, None, None

    # --- THE DEFINITIVE FIX: Use the specific "pvtu" file format hint ---
    # This tells meshio to use its parallel VTK reader.
    mesh = meshio.read(pvtu_files[0], file_format="pvtu")
    x_coords = mesh.points[:, 0]
    
    activation_times = np.full_like(x_coords, np.inf)

    # Loop through each time step (each .pvtu file)
    for i, file_path in enumerate(pvtu_files):
        current_time = i * TIME_STEP
        
        # Also specify the file format here
        mesh = meshio.read(file_path, file_format="pvtu")
        concentration = mesh.point_data["u"]
        
        newly_activated_mask = (concentration > ACTIVATION_THRESHOLD) & (activation_times == np.inf)
        activation_times[newly_activated_mask] = current_time

    # Read the final concentration profile
    final_mesh = meshio.read(pvtu_files[-1], file_format="pvtu")
    final_concentration = final_mesh.point_data["u"]

    return x_coords, activation_times, final_concentration

# --- Main plotting logic (does not need to be changed) ---
if __name__ == "__main__":
    fig, axes = plt.subplots(3, 3, figsize=(12, 10), sharex=True, sharey=True)
    fig.suptitle("Replication of Fisher-Kolmogorov 1D Spreading", fontsize=16)

    alpha_multipliers = [1, 2, 4]
    d_multipliers = [1, 2, 4]
    
    im = None

    for i, d_mult in enumerate(d_multipliers):
        for j, alpha_mult in enumerate(alpha_multipliers):
            ax = axes[i, j]
            dir_name = f"alpha_{alpha_mult}a_d_{d_mult}d"
            full_path = os.path.join(RESULTS_BASE_DIR, dir_name)
            
            x, t_act, c_final = process_one_simulation(full_path)
            
            if x is None: continue

            sort_indices = np.argsort(x)
            x_sorted, t_act_sorted, c_final_sorted = x[sort_indices], t_act[sort_indices], c_final[sort_indices]

            display_grid = np.tile(t_act_sorted, (100, 1))
            y_grid = np.linspace(0, 1, 100)[:, np.newaxis]
            mask = y_grid > c_final_sorted
            display_grid[mask] = np.nan
            
            im = ax.imshow(display_grid, extent=[-1, 1, 0, 1], origin='lower', 
                           aspect='auto', cmap='jet', vmin=0, vmax=TOTAL_TIME)
            
            if i == 0: ax.set_title(f"growth {alpha_mult}xα")
            if j == 2: ax.text(1.1, 0.5, f"spreading {d_mult}xd", rotation=-90,
                               verticalalignment='center', transform=ax.transAxes)

    for ax in axes[-1, :]: ax.set_xlabel("x")
    for ax in axes[:, 0]: ax.set_ylabel("concentration")
    plt.setp(axes, yticks=[0, 1], xticks=[-1, 0, 1])

    fig.subplots_adjust(right=0.85, wspace=0.1, hspace=0.1)
    
    if im:
        cbar_ax = fig.add_axes([0.88, 0.15, 0.03, 0.7])
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.set_label("activation time (t)")
    
    output_filename = "figure3_replication.png"
    plt.savefig(output_filename, dpi=300)
    print(f"\nPlot saved as {output_filename}")
    
    plt.show()