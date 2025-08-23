import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import meshio

# --- Configuration ---
RESULTS_BASE_DIR = "experiment_1d_results"
TOTAL_TIME = 20.0
TIME_STEP = 0.1
ACTIVATION_THRESHOLD = 0.2



# --- Main function to process one simulation run ---
def process_one_simulation(directory_path):
    vtu_files = sorted(glob.glob(os.path.join(directory_path, "solution_[0-9][0-9][0-9].0.vtu")))
    if not vtu_files:
        print(f"Warning: No solution files found in {directory_path}")
        return None, None

    mesh = meshio.read(vtu_files[0])
    x_coords = mesh.points[:, 0]
    activation_times = np.full_like(x_coords, np.inf)

    for i, file_path in enumerate(vtu_files):
        current_time = i * TIME_STEP
        mesh = meshio.read(file_path)
        if "u" not in mesh.point_data:
            print(f"Warning: Data key 'u' not found in {file_path}. Skipping.")
            continue
        concentration = mesh.point_data["u"]
        newly_activated_mask = (concentration > ACTIVATION_THRESHOLD) & (np.isinf(activation_times))
        activation_times[newly_activated_mask] = current_time

    activation_times[np.isinf(activation_times)] = TOTAL_TIME
    return x_coords, activation_times

# --- Main plotting logic ---
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
            
            x, t_act = process_one_simulation(full_path)
            
            if x is None: continue

            sort_indices = np.argsort(x)
            x_sorted, t_act_sorted = x[sort_indices], t_act[sort_indices]
            
            # --- THE FIX: Transpose the Z array to match the shape of X and Y ---
            Y, X = np.meshgrid(np.linspace(0, 1, 2), x_sorted)
            Z = np.vstack([t_act_sorted, t_act_sorted]).T # <-- Added .T here

            levels = np.linspace(0, TOTAL_TIME, 12)
            im = ax.contourf(X, Y, Z, levels=levels, cmap='jet_r', vmin=0, vmax=TOTAL_TIME)
            
            if i == 0: ax.set_title(f"growth {alpha_mult}xα")
            if j == 2: ax.text(1.1, 0.5, f"spreading {d_mult}xd", rotation=-90,
                               verticalalignment='center', transform=ax.transAxes)

    for ax in axes[-1, :]: ax.set_xlabel("x")
    for ax in axes[:, 0]: ax.set_ylabel("concentration")
    plt.setp(axes, yticks=[0, 1], xticks=[-1, 0, 1])

    fig.subplots_adjust(right=0.85, wspace=0.1, hspace=0.1)
    if im:
        c_ax = fig.add_axes([0.88, 0.15, 0.03, 0.7])
        cbar = fig.colorbar(im, cax=c_ax)
        cbar.set_label("activation time (t)")
    
    output_filename = "figure3_replication.png"
    plt.savefig(output_filename, dpi=300)
    print(f"\nPlot saved as {output_filename}")
    
    plt.show()

# import os
# import glob
# import numpy as np
# import matplotlib.pyplot as plt
# import meshio

# # --- Configuration ---
# RESULTS_BASE_DIR = "experiment_1d_results"
# TOTAL_TIME = 20.0
# TIME_STEP = 0.1
# ACTIVATION_THRESHOLD = 0.2

# # --- Main function to process one simulation run ---
# def process_one_simulation(directory_path):
#     vtu_files = sorted(glob.glob(os.path.join(directory_path, "solution_[0-9][0-9][0-9].0.vtu")))
#     if not vtu_files:
#         print(f"Warning: No solution files found in {directory_path}")
#         return None, None, None

#     mesh = meshio.read(vtu_files[0])
#     x_coords = mesh.points[:, 0]
#     activation_times = np.full_like(x_coords, np.inf)

#     for i, file_path in enumerate(vtu_files):
#         current_time = i * TIME_STEP
#         mesh = meshio.read(file_path)
#         if "u" not in mesh.point_data:
#             print(f"Warning: Data key 'u' not found in {file_path}. Skipping.")
#             continue
#         concentration = mesh.point_data["u"]
#         newly_activated_mask = (concentration > ACTIVATION_THRESHOLD) & (np.isinf(activation_times))
#         activation_times[newly_activated_mask] = current_time

#     # We need the final concentration to define the shape of the plot
#     final_concentration = None
#     for file_path in reversed(vtu_files):
#         mesh = meshio.read(file_path)
#         if "u" in mesh.point_data:
#             final_concentration = mesh.point_data["u"]
#             break
            
#     # Handle non-activated areas by setting their time to the max time
#     activation_times[np.isinf(activation_times)] = TOTAL_TIME
            
#     return x_coords, activation_times, final_concentration

# # --- Main plotting logic ---
# if __name__ == "__main__":
#     fig, axes = plt.subplots(3, 3, figsize=(12, 10), sharex=True, sharey=True)
#     fig.suptitle("Replication of Fisher-Kolmogorov 1D Spreading", fontsize=16)

#     alpha_multipliers = [1, 2, 4]
#     d_multipliers = [1, 2, 4]
#     im = None

#     for i, d_mult in enumerate(d_multipliers):
#         for j, alpha_mult in enumerate(alpha_multipliers):
#             ax = axes[i, j]
#             dir_name = f"alpha_{alpha_mult}a_d_{d_mult}d"
#             full_path = os.path.join(RESULTS_BASE_DIR, dir_name)
            
#             x, t_act, c_final = process_one_simulation(full_path)
            
#             if x is None or c_final is None: continue

#             sort_indices = np.argsort(x)
#             x_sorted, t_act_sorted, c_final_sorted = x[sort_indices], t_act[sort_indices], c_final[sort_indices]

#             # --- THE FINAL FIXES ---
#             # 1. Create the 2D image of activation times
#             display_grid = np.tile(t_act_sorted, (100, 1))
            
#             # 2. Create the "wave shape" mask from the final concentration
#             y_grid = np.linspace(0, 1, 100)[:, np.newaxis]
#             mask = y_grid > c_final_sorted
#             display_grid[mask] = np.nan # Make areas outside the shape transparent
            
#             # 3. Create a discrete colormap with 11 distinct colors
#             cmap = plt.get_cmap('jet_r', 11)
            
#             im = ax.imshow(display_grid, extent=[-1, 1, 0, 1], origin='lower', 
#                            aspect='auto', cmap=cmap, vmin=0, vmax=TOTAL_TIME)
            
#             if i == 0: ax.set_title(f"growth {alpha_mult}xα")
#             if j == 2: ax.text(1.1, 0.5, f"spreading {d_mult}xd", rotation=-90,
#                                verticalalignment='center', transform=ax.transAxes)

#     for ax in axes[-1, :]: ax.set_xlabel("x")
#     for ax in axes[:, 0]: ax.set_ylabel("concentration")
#     plt.setp(axes, yticks=[0, 1], xticks=[-1, 0, 1])

#     fig.subplots_adjust(right=0.85, wspace=0.1, hspace=0.1)
#     if im:
#         c_ax = fig.add_axes([0.88, 0.15, 0.03, 0.7])
#         cbar = fig.colorbar(im, cax=c_ax)
#         cbar.set_label("activation time (t)")
    
#     output_filename = "figure3_replication.png"
#     plt.savefig(output_filename, dpi=300)
#     print(f"\nPlot saved as {output_filename}")
    
#     plt.show()