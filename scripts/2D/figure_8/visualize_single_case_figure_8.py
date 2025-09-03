import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import meshio
from stitch_parallel_data import stitch_data # Import our helper function
from matplotlib.patches import Polygon # Import the Polygon tool

# --- Configuration ---
TOTAL_TIME = 48.0
TIME_STEP = 0.24
ACTIVATION_THRESHOLD = 0.25 

# --- Helper function to find the boundary loops of a mesh ---
def find_boundary_loops(mesh):
    print("Finding all boundary loops in the stitched mesh...")
    if 'triangle' not in mesh.cells_dict:
        raise ValueError("Mesh must contain triangles.")
    triangles = mesh.cells_dict['triangle']
    edge_counts = {}
    for tri in triangles:
        for i in range(3):
            edge = tuple(sorted((tri[i], tri[(i + 1) % 3])))
            edge_counts[edge] = edge_counts.get(edge, 0) + 1
    boundary_edges = {edge for edge, count in edge_counts.items() if count == 1}
    if not boundary_edges:
        raise ValueError("No boundary edges found.")
    
    adj = {}
    for v1, v2 in boundary_edges:
        if v1 not in adj: adj[v1] = []
        if v2 not in adj: adj[v2] = []
        adj[v1].append(v2)
        adj[v2].append(v1)

    all_loops, visited_nodes = [], set()
    for start_node in adj:
        if start_node not in visited_nodes:
            path, cn, pn = [start_node], start_node, -1
            while True:
                next_nodes = [n for n in adj[cn] if n != pn]
                if not next_nodes or (len(path) > 2 and path[0] in next_nodes):
                    path.append(path[0])
                    break
                next_n = -1
                for node_opt in next_nodes:
                    if node_opt not in path:
                        next_n = node_opt
                        break
                if next_n == -1: break
                path.append(next_n)
                pn, cn = cn, next_n
            
            loop_nodes = path[:-1]
            if len(loop_nodes) > 2:
                all_loops.append(loop_nodes)
                visited_nodes.update(loop_nodes)

    all_loops.sort(key=len, reverse=True)
    return all_loops


# --- Main function to process one simulation run ---
def process_simulation_data(directory_path):
    pvtu_files = sorted(glob.glob(os.path.join(directory_path, "solution_[0-9][0-9][0-9].pvtu")))
    if not pvtu_files:
        print(f"Error: No solution*.pvtu files found in {directory_path}")
        return None, None, None
        
    print(f"Found {len(pvtu_files)} timesteps to process in {os.path.basename(directory_path)}...")
    first_mesh = stitch_data(directory_path, 0)
    if first_mesh is None: return None, None, None
    
    points = first_mesh.points
    if 'triangle' in first_mesh.cells_dict:
        triangles = first_mesh.cells_dict['triangle']
    else:
        return None, None, None

    activation_times = np.full(len(points), np.inf)
    for t_idx in range(len(pvtu_files)):
        current_time = t_idx * TIME_STEP
        stitched_mesh = stitch_data(directory_path, t_idx)
        if stitched_mesh is None or "u" not in stitched_mesh.point_data:
            continue
        concentration = stitched_mesh.point_data["u"]
        newly_activated_mask = (concentration > ACTIVATION_THRESHOLD) & (np.isinf(activation_times))
        activation_times[newly_activated_mask] = current_time

    activation_times[np.isinf(activation_times)] = TOTAL_TIME
    return points, triangles, activation_times
    
# --- Main plotting logic ---
if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(1)
        
    result_dir = sys.argv[1]
    if not os.path.isdir(result_dir):
        sys.exit(1)

    points, triangles, t_act = process_simulation_data(result_dir)
    
    if points is None:
        sys.exit(1)

    fig, ax = plt.subplots(figsize=(10, 8))
    case_name = os.path.basename(os.path.normpath(result_dir))
    
    # --- THE FIX: Add alpha=1.0 for solid, vibrant colors ---
    im = ax.tricontourf(points[:, 0], points[:, 1], triangles, t_act, 
                       levels=12, cmap='jet_r', vmin=0, vmax=TOTAL_TIME,
                       alpha=1.0, antialiased=True)
    
    stitched_mesh_for_holes = meshio.Mesh(points, [("triangle", triangles)])
    boundary_loops = find_boundary_loops(stitched_mesh_for_holes)
    
    if len(boundary_loops) > 1:
        for hole_loop_indices in boundary_loops[1:]:
            hole_verts_3d = points[hole_loop_indices]
            hole_verts_2d = hole_verts_3d[:, :2]
            
            # Draw the white hole on top of the colors
            hole_patch = Polygon(hole_verts_2d, closed=True, facecolor='white', edgecolor='none', zorder=2)
            ax.add_patch(hole_patch)
            
    ax.set_title(f"Activation Time for: {case_name}", fontsize=16)
    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("activation time (t)")
    cbar.set_ticks([0, TOTAL_TIME])
    cbar.set_ticklabels([r'$t_0$', r'$t_{max}$'])
    
    output_filename = os.path.join(result_dir, f"{case_name}_activation_plot.png")
    plt.savefig(output_filename, dpi=300, bbox_inches='tight', pad_inches=0.1)
    
    print(f"\nPlot saved as {output_filename}")
    plt.show()