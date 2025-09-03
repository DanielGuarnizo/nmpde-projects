
---




# Neurodegenerative Disease Simulation

This project provides a C++ implementation of the reaction-diffusion model for neurodegenerative diseases presented in Weickenmeier et al. (2019), "A physics-based model explains the prion-like features of neurodegeneration...".

The code is designed to be compiled and run in a Linux-based environment (e.g., a Docker container) with the deal.II finite element library. Visualization is handled by a set of Python 3 scripts, which require a Conda environment with ParaView.

### 1. Project Structure

-   `/src`: Contains all C++ source code.
-   `/include`: Contains all C++ header files.
-   `/meshes`: Contains all input mesh files.
-   `/scripts`: Contains all bash and Python scripts for running experiments and visualizations.
-   `/results`: The default location where simulation and visualization outputs are saved.

### 2. Compiling the Code

All compilation should be done inside your configured Linux/Docker environment.

1.  **Load Dependencies:** Before compiling, ensure your environment is set up.
    ```bash
    module load gcc-glibc dealii
    ```
2.  **Run the Build Script:** From the main project directory (`nmpde-neurodegenerative-diseases`), run the provided build script. It will create a `build` directory, configure the project with CMake, and compile all C++ executables.
    ```bash
    ./build.sh
    ```
    This creates the executables `neuro_disease_1D`, `neuro_disease_2D`, etc., inside the `build` folder.

### 3. Running the 1D Experiment (Figure 3 Replication)

This experiment runs a 3x3 parameter sweep to test the effects of growth (`alpha`) and spreading (`d`) on a 1D domain.

#### **Step 1: Run the C++ Simulations**

From the main project directory, execute the experiment script. This will run 9 simulations and save the raw data.

```bash
./scripts/1D/run_1d_experiment.sh
```
-   **Output Location:** `/results/1D_tests/figure_3/`

#### **Step 2: Visualize the 1D Results**

This project includes two different Python scripts to visualize the 1D results.

*   **To create the "Activation Time" plot (like the paper's Figure 3):**
    ```bash
    python3 ./scripts/1D/visualize_1d_experiment.py
    ```
    -   **Output:** `figure3_replication.png` saved in the results directory.

*   **To create the "Space-Time" plot (raw concentration):**
    ```bash
    python3 ./scripts/1D/visualize_spacetime_figure_3.py
    ```
    -   **Output:** `figure3_spacetime_replication.png` saved in the results directory.

### 4. Running the 2D Experiments (Figure 8 & 9)

These experiments run sensitivity analyses on a 2D brain slice mesh. The visualization process requires ParaView.

#### **Step 1: Install ParaView via Conda (One-Time Setup)**

The visualization scripts rely on `pvpython`. The recommended way to install this in the container is with Conda.

1.  **Install Conda/Miniforge:**
    ```bash
    wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
    bash Miniforge3-Linux-x86_64.sh -b -p /opt/miniforge3
    source /opt/miniforge3/bin/activate
    ```
2.  **Create and Activate a ParaView Environment:**
    ```bash
    conda create -n pv-env -c conda-forge paraview -y
    source /opt/miniforge3/bin/activate pv-env
    ```
3.  **Install Python Libraries:** Inside the active environment, install the necessary libraries for plotting.
    ```bash
    pip install matplotlib meshio numpy
    ```
*Note: You must activate the Conda environment (`source /opt/miniforge3/bin/activate pv-env`) in every new terminal session before running the visualization scripts.*


---

#### **4.1 Figure 8: Parameter Sensitivity Analysis**

This runs a 2x2 analysis on the `baseline`, `4x d_ext`, `8x d_axn`, and `2x alpha` cases.

1.  **Run the C++ Simulations:** From the project root, run:
    ```bash
    ./scripts/2D/figure_8/ext_axn_alpha_sensitivity_analysis_figure_8.sh
    ```
    -   **Output Location:** `/results/2D_tests/figure_8/`


2.  **Visualize the Results:** From the project root, activate your Conda environment and run:
    ```bash
    python3 ./scripts/2D/figure_8/visualize_2d_figure_8.py
    ```
    -   **Output:** `figure8_2D_replication.png` saved in the main project directory.

---

#### **4.2 Figure 9: Disease and Fiber Model Analysis**

This runs a full 4x4 analysis on 4 disease types and 4 fiber models.

1.  **Run the C++ Simulations:** From the project root, run:
    ```bash
    ./scripts/2D/figure_9/fiber_seed_sensitivity_analysis_figure_9.sh
    ```
    -   **Output Location:** `/results/2D_tests/figure_9/`

2.  **Visualize the Results:** From the project root, activate your Conda environment and run:
    ```bash
    python3 ./scripts/2D/figure_9/visualize_2d_figure9.py
    ```
    -   **Output:** `figure9_2D_replication.png` saved in the `/results/2D_tests/figure_9/` directory.

---

#### **4.3 Running and Visualizing a Single 2D Case**

To debug or analyze a specific experiment without running the full set, use the single-case scripts.

1.  **Configure and Run a Single Simulation:**
    *   Open the script `./scripts/2D/figure_9/run_single_fig9_test.sh` in a text editor.
    *   In **Step 2**, change the `DISEASE_NAME` and `FIBER_MODEL_NAME` variables to your desired case.
    *   Save the file and run it from the project root:
        ```bash
        ./scripts/2D/figure_9/run_single_fig9_test.sh
        ```
    -   **Output Location:** A new folder like `/results/2D_tests/single_test_fig9/tau_axon_based/`.

2.  **Visualize the Single Result:**
    *   Activate your Conda environment.
    *   Run the single-case visualizer, passing the path to the results folder you just created:
        ```bash
        python3 ./scripts/2D/figure_9/visualize_2d_single_case_figure_9.py results/2D_tests/single_test_fig9/tau_axon_based/
        ```
    -   **Output:** An `activation_plot.png` file will be saved inside that specific results folder.

