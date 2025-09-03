
---

### **README.md**

This project provides a C++ implementation of the reaction-diffusion model for neurodegenerative diseases presented in Weickenmeier et al. (2019), "A physics-based model explains the prion-like features of neurodegeneration...".

The code is designed to be compiled and run in a Linux-based environment (e.g., a Docker container or native Linux) with the deal.II finite element library. Visualization is handled by a set of Python 3 scripts.

### 1. Project Structure

-   `/src`: Contains all C++ source code.
-   `/include`: Contains all C++ header files.
-   `/meshes`: Should contain all input mesh files.
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
    This will create the executables `neuro_disease_1D`, `neuro_disease_2D`, etc., inside the `build` folder.

### 3. Running the 1D Experiment (Figure 3 Replication)

This experiment runs a 3x3 parameter sweep to test the effects of growth (`alpha`) and spreading (`d`) on the 1D reaction-diffusion model.

#### **Step 1: Run the C++ Simulations**

From the main project directory, execute the experiment script. This will run 9 simulations in parallel and save the raw data.

```bash
./scripts/1D/run_1d_experiment.sh
```
-   **Output:** The raw simulation data will be saved in organized subfolders inside `/results/1D_tests/figure_3/`.

#### **Step 2: Visualize the Results**

This project includes two different Python scripts to visualize the 1D results. Make sure you are inside your Conda environment if you have one set up.

*   **To create the "Activation Time" plot (like the paper's Figure 3):**
    This plot shows *when* each point on the line reached a certain concentration. From the main project directory, run:
    ```bash
    python3 ./scripts/1D/visualize_1d_experiment.py
    ```
    -   **Output:** A composite 3x3 image named `figure3_replication.png` will be saved in the `/results/1D_tests/figure_3/` directory.

*   **To create the "Space-Time" plot:**
    This plot shows the raw concentration value at every point in space and time. From the main project directory, run:
    ```bash
    python3 ./scripts/1D/visualize_spacetime_figure_3.py
    ```
    -   **Output:** A composite 3x3 image named `figure3_spacetime_replication.png` will be saved in the `/results/1D_tests/figure_3/` directory.

### 4. Running the 2D Experiment (Figure 8 Replication)

This experiment runs a 2x2 sensitivity analysis on a 2D brain slice mesh. It uses a more complex, two-stage visualization process that requires **ParaView**.

#### **Step 1: Install ParaView (One-Time Setup)**

The visualization scripts rely on `pvpython`, the Python interpreter that comes with ParaView. The recommended way to install this in your container is with Conda.

1.  **Install Conda/Miniforge:**
    ```bash
    wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
    bash Miniforge3-Linux-x86_64.sh -b -p /opt/miniforge3
    source /opt/miniforge3/bin/activate
    ```
2.  **Create a ParaView Environment:**
    ```bash
    conda create -n pv-env -c conda-forge paraview -y
    ```
3.  **Install Python Libraries:** After creating the environment, activate it and install the necessary libraries for the final plotting stage.
    ```bash
    source /opt/miniforge3/bin/activate pv-env
    pip install matplotlib meshio numpy
    ```

#### **Step 2: Run the C++ Simulations**

From the main project directory, execute the experiment script. This will run 4 parallel simulations.

```bash
./scripts/2D/figure_8/ext_axn_alpha_sensitivity_analysis_figure_8.sh
```
-   **Output:** The raw simulation data will be saved in subfolders inside `/results/2D_tests/figure_8/`.

#### **Step 3: Visualize the Results**

The visualization is a two-stage process orchestrated by a single script.

1.  **Activate your Conda environment:**
    ```bash
    source /opt/miniforge3/bin/activate pv-env
    ```
2.  **Run the main visualization script:** From the main project directory, run:
    ```bash
    python3 ./scripts/2D/figure_8/visualize_2d_figure_8.py
    ```
    -   **What it Does:** This script will automatically call `pvpython` to process the raw data for each of the 4 cases, creating intermediate `.png` files. It will then assemble these images into a final 2x2 composite plot.
    -   **Output:** The final image, `figure8_2D_replication.png`, will be saved in the main project directory.

---