# S2S

### 1. Data I/O (`src/file_io/loader.py`)
The first step in the S2S pipeline is reading the raw STARFORGE data. Since STARFORGE is based on the GIZMO codebase, its outputs are standard HDF5 files. The `loader.py` module isolates all `h5py` operations from the rest of the pipeline. It handles:
* **Metadata Extraction:** Pulling cosmological and physical simulation parameters from the HDF5 `Header`.
* **Particle Extraction:** Loading the Meshless Finite Mass (MFM) gas cells (`PartType0`) and sink particles representing stars (`PartType5`) into memory as Python dictionaries.
* **Data Sanitization:** Automatically stripping out boundary and non-physical tracer particles (identifiable by `ParticleIDs >= 1e7`). This ensures that only the physical, star-forming mass is passed to SKIRT9 for radiative transfer modeling.

### 2. Geometric Transformations (`src/geometry/transform.py`)
To optimize the Monte Carlo spatial grid generation in SKIRT and to standardize the synthetic viewing angles, the `transform.py` module spatially normalizes the extracted data. 
* **Centering:** The cloud's geometric median is calculated and shifted to the origin `(0,0,0)`.
* **Principal Axis Alignment:** The module calculates the mass-weighted covariance matrix of the gas cells. By solving for the eigenvectors, it rotates the entire system so that the major axis of the star-forming cloud aligns perfectly with the cartesian axes.
* **Domain Extraction:** A spherical boundary is applied (defined by a user-supplied percentage of the box size) to crop out low-density background material, dramatically reducing SKIRT runtime and memory overhead.
* **Consistency:** All transformations are calculated based on the gas continuum, and the identical translation and rotation matrices are applied to the stellar sink particles to preserve the physical star-gas relative geometry.

### 3. Physical Formatting and Export (`src/processing/formatter.py`)
Once the simulation data is geometrically aligned, it must be translated into physical formats that the SKIRT9 Monte Carlo solvers understand. 
* **Gas Medium Formatting:** The continuous gas data (Meshless Finite Mass cells) is organized into a column-based text file containing coordinates, smoothing lengths, masses, and temperatures. The pipeline simultaneously extracts the spatial domain bounds required to initialize SKIRT's internal spatial grid (e.g., the Voronoi or Octree mesh).
* **Stellar Source Processing:** Because STARFORGE tracks the total bolometric luminosity (including accretion energy), `formatter.py` isolates the *intrinsic* stellar luminosity by calculating and subtracting the accretion luminosity ($L_{\rm acc} \propto G M \dot{M} / R$). Using the intrinsic luminosity and the protostellar radius, the module applies the Stefan-Boltzmann law to derive the effective surface temperature ($T_{\rm eff}$) for each sink particle. This allows SKIRT9 to accurately assign frequency-dependent Spectral Energy Distributions (SEDs) to each star.

### 4. Configuration Automation (`src/config/writer.py` & `replacements.yaml`)
Writing and modifying SKIRT9 XML files manually is highly prone to syntax errors. The S2S pipeline abstracts this entirely.
* **`template.ski`:** A skeletal XML file containing structural definitions for the Monte Carlo algorithms, instrument setups, and probes. All quantitative values are replaced by text placeholders.
* **`replacements.yaml`:** The "Single Source of Truth." Users define their desired physics (e.g., Dust-to-Gas ratio, Grid Refinement levels, Dust models like Draine & Li or THEMIS, and viewing angles) in a simple, readable YAML format.
* **`writer.py`:** Automatically calculates the bounding box of the STARFORGE cloud, merges it with the YAML settings, and compiles the final, executable `.ski` file, perfectly tailored to the specific snapshot being processed.

### 5. Diagnostics and Data Probing (`template.ski` Probes)
To ensure physical accuracy and assist in debugging, the pipeline injects a comprehensive suite of SKIRT Probes into every simulation:
* **Grid Convergence:** `ConvergenceCutsProbe` and `ImportedMediumDensityProbe` are used to verify that the generated Octree spatial grid successfully resolves the raw Lagrangian gas inputs from STARFORGE without losing mass to discretization errors.
* **Internal State Extraction:** Post-simulation probes (`RadiationFieldProbe`, `TemperatureProbe`, `DensityProbe`) extract the fully realized 3D thermodynamic state of the cloud. This provides direct access to localized dust temperatures and the UV/optical mean intensities ($J_\lambda$), enabling deeper analysis beyond simple 2D synthetic imaging.

### 5. Automated Execution (`src/runner/run_sim.py`)
The `run_sim.py` module acts as the master orchestrator for the entire S2S pipeline, designed specifically to support automated batch-processing on high-performance computing (HPC) clusters.
By simply providing an input `.hdf5` path and a target output directory via command-line arguments, the runner will sequentially trigger the data loader, geometric transformer, physical formatter, and XML writer. Finally, it utilizes the Python Toolkit for SKIRT (`PTS9`) to programmatically invoke the SKIRT9 C++ executable, yielding fully synthesized mock observations and 3D physical state probes without any manual intervention.

## Execution & Usage Guide

The S2S pipeline can be executed either interactively on a local machine/Jupyter environment or as an automated batch job on a Kubernetes High-Performance Computing (HPC) cluster.

### Option 1: Local / Interactive Execution
This method is recommended for testing, debugging, or running a single snapshot.

**Prerequisites:**
1. Clone this repository: `git clone https://github.com/AislingA/S2S.git`
2. Ensure you have the SKIRT9 binary installed and the Python Toolkit for SKIRT (`PTS9`) accessible in your Python path.
3. Download the necessary SKIRT Resources (telescope filters) to your local machine and set the environment variable:
   `export SKIRT_RESOURCES_PATH=/path/to/your/Resources`

**Running the Pipeline:**
Navigate to the root of the `S2S` directory and run the master script, pointing it to your raw HDF5 file and desired output directory:
```bash
python src/runner/run_sim.py \
  --input /path/to/snapshot_150.hdf5 \
  --output-dir /path/to/save/outputs/
```

### Option 2: Automated Kubernetes Batch Job (HPC)
For processing multiple, massive STARFORGE snapshots, the pipeline is designed to run asynchronously on a Kubernetes cluster. This method handles dependency installation, S3 data routing, and memory management automatically.

**Workflow:**
1. Ensure your STARFORGE `.hdf5` snapshots and your `skirt_resources` folder are uploaded to your S3 storage (e.g., `s3:sdsu-rosen/`).
2. Configure the `pipeline-job.yaml` file with your target snapshot numbers (e.g., `SNAPSHOTS=(150 200 250)`).
3. Submit the job to the cluster scheduler:
   ```bash
   kubectl apply -f pipeline-job.yaml
   ```
4. Monitor the pipeline's progress in real-time:
   ```bash
   kubectl logs -f job/skirt-batch-job
   ```
Once completed, the pipeline will automatically push the SKIRT `.ski` files, `.txt` physical arrays, and `.fits` mock observations back to your S3 storage bucket under `skirt_outputs/`.
```