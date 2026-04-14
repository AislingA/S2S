The methodology used in the paper to create the synthetic observations in Figure 1 relies on extracting raw physical fluxes from SKIRT and applying observational effects afterward. 

[cite_start]The authors utilized SKIRT 9 to generate broadband images representing the system at different physical distances[cite: 1586]. [cite_start]To create the color composite seen in Figure 1, they extracted the outputs for the Johnson B, V, and I broadband filters, mapping them to the blue, green, and red image channels, respectively[cite: 1624]. [cite_start]To mimic the observational capabilities of the Hubble Space Telescope (HST), they did not rely on SKIRT's internal pipeline; instead, they applied a post-processing 2D Gaussian Point Spread Function (PSF) with a Full Width at Half Maximum (FWHM) of 2.1 pixels[cite: 1625]. [cite_start]Finally, because the dynamic range of star-forming regions is extreme, they applied a logarithmic stretch to the fluxes to make both the bright central starburst and the faint surrounding structures visible simultaneously[cite: 1626].

### Instrument Viability for RGB Composites

Your `template.ski` file defines three distinct instruments. Here is how they apply to creating synthetic RGB images:

* **`optical_uv` (Highly Recommended for Fig 1 Recreation):** This instrument includes the exact Johnson B, V, and I broadbands used in the paper. Because these mimic the optical filters of HST, this is the perfect instrument to use if you want to perform the Gaussian PSF convolution to recreate the exact aesthetic of Figure 1.
* **`jwst` (Highly Recommended for Infrared Mock Observations):** This instrument includes broadbands for NIRCam and MIRI. You can absolutely create a synthetic RGB image using this data by mapping three NIRCam filters (e.g., F090W to Blue, F150W to Green, and F200W to Red). The post-processing convolution would be the same, but to be scientifically accurate, you would use a JWST-specific PSF (which is heavily diffraction-limited and filter-dependent) rather than a simple Gaussian.
* **`logwaves` (Not Recommended for RGB Images):** This instrument records a continuous, highly resolved spectral data cube rather than integrating the light into specific telescope broadbands. It is designed for extracting detailed Spectral Energy Distributions (SEDs) and analyzing fine spectral features across the entire system. While you could technically slice it to form colors, it bypasses the physical filter transmission curves that give synthetic observations their realism. 

### Step-by-Step Post-Processing Plan

To achieve your goal of a 2x2 grid without touching the underlying S2S pipeline code, we will structure the workflow as follows:

**Step 1: Simulation Configuration and Execution**
Before post-processing, you need the four specific vantage points. You will configure your `replacements.yaml` (or batch job) to run the simulation four times, altering the `DISTANCE` and `INCLINATION` variables:
* 1 kpc Face-on: `DISTANCE: 1 kpc`, `INCLINATION: 0 deg`
* 1 kpc Edge-on: `DISTANCE: 1 kpc`, `INCLINATION: 90 deg`
* 10 kpc Face-on: `DISTANCE: 10 kpc`, `INCLINATION: 0 deg`
* 10 kpc Edge-on: `DISTANCE: 10 kpc`, `INCLINATION: 90 deg`

**Step 2: FITS File Ingestion**
Once the pipeline finishes, you will write a standalone post-processing script. This script will locate and load the multi-extension `.fits` datacubes generated specifically by the `optical_uv` (or `jwst`) instrument for all four simulation runs.

**Step 3: Band Extraction**
The script will index into the loaded FITS datacubes to isolate the specific 2D flux arrays corresponding to your desired filters (e.g., extracting the data arrays for the Johnson B, V, and I bands). 

**Step 4: PSF Kernel Generation**
You will define the mathematical kernel that represents your telescope's optics. For the HST recreation, this will be a 2D Gaussian kernel parameterized by the FWHM of 2.1 pixels. If you process the JWST instrument, you would load a synthetic WebbPSF kernel here instead.

**Step 5: Image Convolution**
Using a signal processing library, you will convolve the independent Red, Green, and Blue 2D data arrays with your chosen PSF kernel. This mathematical operation smears the raw physical flux according to the telescope's optical limitations.

**Step 6: Non-Linear Scaling and Compositing**
Because the convolved arrays contain physical units (e.g., MJy/sr) with massive dynamic ranges, you will apply a non-linear stretch (such as a logarithmic or arcsinh function) to the data. The stretched arrays will then be stacked along the third axis to create standard 3D RGB image arrays.

**Step 7: Grid Assembly and Plotting**
Finally, using a plotting library, you will create a figure with a 2x2 subplot grid. The script will render the 1 kpc images on the top row and the 10 kpc images on the bottom row, applying any necessary scale bars or annotations before saving the final visualization.