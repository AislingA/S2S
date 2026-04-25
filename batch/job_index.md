# SKIRT Batch Jobs Index

This document tracks the different Kubernetes batch job files and the specific GitHub configurations required to run them successfully, separated by dust configurations.

**NOTE:** The `.yaml` files in this batch folder only control the cluster execution (percentage, S3 output paths). Because the job clones the S2S repository fresh every time, **you must ensure the correct physics parameters are pushed to GitHub `src/config/replacements.yaml` before applying the job.**

---

## Draine Dust Mix (2007)

### Equilibrium Heating

#### `draine-100-eqheating.yaml`
**Description:** Runs the pipeline extracting 100% of the box, utilizing the Draine & Li (2007) dust model, with thermal equilibrium heating only (Stochastic PAH heating turned OFF to force convergence with the new accretion source). Outputs are saved to `s3:sdsu-rosen/skirt_outputs/draine_100_eqheating/`.

**Pre-Run GitHub Checklist:**
Before running `kubectl apply -f draine-100-eqheating.yaml`, ensure the following are committed and pushed to `src/config/replacements.yaml`:

1.  **Dust Model:** * `DUST_MIX_CONFIG: '<DraineLiDustMix numSilicateSizes="5" numGraphiteSizes="5" numPAHSizes="5"/>'` is uncommented/active.
2.  **Heating Mechanism:** * `DUST_EMISSION_TYPE: "Equilibrium"`
3.  **Convergence:**
    * `ITERATION_MULTIPLIER: "2"` This runs uses a diffuse box and simple thermal physics, and should converge within 10 iterations.

Notes:
- approx 28 iterations
- approx 4 mins per iteration
- approx 4 hours per snapshot

#### `draine-075-eqheating.yaml`
**Description:** Runs the pipeline extracting 75% of the box, utilizing the Draine & Li (2007) dust model, with thermal equilibrium heating only. Outputs are saved to `s3:sdsu-rosen/skirt_outputs/draine_075_eqheating/`.

**Pre-Run GitHub Checklist:**
Before running `kubectl apply -f draine-075-eqheating.yaml`, ensure the following are committed and pushed to `src/config/replacements.yaml`:

1.  **Dust Model:** * `DUST_MIX_CONFIG: '<DraineLiDustMix numSilicateSizes="5" numGraphiteSizes="5" numPAHSizes="5"/>'` is uncommented/active.
2.  **Heating Mechanism:** * `DUST_EMISSION_TYPE: "Equilibrium"`
3.  **Convergence:**
    * `ITERATION_MULTIPLIER: "2"` This runs uses a diffuse box and simple thermal physics, and should converge within 10 iterations.

#### `draine-050-eqheating.yaml`
**Description:** Runs the pipeline extracting 50% of the box, utilizing the Draine & Li (2007) dust model, with thermal equilibrium heating only. Outputs are saved to `s3:sdsu-rosen/skirt_outputs/draine_050_eqheating/`.

**Pre-Run GitHub Checklist:**
Before running `kubectl apply -f draine-050-eqheating.yaml`, ensure the following are committed and pushed to `src/config/replacements.yaml`:

1.  **Dust Model:** * `DUST_MIX_CONFIG: '<DraineLiDustMix numSilicateSizes="5" numGraphiteSizes="5" numPAHSizes="5"/>'` is uncommented/active.
2.  **Heating Mechanism:** * `DUST_EMISSION_TYPE: "Equilibrium"`
3.  **Convergence:**
    * `ITERATION_MULTIPLIER: "3"` This runs is getting denser. Increasing the iteration multipler ensures the accretion does not cause any wobbles.

#### `draine-025-eqheating.yaml`
**Description:** Runs the pipeline extracting 25% of the box, utilizing the Draine & Li (2007) dust model, with thermal equilibrium heating only. Outputs are saved to `s3:sdsu-rosen/skirt_outputs/draine_025_eqheating/`.

**Pre-Run GitHub Checklist:**
Before running `kubectl apply -f draine-025-eqheating.yaml`, ensure the following are committed and pushed to `src/config/replacements.yaml`:

1.  **Dust Model:** * `DUST_MIX_CONFIG: '<DraineLiDustMix numSilicateSizes="5" numGraphiteSizes="5" numPAHSizes="5"/>'` is uncommented/active.
2.  **Heating Mechanism:** * `DUST_EMISSION_TYPE: "Equilibrium"`
3.  **Convergence:**
    * `ITERATION_MULTIPLIER: "4"` In this run, only the dense cores remain. This needs a slightly higher sample size to balance the dense grid.

### Stochastic Heating

#### `draine-100-stochastic.yaml`
**Description:** Runs the pipeline extracting 100% of the box, utilizing the Draine & Li (2007) dust model, with . Outputs are saved to `s3:sdsu-rosen/skirt_outputs/draine_100_stochastic/`.

**Pre-Run GitHub Checklist:**
Before running `kubectl apply -f draine-100-stochastic.yaml`, ensure the following are committed and pushed to `src/config/replacements.yaml`:

1.  **Dust Model:** * `DUST_MIX_CONFIG: '<DraineLiDustMix numSilicateSizes="5" numGraphiteSizes="5" numPAHSizes="5"/>'` is uncommented/active.
2.  **Heating Mechanism:** * `DUST_EMISSION_TYPE: "Stochastic"`
3.  **Convergence:**
    * `ITERATION_MULTIPLIER: "5"` Stochastic PAH heating is expensive. Using a 5x multipler is the minimum to smooth out the noise.

#### `draine-075-stochastic.yaml`
**Description:** Runs the pipeline extracting 75% of the box, utilizing the Draine & Li (2007) dust model, with . Outputs are saved to `s3:sdsu-rosen/skirt_outputs/draine_075_stochastic/`.

**Pre-Run GitHub Checklist:**
Before running `kubectl apply -f draine-075-stochastic.yaml`, ensure the following are committed and pushed to `src/config/replacements.yaml`:

1.  **Dust Model:** * `DUST_MIX_CONFIG: '<DraineLiDustMix numSilicateSizes="5" numGraphiteSizes="5" numPAHSizes="5"/>'` is uncommented/active.
2.  **Heating Mechanism:** * `DUST_EMISSION_TYPE: "Stochastic"`
3.  **Convergence:**
    * `ITERATION_MULTIPLIER: "6"` Removing outer boundaries forces the grid to focus on the noisy star-forming regions.

#### `draine-050-stochastic.yaml`
**Description:** Runs the pipeline extracting 50% of the box, utilizing the Draine & Li (2007) dust model, with . Outputs are saved to `s3:sdsu-rosen/skirt_outputs/draine_050_stochastic/`.

**Pre-Run GitHub Checklist:**
Before running `kubectl apply -f draine-050-stochastic.yaml`, ensure the following are committed and pushed to `src/config/replacements.yaml`:

1.  **Dust Model:** * `DUST_MIX_CONFIG: '<DraineLiDustMix numSilicateSizes="5" numGraphiteSizes="5" numPAHSizes="5"/>'` is uncommented/active.
2.  **Heating Mechanism:** * `DUST_EMISSION_TYPE: "Stochastic"`
3.  **Convergence:**
    * `ITERATION_MULTIPLIER: "8"` High density, stochastic heating, and accretion needs to be heavily smoothed to converge.

#### `draine-025-stochastic.yaml`
**Description:** Runs the pipeline extracting 25% of the box, utilizing the Draine & Li (2007) dust model, with . Outputs are saved to `s3:sdsu-rosen/skirt_outputs/draine_025_stochastic/`.

**Pre-Run GitHub Checklist:**
Before running `kubectl apply -f draine-025-stochastic.yaml`, ensure the following are committed and pushed to `src/config/replacements.yaml`:

1.  **Dust Model:** * `DUST_MIX_CONFIG: '<DraineLiDustMix numSilicateSizes="5" numGraphiteSizes="5" numPAHSizes="5"/>'` is uncommented/active.
2.  **Heating Mechanism:** * `DUST_EMISSION_TYPE: "Stochastic"`
3.  **Convergence:**
    * `ITERATION_MULTIPLIER: "10"` This is the ultimate stress test. All photons are trapped in the dense cores with small grains. This needs the maximum packets.

## THEMIS Dust Mix (year)

### Equilibrium Heating

#### `themis-100-eqheating.yaml`
**Description:** Runs the pipeline extracting 100% of the box, utilizing the THEMIS (year) dust model, with thermal equilibrium heating only (Stochastic PAH heating turned OFF to force convergence with the new accretion source). Outputs are saved to `s3:sdsu-rosen/skirt_outputs/themis_100_eqheating/`.

**Pre-Run GitHub Checklist:**
Before running `kubectl apply -f themis-100-eqheating.yaml`, ensure the following are committed and pushed to `src/config/replacements.yaml`:

1.  **Dust Model:** * `DUST_MIX_CONFIG: '<ThemisDustMix numSilicateSizes="5" numHydrocarbonSizes="5"/>'` is uncommented/active.
2.  **Heating Mechanism:** * `DUST_EMISSION_TYPE: "Equilibrium"`
3.  **Convergence:**
    * `ITERATION_MULTIPLIER: "2"` This runs uses a diffuse box and simple thermal physics, and should converge within 10 iterations.

Notes:
- approx  iterations
- approx  mins per iteration
- approx  hours per snapshot