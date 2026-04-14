# SKIRT Batch Jobs Index

This document tracks the different Kubernetes batch job files and the specific GitHub configurations required to run them successfully, separated by dust configurations.

**NOTE:** The `.yaml` files in this batch folder only control the cluster execution (percentage, S3 output paths). Because the job clones the S2S repository fresh every time, **you must ensure the correct physics parameters are pushed to GitHub `src/config/replacements.yaml` before applying the job.**

---

## Draine Dust Mix (2007)

### `draine-100-eqheating.yaml`
**Description:** Runs the pipeline extracting 100% of the box, utilizing the Draine & Li (2007) dust model, with thermal equilibrium heating only (Stochastic PAH heating turned OFF to force convergence with the new accretion source). Outputs are saved to `s3:sdsu-rosen/skirt_outputs/draine_100_eqheating/`.

**Pre-Run GitHub Checklist:**
Before running `kubectl apply -f draine-100-eqheating.yaml`, ensure the following are committed and pushed to `src/config/replacements.yaml`:

1.  **Dust Model:** * `DUST_MIX_CONFIG: '<DraineLiDustMix numSilicateSizes="5" numGraphiteSizes="5" numPAHSizes="5"/>'` is uncommented/active.
2.  **Heating Mechanism:** * `DUST_EMISSION_TYPE: "Equilibrium"`
3.  **Convergence Safety Net:**
    * `MAX_ITERATIONS: "40"`
    * `ITERATION_MULTIPLIER: "10"`