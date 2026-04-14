# src/visualization/plot_rgb.py

import os
import argparse
import matplotlib.pyplot as plt
from fits_loader import load_rgb_datacube
from image_processing import apply_psf, apply_asinh_stretch

def build_file_paths(snapshot, batch_job, instrument="optical_uv", base_dir="../../outputs/data"):
    """
    Dynamically constructs the expected file paths based on the S2S directory structure.
    
    Because this script is in S2S/src/visualization, the relative path 
    '../../outputs/data' moves up two directories to reach the data folder.
    """
    # Construct the target directory path: S2S/outputs/data/[BATCH JOB]/[SNAPSHOT]/
    target_dir = os.path.join(base_dir, batch_job, f"snapshot_{snapshot}")
    
    # Map the conceptual grid positions to the actual file names
    # "near" = 1 kpc, "far" = 10 kpc
    paths = {
        '1kpc_face':  os.path.join(target_dir, f"snapshot_{snapshot}_{instrument}_near_faceon_total.fits"),
        '1kpc_edge':  os.path.join(target_dir, f"snapshot_{snapshot}_{instrument}_near_edgeon_total.fits"),
        '10kpc_face': os.path.join(target_dir, f"snapshot_{snapshot}_{instrument}_far_faceon_total.fits"),
        '10kpc_edge': os.path.join(target_dir, f"snapshot_{snapshot}_{instrument}_far_edgeon_total.fits")
    }
    return paths

def create_2x2_mock_observation(file_paths, snapshot, batch_job, instrument="optical_uv", fwhm=2.1, out_dir='../../outputs/figures'):
    """
    Orchestrates the loading, processing, and plotting of 4 SKIRT simulations 
    into a single 2x2 synthetic observation grid.
    """
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 12), constrained_layout=True)
    fig.patch.set_facecolor('black')
    
    grid_mapping = {
        '1kpc_face':  {'ax': axes[0, 0], 'title': '1 kpc | Face-on (0°)'},
        '1kpc_edge':  {'ax': axes[0, 1], 'title': '1 kpc | Edge-on (90°)'},
        '10kpc_face': {'ax': axes[1, 0], 'title': '10 kpc | Face-on (0°)'},
        '10kpc_edge': {'ax': axes[1, 1], 'title': '10 kpc | Edge-on (90°)'}
    }

    print(f"Generating 2x2 grid for Snapshot {snapshot} using {instrument}...")

    for key, config in grid_mapping.items():
        ax = config['ax']
        filepath = file_paths.get(key, "")

        ax.set_title(config['title'], color='white', fontsize=14, pad=10)
        ax.axis('off') 
        
        if not os.path.exists(filepath):
            print(f"Warning: Missing file at {filepath}")
            ax.text(0.5, 0.5, 'FILE NOT FOUND', color='red', ha='center', va='center', transform=ax.transAxes, fontsize=12)
            continue
            
        try:
            print(f"  -> Processing {key}...")
            raw_cube = load_rgb_datacube(filepath, instrument=instrument)
            smoothed_cube = apply_psf(raw_cube, fwhm_pixels=fwhm, instrument=instrument)
            
            # Use upper_pct=98 to exactly match your old code's defaults
            final_rgb_image = apply_asinh_stretch(smoothed_cube, lower_pct=1.0, upper_pct=98.0)            
            
            ax.imshow(final_rgb_image, origin='lower')
            
        except Exception as e:
            print(f"  -> Failed to process {key}: {e}")
            ax.text(0.5, 0.5, 'PROCESSING ERROR', color='red', ha='center', va='center', transform=ax.transAxes)

    os.makedirs(out_dir, exist_ok=True)
    output_filename = f"syn_obs_{batch_job}_snap{snapshot}_{instrument}.png"
    output_path = os.path.join(out_dir, output_filename)
    plt.savefig(output_path, dpi=300, facecolor='black', edgecolor='none')
    print(f"\nSuccess! Grid saved as {output_filename}")
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate 2x2 mock observation grids from SKIRT FITS files.")
    parser.add_argument("--snapshot", required=True, type=str, help="The snapshot number (e.g., 150, 200).")
    parser.add_argument("--job", required=True, type=str, help="The name of the batch job directory (e.g., run_01).")
    parser.add_argument("--instrument", type=str, default="optical_uv", choices=['optical_uv', 'jwst'])
    parser.add_argument("--fwhm", type=float, default=2.1, help="PSF FWHM in pixels (default: 2.1). Set to 0 to disable.")
    
    args = parser.parse_args()
    
    dynamic_paths = build_file_paths(snapshot=args.snapshot, batch_job=args.job, instrument=args.instrument)
    
    create_2x2_mock_observation(
        file_paths=dynamic_paths,
        snapshot=args.snapshot,
        batch_job=args.job,
        instrument=args.instrument,
        fwhm=args.fwhm
    )