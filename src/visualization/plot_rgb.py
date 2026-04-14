# src/visualization/plot_rgb.py

import os
import matplotlib.pyplot as plt
from fits_loader import load_rgb_datacube
from image_processing import apply_psf, apply_log_stretch

def create_2x2_mock_observation(file_paths, instrument="optical_uv", stretch_factor=1000, fwhm=2.1):
    """
    Orchestrates the loading, processing, and plotting of 4 SKIRT simulations 
    into a single 2x2 synthetic observation grid.

    Parameters
    ----------
    file_paths: dict
        A dictionary containing the absolute paths to the 4 FITS files.
        !!! Expected keys: '1kpc_face', '1kpc_edge', '10kpc_face', '10kpc_edge'.
    instrument: str
        The instrument to process ("optical_uv" or "jwst").
    stretch_factor: float
        The non-linear logarithmic stretch factor.
    fwhm: float
        The Full Width at Half Maximum for the Gaussian PSF in pixels.
    """
    # Initialize the Matplotlib Figure
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 12), constrained_layout=True)
    
    # Set the background to black to match astronomical images
    fig.patch.set_facecolor('black')
    
    # !!! Map the grid positions to our expected dictionary keys
    grid_mapping = {
        '1kpc_face':  {'ax': axes[0, 0], 'title': '1 kpc | Face-on (0°)'},
        '1kpc_edge':  {'ax': axes[0, 1], 'title': '1 kpc | Edge-on (90°)'},
        '10kpc_face': {'ax': axes[1, 0], 'title': '10 kpc | Face-on (0°)'},
        '10kpc_edge': {'ax': axes[1, 1], 'title': '10 kpc | Edge-on (90°)'}
    }

    print("Starting 2x2 grid generation...")

    # Process each panel in the grid
    for key, config in grid_mapping.items():
        ax = config['ax']
        title = config['title']
        filepath = file_paths.get(key, "")

        # Formatting the subplot
        ax.set_title(title, color='white', fontsize=14, pad=10)
        ax.axis('off')
        
        # Check if the file was provided and exists
        if not os.path.exists(filepath):
            print(f"Warning: File not found for {key} at {filepath}")
            ax.text(0.5, 0.5, 'FILE NOT FOUND', color='red', 
                    ha='center', va='center', transform=ax.transAxes, fontsize=12)
            continue
            
        try:
            print(f"\nProcessing {key}...")            
            # Load the raw physical fluxes (Color, Y, X)
            raw_cube = load_rgb_datacube(filepath, instrument=instrument)
            
            # Apply Telescope Optics (Gaussian PSF Convolution)
            smoothed_cube = apply_psf(raw_cube, fwhm_pixels=fwhm, instrument=instrument)
            
            # Apply Astronomical Scaling (Log Stretch) -> returns (Y, X, Color)
            final_rgb_image = apply_log_stretch(smoothed_cube, stretch_factor=stretch_factor)
            
            # Render the image onto the subplot
            ax.imshow(final_rgb_image, origin='lower')
            
        except Exception as e:
            print(f"Failed to process {key}: {e}")
            ax.text(0.5, 0.5, 'PROCESSING ERROR', color='red', 
                    ha='center', va='center', transform=ax.transAxes)

    # Save and Show the Result
    output_filename = f"synthetic_observation_{instrument}.png"
    plt.savefig(output_filename, dpi=300, facecolor='black', edgecolor='none')
    print(f"\nSuccess! Grid saved as {output_filename}")
    
    plt.show()

if __name__ == "__main__":
    fits_files = {
        '1kpc_face':  "../../data/outputs/snapshot_150_1kpc_face/snapshot_150_optical_uv_total.fits",
        '1kpc_edge':  "../../data/outputs/snapshot_150_1kpc_edge/snapshot_150_optical_uv_total.fits",
        '10kpc_face': "../../data/outputs/snapshot_150_10kpc_face/snapshot_150_optical_uv_total.fits",
        '10kpc_edge': "../../data/outputs/snapshot_150_10kpc_edge/snapshot_150_optical_uv_total.fits"
    }
    
    create_2x2_mock_observation(
        file_paths=fits_files, 
        instrument="optical_uv", 
        stretch_factor=1000, 
        fwhm=2.1
    )