# src/visualization/fits_loader.py

import os
from astropy.io import fits
import numpy as np

def load_rgb_datacube(fits_path, instrument="optical_uv"):
    """
    Loads a SKIRT _total.fits datacube and extracts the specific slices 
    corresponding to the Red, Green, and Blue channels based on the instrument.

    Parameters
    ----------
    fits_path: str
        The absolute or relative path to the SKIRT output FITS file.
        (e.g., '/data/outputs/snapshot_150/snapshot_150_optical_uv_total.fits')
    instrument: str, optional
        The instrument name as defined in the SKIRT .ski file. 
        Supports "optical_uv" or "jwst". Defaults to "optical_uv".

    Returns
    -------
    rgb_cube: numpy.ndarray
        A 3D numpy array of shape (3, Y, X) where index 0 is Red, 
        index 1 is Green, and index 2 is Blue.
    """
    # Verify the file exists before attempting to open it
    if not os.path.exists(fits_path):
        raise FileNotFoundError(f"Cannot find FITS file at: {fits_path}")

    # Open the FITS file using astropy. 
    with fits.open(fits_path) as hdul:
        raw_data = hdul[0].data  
        
    # Define the band indices for the chosen instrument
    if instrument == "optical_uv":
        # template.ski: JOHNSON_I (6), JOHNSON_V (4), JOHNSON_B (3)
        r_idx, g_idx, b_idx = 6, 4, 3
    elif instrument == "jwst":
        # template.ski: F200W (2), F150W (1), F090W (0)
        r_idx, g_idx, b_idx = 2, 1, 0
    else:
        raise ValueError(f"Instrument '{instrument}' not supported. Use 'optical_uv' or 'jwst'.")

    # Extract the 2D slices for each color channel
    r_channel = raw_data[r_idx, :, :]
    g_channel = raw_data[g_idx, :, :]
    b_channel = raw_data[b_idx, :, :]

    # Stack them together into a new (3, Y, X) array and return
    rgb_cube = np.stack([r_channel, g_channel, b_channel], axis=0)

    return rgb_cube