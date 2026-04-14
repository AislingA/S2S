# src/visualization/image_processing.py

import numpy as np
import warnings
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.visualization import AsymmetricPercentileInterval, ImageNormalize, AsinhStretch

def apply_psf(rgb_cube, fwhm_pixels=2.1, instrument="optical_uv"):
    """
    Convolves the raw physical flux arrays with a Point Spread Function (PSF).
    """
    # Allow bypassing the PSF to check pristine dust structures
    if fwhm_pixels <= 0:
        return rgb_cube

    sigma = fwhm_pixels / 2.35482
    kernel = Gaussian2DKernel(x_stddev=sigma)
    convolved_cube = np.zeros_like(rgb_cube)
    
    print(f"Applying PSF convolution (FWHM={fwhm_pixels}px) for {instrument}...")
    for i in range(3):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            convolved_cube[i] = convolve(rgb_cube[i], kernel, normalize_kernel=True)
            
    return convolved_cube

def aggressive_scale(band_data, lower_pct=1.0, upper_pct=98.0):
    """
    Exact replica of the old scaling logic, adapted to safely handle the interval bounds.
    """
    if not np.any(band_data > 0):
        return np.zeros_like(band_data)

    # Uses Asymmetric to correctly respect both your lower and upper bounds
    interval = AsymmetricPercentileInterval(lower_pct, upper_pct)
    vmin, vmax = interval.get_limits(band_data[band_data > 0])
    
    # Hardcoded a=0.01 as it was in your old script
    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch(a=0.01))
    scaled = norm(band_data)
    
    # No np.clip() here. We let Matplotlib natively handle values > 1.0
    return np.ma.filled(scaled, 0)

def apply_asinh_stretch(rgb_cube, lower_pct=1.0, upper_pct=98.0):
    """
    Applies the aggressive astronomical stretch independently to the cube.
    """
    print(f"Applying independent Asinh stretch (Percentiles: {lower_pct}-{upper_pct}%)...")
    
    # Process Red, Green, and Blue exactly as the old code did
    r_band = aggressive_scale(rgb_cube[0], lower_pct, upper_pct)
    g_band = aggressive_scale(rgb_cube[1], lower_pct, upper_pct)
    b_band = aggressive_scale(rgb_cube[2], lower_pct, upper_pct)
    
    # Turn the three (Y, X) arrays into the (Y, X, 3) image
    rgb_image = np.dstack((r_band, g_band, b_band))
    
    return rgb_image