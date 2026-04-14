# src/visualization/image_processing.py

import numpy as np
import warnings
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.visualization import PercentileInterval, ImageNormalize, AsinhStretch

def apply_psf(rgb_cube, fwhm_pixels=2.1, instrument="optical_uv"):
    """
    Convolves the raw physical flux arrays with a Point Spread Function (PSF).
    """
    sigma = fwhm_pixels / 2.35482
    kernel = Gaussian2DKernel(x_stddev=sigma)
    convolved_cube = np.zeros_like(rgb_cube)
    
    print(f"Applying PSF convolution (FWHM={fwhm_pixels}px) for {instrument}...")
    for i in range(3):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            convolved_cube[i] = convolve(rgb_cube[i], kernel, normalize_kernel=True)
            
    return convolved_cube

def apply_asinh_stretch(rgb_cube, lower_pct=1.0, upper_pct=98.0, a_param=0.01):
    """
    Applies an aggressive Asinh stretch to each color channel independently 
    to reveal faint dust structures by deliberately saturating bright point sources.
    
    Parameters
    ----------
    rgb_cube: numpy.ndarray
        The (3, Y, X) array containing the convolved physical fluxes.
    lower_pct: float
        The lower percentile limit for the noise floor.
    upper_pct: float
        The upper percentile limit to saturate bright pixels.
    a_param: float
        The alpha parameter for the Asinh curve transition.

    Returns
    -------
    rgb_image: numpy.ndarray
        A (Y, X, 3) array normalized between 0.0 and 1.0, ready for plotting.
    """
    print(f"Applying independent Asinh stretch (Percentiles: {lower_pct}-{upper_pct}%)...")
    
    stretched_channels = []
    
    # Process Red, Green, and Blue independently to maximize structural contrast
    for i in range(3):
        band_data = rgb_cube[i]
        
        # Safety check: if the channel is completely empty, return zeros
        if not np.any(band_data > 0):
            stretched_channels.append(np.zeros_like(band_data))
            continue
            
        # Isolate positive values to calculate percentiles accurately, ignoring empty space
        positive_data = band_data[band_data > 0]
        interval = PercentileInterval(lower_pct, upper_pct)
        vmin, vmax = interval.get_limits(positive_data)
        
        # Create the normalization object using your custom Asinh logic
        norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch(a=a_param))
        
        # Apply the normalization to the full band data
        scaled = norm(band_data)
        
        # Astropy normalizations return MaskedArrays. Fill masked values with 0.
        filled_scaled = np.ma.filled(scaled, 0)
        
        # Final safety clip strictly between 0.0 and 1.0
        clipped = np.clip(filled_scaled, 0, 1)
        stretched_channels.append(clipped)
        
    # np.dstack takes a list of 2D arrays and stacks them along the 3rd axis
    # turning the three (Y, X) arrays into the (Y, X, 3) image for Matplotlib
    rgb_image = np.dstack(stretched_channels)
    
    return rgb_image