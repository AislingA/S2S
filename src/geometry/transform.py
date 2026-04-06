# src/geometry/transform.py

"""
Module: transform.py
Description: 
    This module performs geometric transformations on hydrodynamic simulation data. 
    It is responsible for shifting the region of interest to the coordinate origin, 
    rotating the system to align with the principal axes of the mass distribution, 
    and applying a spherical spatial mask. These steps are critical for optimizing 
    the spatial grid boundaries and standardizing viewing angles for SKIRT9.
"""

import numpy as np

def center_on_origin(coords, center=None):
    """
    Translates coordinates so that the region of interest is at (0,0,0).

    The translation is defined as:
    $ \vec{r}' = \vec{r} - \vec{r}_{\rm center} $

    Parameters
    ----------
    coords : np.ndarray
        Nx3 array of particle positions.
    center : np.ndarray, optional
        An external (x,y,z) coordinate to center on. If None, it calculates 
        the median of the provided coordinates.

    Returns
    -------
    centered_coords : np.ndarray
        The translated coordinates.
    center : np.ndarray
        The translation vector used (useful for centering other particle types).
    """
    # If no center is forced, calculate the geometric median of the cloud.
    # Using median rather than mean to prevent massive outliers from skewing the center.
    if center is None:
        center = np.median(coords, axis=0)

    # Shift all particles so the calculated center becomes the new origin
    return coords - center, center

def apply_rotation(coords, masses=None, basis=None):
    """
    Rotates coordinates to align with a specific basis. If no basis is 
    provided, it calculates one using the mass-weighted covariance matrix.

    The mass-weighted covariance matrix $C$ is diagonalized to find the principal axes:
    The coordinates are then transformed into the new basis $R$:

    Parameters
    ----------
    coords : np.ndarray
        Nx3 array of centered particle positions.
    masses : np.ndarray, optional
        Array of particle masses used to weight the covariance matrix.
    basis : np.ndarray, optional
        A pre-calculated 3x3 rotation matrix.

    Returns
    -------
    rotated_coords : np.ndarray
        The coordinates transformed into the principarl axis frame.
    basis : np.ndarray
        The 3x3 rotation matrix used for the transformation.
    """
    if basis is None:
        # calculate the mass-weighted covariance matrix of the spatial distribution
        # coords.T is used because np.cov expects variables as rows and observations as columns.
        cov_pos = np.cov(coords.T, aweights=masses)
        try:
            # Diagonalize the covariance matrix to get eigenvalues (w) and eigenvectors (basis)
            w, basis = np.linalg.eigh(cov_pos)
            # Sort eigenvectors by descending eigenvalue (variance/mass spread)
            # This aligns the x-axis with the major axis of the cloud.
            basis = basis[:, w.argsort()[::-1]]
        except np.linalg.LinAlgError:
            # Fallback to the identity matrix if diagonalization fails.
            print("Warning: Matrix did not converge. Using Identity.")
            basis = np.identity(3)

    # Apply the rotation matrix via dot product.
    return coords @ basis, basis

def apply_radius_cut(data_dict, r_extract):
    """
    Filters out particles located outside a specified spherical radius.

    Parameters
    ----------
    data_dict : dict
        Dictionary of particle arrays (must contain 'Coordinates').
    r_extract : float
        The cutoff radius in parsecs.

    Returns
    -------
    dict
        The dictionary containing only particles within the radial bound.
    """
    coords = data_dict['Coordinates']
    # Calculate the squared distance from the origin for computational efficiency
    radii_sq = np.sum(coords**2, axis=1)
    # Create a boolean mask for particles strictly inside the extraction radius
    mask = radii_sq < r_extract**2

    # Apply the spatial mask to all physical properties in the dictionary.
    return {key: val[mask] for key, val in data_dict.items()}

def finalize_dataset(header, pt0_dict, pt5_dict, percentage):
    """
    Applies centering, principal axis rotation, and radius cuts to both Medium (Gas) and 
    Source (Sink) particles.

    Parameters
    ----------
    header : dict
        Metadata extracted from the HDF5 file.
    pt0_dict : dict
        Raw gas particle data.
    pt5_dict : dict
        Raw sink/star particle data.
    percentage : float
        The fraction of the box half-width to extract.

    Returns
    -------
    pt0_final : dict
        Transformed and filtered gas data.
    pt5_final : dict
        Transformed and filtered star data.
    center : np.ndarray
        The translation vector used for centering.
    basis : np.ndarray
        The rotation matrix used for orientation.
    """
    # Define the spherical extraction boundary based on the simulation box size.
    box_size = header['BoxSize']
    r_extract = percentage * (box_size / 2)

    # Transform the Medium (PartType0) to establish the primary coordinate system.
    centered_gas, center = center_on_origin(pt0_dict['Coordinates'])
    rotated_gas, basis = apply_rotation(centered_gas, masses=pt0_dict['Masses'])

    # Update gas dictionary with transformed coordinates and apply radial mask.
    pt0_dict['Coordinates'] = rotated_gas
    pt0_final = apply_radius_cut(pt0_dict, r_extract)

    # Transform the Sources (PartType5) using the Gas coordinate system.
    centered_sinks, _ = center_on_origin(pt5_dict['Coordinates'], center=center)
    rotated_sinks, _ = apply_rotation(centered_sinks, basis=basis)

    # Update the star disctionary with transformed coordinates and apply radial mask.
    pt5_dict['Coordinates'] = rotated_sinks
    pt5_final = apply_radius_cut(pt5_dict, r_extract)

    return pt0_final, pt5_final, center, basis