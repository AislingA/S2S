# src/geometry/transform.py

# To perform geometric transformations on simulation data. This includes 
# centering the region of interest at the origin (0,0,0), rotating the 
# coordinates to align with the principal axes of the mass distribution, 
# and performing radial cuts to define the final simulation volume.

# imports
import numpy as np

def center_on_origin(coords, center=None):
    """
    Translates coordinates so that the region of interest is at (0,0,0).

    Parameters
    ----------
    coords : np.ndarray
        Nx3 array of particle positions.
    center : np.ndarray, optional
        An external (x,y,z) coordinate to center on. If None, it calculates 
        the median of the provided coords.

    Returns
    -------
    centered_coords : np.ndarray
        The translated coordinates.
    center : np.ndarray
        The translation vector used (useful for centering other particle types).
    """
    if center is None:
        center = np.median(coords, axis=0)

    return coords - center, center

def apply_rotation(coords, masses=None, basis=None):
    """
    Rotates coordinates to align with a specific basis. If no basis is 
    provided, it calculates one using the mass-weighted covariance matrix.

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
        The coordinates transformed into the new basis.
    basis : np.ndarray
        The 3x3 rotation matrix used for the transformation.
    """
    if basis is None:
        # calculating the basis from principal axes
        cov_pos = np.cov(coords.T, aweights=masses)
        try:
            w, basis = np.linalg.eigh(cov_pos)
            basis = basis[:, w.argsort()[::-1]]
        except np.linalg.LinAlgError:
            print("Warning: Matrix did not converge. Using Identity.")
            basis = np.identity(3)

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
    radii_sq = np.sum(coords**2, axis=1)
    mask = radii_sq < r_extract**2

    return {key: val[mask] for key, val in data_dict.items()}

def finalize_dataset(header, pt0_dict, pt5_dict, percentage):
    """
    Applies centering, rotation, and radius cuts to both Medium (Gas) and 
    Source (Sink) particles.

    Parameters
    ----------
    header : dict
        Metadata from the snapshot.
    pt0_dict : dict
        Raw gas particle data.
    pt5_dict : dict
        Raw sink/star particle data.
    percentage : float
        The fraction of the box half-width (BoxSize/2) to extract.

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
    # setting up dimensions
    box_size = header['BoxSize']
    r_extract = percentage * (box_size / 2)

    # transforming mediums (PartType0) and establishing the coordinate system
    centered_gas, center = center_on_origin(pt0_dict['Coordinates'])
    rotated_gas, basis = apply_rotation(centered_gas, masses=pt0_dict['Masses'])

    # updating dict with transformed coords for the mask
    pt0_dict['Coordinates'] = rotated_gas
    pt0_final = apply_radius_cut(pt0_dict, r_extract)

    # transforming sources (PartType5) using the same system
    centered_sinks, _ = center_on_origin(pt5_dict['Coordinates'], center=center)
    rotated_sinks, _ = apply_rotation(centered_sinks, basis=basis)

    pt5_dict['Coordinates'] = rotated_sinks
    pt5_final = apply_radius_cut(pt5_dict, r_extract)

    return pt0_final, pt5_final, center, basis