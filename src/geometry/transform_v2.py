"""
transform_v2.py

Geometric transforms for the S2S pipeline.
"""
from __future__ import annotations

import numpy as np

def center_on_origin(coords: np.ndarray, center: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """
    Shift coordinates so that the chosen center lies at the origin.

    Parameters
    ----------
    coords: np.ndarray
        Array of particle coordinates with shape (N, 3).
    center: np.ndarray | None, optional
        Coordinate to center on. If None, use the median position.

    Returns
    -------
    tuple
        centered_coords, center_used
    """
    if center is None:
            center = np.median(coords, axis=0)

    centered_coords = coords - center
    return centered_coords, center

def compute_rotation_basis(coords: np.ndarray, masses: np.ndarray | None = None) -> np.ndarray:
    """
    Compute the principal-axis rotation basis from the coordinate distribution.

    Parameters
    ----------
    coords: np.ndarray
        Centered coordinate array with shape (N, 3).
    masses: np.ndarray | None, optional
        Particle masses used for mass-weighted covariance.

    Returns
    -------
    np.ndarray
        Rotation basis matrix with shape (3, 3).
    """
    cov_pos = np.cov(coords.T, aweights=masses)

    try:
        eigenvalues, basis = np.linalg.eigh(cov_pos)
    except np.linalg.LinAlgError:
         return np.identity(3)
    
    basis = basis[:, eigenvalues.argsort()[::-1]]
    return basis

def apply_rotation(coords: np.ndarray, basis: np.ndarray) -> np.ndarray:
    """
    Rotate coordinates using a precomputed basis matrix.

    Parameters
    ----------
    coords: np.ndarray
        Centered coordinate array with shape (N, 3).
    basis: np.ndarray
        Rotation basis matrix.
    
    Returns
    -------
    np.ndarray
        Rotated coordinates.
    """
    return coords @ basis

def apply_radius_cut(particle_data:dict, r_extract: float) -> dict:
    """
    Keep only particles inside the spherical extraction radius.

    Parameters
    ----------
    particle_data: dict
        Dictionary of particle arrays. Must include 'Coordinates'.
    r_extract: float
        Spherical extraction radius in parsecs.

    Returns
    -------
    dict
        New particle-data dictionary containing only particles inside the radius.
    """
    coords = particle_data["Coordinates"]
    radii_sq = np.sum(coords**2, axis=1)
    mask = radii_sq < r_extract**2

    filtered_particle_data = {
         key: value[mask]
         for key, value in particle_data.items()
    }

    return filtered_particle_data

def transform_particle_data(
        particle_data: dict,
        center: np.ndarray,
        basis: np.ndarray,
        r_extract: float
) -> dict:
    """
    Apply centering, rotation, and spherical extraction to one particle dataset.

    Parameters
    ----------
    particle_data: dict
        Dictionary of particle arrays. Must include 'Coordinates'.
    center: np.ndarray
        Center coordinate to subtract.
    basis: np.ndarray
        Rotation basis matrix.
    r_extract: float
        Spherical extraction radius in parsecs.

    Returns
    -------
    dict
        Transformed and filtered particle-data dictionary.
    """
    centered_coords, _ = center_on_origin(particle_data["Coordinates"], center=center)
    rotated_coords = apply_rotation(centered_coords, basis)

    transformed_particle_data = {
         key: value.copy() if hasattr(value, "copy") else value
         for key, value in particle_data.items()
    }
    transformed_particle_data["Coordinates"] = rotated_coords

    return apply_radius_cut(transformed_particle_data, r_extract)

def finalize_dataset(
    header: dict,
    pt0_data: dict,
    pt5_data: dict,
    percentage: float
) -> tuple[dict, dict, np.ndarray, np.ndarray]:
    """
    Apply the full geometric preprocessing workflow to gas and sink particles.

    Parameters
    ----------
    header: dict
        Snapshot metadata dictionary.
    pt0_data: dict
        Gas particle data.
    pt5_data: dict
        Sink/star particle data.
    percentage: float
        Fraction of the box half-width used for the extraction radius.

    Returns
    -------
    tuple
        pt0_final, pt5_final, center, basis
    """
    box_size = header["BoxSize"]
    r_extract = percentage * (box_size / 2)

    centered_gas, center = center_on_origin(pt0_data["Coordinates"])
    basis = compute_rotation_basis(centered_gas, masses=pt0_data["Masses"])

    pt0_final = transform_particle_data(
         pt0_data,
         center=center,
         basis=basis,
         r_extract=r_extract
    )

    pt5_final = transform_particle_data(
         pt5_data,
         center=center,
         basis=basis,
         r_extract=r_extract
    )

    return pt0_final, pt5_final, center, basis