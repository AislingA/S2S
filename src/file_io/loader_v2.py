"""
loader_v2.py

Snapshot I/O utilities for the S2S pipeline.
"""

from __future__ import annotations
from typing import Any

def load_snapshot(snapshot_path: str):
    """
    Open a STARFORGE snapshot in read-only mode.

    Parameters
    ----------
    snapshot_path: str
        Absolute or relative path to the STARFORGE .hdf5 snapshot.

    Returns
    -------
    h5py.File
        Open HDF5 file object.
    """
    import h5py

    return h5py.File(snapshot_path, "r")

def get_header_data(file_obj) -> dict:
    """
    Extract the snapshot header attributes into a standard Python dictionary.

    Parameters
    ----------
    file_obj: h5py.File
        Open HDF5 snapshot file,
    
    Returns
    -------
    dict
        Header attributes.
    """
    return dict(file_obj["Header"].attrs.items())

def validate_particle_data(
        particle_data: dict,
        required_fields: list[str] | None = None,
) -> dict:
    """
    Validate that particle arrays are internally consistent.

    Parameters
    ----------
    particle_data: dict
        Dictionary of particle property arrays.
    required_fields: list[str] | None, optional
        Field names that must be present.

    Returns
    -------
    dict
        The original particle data if validation passes.

    Raises
    ------
    KeyError
        If a required field is missing.
    ValueError
        If array lengths are inconsistent.
    """
    if required_fields is not None:
        missing_fields = [field for field in required_fields if field not in particle_data]
        if missing_fields:
            raise KeyError(f"Missing required particle fields: {', '.join(missing_fields)}")
        
    if not particle_data:
        return particle_data
    
    first_key = next(iter(particle_data))
    expected_length = len(particle_data[first_key])

    for key, value in particle_data.items():
        if len(value) != expected_length:
            raise ValueError(f" Particle field '{key}' has length {len(value)} but expected {expected_length}.")
        
    return particle_data

def get_particle_data(file_obj, part_type: int) -> dict:
    """
    Read all datasets for a given GIZMO particle type into memory.

    Parameters
    ----------
    file_obj: h5py.File
        Open HDF5 snapshot file.
    part_type: int
        GIZMO particle type index.

    Returns
    -------
    dict
        Dictionary of particle arrays. Returns an empty dict if the particle type is not present.
    """
    group_name = f"PartType{part_type}"

    if group_name not in file_obj:
        return {}
    
    particle_data = {
        key: file_obj[group_name][key][()] for key in file_obj[group_name].keys()
    }

    return validate_particle_data(particle_data)

def filter_by_id(
        particle_data: dict,
        id_threshold: int = int(1e7)        
) -> tuple[dict, dict[str, Any]]:
    """
    Remove particles with IDs greater than or equal to the threshold.

    Parameters
    ----------
    particle_data: dict
        Dictionary of particle property arrays.
    id_threshold: int, optional
        Maximum allowed particle ID.
    
    Returns
    -------
    tuple
        filtered particle data, filter stats
    """
    if not particle_data:
        return {}, {
            "original_count": 0,
            "filtered_count": 0,
            "id_threshold": id_threshold,
            "filter_applied": False
        }
    
    if "ParticleIDs" not in particle_data:
        return particle_data, {
            "original_count": len(next(iter(particle_data.values()))),
            "filtered_count": len(next(iter(particle_data.values()))),
            "id_threshold": id_threshold,
            "filter_applied": False
        }
    
    particle_data = validate_particle_data(
        particle_data,
        required_fields=["ParticleIDs"]
    )

    mask = particle_data["ParticleIDs"] < id_threshold

    filtered_particle_data = {
        key: value[mask] for key, value in particle_data.items()
    }

    filter_stats = {
        "original_count": len(particle_data["ParticleIDs"]),
        "filtered_count": len(filtered_particle_data["ParticleIDs"]),
        "id_threshold": id_threshold,
        "filter_applied": True
    }

    return filtered_particle_data, filter_stats