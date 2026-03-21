# loader.py

# This module serves as the primary interface for raw STARFORGE snapshot 
# data. Its goal is to handle the low-level HDF5 file I/O, allowing the rest of 
# the pipeline to work with clean Python dictionaries. It extracts metadata 
# (Header), raw particle data, and provides initial ID-based cleaning to 
# remove boundary or non-physical particles.

# imports
import h5py
import numpy as np

def load_snapshot(snapshot_path):
    """
    Opens an HDF5 snapshot file in read-only mode using h5py.

    Parameters
    ----------
    snapshot_path : str
        The absolute or relative file path to the STARFORGE .hdf5 snapshot.

    Returns
    -------
    f : h5py.File
        The opened HDF5 file object.
    """
    return h5py.File(snapshot_path, 'r')

def get_header_data(file_obj):
    """
    Extracts simulation metadata from the snapshot 'Header' group.

    Parameters
    ----------
    file_obj : h5py.File
        The open HDF5 file object returned by load_snapshot.

    Returns
    -------
    dict
        A dictionary containing header attributes.
    """
    return dict(file_obj['Header'].attrs.items())

def get_particle_data(file_obj, part_type):
    """
    Reads all available datasets for a specific GIZMO/STARFORGE particle type.

    Parameters
    ----------
    file_obj : h5py.File
        The open HDF5 file object.
    part_type : int
        The GIZMO particle type index (0 for Gas/Medium, 5 for Sinks/Stars).

    Returns
    -------
    dict
        A dictionary where keys are dataset names ('Coordinates', 'Masses') 
        and values are the corresponding NumPy arrays. Returns an empty dict 
        if the part_type is not present in the file.
    """
    group_name = f'PartType{part_type}'
    if group_name not in file_obj:
        return {}
    
    return {key: file_obj[group_name][key][()] for key in file_obj[group_name].keys()}

def filter_by_id(pt_data, id_threshold=int(1e7)):
    """
    Filters out particles with IDs exceeding a specified threshold.

    Parameters
    ----------
    pt_data : dict
        A dictionary of particle arrays (must include 'ParticleIDs').
    id_threshold : int, optional
        The maximum allowable ParticleID. Particles with IDs >= this value 
        are discarded. Default is 1e7.

    Returns
    -------
    dict
        A dictionary containing the filtered arrays, maintaining the original 
        dictionary structure.
    """
    if 'ParticleIDs' not in pt_data:
        return pt_data
    
    mask = pt_data['ParticleIDs'] < id_threshold

    # applying mask to every dataset in this particle type
    filtered_data = {key: val[mask] for key, val in pt_data.items()}

    print(f"Filtered {len(pt_data['ParticleIDs'])} down to {len(filtered_data['ParticleIDs'])} particles.")
    return filtered_data