# src/file_io/loader.py

"""
Module: loader.py
Description: 
    This module serves as the primary I/O interface for raw STARFORGE (GIZMO) 
    snapshot data. It handles the low-level HDF5 file operations, converting 
    binary datasets into clean, accessible Python dictionaries. It extracts 
    global metadata, raw particle arrays (gas and sink particles), and 
    provides initial ID-based filtering to eliminate non-physical boundary 
    particles before radiative transfer post-processing.
"""

import h5py

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
    # Open the HDF5 file in read-only mode ('r') to prevent overwrites
    # of the raw simulation data.
    return h5py.File(snapshot_path, 'r')

def get_header_data(file_obj):
    """
    Extracts simulation metadata from the snapshot's 'Header' group.

    Parameters
    ----------
    file_obj : h5py.File
        The open HDF5 file object returned by load_snapshot.

    Returns
    -------
    dict
        A dictionary containing header attributes.
    """
    # Access the 'Header' group and extract its attributes (.attrs).
    # Cast the items to a standard Python dictionary for easier access.
    return dict(file_obj['Header'].attrs.items())

def get_particle_data(file_obj, part_type):
    """
    Reads all available datasets for a specific GIZMO/STARFORGE particle type into memory.

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
    # Construct the standard GIZMO group name for the requested particle type
    group_name = f'PartType{part_type}'
    # Check if the particel type exists in this specific snapshot to avoid KeyErrors
    if group_name not in file_obj:
        return {}
    # Iterate through all datasets within the group
    # The [()] syntax instructs the h5py to read the entire dataset from disk
    # into memory as an array.
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
    # If the dataset does not have ParticleIDs, we cannot filter it. Return as is.
    if 'ParticleIDs' not in pt_data:
        return pt_data
    
    # Create a boolean array where True indicates a valid, physical particle.
    mask = pt_data['ParticleIDs'] < id_threshold

    # Applying mask to every dataset in this particle type
    filtered_data = {key: val[mask] for key, val in pt_data.items()}

    # Ouput a brief log of the filteration results for diagnostic tracking.
    print(f"Filtered {len(pt_data['ParticleIDs'])} down to {len(filtered_data['ParticleIDs'])} particles.")
    return filtered_data