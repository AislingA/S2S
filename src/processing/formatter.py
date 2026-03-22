# src/processing/formatter.py

# imports
import numpy as np
import os
import constants as c

def get_intrinsic_luminosity(pt5_data, f_acc=0.75):
    """
    Calculates the intrinsic stellar luminosity by subtracting the 
    accretion luminosity from the total bolometric luminosity.

    source where eq and f acc came from 
    """
    M_star_kg = pt5_data['BH_Mass'] * c.M_sun
    M_dot_kg_s = pt5_data['BH_Mdot'] * c.M_sun / c.pc_to_m
    R_star_m = pt5_data['ProtoStellarRadius_inSolar'] * c.R_sun
    L_tot = pt5_data['StarLuminosity_Solar']
    
    L_acc_Lsun = np.zeros_like(L_tot)
    valid = pt5_data['ProtoStellarRadius_inSolar'] > 0
    
    L_acc_W = f_acc * c.G * M_star_kg[valid] * M_dot_kg_s[valid] / R_star_m[valid]
    L_acc_Lsun[valid] = L_acc_W / c.L_sun
    
    return np.clip(L_tot - L_acc_Lsun, 0, None)

def compute_stellar_temperature(luminosity, star_radius):
    """
    Computes the surface effective temperature of a star using the 
    Stefan-Boltzmann Law.

    Parameters
    ----------
    luminosity : np.ndarray
        Stellar luminosity in units of Lsun.
    star_radius : np.ndarray
        Stellar radius in units of Rsun.

    Returns
    -------
    np.ndarray
        Effective temperature in Kelvin (K).
    """
    return c.T_sun * (luminosity / star_radius**2)**0.25

def print_stats(name, data_dict, verbose=False):
    """
    Utility function to print the dynamic range of datasets for verification.

    Parameters
    ----------
    name : str
        Identifier for the data group ("Source" or "Medium").
    data_dict : dict
        The dictionary of arrays to summarize.
    verbose : bool, optional
        If False, the function returns without printing. Default is False.
    """
    if not verbose:
        return
    print(f'--- {name} Data min and max values ---')
    for key, val in data_dict.items():
        if isinstance(val, np.ndarray) and val.ndim > 0:
            print(f'{key}: Min {np.min(val):.2e}, Max {np.max(val):.2e}')

def format_source_file(pt5_data, output_path, verbose=False):
    """
    Generates a SKIRT-compatible text file for stellar sources (PartType5).

    Parameters
    ----------
    pt5_data : dict
        Cleaned and transformed sink/star particle data.
    output_path : str
        The file path where the .txt source file will be saved.
    verbose : bool, optional
        Whether to print statistical summaries of the exported data.

    Returns
    -------
    str
        The path to the created source input file.
    """
    intrinsic_lums = get_intrinsic_luminosity(pt5_data)
    valid_mask = (pt5_data['ProtoStellarRadius_inSolar'] > 0) & (intrinsic_lums > 0)
    
    clean_lums = intrinsic_lums[valid_mask]
    clean_radii = pt5_data['ProtoStellarRadius_inSolar'][valid_mask]
    clean_coords = pt5_data['Coordinates'][valid_mask]
    clean_h = pt5_data['BH_AccretionLength'][valid_mask]

    if verbose and len(clean_lums) < len(pt5_data['StarLuminosity_Solar']):
        print(f"Warning: Dropped {len(pt5_data['StarLuminosity_Solar']) - len(clean_lums)} sinks with 0 radius/luminosity to prevent SKIRT crash.")

    # Calculating T_eff
    temp = compute_stellar_temperature(clean_lums, clean_radii)

    # Structuring for SKIRT: x, y, z, h, R, T
    # Converting solar radii to kilometers
    radius_km = clean_radii * 695700

    export_data = np.column_stack([
        clean_coords[:, 0],
        clean_coords[:, 1],
        clean_coords[:, 2],
        clean_h,
        radius_km,
        temp
    ])

    print_stats("Source", pt5_data, verbose)

    header = "# x(pc) y(pc) z(pc) h(pc) R(km) T(K)"
    np.savetxt(output_path, export_data, fmt='%.6e', header=header, comments='')
    return output_path

def format_gas_file(pt0_data, output_path, verbose=False):
    """
    Generates a SKIRT-compatible text file for the gas medium (PartType0).

    Parameters
    ----------
    pt0_data : dict
        Cleaned and transformed gas particle data.
    output_path : str
        The file path where the .txt gas file will be saved.
    verbose : bool, optional
        Whether to print statistical summaries of the exported data.

    Returns
    -------
    tuple
        (str, float, float, float, float, float, float)
        The output path followed by the (xmin, xmax, ymin, ymax, zmin, zmax) 
        bounds of the gas distribution in parsecs.
    """
    # Structuring for SKIRT: x, y, z, h, M, T
    export_data = np.column_stack([
        pt0_data['Coordinates'][:, 0],
        pt0_data['Coordinates'][:, 1],
        pt0_data['Coordinates'][:, 2],
        pt0_data['SmoothingLength'],
        pt0_data['Masses'],
        pt0_data['Temperature']
    ])

    print_stats("Medium", pt0_data, verbose)

    # Calculating bounds for the SKIRT geometry setup (Fixed syntax errors)
    coords = pt0_data['Coordinates']
    bounds = (coords[:, 0].min(), coords[:, 0].max(),
              coords[:, 1].min(), coords[:, 1].max(),
              coords[:, 2].min(), coords[:, 2].max())

    header = "# x(pc) y(pc) z(pc) h(pc) M(Msun) T(K)"
    np.savetxt(output_path, export_data, fmt='%.6e', header=header, comments='')
    
    # Return output path and unpack the bounds tuple
    return (output_path, *bounds)