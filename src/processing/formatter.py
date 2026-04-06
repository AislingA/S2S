# src/processing/formatter.py

"""
Module: formatter.py
Description: 
    This module handles the extraction, physical derivation, and formatting of 
    STARFORGE simulation data into SKIRT9-compatible text files. It computes 
    intrinsic stellar luminosities by subtracting accretion energy, calculates 
    effective surface temperatures using the Stefan-Boltzmann law, and exports 
    both the continuous gas medium and discrete stellar sources into structured 
    column formats.
"""

import numpy as np
import constants as c

def get_intrinsic_luminosity(pt5_data, f_acc=0.75):
    """
    Calculates the intrinsic stellar luminosity by subtracting the 
    accretion luminosity from the total bolometric luminosity tracked by STARFORGE.

    The accretion luminosity is derived from the gravitational potential energy 
    released by infalling matter:
    $ L_{\\rm acc} = f_{\\rm acc} \\frac{G M_* \\dot{M}}{R_*} $
    
    Reference: Offner et al. (2009); Grudić et al. (2021). The parameter f_acc 
    represents the fraction of accretion energy radiated away, assuming the 
    remaining 25% drives protostellar outflows/jets.

    Parameters
    ----------
    pt5_data: dict
        Dictionary containing the sink particle properties.
    f_acc: float, optional
        Fraction of accretion energy radiated. Default is 0.75.

    Returns
    -------
    np.ndarray
        Array of intrinsic stellar luminosities in units of Solar Luminosity (Lsun).
        Values are clipped to a minimum of 0 to precent unphysical negative luminosities.
    """
    M_star_kg = pt5_data['BH_Mass'].astype(np.float64) * c.M_sun
    # Mass accretion rate must be Mass/Time
    M_dot_kg_s = pt5_data['BH_Mdot'].astype(np.float64) * c.M_sun / c.yr_to_s
    R_star_m = pt5_data['ProtoStellarRadius_inSolar'] * c.R_sun
    L_tot = pt5_data['StarLuminosity_Solar']
    
    L_acc_Lsun = np.zeros_like(L_tot)
    # Create a mask to prevent division by zero for zero-radius sinks
    valid = pt5_data['ProtoStellarRadius_inSolar'] > 0
    
    # Calculate L_acc in Watts, then convert to Solar Luminosities
    L_acc_W = f_acc * c.G * M_star_kg[valid] * M_dot_kg_s[valid] / R_star_m[valid]
    L_acc_Lsun[valid] = L_acc_W / c.L_sun
    
    # Subtract accretion luminosity from total to get intrinsic, floor at 0
    return np.clip(L_tot - L_acc_Lsun, 0, None)

def compute_stellar_temperature(luminosity, star_radius):
    """
    Computes the surface effective temperature of a star using the 
    Stefan-Boltzmann Law, scaled relative to the Sun:
    $ \\frac{L}{L_\\odot} = \\left(\\frac{R}{R_\\odot}\\right)^2 \\left(\\frac{T}{T_\\odot}\\right)^4 $

    Parameters
    ----------
    luminosity: np.ndarray
        Stellar luminosity in units of Lsun.
    star_radius: np.ndarray
        Stellar radius in units of Rsun.

    Returns
    -------
    np.ndarray
        Effective surface temperature in Kelvin (K).
    """
    # Isolate T and solve using standard Solar scaling
    return c.T_sun * (luminosity / star_radius**2)**0.25

def print_stats(name, data_dict, verbose=False):
    """
    Utility function to print the dynamic range of datasets for verification.

    Parameters
    ----------
    name: str
        Identifier for the data group ("Source" or "Medium").
    data_dict: dict
        The dictionary of arrays to summarize.
    verbose: bool, optional
        If False, the function returns without printing. Default is False.
    """
    if not verbose:
        return
    print(f'--- {name} Data min and max values ---')
    for key, val in data_dict.items():
        # Only print stats for populated NumPy arrays
        if isinstance(val, np.ndarray) and val.ndim > 0:
            print(f'{key}: Min {np.min(val):.2e}, Max {np.max(val):.2e}')

def format_source_file(pt5_data, output_path, verbose=False):
    """
    Generates a SKIRT-compatible text file for stellar sources (PartType5).

    Extracts spatial coordinates, smoothing lengths, radii, and derives the
    effective temperature to assign proper SEDs in SKIRT.

    Parameters
    ----------
    pt5_data: dict
        Cleaned and transformed sink/star particle data.
    output_path: str
        The file path where the .txt source file will be saved.
    verbose: bool, optional
        Whether to print statistical summaries of the exported data.

    Returns
    -------
    str
        The path to the created source input file.
    """
    intrinsic_lums = get_intrinsic_luminosity(pt5_data)
    # Filter out sinks with zero radius or zero intrinsic luminosity
    valid_mask = (pt5_data['ProtoStellarRadius_inSolar'] > 0) & (intrinsic_lums > 0)
    
    clean_lums = intrinsic_lums[valid_mask]
    clean_radii = pt5_data['ProtoStellarRadius_inSolar'][valid_mask]
    clean_coords = pt5_data['Coordinates'][valid_mask]
    clean_h = pt5_data['BH_AccretionLength'][valid_mask]

    if verbose and len(clean_lums) < len(pt5_data['StarLuminosity_Solar']):
        print(f"Warning: Dropped {len(pt5_data['StarLuminosity_Solar']) - len(clean_lums)} sinks with 0 radius/luminosity to prevent SKIRT crash.")

    # Deriving effective temperature for SED assignment.
    temp = compute_stellar_temperature(clean_lums, clean_radii)

    # Structure data for SKIRT: x, y, z, h, R, T
    # Converting solar radii to kilometers (SKIRT standard unit for radius)
    radius_km = clean_radii * c.R_sun * c.m_to_km

    export_data = np.column_stack([
        clean_coords[:, 0],
        clean_coords[:, 1],
        clean_coords[:, 2],
        clean_h,
        radius_km,
        temp
    ])

    print_stats("Source", pt5_data, verbose)

    # Write to file with SKIRT-compatible header
    header = "# x(pc) y(pc) z(pc) h(pc) R(km) T(K)"
    np.savetxt(output_path, export_data, fmt='%.6e', header=header, comments='')
    
    return output_path

def format_gas_file(pt0_data, output_path, verbose=False):
    """
    Generates a SKIRT-compatible text file for the gas medium (PartType0).

    Parameters
    ----------
    pt0_data: dict
        Cleaned and transformed gas particle data.
    output_path: str
        The file path where the .txt gas file will be saved.
    verbose: bool, optional
        Whether to print statistical summaries of the exported data.

    Returns
    -------
    tuple
        (output_path, xmin, xmax, ymin, ymax, zmin, zmax)
        The path to the created file and the spatial bounds of the gas 
        distribution in parsecs (required for SKIRT spatial grid setup).
    """
    # Structure data for SKIRT: x, y, z, h, M, T
    export_data = np.column_stack([
        pt0_data['Coordinates'][:, 0],
        pt0_data['Coordinates'][:, 1],
        pt0_data['Coordinates'][:, 2],
        pt0_data['SmoothingLength'],
        pt0_data['Masses'],
        pt0_data['Temperature']
    ])

    print_stats("Medium", pt0_data, verbose)

    # Calculate spatial bounds to dynamically define the SKIRT domain.
    coords = pt0_data['Coordinates']
    bounds = (
        coords[:, 0].min(), coords[:, 0].max(),
        coords[:, 1].min(), coords[:, 1].max(),
        coords[:, 2].min(), coords[:, 2].max()
        )

    # Write to file with SKIRT-compatible header
    header = "# x(pc) y(pc) z(pc) h(pc) M(Msun) T(K)"
    np.savetxt(output_path, export_data, fmt='%.6e', header=header, comments='')
    
    return (output_path, *bounds)