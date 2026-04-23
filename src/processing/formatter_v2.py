"""
formatter_v2.py

Physical derivation and SKIRT table formatting utilities for the S2S pipeline.
"""
from __future__ import annotations

import numpy as np

import constants as c

def get_intrinsic_luminosity(
    pt5_data: dict,
    f_acc: float = 0.75
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute intrinsic and accretion luminosities for sink particles.

    Parameters
    ----------
    pt5_data: dict
        Sink particle data dictionary.
    f_acc: float, optional
        Fraction of accretion energy radiated away.

    Returns
    -------
    tuple
        intrinsic_luminosity, accretion_luminosity
    """
    stellar_mass_kg = pt5_data["BH_Mass"].astype(np.float64) * c.M_sun
    accretion_rate_kg_s = pt5_data["BH_Mdot"].astype(np.float64) * c.M_sun / c.yr_to_s
    stellar_radius_m = pt5_data["ProtoStellarRadius_inSolar"].astype(np.float64) * c.R_sun
    total_luminosity = pt5_data["StarLuminosity_Solar"].astype(np.float64)

    accretion_luminosity = np.zeros_like(total_luminosity)

    valid_radius = stellar_radius_m > 0
    if np.any(valid_radius):
        accretion_luminosity_watts = (
            f_acc * c.G * stellar_mass_kg[valid_radius] * accretion_rate_kg_s[valid_radius] / stellar_radius_m[valid_radius]
        )
        accretion_luminosity[valid_radius] = accretion_luminosity_watts / c.L_sun

    intrinsic_luminosity = np.clip(total_luminosity - accretion_luminosity, 0, None)

    return intrinsic_luminosity, accretion_luminosity

def compute_stellar_temperature(
    luminosity: np.ndarray,
    star_radius: np.ndarray
) -> np.ndarray:
    """
    Compute effective stellar temperatures from luminosity and radius.

    Parameters
    ----------
    luminosity: np.ndarray
        Stellar luminosity in solar luminosity units.
    star_radius: np.ndarray
        Stellar radius in solar radius units.

    Returns
    -------
    np.ndarray
        Effective temperature in Kelvin.
    """
    temperature = np.zeros_like(luminosity, dtype=np.float64)

    valid_radius = star_radius > 0
    if np.any(valid_radius):
        temperature[valid_radius] = c.T_sun * (luminosity[valid_radius] / star_radius[valid_radius]**2)**0.25

    return temperature


def build_source_table(
    pt5_data: dict,
    luminosity: np.ndarray,
    label: str
) -> tuple[np.ndarray, str, dict]:
    """
    Build a SKIRT-compatible source table for either intrinsic or accretion
    emission.

    Parameters
    ----------
    pt5_data: dict
        Sink particle data dictionary.
    luminosity: np.ndarray
        Luminosity array to export.
    label: str
        Descriptive label for the table type.

    Returns
    -------
    tuple
        export_data, header, source_stats
    """
    valid_mask = (
        (pt5_data["ProtoStellarRadius_inSolar"] > 0) & 
        (luminosity > 0)
    )

    clean_luminosity = luminosity[valid_mask]
    clean_radii = pt5_data["ProtoStellarRadius_inSolar"][valid_mask]
    clean_coords = pt5_data["Coordinates"][valid_mask]
    clean_smoothing_length = pt5_data["BH_AccretionLength"][valid_mask]

    temperature = compute_stellar_temperature(clean_luminosity, clean_radii)
    radius_km = clean_radii * c.R_sun * c.m_to_km

    export_data = np.column_stack([
        clean_coords[:, 0],
        clean_coords[:, 1],
        clean_coords[:, 2],
        clean_smoothing_length,
        radius_km,
        temperature
    ])

    header = "# x(pc) y(pc) z(pc) h(pc) R(km) T(K)"

    source_stats = {
        "label": label,
        "original_count": len(pt5_data["ProtoStellarRadius_inSolar"]),
        "exported_count": len(clean_luminosity),
        "dropped_count": len(pt5_data["ProtoStellarRadius_inSolar"]) - len(clean_luminosity)
    }

    return export_data, header, source_stats

def build_gas_table(pt0_data: dict) -> tuple[np.ndarray, str, tuple[float, float, float, float, float, float]]:
    """
    Build a SKIRT-compatible gas table and derive the spatial bounds.

    Parameters
    ----------
    pt0_data: dict
        Gas particle data dictionary.

    Returns
    -------
    tuple
        export_data, header, bounds
    """
    coords = pt0_data["Coordinates"]

    export_data = np.column_stack([
        coords[:, 0],
        coords[:, 1],
        coords[:, 2],
        pt0_data["SmoothingLength"],
        pt0_data["Masses"],
        pt0_data["Temperature"]
    ])

    bounds = (
        float(coords[:, 0].min()),
        float(coords[:, 0].max()),
        float(coords[:, 1].min()),
        float(coords[:, 1].max()),
        float(coords[:, 2].min()),
        float(coords[:, 2].max())
    )

    header = "# x(pc) y(pc) z(pc) h(pc) M(Msun) T(K)"

    return export_data, header, bounds


def write_table_file(
    output_path: str,
    export_data: np.ndarray,
    header: str
) -> str:
    """
    Write a formatted data table to disk.

    Parameters
    ----------
    output_path: str
        Path to the output text file.
    export_data: np.ndarray
        Table data to write.
    header: str
        Header string for the output file.

    Returns
    -------
    str
        Output path.
    """
    np.savetxt(output_path, export_data, fmt="%.6e", header=header, comments="")
    return output_path