# src/constants.py

"""
Module: constants.py
Description: 
    A centralized repository for physical and astronomical constants used 
    throughout the S2S pipeline. All physical constants are strictly defined 
    in SI units (MKS system: meters, kilograms, seconds, Watts) to ensure 
    dimensional consistency during thermodynamic and radiative calculations.
    Solar nominal values follow standard astrophysical conventions (IAU 2015).
"""

# --- Fundamental Physics Constants ---
G = 6.67430e-11        # Gravitational constant [m^3 kg^-1 s^-2]

# --- Solar Nominal Values (IAU 2015 Standard) ---
M_sun = 1.9884e30      # Mass of the Sun [kg]
R_sun = 6.957e8        # Radius of the Sun [m]
L_sun = 3.828e26       # Nominal Solar Luminosity [W]
T_sun = 5777.0         # Effective Surface Temperature of the Sun [K]

# --- Unit Conversions ---
pc_to_m = 3.085677e16  # Parsecs to meters [m/pc]
yr_to_s = 3.15576e7    # Julian year to seconds [s/yr] (365.25 days)
m_to_km = 1.0e-3       # Meters to kilometers [km/m]