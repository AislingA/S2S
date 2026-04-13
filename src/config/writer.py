# src/config/writer.py

"""
Module: writer.py
Description:
    This module programmatically generates SKIRT9 XML configuration files (.ski).
    It merges dynamic spatial boundaries (derived from the HDF5 data) with 
    user-defined physical and numerical parameters stored in `replacements.yaml`.
    This abstraction ensures the user never has to manually edit the XML file.
"""

import yaml

def load_config(yaml_path):
    """
    Loads the user-defined configuration from a YAML file.

    Parameters
    ----------
    yaml_path: str
        The path to the replacements.yaml file.

    Returns
    -------
    dict
        A nested dictionary containing the configuration settings.
    """
    with open(yaml_path, 'r') as f:
        return yaml.safe_load(f)

def get_default_replacements(src_file, acc_file, gas_file, bounds):
    """
    Calculates the spatial bounding box and field of view for SKIRT instruments 
    based on the transformed STARFORGE particle data.

    Parameters
    ----------
    src_file: str
        Path to the formatted stellar source text file.
    acc_file: str
        Path to the formatted accretion source text file.
    gas_file: str
        Path to the formatted gas medium text file.
    bounds: tuple
        A tuple of floats representing (xmin, xmax, ymin, ymax, zmin, zmax) in parsecs.

    Returns
    -------
    dict
        A dictionary mapping spatial placeholders to their computed values.
    """
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    
    # Calculate the total physical extent in each dimension
    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin
    
    # The instrument Field of View (FOV) must encompass the largest dimension 
    # to ensure no part of the cloud is cut off in the mock images.
    max_fov = max(dx, dy, dz)
    
    return {
        "SOURCEFILE": src_file,
        "ACCRETIONFILE": acc_file,
        "GASFILE": gas_file,
        "XMIN": f"{xmin} pc",
        "XMAX": f"{xmax} pc",
        "YMIN": f"{ymin} pc",
        "YMAX": f"{ymax} pc",
        "ZMIN": f"{zmin} pc",
        "ZMAX": f"{zmax} pc",
        "FOVX": f"{max_fov} pc",
        "FOVY": f"{max_fov} pc"
    }

def apply_yaml_replacements(template_path, output_path, runtime_replacements, yaml_path='replacements.yaml'):
    """
    Injects dynamic data and user YAML configurations into the SKIRT XML template.

    Parameters
    ----------
    template_path: str
        Path to the base template.ski file containing placeholders.
    output_path: str
        Path where the final, executable .ski file will be saved.
    runtime_replacements: dict
        Dynamically generated spatial variables from `get_default_replacements`.
    yaml_path: str, optional
        Path to the replacements.yaml file. Default is 'replacements.yaml'.

    Returns
    -------
    str
        The file path to the completed SKIRT configuration file.
    """
    config_data = load_config(yaml_path)
    final_replacements = runtime_replacements.copy()

    # Flatten the nested YAML structure so all placeholder keys are at the top level
    if config_data:
        for section, settings in config_data.items():
            if isinstance(settings, dict):
                final_replacements.update(settings)

    # Read the raw XML template
    with open(template_path, 'r') as f:
        content = f.read()

    # Find and replace all placeholders with their respective values
    for placeholder, value in final_replacements.items():
        content = content.replace(placeholder, str(value))

    # Write the ready-to-run simulation file
    with open(output_path, 'w') as f:
        f.write(content)

    print(f"Config written successfully to: {output_path}")
    return output_path