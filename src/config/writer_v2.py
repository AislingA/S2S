"""
writer_v2.py

SKIRT configuration utilities for the S2S pipeline.
"""

from __future__ import annotations

import re


def load_config(yaml_path: str) -> dict:
    """
    Load the YAML configuration used to fill SKIRT template placeholders.

    Parameters
    ----------
    yaml_path: str
        Path to the YAML configuration file.

    Returns
    -------
    dict
        Parsed YAML configuration dictionary.
    """
    import yaml

    with open(yaml_path, "r") as f:
        config_data = yaml.safe_load(f)

    if config_data is None:
        return {}
    
    if not isinstance(config_data, dict):
        raise ValueError("YAML Configuration must contain a top-level dictionary.")

    return config_data


def flatten_run_config(config_data: dict) -> dict[str, object]:
    flattened: dict[str, object] = {}

    for section, settings in config_data.items():
        if not isinstance(settings, dict):
            raise ValueError(
                f"Configuration section '{section}' must contain a dictionary."
            )
        overlapping_keys = set(flattened).intersection(settings)
        if overlapping_keys:
            raise ValueError(
                "Duplicate configuration keys found across sections: "
                f"{', '.join(sorted(overlapping_keys))}"
            )
        flattened.update(settings)

    return flattened


def get_default_replacements(
    src_file: str,
    acc_file: str,
    gas_file: str,
    bounds: tuple[float, float, float, float, float, float]
) -> dict[str, str]:
    """
    Build runtime replacement values derived from the exported particle tables.

    Parameters
    ----------
    src_file: str
        Path to the intrinsic source file.
    acc_file: str
        Path to the accretion source file.
    gas_file: str
        Path to the gas medium file.
    bounds: tuple
        Spatial bounds as (xmin, xmax, ymin, ymax, zmin, zmax).

    Returns
    -------
    dict
        Replacement dictionary for dynamic placeholders.
    """
    xmin, xmax, ymin, ymax, zmin, zmax = bounds

    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin

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


def collect_template_placeholders(template_text: str) -> set[str]:
    """
    Collect uppercase placeholder tokens from the SKIRT template text.

    Parameters
    ----------
    template_text: str
        Raw template text.

    Returns
    -------
    set[str]
        Set of placeholder names found in the template.
    """
    placeholder_pattern = re.compile(r"\b[A-Z][A-Z0-9_]+\b")
    return set(placeholder_pattern.findall(template_text))


def render_template(
    template_path: str,
    output_path: str,
    runtime_replacements: dict[str, str],
    yaml_path: str
) -> str:
    """
    Render a final SKIRT .ski file from the template and replacement values.

    Parameters
    ----------
    template_path: str
        Path to the SKIRT template file.
    output_path: str
        Path where the rendered .ski file should be written.
    runtime_replacements: dict[str, str]
        Dynamically generated replacements from the pipeline run.
    yaml_path: str
        Path to the YAML config file.

    Returns
    -------
    str
        Output path of the rendered .ski file.
    """
    config_data = load_config(yaml_path)
    final_replacements = runtime_replacements.copy()

    for section, settings in config_data.items():
        if not isinstance(settings, dict):
            raise ValueError(f"Configuration section '{section}' must contain a dictionary of placeholder settings.")
        duplicate_keys = set(final_replacements).intersection(settings)
        if duplicate_keys:
            raise ValueError(
                "Duplicate replacement keys found while rendering template: "
                f"{', '.join(sorted(duplicate_keys))}"
            )
        final_replacements.update(settings)

    with open(template_path, "r") as f:
        template_text = f.read()

    template_placeholders = collect_template_placeholders(template_text)
    replacement_keys = set(final_replacements.keys())

    missing_placeholders = template_placeholders - replacement_keys
    unused_replacements = replacement_keys - template_placeholders

    if missing_placeholders:
        raise ValueError(f"Missing replacements for template placeholders: {', '.join(sorted(missing_placeholders))}")

    if unused_replacements:
        raise ValueError(f"Unused replacement keys not found in template: {', '.join(sorted(unused_replacements))}")

    rendered_text = template_text
    for placeholder, value in final_replacements.items():
        rendered_text = rendered_text.replace(placeholder, str(value))

    with open(output_path, "w") as f:
        f.write(rendered_text)

    return output_path
