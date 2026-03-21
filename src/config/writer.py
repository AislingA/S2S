# src/config/writer

# imports
import yaml
import os

def load_config(yaml_path):
    with open(yaml_path, 'r') as f:
        return yaml.safe_load(f)

def get_default_replacements(src_file, gas_file, bounds):
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin
    max_fov = max(dx, dy, dz)
    return {
        "SOURCEFILE": src_file,
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
    config_data = load_config(yaml_path)
    final_replacements = runtime_replacements.copy()

    if config_data:
        for section, settings in config_data.items():
            if isinstance(settings, dict):
                final_replacements.update(settings)

    with open(template_path, 'r') as f:
        content = f.read()

    for placeholder, value in final_replacements.items():
        content = content.replace(placeholder, str(value))

    with open(output_path, 'w') as f:
        f.write(content)

    print(f"Config written to: {output_path}")
    return output_path