# src/runner/run_sim.py

"""
Module: run_sim.py
Description: 
    The master execution script for the S2S (STARFORGE-to-SKIRT) pipeline. 
    It orchestrates the entire workflow: loading raw HDF5 data, filtering 
    unphysical particles, applying principal axis geometric transformations, 
    deriving physical stellar/dust properties, generating the SKIRT XML 
    configuration, and finally executing the SKIRT9 Monte Carlo binary via PTS9.
"""

import time
import os
import argparse
import PTS9.simulation as sm

# Internal pipeline imports
from file_io.loader import load_snapshot, get_header_data, get_particle_data, filter_by_id
from geometry.transform import finalize_dataset
from processing.formatter import format_source_file, format_gas_file
from config.writer import get_default_replacements, apply_yaml_replacements

def run_pipeline(snapshot_path, out_dir, percentage=1.0, verbose=True):
    """
    Executes the full end-to-end pipeline from raw snapshot to SKIRT output.

    Parameters
    ----------
    snapshot_path: str
        Absolute or relative path to the raw STARFORGE .hdf5 snapshot.
    out_dir: str
        Directory path where the formatted text files, the .ski config, 
        and the final SKIRT outputs will be saved.
    percentage: float, optional
        Fraction of the simulation box half-width to extract (defaults to 1.0).
    verbose: bool, optional
        If True, prints statistical diagnostics at each pipeline stage.

    Returns
    -------
    sim: PTS9.simulation.Simulation
        The completed PTS9 simulation object, or None if an error occurred.
    """
    start_time = time.time()

    # Initialize Output Environment
    os.makedirs(out_dir, exist_ok=True)
    if not os.path.exists(snapshot_path):
        print(f"Error: Snapshot {snapshot_path} not found.")
        return None
    
    try:
        # Data Ingestion (loader.py)
        # Using a context manager ('with') ensures the HDF5 file safely closes 
        # after data extraction, preventing memory leaks or file corruption.
        print(f"Loading snapshot: {snapshot_path}")
        with load_snapshot(snapshot_path) as f:
            header = get_header_data(f)
            raw_gas = get_particle_data(f, part_type=0)
            raw_sinks = get_particle_data(f, part_type=5)

        # Data Sanitization (loader.py)
        print("Filtering boundary particles...")
        gas_clean = filter_by_id(raw_gas)
        sinks_clean = filter_by_id(raw_sinks)

        # Geometric Normalization (transform.py)
        print("Applying geometric transformations...")
        pt0, pt5, center, basis = finalize_dataset(header, gas_clean, sinks_clean, percentage)

        # Physical Formatting & Export (formatter.py)
        print("Formatting text inputs for SKIRT...")
        base_name = os.path.basename(snapshot_path).replace('.hdf5', '')
        
        src_path = os.path.join(out_dir, f"{base_name}_src.txt")
        gas_path = os.path.join(out_dir, f"{base_name}_gas.txt")
        ski_output = os.path.join(out_dir, f"{base_name}.ski")

        # Export stellar sources with derived intrinsic luminosities/temperatures
        format_source_file(pt5, src_path, verbose=verbose)
        
        # Export gas medium and capture the spatial bounding box
        gas_info = format_gas_file(pt0, gas_path, verbose=verbose)
        
        # Unpack the returned tuple: first element is path, the rest are bounds
        gas_file, *bounds = gas_info

        # Configuration Generation (writer.py)
        print("Generating SKIRT configuration file...")
        # Dynamically resolve the absolute paths to the config directory
        base_dir = os.path.dirname(os.path.abspath(__file__)) 
        src_dir = os.path.dirname(base_dir) 
        
        template_path = os.path.join(src_dir, 'config', 'template.ski')
        yaml_path = os.path.join(src_dir, 'config', 'replacements.yaml')

        # Map dynamic bounds and apply user-defined physics from YAML
        replacements = get_default_replacements(src_path, gas_path, bounds)
        apply_yaml_replacements(template_path, ski_output, replacements, yaml_path)
        
        # Radiative Transfer Execution (PTS9)
        print('\n--- Executing SKIRT Simulation ---')
        skirt = sm.Skirt()
        # 'brief' console mode prevents the terminal from being overwhelmed by photon tracking logs
        sim = skirt.execute(ski_output, console='brief')

        elapsed = time.time() - start_time
        print(f"\nPipeline complete. Execution time: {elapsed:.2f} seconds")
        return sim

    except Exception as e:
        print(f"Pipeline failed with error: {e}")
        return None
    
if __name__ == "__main__":
    # Command-line interface for HPC batch scheduling
    parser = argparse.ArgumentParser(description="Run the S2S pipeline on a STARFORGE snapshot.")
    
    parser.add_argument("--input", required=True, help="Path to the raw STARFORGE .hdf5 snapshot.")
    parser.add_argument("--output-dir", required=True, help="Directory to save all output files.")
    
    args = parser.parse_args()

    run_pipeline(snapshot_path=args.input, out_dir=args.output_dir, percentage=1.0)