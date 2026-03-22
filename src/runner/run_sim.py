# src/runner/run_sim.py

import time
import os
import argparse  # Added missing import
import PTS9.simulation as sm

from file_io.loader import load_snapshot, get_header_data, get_particle_data, filter_by_id
from geometry.transform import finalize_dataset
from processing.formatter import format_source_file, format_gas_file
from config.writer import get_default_replacements, apply_yaml_replacements

def run_pipeline(snapshot_path, out_dir, percentage=1.0, verbose=True):
    """
    Executes the full end-to-end pipeline from raw snapshot to SKIRT output.
    """
    start_time = time.time()

    # Create the output directory if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)

    if not os.path.exists(snapshot_path):
        print(f"Error: Snapshot {snapshot_path} not found.")
        return None
    
    try:
        print(f"Loading snapshot: {snapshot_path}")
        with load_snapshot(snapshot_path) as f:
            header = get_header_data(f)
            raw_gas = get_particle_data(f, part_type=0)
            raw_sinks = get_particle_data(f, part_type=5)

        print("Filtering boundary particles...")
        gas_clean = filter_by_id(raw_gas)
        sinks_clean = filter_by_id(raw_sinks)

        print("Applying geometric transformations...")
        pt0, pt5, center, basis = finalize_dataset(header, gas_clean, sinks_clean, percentage)

        print("Formatting text inputs for SKIRT...")
        base_name = os.path.basename(snapshot_path).replace('.hdf5', '')
        
        # out_dir is now properly scoped
        src_path = os.path.join(out_dir, f"{base_name}_src.txt")
        gas_path = os.path.join(out_dir, f"{base_name}_gas.txt")
        ski_output = os.path.join(out_dir, f"{base_name}.ski")

        format_source_file(pt5, src_path, verbose=verbose)
        
        gas_info = format_gas_file(pt0, gas_path, verbose=verbose)
        gas_file, *bounds = gas_info

        print("Generating SKIRT configuration file...")
        base_dir = os.path.dirname(os.path.abspath(__file__)) 
        src_dir = os.path.dirname(base_dir) 
        
        template_path = os.path.join(src_dir, 'config', 'template.ski')
        yaml_path = os.path.join(src_dir, 'config', 'replacements.yaml')

        replacements = get_default_replacements(src_path, gas_path, bounds)
        apply_yaml_replacements(template_path, ski_output, replacements, yaml_path)
        
        print('\n--- Executing SKIRT Simulation ---')
        skirt = sm.Skirt()
        sim = skirt.execute(ski_output, console='brief')

        elapsed = time.time() - start_time
        print(f"\nPipeline complete. Execution time: {elapsed:.2f} seconds")
        return sim

    except Exception as e:
        print(f"Pipeline failed with error: {e}")
        return None
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run SKIRT pipeline on a snapshot.")
    
    # Updated to match the YAML flags exactly
    parser.add_argument("--input", required=True, help="Path to the raw STARFORGE .hdf5 snapshot.")
    parser.add_argument("--output-dir", required=True, help="Directory to save all output files.")
    
    args = parser.parse_args()

    # Pass the arguments in the correct order based on the new function signature
    run_pipeline(snapshot_path=args.input, out_dir=args.output_dir, percentage=1.0)