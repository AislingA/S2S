"""
run_sim_v2.py

This version keeps the workflow readable for research use while separating the
main pipeline stages into small functions that are easier to test, debug, and
reuse.
"""

from __future__ import annotations

import argparse
import os
import pathlib
import sys
import time
from dataclasses import dataclass
from typing import Any

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from file_io.loader_v2 import load_snapshot, get_header_data, get_particle_data, filter_by_id
from geometry.transform_v2 import finalize_dataset
from processing.formatter_v2 import (
    get_intrinsic_luminosity,
    build_source_table,
    build_gas_table,
    write_table_file
)
from config.writer_v2 import get_default_replacements, render_template

@dataclass
class PipelineResult:
    """
    Stores the outcome of one pipeline run.

    This gives the caller more information than returning only a SKIRT
    simulation objects or None.
    """

    success: bool
    message: str
    elapsed_seconds: float
    simulation: Any = None
    output_files: dict[str, str] | None = None
    bounds: tuple[float, float, float, float, float, float] | None = None

def prepare_snapshot_data(snapshot_path: str, percentage: float) -> tuple[dict, dict, dict]:
    """
    Load the snapshot, filter by non-physical particles, and apply the geometric
    transformations needed before SKIRT input formatting.

    Returns
    -------
    tuple
        header, transformed_gas_data, transformed_sink_data
    """
    with load_snapshot(snapshot_path) as f:
        header = get_header_data(f)
        raw_gas = get_particle_data(f, part_type=0)
        raw_sinks = get_particle_data(f, part_type=5)

    gas_clean, _ = filter_by_id(raw_gas)
    sinks_clean, _ = filter_by_id(raw_sinks)

    pt0_final, pt5_final, _, _ = finalize_dataset(
        header,
        gas_clean,
        sinks_clean,
        percentage
    )

    return header, pt0_final, pt5_final

def build_skirt_inputs(
        snapshot_path: str,
        out_dir: str,
        pt0_data: dict,
        pt5_data: dict
) -> tuple[dict[str, str], tuple[float, float, float, float, float, float]]:
    """
    Write the SKIRT source/gas text files and generate the final .ski file.

    Returns
    -------
    tuple
        output_files, bounds
    """
    base_name = os.path.basename(snapshot_path).replace(".hdf5", "")

    src_path = os.path.join(out_dir, f"{base_name}_src.txt")
    acc_path = os.path.join(out_dir, f"{base_name}_acc.txt")
    gas_path = os.path.join(out_dir, f"{base_name}_gas.txt")
    ski_path = os.path.join(out_dir, f"{base_name}.ski")

    intrinsic_luminosity, accretion_luminosity = get_intrinsic_luminosity(pt5_data)

    src_data, src_header, _ = build_source_table(
        pt5_data,
        intrinsic_luminosity,
        label="intrinsic",
    )

    acc_data, acc_header, _ = build_source_table(
        pt5_data,
        accretion_luminosity,
        label="accretion",
    )

    gas_data, gas_header, bounds = build_gas_table(pt0_data)

    write_table_file(src_path, src_data, src_header)
    write_table_file(acc_path, acc_data, acc_header)
    write_table_file(gas_path, gas_data, gas_header)

    base_dir = os.path.dirname(os.path.abspath(__file__))
    src_dir = os.path.dirname(base_dir)

    template_path = os.path.join(src_dir, "config", "template.ski")
    yaml_path = os.path.join(src_dir, "config", "replacements.yaml")

    replacements = get_default_replacements(src_path, acc_path, gas_path, bounds, yaml_path)
    render_template(template_path, ski_path, replacements, yaml_path)

    output_files = {
        "source_file": src_path,
        "accretion_file": acc_path,
        "gas_file": gas_path,
        "ski_file": ski_path
    }

    return output_files, bounds

def run_skirt(ski_path: str, out_dir: str):
    """
    Execute the SKIRT simulation using PTS9 and return the simulation object.
    """
    import PTS9.simulation as sm

    skirt = sm.Skirt()

    original_dir = os.getcwd()
    os.chdir(out_dir)

    try:
        sim = skirt.execute(ski_path, console="brief")
    finally:
        os.chdir(original_dir)

    return sim


def run_pipeline(
        snapshot_path: str,
        out_dir: str,
        percentage: float = 1.0
) -> PipelineResult:
    """
    Run the full S2S pipeline from raw STARFORGE snapshot to SKIRT execution.
    """
    start_time = time.time()

    os.makedirs(out_dir, exist_ok=True)

    if not os.path.exists(snapshot_path):
        elapsed = time.time() - start_time
        return PipelineResult(
            success=False,
            message=f"Snapshot not found: {snapshot_path}",
            elapsed_seconds=elapsed
        )
    
    try:
        _, pt0_data, pt5_data = prepare_snapshot_data(snapshot_path, percentage)

        output_files, bounds = build_skirt_inputs(
            snapshot_path=snapshot_path,
            out_dir=out_dir,
            pt0_data=pt0_data,
            pt5_data=pt5_data
        )

        sim = run_skirt(output_files["ski_file"], out_dir)

        elapsed = time.time() - start_time
        return PipelineResult(
            success=True,
            message="Pipeline completed successfully.",
            elapsed_seconds=elapsed,
            simulation=sim,
            output_files=output_files,
            bounds=bounds
        )
    
    except Exception as e:
        elapsed = time.time() - start_time
        return PipelineResult(
            success=False,
            message=str(e),
            elapsed_seconds=elapsed
        )

def build_parser() -> argparse.ArgumentParser:
    """
    Build the command-line argument parser for the pipeline runner.
    """
    parser = argparse.ArgumentParser(
        description="Run the S2S pipeline on a STARFORGE snapshot."
    )
    parser.add_argument("--input", required=True, help="Path to the raw STARFORGE .hdf5 snapshot.")
    parser.add_argument("--output-dir", required=True, help="Directory to save all output files.")
    parser.add_argument(
        "--percentage",
        type=float,
        default=1.0,
        help="Fraction of the simulation box half-width to extract."
    )
    return parser

def main() -> None:
    """
    Parse command-line arguments and launch the pipeline.
    """
    parser = build_parser()
    args = parser.parse_args()

    result = run_pipeline(
        snapshot_path=args.input,
        out_dir=args.output_dir,
        percentage=args.percentage
    )

    if result.success:
        print(result.message)
    else:
        print(f"Pipeline failed: {result.message}")

if __name__ == "__main__":
    main()