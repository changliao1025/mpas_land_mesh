#!/usr/bin/env python3
"""Yellow River (Huanghe) MPAS mesh workflow.

This example demonstrates using mpas_land_mesh's ULRMWorkflow to orchestrate
preprocessing and mesh generation from a single YAML/JSON configuration.
"""

import argparse
import logging
import os
from pathlib import Path

from mpas_land_mesh import ULRMWorkflow

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        "-c",
        default=None,
        help="Path to workflow YAML/JSON configuration file (default: workflow_config.yaml)",
    )
    parser.add_argument(
        "--preprocess-only",
        action="store_true",
        help="Run only the preprocessing step and compile mesh_ready_workflow.yaml",
    )
    args = parser.parse_args()

    # ── Step 1: Resolve configuration path ────────────────────────────────────
    if args.config:
        sFilename_workflow = Path(args.config).expanduser().resolve()
    else:
        sFilename_workflow = Path(__file__).with_name("workflow_config.yaml")
        if not sFilename_workflow.exists():
            sFilename_workflow = Path(__file__).with_name("workflow_config.json")

    logger.info("Step 1: Initializing ULRM workflow from: %s", sFilename_workflow)
    workflow = ULRMWorkflow(sFilename_workflow)

    # ── Step 2: Preprocessing & Mesh-Ready Compilation ───────────────────────
    logger.info("Step 2: Preprocessing datasets...")
    preprocessed = workflow.preprocess()
    logger.info("   * Output workspace: %s", workflow.output_workspace)

    logger.info("Step 3: Compiling mesh-ready configuration...")
    mesh_ready_cfg = workflow.compile_mesh_ready_config()
    mesh_ready_file = os.path.join(
        workflow.output_workspace, "mesh_ready_workflow.yaml"
    )
    logger.info("   * Mesh-ready specification written to: %s", mesh_ready_file)

    if args.preprocess_only:
        logger.info("Preprocess-only flag specified; halting before mesh generation.")
        return

    # ── Step 3: Run Mesh Generation or HPC Batch Script Generation ───────────
    logger.info("Step 4: Executing mesh workflow via pyearthmesh...")
    mesh_workflow = workflow.create_mesh_workflow()
    if workflow.config["workflow"].get("run_on_hpc"):
        hours = workflow.config.get("jigsaw", {}).get("hours", 5)
        logger.info("Generating HPC SLURM job scripts (walltime=%dh)...", hours)
        mesh_workflow.create_hpc_job(hours=hours)
        logger.info("SLURM job scripts created in: %s", mesh_workflow.output_workspace)
    else:
        logger.info("Setting up MPAS workspace...")
        mesh_workflow.setup()
        logger.info("Running MPAS mesh generation...")
        mesh_workflow.run()

    logger.info("Done.")


if __name__ == "__main__":
    main()
