# YAML Workflow Configuration & Two-Tier Compiler Architecture

## Summary

The workflow configuration for `mpas_land_mesh` adopts the unified **Two-Tier Compiler Architecture** established in [`pyearthmesh`](https://github.com/e3sm-project/pyearthmesh) (see reference specification: `pyearthmesh/plans/two_tier_compiler_design.md`).

A **single, human-readable, comment-annotated YAML file** (`workflow_config.yaml`) serves as the **primary source of truth** for an entire run. It completely replaces hardcoded script parameters, repetitive on-disk edits, and intermediate flat configuration decks (`mpas_configuration.json`, `jigsaw_configuration.json`), while maintaining full backward-compatibility with JSON (`workflow_config.json`).

The compiler boundary translates the high-level domain specification directly into engine-native parameters in memory, derives grid spacings and boundary geometries, and orchestrates either local mesh generation or SLURM HPC batch job generation.

---

## Architectural Layers

```text
User Intent (YAML with # comments / JSON)
        │
        ▼
Tier 1: Domain Specification (workflow_config.yaml / workflow_config.json)
  ├── workflow: run metadata, preprocessing flags, execution targets
  ├── resolutions: physical cell resolutions in kilometres
  ├── paths: multi-platform dataset and output workspaces
  └── jigsaw: engine-specific numerical knobs, spacing grid dimensions, SLURM hours/queues
        │
        ▼  [load_workflow_config() / WorkflowSpec]
Compiler Boundary: compile_workflow(spec_or_config)
  ├── Strips / parses comments seamlessly (via PyYAML)
  ├── Merges with DEFAULT_CONFIG defaults
  ├── Resolves relative paths against config file location
  ├── Resolves OS-specific paths (Darwin / Linux / Windows / default)
  ├── Derives pBoundary_wkt from GeoJSON boundary mask
  ├── Computes spherical grid spacing (nrow_space, ncolumn_space)
  └── Compiles native JIGSAW deck in memory
        │
        ▼  [create_workflow()]
Tier 2: Runtime Workflow Adapter (MPASWorkflow / BaseCase)
  ├── case.output_workspace -> resolved output directory
  ├── case.mesh_file        -> target GeoJSON/NetCDF mesh path
  ├── case.setup()          -> prepares directory structure and inputs
  ├── case.run()            -> executes generator, returns MPAS cells
  └── case.create_hpc_job() -> emits run script and SLURM submit.job
```

---

## Why YAML is the Primary Tier 1 Format

While standard JSON is effective for machine-to-machine exchange, it possesses major shortcomings for scientific modeling:
- **No Comments Allowed:** Standard JSON strictly prohibits comments (`#` or `//`). Modelers cannot explain parameter choices, document data sources, or explain resolution tradeoffs.
- **Syntax Rigidity:** Quoted keys, trailing commas, and bracket matching make manual editing error-prone.

By adopting **YAML** as the primary human-facing format:
1. **Inline & Block Comments (`# ...`):** Users can thoroughly document grid parameters, data provenance, and HPC queue selections.
2. **Cleaner Syntax:** Whitespace indentation eliminates bracket clutter and improves readability.
3. **Dual-Format Parity:** Standard JSON is a syntactic subset of YAML; `pyearthmesh` and `mpas_land_mesh` load both `.yaml`/`.yml` and `.json` identically.
4. **Directory & Extensionless Stem Auto-Discovery:** `create_workflow(dir_path)` and `create_workflow(stem_without_ext)` seamlessly resolve `workflow_config.yaml`, `workflow_config.yml`, or `workflow_config.json`.

---

## Configuration Schema

A workflow configuration file contains four primary sections:

1. **`workflow`**: Run metadata, preprocessing toggles, and execution mode:
   - `mesh_type`: Target mesh family (`"mpas"`).
   - `case_index`: Integer run identifier (e.g. `1`, `3`).
   - `date`: Run date string (`YYYYMMDD`) or `null` for today's date.
   - `model`: Meshing engine (`"jigsaw"`).
   - `simplify_hydrosheds_river_network`: Boolean flag for river network simplification.
   - `process_coastline`: Boolean flag for coastline buffering and island removal.
   - `debug`: Diagnostic logging verbosity (`true`/`false`).
   - `largest_outlets`: Number of largest river basins/outlets to process.
   - `run_on_hpc`: Boolean flag. When `true`, creates SLURM batch scripts; when `false`, runs mesh generation locally.

2. **`resolutions`**: Physical target resolutions in **kilometres**:
   - `ocean_km`: Background ocean mesh resolution.
   - `land_km`: Background land resolution.
   - `river_network_km`: Refinement resolution along simplified river networks.
   - `coastline_km`: Refinement resolution along coastline buffer zones.
   - `lake_boundary_km`: Refinement resolution along lake shorelines.
   - `watershed_boundary_km`: Refinement resolution along watershed divides.
   - `dam_km`: Optional dam point-feature refinement.
   - `city_km`: Optional urban point-feature refinement.

3. **`paths`**: Input datasets, output workspaces, and intermediate vectors/rasters with multi-platform support.

4. **`jigsaw`**: Engine-specific numerical knobs, spacing grid definitions, and HPC submission parameters:
   - `standalone`: Standalone directory naming convention toggle (`true`/`false`).
   - `create_directory`: Auto-create output workspace directories (`true`/`false`).
   - `hours`: Wall-clock limit for SLURM jobs.
   - `slurm`: SLURM partition name (e.g. `"slurm"`, `"cpu_single"`).
   - Spacing grid dimensions (`ncolumn_space`, `nrow_space`) and feature flags (`iFlag_geom_*`, `iFlag_spac_*`).

---

## Fully Commented Reference Workflow (`workflow_config.yaml`)

```yaml
# ==============================================================================
# Tier 1 Specification: Variable-Resolution MPAS Unstructured Mesh Workflow
# ==============================================================================

# High-level workflow orchestration and execution targets
workflow:
  date: "20260901"                                # Run date (YYYYMMDD); null defaults to today
  mesh_type: mpas                                 # Target mesh family (MPAS / JIGSAW)
  case_index: 1                                   # Case index number
  model: jigsaw                                   # Core mesh generation backend
  simplify_hydrosheds_river_network: false        # Optional line simplification for rivers
  process_coastline: false                        # Coastline masking and buffering
  debug: true                                     # Enable diagnostic logging
  largest_outlets: 10                             # Retain top N largest basin outlets
  run_on_hpc: false                               # Set true to generate SLURM batch job scripts

# Multi-scale target cell resolutions (in kilometers)
resolutions:
  ocean_km: 100                                   # Coarse open ocean resolution (100 km)
  land_km: 50.0                                   # Continental land background resolution (50 km)
  river_network_km: 5                             # River corridor refinement (5 km)
  coastline_km: 10                                # Near-shore transition resolution (10 km)
  lake_boundary_km: 5                             # Lake perimeter refinement (5 km)
  watershed_boundary_km: 5                        # Watershed divide refinement (5 km)

# Cross-platform dataset and workspace paths
paths:
  input_workspace:
    Darwin: /Users/me/data/hexwatershed
    Linux: /qfs/people/user/data/input
    Windows: D:/data/hexwatershed
    default: ./input
  output_workspace:
    Darwin: /Users/me/scratch/earthmesh/mpas
    Linux: /scratch/user/output
    Windows: D:/scratch/output
    default: ./output
  river_network_vector:
    Darwin: /Users/me/data/vector/rivers.geojson
    Linux: /data/vector/rivers.geojson
    default: ./data/rivers.geojson
  river_network_raster:
    Darwin: /Users/me/data/raster/rivers.tif
    Linux: /data/raster/rivers.tif
    default: ./data/rivers.tif
  lake_boundary_vector: ./data/lake_boundary.geojson
  lake_boundary_raster: ./data/lake_boundary.tif
  watershed_boundary_raster: ./data/basin_boundary.tif
  coastline_raster: ./data/land_ocean_mask.tif
  region_geometry: ./data/region.geojson
  dam_vector: ""

# Backend tuning parameters and spacing field flags for JIGSAW
jigsaw:
  standalone: true                                # Run standalone generator without external wrapper
  create_directory: true                          # Auto-create case subfolder in workspace
  hours: 3                                        # SLURM walltime allocation (hours)
  slurm: cpu_single                               # SLURM partition or queue name
  ncolumn_space: 3600                             # Spacing raster longitudinal grid cells (0.1 deg)
  nrow_space: 1800                                # Spacing raster latitudinal grid cells (0.1 deg)
  iFlag_geom: true                                # Build boundary geometry file (.msh)
  iFlag_spac: true                                # Build spacing field raster (.msh)
  iFlag_spac_ocean: true                          # Apply ocean spacing constraint
  iFlag_spac_land: true                           # Apply land spacing constraint
  iFlag_geom_river_network: true                  # Enforce river network feature constraints
  iFlag_spac_river_network: true                  # Enforce river network resolution refinement
  iFlag_geom_lake_boundary: true                  # Enforce lake boundary feature constraints
  iFlag_spac_lake_boundary: true                  # Enforce lake boundary resolution refinement
  iFlag_spac_watershed_boundary: true             # Refine watershed divide boundaries
  iFlag_spac_coastline: true                      # Refine coastline transition zone
```

---

## Multi-Platform Path Resolution

Every path entry can either be a simple string or a dictionary keyed by operating system (`platform.system()`):

```yaml
output_workspace:
  Darwin: /Users/me/scratch/earthmesh/mpas
  Linux: /compyfs/user/model/output
  Windows: D:/scratch/output
  default: ./output
```

- `Darwin` is automatically selected on macOS.
- `Linux` is automatically selected on Linux / HPC clusters.
- `Windows` is automatically selected on Windows.
- `default` is used if the current operating system is not listed.
- Forward slashes (`/`) should be used for all paths across all operating systems.
- Relative paths are resolved relative to the directory containing the configuration file.

---

## Python Workflow Driver

Workflow runner scripts are concise and format-agnostic:

```python
#!/usr/bin/env python3
import argparse
from pathlib import Path
from pyearthmesh.config.config_manager import create_workflow

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        "-c",
        default="workflow_config.yaml",
        help="Path to workflow YAML or JSON configuration file",
    )
    args = parser.parse_args()

    sFilename_workflow = Path(args.config).expanduser().resolve()

    # 1. Load configuration and compile workflow object in memory
    workflow = create_workflow(sFilename_workflow)

    # 2. Execute according to configuration settings
    if workflow.workflow.get("run_on_hpc"):
        workflow.create_hpc_job()
    else:
        workflow.setup()
        mesh = workflow.run()

if __name__ == "__main__":
    main()
```

---

## Standardized Examples

All workflow examples in [`examples/`](../examples/) follow this convention with dual-format parity:

- [**`examples/minimal_workflow/`**](../examples/minimal_workflow/): Minimal global test case.
- [**`examples/ocn30coast10lnd10/`**](../examples/ocn30coast10lnd10/): Global 30km ocean, 10km land, 10km coastline mesh.
- [**`examples/qinghaihu/`**](../examples/qinghaihu/): Regional Qinghai Lake refined mesh (river, lake, and watershed features).
- [**`examples/huanghe/`**](../examples/huanghe/): Yellow River basin workflow template.
