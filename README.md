# Unified Land-River Mesh (ULRM) Workflow Module

A focused package for running unified land-river mesh workflows for MPAS, built on top of the `pyearth` geospatial ecosystem.

## Architecture & Dependencies

This module is part of the layered geospatial mesh architecture:
```
mpas_land_mesh (land-river mesh workflows)
  ├── pyearthmesh (unstructured & structured mesh generation, JIGSAW, workflow compiler)
  │     └── pyearth (foundational GIS, geometry, raster, vector, system)
  ├── pyearthbuffer (geodesic GCS buffer operations)
  ├── pyearthriver (river network topology, flowlines, rivergraphs)
  ├── hexwatershed_utility (hydrological river network preprocessing & simplification)
  └── pyearth (shared GIS and scientific utilities)
```

- **`pyearth`**: Foundation layer providing spatial references, geometry, GDAL raster/vector I/O, format conversion, and system environment utilities.
- **`pyearthmesh`**: Unstructured (MPAS, JIGSAW) and structured mesh generation, workflow configuration schemas, and compilation.
- **`pyearthriver`**: River network processing, flowline graphs, and hydrological features.
- **`pyearthbuffer`**: Geodesic buffer calculations on the WGS84 ellipsoid for GCS vector data.
- **`hexwatershed_utility`**: River network preprocessing, topology building, and line simplification.
- **`mpas_land_mesh`**: High-level unified land-river mesh workflow execution, coastline/ocean masking, and basin/river network simplification.

## Structure

```
mpas_land_mesh/
├── __init__.py                     # Module initialization
├── README.md                       # This file
├── INSTALL.md                      # Installation instructions
├── USAGE.md                        # Detailed usage guide
├── MIGRATION.md                    # Migration notes
├── requirements.txt                # Minimal dependencies
├── setup.py                        # Package setup
├── preprocessing/                  # Preprocessing re-exports from hexwatershed_utility
│   └── __init__.py
└── utilities/                      # Domain-specific utilities
    ├── __init__.py
    ├── constants.py                # Physical and geometric constants
    └── mesh_stats.py               # MPAS mesh statistics
```

> **Note on imports:** Data classes and core capabilities are shared across packages:
> - Base geometry and topology (`pyvertex`, `pyedge`, `pynvector`), GIS, vector, raster, and spatial references are provided by `pyearth`.
> - River network topology (`pyrivergraph`, `pyflowline`, `pyconfluence`) and export functions are provided by `pyearthriver`.
> - Geodesic buffer operations on GCS WGS84 are provided by `pyearthbuffer`.
> - Mesh cell definitions (`pympas`, `pycelllink`, `pymeshcell`), mesh algorithms, and workflow configuration are provided by `pyearthmesh`.

## Key Dependencies

Minimal required packages:
- `pyearth` — foundational GIS, geometry, raster, vector, and spatial reference
- `pyearthmesh` — unstructured/structured mesh generation and workflow configuration
- `pyearthriver` — river network structures and flowlines
- `pyearthbuffer` — geodesic buffer calculations
- `hexwatershed_utility` — river network preprocessing and simplification
- `gdal` / `osgeo` — geospatial data processing
- `numpy>=1.20.0` — array operations
- `scipy>=1.7.0` — spatial operations
- `netCDF4` — MPAS NetCDF mesh handling
- `rtree` — spatial indexing for river network processing

## Core Functions

### Mesh Statistics (`mpas_land_mesh.utilities`)
- `get_mesh_stats()` — compute mesh cell count, resolution, and area metrics from MPAS NetCDF

### River Network Processing (`hexwatershed_utility.preprocess.feature.river`)
- `simplify_river_network()` (aliased as `simplify_hydrorivers_network`) — simplify river networks by drainage area and distance tolerance
- `tag_river_outlet()` — tag river outlets in a flowline dataset
- `get_outlet_location()` — get outlet coordinates
- `crop_river_network_by_boundary()` — crop flowlines to region boundary
- `convert_geometry_flowline()` — convert OGR geometry to a flowline object
- `precompute_flowline_geometries()` — precompute spatial bounds for flowlines
- `precompute_flowline_geometries_by_segment()` — precompute bounds grouped by segment

### Coastline & Mask Processing (`hexwatershed_utility.preprocess.feature.coastline`)
- `create_land_ocean_mask()` / `create_land_ocean_mask_from_naturalearth()` — create land/ocean mask from Natural Earth data
- `create_land_ocean_mask_from_shoreline()` — create land/ocean mask from custom shoreline vector
- `create_land_ocean_mask_from_basins()` — create land/ocean mask from regional basin boundaries
- `extend_land_mask_to_rivers()` (aliased as `fix_naturalearth_hydrosheds_incompatibility`) — reconcile land polygons and coastal flowlines
- `geometries_bbox_overlap()` — test bounding-box overlap between geometries
- `remove_internal_polygon()` — remove interior rings (lakes/holes) from land polygons

### Configuration & Workflow Management (`pyearthmesh.config.config_manager`)
- `load_workflow_config()` — load single sectioned workflow YAML/JSON configuration
- `dump_workflow_config()` — export workflow configuration to YAML or JSON
- `create_workflow()` — instantiate mesh workflow object and case adapter
- `compile_workflow()` — compile user-tier workflow spec into backend mesh parameters
- `DEFAULT_CONFIG` — reference default dictionary schema for all workflows

## Installation

```bash
# Recommended: conda environment with core dependencies
conda install -c conda-forge gdal numpy scipy rtree pyyaml

# Install pyearth, pyearthriver, and pyearthmesh
pip install pyearth pyearthriver pyearthmesh pyyaml
```

## Usage

```python
from pyearthmesh.config.config_manager import create_workflow, load_workflow_config
from mpas_land_mesh.preprocessing import simplify_hydrorivers_networks

# --- Step 1: Preprocess river network ---
simplify_hydrorivers_networks(
    sFilename_flowline_hydroshed_in='HydroRIVERS_v10.shp',
    sFilename_flowline_hydroshed_out='river_simplified.geojson',
    dDistance_tolerance_in=10_000,       # metres
    dDrainage_area_threshold_in=1e10,    # m²
    nOutlet_largest=10,
)

# --- Step 2: Run mesh generation via unified YAML workflow (with # comments) ---
workflow = create_workflow('workflow_config.yaml')

# Run locally or create HPC job scripts
if workflow.run_on_hpc:
    workflow.create_hpc_job(hours=5)
else:
    workflow.setup()
    workflow.run()
```

See [`examples/run_minimal_workflow.py`](../examples/run_minimal_workflow.py) for a complete end-to-end example.

## Benefits

1. **Lightweight** — only essential functions, no bloated dependencies
2. **Standalone** — distributable as a single package
3. **Fast installation** — minimal dependency resolution time
4. **Clear API** — simplified interface for common workflows
5. **Maintainable** — easier to debug and modify

## Layered Architecture

| Package | Role |
|---|---|
| `pyearth` | Foundational GIS, geometry, coordinate transforms, raster/vector operations |
| `pyearthbuffer` | Geodesic buffer calculations on WGS84 for GCS vector data |
| `pyearthriver` | Hydrological flowline data structures, river networks, vector exports |
| `hexwatershed_utility` | River network & coastline preprocessing, flowline simplification, land-ocean masking |
| `pyearthmesh` | Unstructured and structured mesh generation, JIGSAW wrapper, config schema |
| `mpas_land_mesh` | Unified Land-River Mesh (ULRM) preprocessing, masking, workflow execution |

## License

Same as parent project (see `LICENSE` file in root directory).
