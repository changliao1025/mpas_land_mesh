# Unified Land-River Mesh (ULRM) Workflow Module

A lightweight, standalone module for running the land-river mesh workflow without installing full dependencies.

## Overview

This module extracts the essential functions from the full workflow and packages them into a minimal, self-contained library. It eliminates the need to install:
- `pyearth`
- `hexwatershed_utility`
- `pyflowline` (only config utilities are kept)
- Other heavy dependencies

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
├── classes/                        # Core data-model classes
│   ├── confluence.py
│   ├── edge.py
│   ├── flowline.py
│   ├── jigsawcase.py               # JIGSAW case object
│   ├── link.py
│   ├── mpas.py
│   ├── nvector.py
│   ├── rivergraph.py
│   └── vertex.py
├── mesh/                           # Mesh generation
│   ├── jigsaw/                     # JIGSAW mesh utilities
│   │   ├── inpoly2.py
│   │   ├── loadgeo.py
│   │   ├── run_jigsaw.py
│   │   ├── saveesm.py
│   │   └── utility.py
│   └── mpas/                       # MPAS mesh utilities
│       ├── convert_attributes.py
│       ├── create_mpas_mesh.py
│       └── utilities/
│           ├── compose.py
│           ├── mpasmsh.py
│           └── spacing.py
├── preprocessing/                  # Data preprocessing
│   ├── __init__.py
│   ├── river_networks.py           # River network processing
│   └── coastlines.py               # Coastline processing (optional deps)
└── utilities/                      # Core utility functions
    ├── __init__.py
    ├── change_json_key_value.py    # JSON configuration helpers
    ├── config_manager.py           # JIGSAW configuration manager
    ├── constants.py                # Physical constants
    ├── gcsbuffer.py                # Geographic buffer utilities
    ├── geometry.py                 # Spherical geometry calculations
    ├── io.py                       # I/O helpers (direct import only)
    ├── object.py                   # Object helpers (direct import only)
    ├── raster.py                   # Raster processing
    ├── spatial_reference.py        # Coordinate reprojection
    ├── system.py                   # System/path utilities
    └── vector.py                   # Vector file operations
```

> **Note on circular imports:** `utilities/object.py` and `utilities/io.py` depend on
> `mpas_land_mesh.classes`, which in turn imports from `mpas_land_mesh.utilities`.
> These two modules are **not** re-exported from `mpas_land_mesh.utilities` at
> package-init time to avoid circular imports.  Import them directly:
>
> ```python
> from mpas_land_mesh.utilities.object import find_vertex_on_edge, slerp
> from mpas_land_mesh.utilities.io import export_vertex_to_geojson
> ```

## Key Dependencies

Minimal required packages:
- `gdal` / `osgeo` — geospatial data processing (vector & raster)
- `numpy>=1.20.0` — array operations
- `scipy>=1.7.0` — spatial operations (e.g. `ndimage` for raster buffering)
- `rtree` — spatial indexing for river network processing
- `json`, `os`, `glob`, `shutil` — standard library

Optional (for coastline processing):
- `pyearth` — required by `mpas_land_mesh.preprocessing.coastlines`
- `pyearthbuffer` — required by `mpas_land_mesh.preprocessing.coastlines`

## Core Functions

### Vector Operations (`mpas_land_mesh.utilities`)
- [`get_field_and_value()`](utilities/vector.py) — extract field attributes from a vector file
- [`add_field_to_vector_file()`](utilities/vector.py) — add a field to a vector file
- [`merge_features()`](utilities/vector.py) — merge vector features into one file
- [`write_wkt_to_vector_file()`](utilities/vector.py) — write WKT geometry to a vector file
- [`remove_small_polygon()`](utilities/vector.py) — remove polygons below an area threshold
- [`get_vector_format_from_filename()`](utilities/vector.py) — resolve OGR driver from filename
- [`get_available_vector_formats()`](utilities/vector.py) — list GDAL-available vector formats

### Raster Operations (`mpas_land_mesh.utilities`)
- [`convert_vector_to_global_raster()`](utilities/raster.py) — rasterize a vector file to a global GeoTIFF
- [`create_raster_buffer_zone()`](utilities/raster.py) — dilate non-zero pixels in a raster
- [`fix_raster_antimeridian_issue()`](utilities/raster.py) — fix antimeridian discontinuities
- [`gdal_read_geotiff_file()`](utilities/raster.py) — read a GeoTIFF and return data + metadata

### Geographic Buffer Utilities (`mpas_land_mesh.utilities`)
- [`meters_to_degrees()`](utilities/gcsbuffer.py) — convert metres to degrees at a given latitude
- [`create_geometry_buffer_degrees()`](utilities/gcsbuffer.py) — buffer an OGR geometry in degrees
- [`create_wkt_buffer_degrees()`](utilities/gcsbuffer.py) — buffer a WKT geometry in degrees
- [`create_wkt_buffer_distance()`](utilities/gcsbuffer.py) — buffer a WKT geometry by a metric distance
- [`create_file_buffer_degrees()`](utilities/gcsbuffer.py) — buffer all features in a vector file
- [`create_buffer_with_meter_distance()`](utilities/gcsbuffer.py) — buffer a vector file in metres

### Spherical Geometry (`mpas_land_mesh.utilities`)
- [`calculate_distance_based_on_longitude_latitude()`](utilities/geometry.py) — great-circle distance
- [`haversine()`](utilities/geometry.py) — haversine formula
- [`calculate_polygon_area()`](utilities/geometry.py) — polygon area on the sphere
- [`calculate_spherical_triangle_area()`](utilities/geometry.py) — spherical triangle area
- [`convert_longitude_latitude_to_sphere_3d()`](utilities/geometry.py) — lon/lat → 3-D Cartesian
- [`find_minimal_enclosing_polygon()`](utilities/geometry.py) — minimal enclosing polygon
- [`split_international_date_line_polygon_coordinates()`](utilities/geometry.py) — IDL split

### Spatial Reference (`mpas_land_mesh.utilities`)
- [`reproject_coordinates()`](utilities/spatial_reference.py) — reproject a single coordinate pair
- [`reproject_coordinates_batch()`](utilities/spatial_reference.py) — reproject an array of coordinates

### River Network Processing (`mpas_land_mesh.preprocessing`)
- [`simplify_hydrorivers_networks()`](preprocessing/river_networks.py) — simplify HydroRIVERS networks by drainage area and distance tolerance
- [`tag_river_outlet()`](preprocessing/river_networks.py) — tag river outlets in a flowline dataset
- [`get_outlet_location()`](preprocessing/river_networks.py) — get outlet coordinates
- [`convert_geometry_flowline()`](preprocessing/river_networks.py) — convert OGR geometry to a flowline object
- [`precompute_flowline_geometries()`](preprocessing/river_networks.py) — precompute spatial bounds for flowlines
- [`precompute_flowline_geometries_by_segment()`](preprocessing/river_networks.py) — precompute bounds grouped by segment

### Coastline Processing (`mpas_land_mesh.preprocessing`)

> **Optional:** these functions require `pyearth` and `pyearthbuffer`.  If those
> packages are not installed the symbols are simply absent from the package
> namespace; import them directly from the module when needed.

- [`create_land_ocean_mask()`](preprocessing/coastlines.py) — create a land/ocean mask
- [`create_land_ocean_mask_from_naturalearth()`](preprocessing/coastlines.py) — create mask from Natural Earth data
- [`fix_naturalearth_hydrosheds_incompatibility()`](preprocessing/coastlines.py) — reconcile Natural Earth and HydroSHEDS coastlines
- [`geometries_bbox_overlap()`](preprocessing/coastlines.py) — test bounding-box overlap between geometries

### Configuration Management (`mpas_land_mesh.utilities`)
- [`create_jigsaw_template_configuration_file()`](utilities/config_manager.py) — write a JIGSAW JSON config template
- [`read_jigsaw_configuration_file()`](utilities/config_manager.py) — read a JIGSAW JSON config and return a `jigsawcase` object
- [`JigsawConfigManager`](utilities/config_manager.py) — class with `get_default_config()`, `create_template_config()`, `load_config()`

### JSON Utilities (`mpas_land_mesh.utilities`)
- [`change_json_key_value()`](utilities/change_json_key_value.py) — safely update a single key in a JSON file
- [`change_json_keys_values()`](utilities/change_json_key_value.py) — safely update multiple keys in a JSON file

## Installation

```bash
# Recommended: conda for easier GDAL installation
conda install -c conda-forge gdal numpy scipy rtree

# Or pip (GDAL must be installed separately)
pip install numpy scipy rtree
# Linux: apt-get install python3-gdal
# macOS: brew install gdal && pip install gdal
```

For optional coastline processing:
```bash
pip install pyearth pyearthbuffer
```

## Usage

```python
from mpas_land_mesh.utilities import (
    get_field_and_value,
    merge_features,
    add_field_to_vector_file,
    convert_vector_to_global_raster,
    create_jigsaw_template_configuration_file,
    read_jigsaw_configuration_file,
    change_json_key_value,
)
from mpas_land_mesh.preprocessing import (
    simplify_hydrorivers_networks,
)

# --- Step 1: simplify river network ---
simplify_hydrorivers_networks(
    sFilename_flowline_hydroshed_in='HydroRIVERS_v10.shp',
    sFilename_flowline_hydroshed_out='river_simplified.geojson',
    dDistance_tolerance_in=10_000,       # metres
    dDrainage_area_threshold_in=1e10,    # m²
    nOutlet_largest=10,
)

# Rasterize the simplified network
convert_vector_to_global_raster(
    'river_simplified.geojson',
    'river_network.tif',
    dResolution_x_in=30/3600 * 10,
    dResolution_y_in=30/3600 * 10,
)

# --- Step 2: create JIGSAW configuration ---
create_jigsaw_template_configuration_file('jigsaw_config.json')

change_json_key_value('jigsaw_config.json', 'sWorkspace_output', '/path/to/output')
change_json_key_value('jigsaw_config.json', 'iFlag_geom', 'true')
change_json_key_value('jigsaw_config.json', 'iFlag_spac', 'true')
change_json_key_value('jigsaw_config.json', 'dResolution_ocean', 100)
change_json_key_value('jigsaw_config.json', 'dResolution_land',  100)
change_json_key_value('jigsaw_config.json', 'dResolution_river_network', 10)
change_json_key_value('jigsaw_config.json', 'sFilename_river_network_raster', 'river_network.tif')

# --- Step 3: read config and run JIGSAW ---
oJigsaw = read_jigsaw_configuration_file(
    'jigsaw_config.json',
    iCase_index_in=1,
    sDate_in='20260401',
    iFlag_create_directory_in=1,
)
oJigsaw._jigsaw_create_hpc_job(sSlurm_in='slurm', hours_in=5)
# Submit the generated SLURM script manually
```

See [`examples/run_minimal_workflow.py`](../examples/run_minimal_workflow.py) for a complete end-to-end example.

## Benefits

1. **Lightweight** — only essential functions, no bloated dependencies
2. **Standalone** — distributable as a single package
3. **Fast installation** — minimal dependency resolution time
4. **Clear API** — simplified interface for common workflows
5. **Maintainable** — easier to debug and modify

## Original Workflow Dependencies Eliminated

| Original | Replacement |
|---|---|
| `pyearth.toolbox.*` | Local implementations in `utilities/` |
| `hexwatershed_utility.*` | Core functions extracted to `utilities/` and `preprocessing/` |
| `pyflowline.*` (most) | Only config utilities kept; replaced by `JigsawConfigManager` |
| Heavy visualization libraries | Not needed for core workflow |

## License

Same as parent project (see `LICENSE` file in root directory).
