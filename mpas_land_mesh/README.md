# Unified Land-River Mesh (ULRM) Workflow Module

A lightweight, standalone module for running the land-river mesh workflow without installing full dependencies.

## Overview

This module extracts the essential functions from the full workflow and packages them into a minimal, self-contained library. It eliminates the need to install:
- `pyearth`
- `hexwatershed_utility`
- `pyflowline` (only config utilities)
- Other heavy dependencies

## Structure

```
mpas_land_mesh/
├── __init__.py                 # Module initialization
├── README.md                   # This file
├── core/                       # Core functionality
│   ├── __init__.py
│   ├── vector_utils.py         # Vector file operations
│   ├── raster_utils.py         # Raster processing
│   └── config_utils.py         # Configuration management
├── preprocessing/              # Data preprocessing
│   ├── __init__.py
│   ├── river_network.py        # River network processing
│   ├── coastline.py            # Coastline processing
│   └── watershed.py            # Watershed boundary processing
└── examples/                   # Usage examples
    └── run_minimal_workflow.py
```

## Key Dependencies

Minimal required packages:
- `gdal/osgeo` - Geospatial data processing
- `json` - Configuration management
- `os`, `glob`, `shutil` - File operations
- `datetime` - Timestamps

## Core Functions

### Vector Operations
- `get_field_and_value()` - Extract field attributes
- `add_field_to_vector_file()` - Add fields to vector
- `merge_features()` - Merge vector features
- `gdal_write_wkt_to_vector_file()` - Write WKT to file

### Raster Operations
- `convert_vector_to_global_raster()` - Vector to raster conversion
- `create_raster_buffer_zone()` - Create buffer zones
- `fix_raster_antimeridian_issue()` - Fix antimeridian issues

### River Network Processing
- `simplify_hydrorivers_networks()` - Simplify river networks
- `tag_river_outlet()` - Tag river outlets
- `get_outlet_location()` - Get outlet coordinates

### Coastline Processing
- `create_land_ocean_mask_from_naturalearth()` - Create land-ocean mask
- `fix_naturalearth_hydrosheds_incompatibility()` - Fix data incompatibility

### Watershed Processing
- `find_minimal_hydrobasins_watershed_boundary()` - Find watershed boundaries

### Configuration Management
- `create_pyflowline_template_configuration_file()` - Create config templates
- `pyflowline_read_configuration_file()` - Read config files
- `change_json_key_value()` - Modify JSON configurations

## Installation

```bash
# Install minimal dependencies
pip install gdal numpy

# Or use conda for easier GDAL installation
conda install -c conda-forge gdal numpy
```

## Usage

```python
from mpas_land_mesh.preprocessing import process_river_network
from mpas_land_mesh.preprocessing import process_coastline
from mpas_land_mesh.core.config_utils import create_workflow_config

# Set up configuration
config = create_workflow_config(
    output_dir='/path/to/output',
    resolution_ocean=18,
    resolution_land=12,
    resolution_coastline=6
)

# Process river network
river_output = process_river_network(config)

# Process coastline
coastline_output = process_coastline(config)
```

## Benefits

1. **Lightweight**: Only essential functions, no bloated dependencies
2. **Standalone**: Can be distributed as a single package
3. **Fast Installation**: Minimal dependency resolution time
4. **Clear API**: Simplified interface for common workflows
5. **Maintainable**: Easier to debug and modify

## Original Workflow Dependencies Eliminated

- `pyearth.toolbox.*` - Replaced with local implementations
- `hexwatershed_utility.*` - Core functions extracted
- `pyflowline.*` (most) - Only config utils kept
- Heavy visualization libraries - Not needed for core workflow

## License

Same as parent project (see LICENSE file in root directory)
