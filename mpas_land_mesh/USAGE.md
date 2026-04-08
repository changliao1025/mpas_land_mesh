# Usage Guide

## Overview

The ULRM Workflow module provides a streamlined interface for creating unified land-river meshes without requiring the full dependency stack.

## Basic Usage

### 1. Import the Module

```python
from mpas_land_mesh.core import (
    get_field_and_value,
    add_field_to_vector_file,
    merge_features,
    convert_vector_to_global_raster,
    create_config_template
)

from mpas_land_mesh.preprocessing import (
    get_outlet_location,
    simplify_river_network
)
```

### 2. Configuration Management

#### Create Configuration File

```python
from mpas_land_mesh.core.config_utils import create_config_template

config = create_config_template(
    'workflow_config.json',
    iCase_index=1,
    sWorkspace_output='/path/to/output',
    dResolution_ocean=18.0,
    dResolution_land=12.0,
    dResolution_coastline=6.0
)
```

#### Modify Configuration

```python
from mpas_land_mesh.core.config_utils import change_json_key_value

change_json_key_value(
    'workflow_config.json',
    'dResolution_ocean',
    30.0
)
```

### 3. Vector Operations

#### Extract Fields

```python
from mpas_land_mesh.core.vector_utils import get_field_and_value

fields, values = get_field_and_value('input.geojson')
print(f"Fields: {fields}")
print(f"Values: {values}")
```

#### Merge Features

```python
from mpas_land_mesh.core.vector_utils import merge_features

merge_features(
    'input.geojson',
    'merged_output.geojson',
    iFlag_force=True
)
```

### 4. Raster Operations

#### Convert Vector to Raster

```python
from mpas_land_mesh.core.raster_utils import convert_vector_to_global_raster

convert_vector_to_global_raster(
    'coastline.geojson',
    'coastline_raster.tif',
    dResolution_x_in=0.1,  # degrees
    dResolution_y_in=0.1,
    iFlag_boundary_only_in=1,
    dFill_value_in=1
)
```

#### Create Buffer Zone

```python
from mpas_land_mesh.core.raster_utils import create_raster_buffer_zone

create_raster_buffer_zone(
    'input.tif',
    'buffered.tif',
    iBuffer_value_in=2,
    iFill_value_in=1
)
```

### 5. River Network Processing

#### Get Outlet Location

```python
from mpas_land_mesh.preprocessing.river_network import get_outlet_location

longitude, latitude = get_outlet_location('river_network.geojson')
print(f"Outlet at: {longitude}, {latitude}")
```

#### Simplify River Network

```python
from mpas_land_mesh.preprocessing.river_network import simplify_river_network

simplify_river_network(
    'HydroRIVERS.shp',
    'simplified.geojson',
    dTolerance=1000.0,  # meters
    dDrainage_threshold=1e8  # m²
)
```

## Complete Workflow Example

```python
import os
from datetime import datetime
from mpas_land_mesh.core.config_utils import (
    create_config_template,
    create_jigsaw_config_template
)
from mpas_land_mesh.core.vector_utils import merge_features
from mpas_land_mesh.core.raster_utils import convert_vector_to_global_raster

# Setup
output_dir = './workflow_output'
os.makedirs(output_dir, exist_ok=True)

# Create configuration
config = create_config_template(
    os.path.join(output_dir, 'config.json'),
    sDate=datetime.now().strftime('%Y%m%d'),
    sWorkspace_output=output_dir,
    dResolution_ocean=18.0,
    dResolution_land=12.0
)

jigsaw_config = create_jigsaw_config_template(
    os.path.join(output_dir, 'jigsaw_config.json'),
    dResolution_ocean=18.0,
    dResolution_land=12.0
)

print(f"Configuration created in: {output_dir}")

# Process data (if you have input files)
# merge_features('input1.geojson', 'merged.geojson')
# convert_vector_to_global_raster('merged.geojson', 'raster.tif', 0.1, 0.1)
```

## Limitations of Minimal Module

This minimal module provides **configuration and basic utilities** but has placeholders for:

1. **River network simplification** - Basic implementation provided; use `hexwatershed_utility` for production
2. **Coastline processing** - Requires Natural Earth data and full `hexwatershed_utility`
3. **Watershed boundary extraction** - Requires HydroBasins data and full `hexwatershed_utility`
4. **Mesh generation** - Requires `pyflowline` and JIGSAW

## Migration to Full Workflow

When ready for production, install full dependencies:

```bash
pip install hexwatershed_utility pyflowline pyhexwatershed
```

Then use the original workflow scripts in [`examples/global/`](../../examples/global/).

## API Reference

### Core Utilities

- [`vector_utils.py`](core/vector_utils.py) - Vector file operations
- [`raster_utils.py`](core/raster_utils.py) - Raster processing
- [`config_utils.py`](core/config_utils.py) - Configuration management

### Preprocessing

- [`river_network.py`](preprocessing/river_network.py) - River network processing
- [`coastline.py`](preprocessing/coastline.py) - Coastline processing
- [`watershed.py`](preprocessing/watershed.py) - Watershed processing

## Support

For issues or questions:
1. Check the main project documentation
2. Review example workflows in [`examples/`](examples/)
3. Consult the original workflow files for comparison
