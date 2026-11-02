# Usage Guide

## Overview

The `mpas_land_mesh` package runs unified land-river mesh workflows for MPAS by integrating the specialized packages in the `pyearth` ecosystem:
- **`pyearth`**: Foundational GIS, vector, raster, spatial reference, and geometry operations.
- **`pyearthbuffer`**: Geodesic buffer operations on the WGS84 ellipsoid.
- **`hexwatershed_utility`**: Preprocessing for river networks, flowline simplification, and coastline/mask generation.
- **`pyearthmesh`**: Mesh generation (MPAS, JIGSAW) and unified YAML/JSON workflow configuration.
- **`mpas_land_mesh`**: Domain constants, MPAS mesh statistics, and workflow execution.

---

## Key Modules & Imports

### 1. Vector Operations (`pyearth`)

Extract field attributes or merge vector features using `pyearth`:

```python
from pyearth.toolbox.management.vector.fields import (
    add_field_to_vector_file,
    get_field_and_value,
    get_field_value,
)
from pyearth.toolbox.management.vector.merge_features import merge_features

# Extract fields and attribute dictionaries
fields, values = get_field_and_value("boundary.geojson")
print(f"Fields: {fields}")
print(f"First feature attributes: {values[0] if values else None}")

# Merge vector features into a single geometry
merge_features("input.geojson", "merged_output.geojson")
```

### 2. Raster Operations (`pyearth`)

Rasterize vectors and buffer raster masks:

```python
from pyearth.toolbox.analysis.image.raster_process import (
    create_raster_buffer_zone,
    fix_raster_antimeridian_issue,
)
from pyearth.toolbox.conversion.convert_vector_to_global_raster import (
    convert_vector_to_global_raster,
)

# Convert vector to global raster
convert_vector_to_global_raster(
    "coastline.geojson",
    "coastline_raster.tif",
    dResolution_x=0.1,  # degrees
    dResolution_y=0.1,
    iFlag_boundary_only_in=0,
    dFill_value_in=2,
)

# Dilate/buffer raster mask
create_raster_buffer_zone(
    "coastline_raster.tif",
    "coastline_buffered.tif",
    iBuffer_value_in=2,
    iBuffer_pixel_in=1,
)
```

### 3. River Network Preprocessing (`hexwatershed_utility`)

Simplify river networks and identify basin outlets:

```python
from hexwatershed_utility.preprocess.feature.river import (
    get_outlet_location,
    simplify_river_network,
    tag_river_outlet,
)

# Simplify river network by distance tolerance and drainage area
simplify_river_network(
    sFilename_flowline_hydroshed_in="HydroRIVERS.shp",
    sFilename_flowline_hydroshed_out="simplified_rivers.geojson",
    dDistance_tolerance_in=10_000.0,  # meters
    dDrainage_area_threshold_in=1e8,  # m²
    nOutlet_largest=10,
)

# Get coordinates of river outlet
lon, lat = get_outlet_location("simplified_rivers.geojson")
print(f"River outlet at: {lon}, {lat}")
```

### 4. Coastline & Land-Ocean Masking (`hexwatershed_utility`)

Create land-ocean masks and reconcile river mouths:

```python
from hexwatershed_utility.preprocess.feature.coastline import (
    create_land_ocean_mask_from_naturalearth,
    create_land_ocean_mask_from_shoreline,
    extend_land_mask_to_rivers,
)

# Generate land-ocean mask from Natural Earth
sRaster_mask, sVector_mask = create_land_ocean_mask_from_naturalearth(
    sWorkspace_coastline_output="./output_coastline",
    dResolution_x_in=30 / 3600.0,
    dResolution_y_in=30 / 3600.0,
    dThreshold_area_island=1000.0,  # km²
    iRaster_buffer_pixel=1,
)

# Reconcile protruding river mouths with land boundary
extend_land_mask_to_rivers(
    aFilename_flowlines=["simplified_rivers.geojson"],
    sFilename_vector_land=sVector_mask,
    sFilename_vector_land_updated="./output_coastline/land_ocean_mask_updated.geojson",
)
```

### 5. Mesh Statistics & Constants (`mpas_land_mesh`)

Inspect generated MPAS NetCDF mesh files and use standard constants:

```python
from mpas_land_mesh.utilities.constants import (
    DRAINAGE_AREA_MULTIPLIER,
    ISLAND_AREA_MULTIPLIER,
    KM2_TO_M2,
    earth_radius,
)
from mpas_land_mesh.utilities.mesh_stats import get_mesh_stats

# Read statistics from an MPAS NetCDF mesh
stats = get_mesh_stats("mpas_mesh.nc")
print(f"Cells: {stats['nCells']}")
print(f"Min cell resolution: {stats['dLength_min']:.2f} m")
print(f"Max cell resolution: {stats['dLength_max']:.2f} m")
```

### 6. Workflow Configuration & Execution (`pyearthmesh`)

Load, mutate, and compile workflow configurations:

```python
from pyearthmesh.config.config_manager import (
    DEFAULT_CONFIG,
    compile_workflow,
    create_workflow,
    dump_workflow_config,
    load_workflow_config,
)

# Load YAML/JSON workflow
config = load_workflow_config("workflow_config.yaml")

# Instantiate workflow engine
workflow = create_workflow(config)
if workflow.run_on_hpc:
    workflow.create_hpc_job(hours=5)
else:
    workflow.setup()
    workflow.run()
```
