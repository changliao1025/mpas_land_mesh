# JSON Workflow Configuration

## Summary

The workflow configuration was moved from hard-coded Python variables into JSON files. The Python workflow remains responsible for processing, while the JSON file contains run-specific settings such as flags, resolutions, paths, and JIGSAW job options.

This separates reusable workflow logic from machine- and run-specific configuration.

## Before

The Python example contained values such as:

```python
iFlag_process_coastline = 1
dResolution_ocean = 30
sWorkspace_output = '/path/to/output'
sFilename_geojson_geometery_feature = '/path/to/region.geojson'
```

Changing a run required editing the Python source code. Paths were also tied to one operating system or compute environment.

## After

A workflow JSON file contains four sections:

- `workflow`: run metadata and processing flags
- `resolutions`: mesh and preprocessing resolutions in kilometres
- `paths`: input and output locations
- `jigsaw`: JIGSAW and HPC submission settings

The Python entry point loads the configuration:

```python
from mpas_land_mesh.utilities.workflow_config import load_workflow_config

config = load_workflow_config(args.config)
workflow = config['workflow']
resolutions = config['resolutions']
paths = config['paths']
jigsaw = config['jigsaw']
```

The loader is implemented in [`mpas_land_mesh/utilities/workflow_config.py`](../mpas_land_mesh/utilities/workflow_config.py).

## Configuration Example

```json
{
  "workflow": {
    "date": null,
    "mesh_type": "mpas",
    "case_index": 1,
    "model": "jigsaw",
    "simplify_hydrosheds_river_network": false,
    "process_coastline": false,
    "debug": true,
    "largest_outlets": 10
  },
  "resolutions": {
    "ocean_km": 100,
    "land_km": 100,
    "river_network_km": 10,
    "coastline_km": 10
  },
  "paths": {
    "input_workspace": {
      "Windows": "./input",
      "Linux": "/data/project/input",
      "Darwin": "./input",
      "default": "./input"
    },
    "output_workspace": {
      "Windows": "./output",
      "Linux": "/scratch/project/output",
      "Darwin": "./output",
      "default": "./output"
    },
    "hydrosheds_rivers": {
      "Windows": "./input/HydroRIVERS_v10.gpkg",
      "Linux": "/data/project/input/HydroRIVERS_v10.gpkg",
      "Darwin": "./input/HydroRIVERS_v10.gpkg",
      "default": "./input/HydroRIVERS_v10.gpkg"
    },
    "region_geometry": {
      "Windows": "./input/region.geojson",
      "Linux": "/data/project/input/region.geojson",
      "Darwin": "./input/region.geojson",
      "default": "./input/region.geojson"
    },
    "dam_vector": {
      "Windows": "",
      "Linux": "",
      "Darwin": "",
      "default": ""
    }
  },
  "jigsaw": {
    "standalone": true,
    "create_directory": true,
    "hours": 5,
    "slurm": "slurm"
  }
}
```

## Platform-Specific Paths

Each path can be either a normal string or an object keyed by platform:

```json
"output_workspace": {
  "Windows": "D:/scratch/mpas/output",
  "Linux": "/scratch/mpas/output",
  "Darwin": "/Users/me/mpas/output",
  "default": "./output"
}
```

The loader uses Python's `platform.system()` values:

- `Windows` for Windows
- `Linux` for Linux
- `Darwin` for macOS
- `default` when no matching platform entry exists

Use `/` in JSON paths, including Windows paths. For example, use `D:/data/input` rather than `D:\\data\\input`.

Relative paths are resolved relative to the JSON file, not the current terminal directory. This makes a configuration portable when its input and output directories are kept beside it.

An empty string is preserved as an empty path and is appropriate for optional inputs such as an unused dam vector.

## Defaults and Validation

The loader provides defaults for missing settings. It also rejects unknown keys within the supported sections, which helps catch spelling mistakes in configuration files.

If `workflow.date` is `null`, the loader fills it with the current date in `YYYYMMDD` format.

The loader resolves `~` in paths. Environment variables such as `$HOME` and `%USERPROFILE%` are not currently expanded.

## Running a Workflow

From the repository root, pass the JSON file as the positional argument:

```text
python examples/run_workflow.py examples/minimal_workflow/minimal_workflow.json
```

For the 30 km ocean example:

```text
python examples/run_workflow.py examples/ocn30coast10lnd10/ocn30coast10lnd10.json
```

The configuration files created during this migration are:

- [`examples/minimal_workflow/minimal_workflow.json`](minimal_workflow/minimal_workflow.json)
- [`examples/ocn30coast10lnd10/ocn30coast10lnd10.json`](ocn30coast10lnd10/ocn30coast10lnd10.json)
- [`examples/qinghaihu/ocn100coast10lnd5riv4lak3.json`](qinghaihu/ocn100coast10lnd5riv4lak3.json)

## JSON Instead of YAML

JSON was retained because it is supported by Python's standard library, is strict and predictable, and is already used by the JIGSAW configuration workflow. YAML could be more convenient for comments and manual editing, but it would add a dependency and introduce additional parsing rules.

For this project, schema validation and clear platform-aware paths provide more value than changing the file format.
