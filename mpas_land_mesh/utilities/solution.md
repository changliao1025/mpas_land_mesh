Viewed ocn30coast10lnd10.py:1-100
Viewed jigsawcase.py:1-100
Viewed jigsawcase.py:100-200

### 1. Is it possible and reasonable to merge them?

**Yes, it is both possible and highly reasonable.**

Currently, the two files represent two halves of the same problem:
* [`workflow_config.py`](file:///Users/changliao/workspace/python/mpas_land_mesh/mpas_land_mesh/utilities/workflow_config.py) handles the **user-facing input** (paths, flags, resolutions in km, HPC options).
* [`config_manager.py`](file:///Users/changliao/workspace/python/mpas_land_mesh/mpas_land_mesh/utilities/config_manager.py) handles the **engine-facing setup** (50+ numerical and geometric parameters for JIGSAW).

Because they are separated without a clean bridge between them, the user scripts (such as [`examples/run_workflow.py#L168-L200`](file:///Users/changliao/workspace/python/mpas_land_mesh/examples/run_workflow.py#L168-L200)) currently have to do this:

```python
# Current awkward glue code in run_workflow.py:
create_jigsaw_template_configuration_file(sFilename_jigsaw_configuration_json)
change_json_key_value(sFilename_jigsaw_configuration_json, "sWorkspace_output", sWorkspace_output)
oJigsaw = read_jigsaw_configuration_file(...)
copy2(sFilename_jigsaw_configuration_json, sFilename_jigsaw_configuration_copy)

# 20+ calls modifying JSON on disk!
change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_geom", "true")
change_json_key_value(sFilename_jigsaw_configuration_copy, "dResolution_ocean", dResolution_ocean)
...
```

Merging them into a unified configuration system (e.g. under `mpas_land_mesh/configuration/` or a single `mpas_land_mesh/utilities/config.py`) eliminates this fragile on-disk editing pattern and brings configuration management under a single, consistent API.

---

### 2. What is the best way to manage configuration for both User Input and Model Setup?

In scientific modeling and mesh generation pipelines, the gold standard is a **two-tier "compiler" pattern**:

```
 ┌────────────────────────────────────────────────────────┐
 │ 1. High-Level User Config (workflow_config.json)       │
 │    - Domain concepts: resolutions in km, datasets      │
 │    - Cross-platform paths                              │
 │    - Feature toggles (e.g. process_coastline: true)    │
 │    - Optional expert overrides                         │
 └───────────────────────────┬────────────────────────────┘
                             │
                             │ Translated / Derived automatically
                             ▼
 ┌────────────────────────────────────────────────────────┐
 │ 2. Low-Level Model Config (JIGSAW deck)                │
 │    - 50+ JIGSAW flags, sizing functions, kernels       │
 │    - Derived grid dimensions (ncolumn, nrow)           │
 │    - Exact file paths to generated rasters/vectors     │
 └───────────────────────────┬────────────────────────────┘
                             │
                             ▼
 ┌────────────────────────────────────────────────────────┐
 │ 3. Execution Engine (jigsawcase)                       │
 └────────────────────────────────────────────────────────┘
```

#### Why keep the two tiers?
* **Users shouldn't need to know JIGSAW internals**: Most users only care about *“generate an MPAS mesh with 30 km ocean, 10 km land, using HydroSHEDS river data.”* They should not have to manually set `mesh_rad2: 1.5`, `optm_iter: 32`, `dhdx_lim: 0.25`, or calculate `ncolumn_space = int(360 / (30.0 / 3600 * res))`.
* **Engine parameters still need defaults and overrides**: Advanced users occasionally *do* need to tune JIGSAW parameters (like optimization tolerances or mesh kernels).

---

### 3. Recommended Architecture: A Unified Configuration System

Instead of two disconnected modules, a unified configuration module (e.g., [`mpas_land_mesh/utilities/config.py`](file:///Users/changliao/workspace/python/mpas_land_mesh/mpas_land_mesh/utilities) or `mpas_land_mesh/configuration/`) can provide:

#### A. Single High-Level User JSON with an optional `jigsaw_parameters` block
The user provides a clean JSON:
```json
{
  "workflow": {
    "date": null,
    "case_index": 1,
    "process_coastline": true,
    "simplify_hydrosheds_river_network": true
  },
  "resolutions": {
    "ocean_km": 30,
    "land_km": 10,
    "river_network_km": 10,
    "coastline_km": 10
  },
  "paths": {
    "output_workspace": "./output",
    "hydrosheds_rivers": {
      "Darwin": "~/data/hydrosheds/rivers.shp",
      "Linux": "/compyfs/data/hydrosheds/rivers.shp"
    }
  },
  "jigsaw": {
    "standalone": true,
    "hours": 5,
    "slurm": "slurm",
    "overrides": {
      "optm_iter": 64
    }
  }
}
```

#### B. In-Memory Compilation (No `change_json_key_value` disk churn)
The unified module automatically translates high-level workflow settings into the JIGSAW configuration dictionary:

```python
class WorkflowConfig:
    """Parses, validates, and resolves platform paths for workflow JSON."""
    @classmethod
    def load(cls, filepath: str | Path) -> "WorkflowConfig":
        ...

    def build_jigsaw_config(self, generated_files: dict = None) -> dict:
        """Derives the complete JIGSAW input deck in-memory.

        1. Starts from JIGSAW default values (mesh_rad2, optm_iter, etc.).
        2. Translates km resolutions into grid spacing and degrees.
        3. Sets feature flags based on self.workflow settings.
        4. Populates paths to vectors/rasters produced during preprocessing.
        5. Applies any expert overrides specified in config['jigsaw']['overrides'].
        """
        ...
```

#### C. Clean separation of concerns with `jigsawcase`
Currently, [`read_jigsaw_configuration_file()`](file:///Users/changliao/workspace/python/mpas_land_mesh/mpas_land_mesh/utilities/config_manager.py#L169) does three unrelated things:
1. Reads a JSON file from disk.
2. Creates directories via `Path.mkdir`.
3. Instantiates and returns a [`jigsawcase`](file:///Users/changliao/workspace/python/mpas_land_mesh/mpas_land_mesh/classes/jigsawcase.py) object.

The cleaner pattern is:
* The config manager produces the **configuration dictionary** (validated and complete).
* [`jigsawcase`](file:///Users/changliao/workspace/python/mpas_land_mesh/mpas_land_mesh/classes/jigsawcase.py) receives that dictionary directly:
  ```python
  jigsaw_config = workflow_config.build_jigsaw_config(...)
  case = jigsawcase(jigsaw_config)
  case.jigsaw_setup()
  case.jigsaw_run()
  ```

---

### Summary Recommendation

| Action | Verdict | Rationale |
| :--- | :--- | :--- |
| **Merge into one module / package** | **Yes (Recommended)** | Eliminates code duplication, unifies path resolution, and avoids maintaining two separate configuration loaders. |
| **Maintain the two levels of abstraction** | **Yes (Recommended)** | Keep the clean separation between **User Intent** (high-level, km units, cross-platform) and **Engine Deck** (low-level JIGSAW knobs), but automate the translation between them. |
| **Deprecate `change_json_key_value` disk patching** | **Yes (Recommended)** | Build and mutate the engine config in Python memory; only serialize to JSON once when saving run artifacts. |