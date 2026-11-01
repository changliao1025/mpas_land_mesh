# Examples & Workflow Configuration Guide

## Summary

The workflow configuration for `mpas_land_mesh` has been modernized to adopt the unified **Two-Tier Compiler Architecture** established in [`pyearthmesh`](https://github.com/e3sm-project/pyearthmesh) (see reference specification: [`plan/YAML_WORKFLOW_CONFIGURATION.md`](../plan/YAML_WORKFLOW_CONFIGURATION.md)).

A **single, comment-annotated YAML file** (`workflow_config.yaml`) is now the **primary source of truth** for an entire run. Standard JSON (`workflow_config.json`) is fully supported for backward compatibility and programmatic workflows.

The unified workflow eliminates:
- Hardcoded Python script parameters
- Intermediate generated flat configuration files (`mpas_configuration.json`, `jigsaw_configuration.json`)
- Repetitive on-disk JSON modifications (`change_json_key_value`)

The workflow class (`MPASWorkflow`) translates the sectioned configuration in memory, derives grid spacings and boundary geometries, and orchestrates either local execution or SLURM HPC job script generation.

---

## Architectural Layers

```text
                     ┌──────────────────────────────────────────────┐
                     │            workflow_config.yaml              │
                     │    • Comments explaining parameters (#)      │
                     │    • Domain intent (resolutions in km)       │
                     │    • Cross-platform dataset paths            │
                     │    • Preprocessing & feature toggles         │
                     │    • Execution mode (run_on_hpc: true/false) │
                     └──────────────────────┬───────────────────────┘
                                            │
                                            ▼
                     ┌──────────────────────────────────────────────┐
                     │          load_workflow_config()              │
                     │    • Parses YAML (or JSON) seamlessly        │
                     │    • Validates schema against defaults       │
                     │    • Resolves OS-specific paths & relative   │
                     │    • Normalizes run date (YYYYMMDD)          │
                     └──────────────────────┬───────────────────────┘
                                            │
                                            ▼
                     ┌──────────────────────────────────────────────┐
                     │       create_workflow(config_or_path)        │
                     │    • Factory selects MPASWorkflow            │
                     │    • Derives pBoundary_wkt from GeoJSON      │
                     │    • Calculates grid spacing (ncolumn, nrow) │
                     │    • Compiles native JIGSAW deck in memory   │
                     │    • Builds case adapter (mpascase)          │
                     └──────────────────────┬───────────────────────┘
                                            │
                            ┌───────────────┴───────────────┐
                            ▼                               ▼
                  [workflow.run_on_hpc: true]   [workflow.run_on_hpc: false]
                            │                               │
                            ▼                               ▼
                    workflow.create_hpc_job()        workflow.setup()
                    • Generates run_mpas.py          • Prepares out/tmp dirs
                    • Generates submit.job           workflow.run()
                                                     • Generates MPAS mesh
```

---

## Standardized Examples

All workflow examples in this directory follow this unified convention:

- [**`minimal_workflow/`**](./minimal_workflow/): Minimal global test case.
- [**`ocn30coast10lnd10/`**](./ocn30coast10lnd10/): Global 30km ocean, 10km land, 10km coastline mesh.
- [**`qinghaihu/`**](./qinghaihu/): Regional Qinghai Lake refined mesh (river, lake, and watershed features).
- [**`huanghe/`**](./huanghe/): Yellow River basin workflow template.

Each example directory contains:
- `workflow_config.yaml`: Comment-annotated primary configuration.
- `workflow_config.json`: Backward-compatible programmatic JSON configuration.
- `run_workflow.py`: Unified Python workflow driver.

---

## Running an Example

```bash
cd examples/qinghaihu

# Run using the default workflow_config.yaml
python run_workflow.py

# Or specify a custom YAML or JSON configuration file
python run_workflow.py --config workflow_config.yaml
```
