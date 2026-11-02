"""Configuration manager for ULRM (Unified Land-River Mesh) workflows.

Manages loading, validation, multi-platform path resolution, and serialization
of user-facing ULRM configuration files (YAML or JSON).
"""

from __future__ import annotations

import copy
import json
import platform
from datetime import datetime
from pathlib import Path
from typing import Any, Union

import yaml

DEFAULT_ULRM_CONFIG: dict[str, Any] = {
    "workflow": {
        "date": None,
        "mesh_type": "mpas",
        "case_index": 1,
        "model": "jigsaw",
        "debug": True,
        "run_on_hpc": False,
    },
    "preprocessing": {
        "simplify_rivers": False,
        "largest_outlets": 10,
        "drainage_area_multiplier": 100.0,
        "process_coastline": False,
        "island_area_multiplier": 10.0,
        "raster_buffer_pixel": 2,
        "source": "naturalearth",  # 'naturalearth', 'shoreline', or 'basins'
    },
    "resolutions": {
        "ocean_km": 100.0,
        "land_km": 50.0,
        "river_network_km": 10.0,
        "coastline_km": 10.0,
        "lake_boundary_km": None,
        "watershed_boundary_km": None,
        "dam_km": None,
        "city_km": None,
    },
    "paths": {
        "input_workspace": "./input",
        "output_workspace": "./output",
        "hydrosheds_rivers": "",
        "region_geometry": "",
        "dam_vector": "",
        "naturalearth_land": "",
        "basins_boundary": "",
        # Direct mesh-ready inputs if skipping preprocessing
        "river_network_vector": "",
        "river_network_raster": "",
        "coastline_raster": "",
        "coastline_vector": "",
    },
    "jigsaw": {
        "standalone": True,
        "create_directory": True,
        "hours": 5,
        "slurm": "slurm",
        "ncolumn_space": 360,
        "nrow_space": 180,
        "iFlag_geom": True,
        "iFlag_spac": True,
        "iFlag_spac_ocean": True,
        "iFlag_spac_land": True,
        "iFlag_spac_coastline": True,
        "iFlag_geom_river_network": True,
        "iFlag_spac_river_network": True,
    },
}


def _select_platform_path(value: Any) -> str:
    """Select OS-specific path from a string or mapping."""
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if not isinstance(value, dict):
        raise ValueError(f"Path value must be a string or dict, got: {type(value)}")

    current_os = platform.system()
    selected = value.get(current_os, value.get("default", ""))
    if not isinstance(selected, str):
        raise ValueError(
            f"Path for {current_os} must be a string, got: {type(selected)}"
        )
    return selected


def _resolve_paths(paths: dict[str, Any], base_dir: Path) -> dict[str, str]:
    """Resolve multi-platform path dictionaries and normalize relative paths."""
    resolved: dict[str, str] = {}
    for key, val in paths.items():
        chosen = _select_platform_path(val)
        if not chosen:
            resolved[key] = ""
            continue
        p = Path(chosen)
        if not p.is_absolute():
            resolved[key] = str((base_dir / p).resolve())
        else:
            resolved[key] = str(p.resolve())
    return resolved


def _deep_merge(defaults: dict[str, Any], overrides: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge overrides into defaults dictionary."""
    result = copy.deepcopy(defaults)
    for k, v in overrides.items():
        if k in result and isinstance(result[k], dict) and isinstance(v, dict):
            result[k] = _deep_merge(result[k], v)
        else:
            result[k] = copy.deepcopy(v)
    return result


def load_ulrm_config(
    config_or_path: str | Path | dict[str, Any],
) -> dict[str, Any]:
    """Load, resolve, and normalize a ULRM workflow configuration.

    Parameters
    ----------
    config_or_path : str, Path, or dict
        Path to YAML or JSON config, or existing config dictionary.

    Returns
    -------
    dict[str, Any]
        Fully resolved and validated configuration dictionary.
    """
    base_dir = Path.cwd()
    if isinstance(config_or_path, (str, Path)):
        p = Path(config_or_path).expanduser().resolve()
        if not p.is_file():
            raise FileNotFoundError(f"Configuration file not found: {p}")
        base_dir = p.parent
        with open(p, encoding="utf-8") as f:
            raw = yaml.safe_load(f) or {}
    elif isinstance(config_or_path, dict):
        raw = copy.deepcopy(config_or_path)
    else:
        raise TypeError(f"Expected str, Path, or dict, got: {type(config_or_path)}")

    # Merge with defaults
    config = _deep_merge(DEFAULT_ULRM_CONFIG, raw)

    # Legacy backward-compatibility mappings from workflow -> preprocessing
    wf = config.get("workflow", {})
    prep = config.setdefault("preprocessing", {})
    if "simplify_hydrosheds_river_network" in wf and "simplify_rivers" not in raw.get(
        "preprocessing", {}
    ):
        prep["simplify_rivers"] = bool(wf["simplify_hydrosheds_river_network"])
    if "process_coastline" in wf and "process_coastline" not in raw.get(
        "preprocessing", {}
    ):
        prep["process_coastline"] = bool(wf["process_coastline"])
    if "largest_outlets" in wf and "largest_outlets" not in raw.get(
        "preprocessing", {}
    ):
        prep["largest_outlets"] = int(wf["largest_outlets"])

    # Normalize date
    if not wf.get("date"):
        wf["date"] = datetime.now().strftime("%Y%m%d")

    # Resolve paths (multi-platform OS dictionaries & relative paths)
    if "paths" in config:
        config["paths"] = _resolve_paths(config["paths"], base_dir)

    return config


def dump_ulrm_config(config: dict[str, Any], path: str | Path) -> Path:
    """Save configuration dictionary to a YAML or JSON file.

    Parameters
    ----------
    config : dict[str, Any]
        Configuration dictionary to serialize.
    path : str or Path
        Target file path (.yaml, .yml, or .json).

    Returns
    -------
    Path
        Path to the written file.
    """
    p = Path(path).expanduser().resolve()
    p.parent.mkdir(parents=True, exist_ok=True)

    if p.suffix.lower() == ".json":
        with open(p, "w", encoding="utf-8") as f:
            json.dump(config, f, indent=2)
    else:
        with open(p, "w", encoding="utf-8") as f:
            yaml.safe_dump(config, f, sort_keys=False, default_flow_style=False)

    return p
