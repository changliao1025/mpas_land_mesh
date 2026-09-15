"""Load and validate JSON configuration for the example workflow."""

import json
import platform
from pathlib import Path
from typing import Any


DEFAULT_CONFIG = {
    "workflow": {"date": None, "mesh_type": "mpas", "case_index": 1, "model": "jigsaw",
                 "simplify_hydrosheds_river_network": False, "process_coastline": False,
                 "debug": True, "largest_outlets": 10},
    "resolutions": {"ocean_km": 100, "land_km": 100, "river_network_km": 10, "coastline_km": 10},
    "paths": {"input_workspace": "", "output_workspace": "./output", "hydrosheds_rivers": "",
              "region_geometry": "", "dam_vector": ""},
    "jigsaw": {"standalone": True, "create_directory": True, "hours": 5, "slurm": "slurm"},
}


def _merge(defaults: dict[str, Any], values: dict[str, Any], section: str) -> dict[str, Any]:
    unknown = set(values) - set(defaults)
    if unknown:
        raise ValueError(f"Unknown {section} setting(s): {', '.join(sorted(unknown))}")
    result = defaults.copy()
    result.update(values)
    return result


def _select_platform_path(value: Any) -> str:
    if isinstance(value, str):
        return value
    if not isinstance(value, dict):
        raise ValueError("Path values must be strings or platform-specific objects")

    current_platform = platform.system()
    selected = value.get(current_platform, value.get("default", ""))
    if not isinstance(selected, str):
        raise ValueError(f"Path for {current_platform} must be a string")
    return selected


def load_workflow_config(filename: str | Path) -> dict[str, Any]:
    """Load a workflow JSON file, apply defaults, and resolve relative paths."""
    filename = Path(filename).expanduser().resolve()
    with filename.open(encoding="utf-8") as config_file:
        raw = json.load(config_file)
    if not isinstance(raw, dict):
        raise ValueError("Workflow configuration must contain a JSON object")
    config = {section: _merge(DEFAULT_CONFIG[section], raw.get(section, {}), section)
              for section in DEFAULT_CONFIG}
    config["config_file"] = filename
    for key, value in config["paths"].items():
        value = _select_platform_path(value)
        if value:
            path = Path(value).expanduser()
            if not path.is_absolute():
                path = filename.parent / path
            config["paths"][key] = str(path.resolve())
        else:
            config["paths"][key] = ""
    if not config["workflow"]["date"]:
        from datetime import datetime
        config["workflow"]["date"] = datetime.now().strftime("%Y%m%d")
    return config