"""
Core utilities for ULRM workflow
"""

# --- config (pyearthmesh unified) ---
from pyearthmesh.config.config_manager import (
    DEFAULT_CONFIG,
    compile_workflow,
    create_workflow,
    dump_workflow_config,
    load_workflow_config,
)

# --- change_json_key_value (pyearthmesh unified) ---
from pyearthmesh.utility.change_json_key_value import (
    change_json_key_value,
    change_json_keys_values,
)

# --- constants ---
from .constants import (
    DRAINAGE_AREA_MULTIPLIER,
    IDL_OFFSET,
    IDL_TOLERANCE,
    ISLAND_AREA_MULTIPLIER,
    KM2_TO_M2,
    IDL_threshold,
    earth_radius,
)

# --- mesh_stats ---
from .mesh_stats import get_mesh_stats

__all__ = [
    # constants
    "DRAINAGE_AREA_MULTIPLIER",
    "IDL_OFFSET",
    "IDL_TOLERANCE",
    "IDL_threshold",
    "ISLAND_AREA_MULTIPLIER",
    "KM2_TO_M2",
    "earth_radius",
    # mesh_stats
    "get_mesh_stats",
    # config (pyearthmesh unified)
    "DEFAULT_CONFIG",
    "load_workflow_config",
    "dump_workflow_config",
    "create_workflow",
    "compile_workflow",
    # change_json_key_value (pyearthmesh unified)
    "change_json_key_value",
    "change_json_keys_values",
]
