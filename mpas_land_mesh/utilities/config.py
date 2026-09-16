"""Unified configuration manager for mpas_land_mesh.

This module unifies high-level user workflow configurations (workflow options,
preprocessing flags, resolutions in km, cross-platform paths, HPC job settings)
and low-level JIGSAW mesh generator decks into a cohesive two-tier architecture:

1. User Tier (WorkflowConfig):
   Loads human-readable JSON files, validates options against defaults, resolves
   platform-specific paths, and provides dictionary-compatible access.
2. Engine Tier (JigsawConfigManager & build_jigsaw_config):
   Derives and compiles the detailed 50+ parameter JIGSAW input deck in-memory,
   translating physical kilometer resolutions into grid dimensions and degrees,
   setting feature flags, mapping preprocessing output paths, and applying expert
   overrides without requiring repetitive on-disk JSON file modifications.
"""

from __future__ import annotations

import copy
import json
import os
import platform
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping


# ---------------------------------------------------------------------------
# Default Configurations
# ---------------------------------------------------------------------------

DEFAULT_WORKFLOW_CONFIG: dict[str, Any] = {
    "workflow": {
        "date": None,
        "mesh_type": "mpas",
        "case_index": 1,
        "model": "jigsaw",
        "simplify_hydrosheds_river_network": False,
        "process_coastline": False,
        "debug": True,
        "largest_outlets": 10,
    },
    "resolutions": {
        "ocean_km": 100,
        "land_km": 100,
        "river_network_km": 10,
        "coastline_km": 10,
        "lake_boundary_km": None,
        "watershed_boundary_km": None,
    },
    "paths": {
        "input_workspace": "",
        "output_workspace": "./output",
        "hydrosheds_rivers": "",
        "region_geometry": "",
        "dam_vector": "",
        "lake_boundary": "",
        "watershed_boundary": "",
    },
    "jigsaw": {
        "standalone": True,
        "create_directory": True,
        "hours": 5,
        "slurm": "slurm",
        "overrides": {},
    },
}

# Alias for backward compatibility
DEFAULT_CONFIG = DEFAULT_WORKFLOW_CONFIG

DEFAULT_JIGSAW_CONFIG: dict[str, Any] = {
    # Grid resolution parameters
    "ncolumn_space": 360,  # Number of columns in spacing grid (longitude)
    "nrow_space": 180,  # Number of rows in spacing grid (latitude)
    "dSpac_value": 100.0,  # Default spacing value
    # Feature flags for geometry generation
    "iFlag_geom": False,  # Enable geometry generation
    "iFlag_spac": False,  # Enable spacing function generation
    "iFlag_init": False,  # Enable initialization mesh
    "iFlag_opts": False,  # Enable custom options
    # Environment type flags
    "iFlag_spac_ocean": False,  # Apply ocean-specific spacing
    "iFlag_spac_land": True,  # Apply land-specific spacing
    # Point feature geometry and spacing flags
    "iFlag_geom_dam": False,  # Include dam geometries
    "iFlag_spac_dam": False,  # Apply dam-specific spacing
    "iFlag_geom_city": False,  # Include city geometries
    "iFlag_spac_city": False,  # Apply city-specific spacing
    # Line feature geometry and spacing flags
    "iFlag_geom_river_network": False,  # Include river network geometries
    "iFlag_spac_river_network": False,  # Apply river network-specific spacing
    "iFlag_geom_coastline": False,  # Include coastline geometries
    "iFlag_spac_coastline": False,  # Apply coastline-specific spacing
    # Polygon feature geometry and spacing flags
    "iFlag_geom_watershed_boundary": False,  # Include watershed boundary geometries
    "iFlag_spac_watershed_boundary": False,  # Apply watershed boundary-specific spacing
    "iFlag_geom_lake_boundary": False,  # Include lake boundary geometries
    "iFlag_spac_lake_boundary": False,  # Apply lake boundary-specific spacing
    # Resolution parameters for different features (in degrees)
    "dResolution_land": 45.0,  # Land feature spacing
    "dResolution_dam": 4.0,  # Dam feature spacing
    "dResolution_city": 4.0,  # City feature spacing
    "dResolution_river_network": 4.0,  # River network feature spacing
    "dResolution_coastline": 4.0,  # Coastline feature spacing
    "dResolution_watershed_boundary": 4.0,  # Watershed boundary feature spacing
    "dResolution_lake_boundary": 4.0,  # Lake boundary feature spacing
    # Mesh type identifiers
    "geom_mshID": "ellipsoid-mesh",  # Geometry mesh type
    "spac_mshID": "ellipsoid-grid",  # Spacing grid type
    # Earth/sphere parameters
    "FULL_SPHERE_RADIUS": 6371.0,  # Earth radius in km
    # Gradient limiting
    "dhdx_lim": 0.25,  # Gradient limit for mesh sizing
    # File paths for features
    "sFilename_dam_vector": None,  # Dam vector file
    "sFilename_dam_raster": None,  # Dam raster file
    "sFilename_city_vector": None,  # City vector file
    "sFilename_city_raster": None,  # City raster file
    "sFilename_river_network_vector": None,  # River network vector file
    "sFilename_river_network_raster": None,  # River network raster file
    "sFilename_coastline_vector": None,  # Coastline vector file
    "sFilename_coastline_raster": None,  # Coastline raster file
    "sFilename_watershed_boundary_vector": None,  # Watershed boundary vector file
    "sFilename_watershed_boundary_raster": None,  # Watershed boundary raster file
    "sFilename_lake_boundary_vector": None,  # Lake boundary vector file
    "sFilename_lake_boundary_raster": None,  # Lake boundary raster file
    # Mesh sizing parameters
    "hfun_hmax": "inf",  # Max. refinement function value
    "hfun_hmin": 0.0,  # Min. refinement function value
    "hfun_scal": "absolute",  # Scaling type: "relative" or "absolute"
    "mesh_dims": 2,  # Mesh dimension (2 for surface)
    "bisection": -1,  # Bisection method (-1 for heuristic)
    # Optimization parameters
    "optm_qlim": 0.95,  # Quality limit for optimization
    "optm_iter": 32,  # Number of optimization iterations
    "optm_qtol": 1.0e-05,  # Quality tolerance
    # Core mesh sizing and quality parameters
    "mesh_rad2": 1.5,  # Max. radius-edge ratio
    "mesh_rad3": 2.0,  # Max. radius-circumsphere ratio for tetras
    "mesh_eps1": 0.333,  # Min. mesh quality threshold
    "mesh_eps2": 0.333,  # Min. mesh quality threshold for tetra
    "mesh_top": 1,  # Mesh topology (1 for manifold surface)
    "mesh_iter": 3,  # Mesh iteration limit
    # Verbosity and iterations
    "verbosity": 0,  # Verbosity level (0-3)
    # File paths (will be populated based on workspace)
    "geom_file": None,  # Input geometry file
    "hfun_file": None,  # Input mesh-size file
    "mesh_file": None,  # Output mesh file
    # Algorithm selection
    "mesh_kern": "delfront",  # Meshing kernel: "delfront" or "delaunay"
    "optm_kern": "odt+dqdx",  # Optimisation kernel
    # Region boundary
    "geom_feat": True,  # Detect sharp features in geometry
    # Output options
    "mesh_type": "euclidean-mesh",  # Mesh type (euclidean-mesh or ellipsoid-mesh)
    "output_formats": ["vtk", "gmsh"],  # Output formats to generate
}


# ---------------------------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------------------------

def _merge(defaults: dict[str, Any], values: dict[str, Any], section: str) -> dict[str, Any]:
    """Merge user values with default settings, validating against known keys."""
    # Allow arbitrary expert settings in the jigsaw section
    if section != "jigsaw":
        unknown = set(values) - set(defaults)
        if unknown:
            raise ValueError(f"Unknown {section} setting(s): {', '.join(sorted(unknown))}")
    result = defaults.copy()
    result.update(values)
    return result


def _select_platform_path(value: Any) -> str:
    """Resolve a path value that may be a string or a platform-keyed mapping."""
    if isinstance(value, (str, Path)):
        return str(value)
    if not isinstance(value, dict):
        raise ValueError("Path values must be strings or platform-specific objects")

    current_platform = platform.system()
    selected = value.get(current_platform, value.get("default", ""))
    if not isinstance(selected, (str, Path)):
        raise ValueError(f"Path for {current_platform} must be a string")
    return str(selected)


# ---------------------------------------------------------------------------
# WorkflowConfig Class
# ---------------------------------------------------------------------------

class WorkflowConfig(dict):
    """High-level workflow configuration supporting validation, path resolution,
    and automatic JIGSAW deck compilation.

    Inherits from dict for 100% backward compatibility with dictionary indexing
    (e.g., config['workflow']), while providing property accessors and high-level
    utility methods.
    """

    @classmethod
    def load(cls, filename: str | Path) -> WorkflowConfig:
        """Load a workflow JSON file, apply defaults, and resolve relative paths.

        Args:
            filename: Path to the workflow JSON configuration file.

        Returns:
            WorkflowConfig instance populated with resolved settings.
        """
        filename = Path(filename).expanduser().resolve()
        with filename.open(encoding="utf-8") as config_file:
            raw = json.load(config_file)
        if not isinstance(raw, dict):
            raise ValueError("Workflow configuration must contain a JSON object")

        return cls.from_dict(raw, config_file=filename)

    @classmethod
    def from_dict(
        cls,
        raw: dict[str, Any],
        config_file: str | Path | None = None,
    ) -> WorkflowConfig:
        """Create a WorkflowConfig from a raw dictionary.

        Args:
            raw: Dictionary containing workflow settings.
            config_file: Optional path to the configuration file for resolving relative paths.

        Returns:
            WorkflowConfig instance populated with defaults and resolved paths.
        """
        cfg_file_path = Path(config_file).expanduser().resolve() if config_file else None
        base_dir = cfg_file_path.parent if cfg_file_path else None

        merged: dict[str, Any] = {}
        for section in DEFAULT_WORKFLOW_CONFIG:
            merged[section] = _merge(
                DEFAULT_WORKFLOW_CONFIG[section],
                raw.get(section, {}),
                section,
            )

        if cfg_file_path:
            merged["config_file"] = cfg_file_path

        # Resolve paths
        for key, value in merged["paths"].items():
            if value is None:
                merged["paths"][key] = ""
                continue
            resolved_str = _select_platform_path(value)
            if resolved_str:
                p = Path(resolved_str).expanduser()
                if not p.is_absolute() and base_dir is not None:
                    p = base_dir / p
                merged["paths"][key] = str(p.resolve())
            else:
                merged["paths"][key] = ""

        # Ensure date is set
        if not merged["workflow"]["date"]:
            merged["workflow"]["date"] = datetime.now().strftime("%Y%m%d")

        return cls(merged)

    # -----------------------------------------------------------------------
    # Property Accessors
    # -----------------------------------------------------------------------

    @property
    def workflow(self) -> dict[str, Any]:
        """Workflow metadata and preprocessing toggles."""
        return self.setdefault("workflow", {})

    @property
    def resolutions(self) -> dict[str, Any]:
        """Target resolutions in physical units (km)."""
        return self.setdefault("resolutions", {})

    @property
    def paths(self) -> dict[str, Any]:
        """Input datasets and output directory paths."""
        return self.setdefault("paths", {})

    @property
    def jigsaw(self) -> dict[str, Any]:
        """JIGSAW run and HPC scheduling options."""
        return self.setdefault("jigsaw", {})

    @property
    def config_file(self) -> Path | None:
        """Path to the loaded configuration file, if loaded from disk."""
        val = self.get("config_file")
        return Path(val) if val else None

    # -----------------------------------------------------------------------
    # JIGSAW Deck Compilation
    # -----------------------------------------------------------------------

    def build_jigsaw_config(
        self,
        generated_files: Mapping[str, str | Path] | None = None,
        output_workspace: str | Path | None = None,
        **overrides: Any,
    ) -> dict[str, Any]:
        """Compile high-level workflow settings into a complete JIGSAW configuration deck.

        Translates physical kilometer resolutions into grid dimensions and degrees,
        configures geometry and spacing feature flags, maps preprocessed intermediate
        dataset paths, and applies any expert overrides.

        Args:
            generated_files: Mapping of preprocessed file keys to paths. Supported keys:
                - "river_network_vector": Path to simplified river network vector
                - "river_network_raster": Path to river network raster TIFF
                - "coastline_raster": Path to coastline / land-ocean mask raster TIFF
                - "land_ocean_mask": Path to land-ocean mask vector
                - "dam_vector": Override path to dam vector file
                - "lake_boundary_vector": Path to lake boundary vector
                - "lake_boundary_raster": Path to lake boundary raster
                - "watershed_boundary_vector": Path to watershed boundary vector
                - "watershed_boundary_raster": Path to watershed boundary raster
            output_workspace: Optional override for the case output directory.
                Defaults to paths["output_workspace"].
            **overrides: Arbitrary JIGSAW parameters to override (e.g. optm_iter=64).

        Returns:
            Dictionary containing full JIGSAW configuration ready for jigsawcase / run_jigsaw.
        """
        deck = copy.deepcopy(DEFAULT_JIGSAW_CONFIG)
        gen = {k: str(v) for k, v in (generated_files or {}).items()}

        # 1. Grid resolution and spacing calculations
        coastline_km = float(self.resolutions.get("coastline_km", 10.0))
        dResolution_x_in = 30.0 / 3600.0 * coastline_km
        dResolution_y_in = dResolution_x_in
        nrow = int(180.0 / dResolution_y_in)
        ncolumn = int(360.0 / dResolution_x_in)

        deck["ncolumn_space"] = ncolumn
        deck["nrow_space"] = nrow
        deck["dResolution_coastline"] = coastline_km

        if "ocean_km" in self.resolutions and self.resolutions["ocean_km"] is not None:
            deck["dResolution_ocean"] = float(self.resolutions["ocean_km"])
        if "land_km" in self.resolutions and self.resolutions["land_km"] is not None:
            deck["dResolution_land"] = float(self.resolutions["land_km"])
        if "river_network_km" in self.resolutions and self.resolutions["river_network_km"] is not None:
            deck["dResolution_river_network"] = float(self.resolutions["river_network_km"])
        if self.resolutions.get("lake_boundary_km") is not None:
            deck["dResolution_lake_boundary"] = float(self.resolutions["lake_boundary_km"])
        if self.resolutions.get("watershed_boundary_km") is not None:
            deck["dResolution_watershed_boundary"] = float(self.resolutions["watershed_boundary_km"])

        # 2. General meshing flags
        deck["iFlag_geom"] = True
        deck["iFlag_spac"] = True
        deck["iFlag_spac_ocean"] = True
        deck["iFlag_spac_land"] = True

        # 3. Coastline configuration
        process_coastline = bool(self.workflow.get("process_coastline", False))
        if process_coastline or "coastline_raster" in gen:
            deck["iFlag_spac_coastline"] = True
        if "coastline_raster" in gen:
            deck["sFilename_coastline_raster"] = gen["coastline_raster"]

        # 4. River network configuration
        simplify_rivers = bool(self.workflow.get("simplify_hydrosheds_river_network", False))
        has_rivers = (
            simplify_rivers
            or "river_network_vector" in gen
            or "river_network_raster" in gen
        )
        if has_rivers:
            deck["iFlag_geom_river_network"] = True
            deck["iFlag_spac_river_network"] = True
        if "river_network_vector" in gen:
            deck["sFilename_river_network_vector"] = gen["river_network_vector"]
        if "river_network_raster" in gen:
            deck["sFilename_river_network_raster"] = gen["river_network_raster"]

        # 5. Dam configuration
        dam_path = gen.get("dam_vector") or self.paths.get("dam_vector")
        if dam_path and os.path.isfile(str(dam_path)):
            deck["iFlag_geom_dam"] = True
            deck["sFilename_dam_vector"] = str(dam_path)

        # 6. Land ocean mask & other features
        if "land_ocean_mask" in gen:
            deck["sFilename_land_ocean_mask"] = gen["land_ocean_mask"]

        if "lake_boundary_vector" in gen:
            deck["iFlag_geom_lake_boundary"] = True
            deck["sFilename_lake_boundary_vector"] = gen["lake_boundary_vector"]
        if "lake_boundary_raster" in gen:
            deck["iFlag_spac_lake_boundary"] = True
            deck["sFilename_lake_boundary_raster"] = gen["lake_boundary_raster"]

        if "watershed_boundary_vector" in gen:
            deck["iFlag_geom_watershed_boundary"] = True
            deck["sFilename_watershed_boundary_vector"] = gen["watershed_boundary_vector"]
        if "watershed_boundary_raster" in gen:
            deck["iFlag_spac_watershed_boundary"] = True
            deck["sFilename_watershed_boundary_raster"] = gen["watershed_boundary_raster"]

        # 7. Metadata and workspace paths
        out_ws = str(output_workspace) if output_workspace else self.paths.get("output_workspace", "./output")
        deck["sWorkspace_output"] = str(out_ws)
        deck["iCase_index"] = int(self.workflow.get("case_index", 1))
        deck["sDate"] = str(self.workflow.get("date", datetime.now().strftime("%Y%m%d")))
        deck["sModel"] = str(self.workflow.get("model", "jigsaw"))
        deck["iFlag_standalone"] = 1 if self.jigsaw.get("standalone", True) else 0

        # 8. Expert overrides from jigsaw.overrides or top-level jigsaw options
        jigsaw_cfg = self.jigsaw
        # If user put an explicit "overrides" mapping in the jigsaw section
        if "overrides" in jigsaw_cfg and isinstance(jigsaw_cfg["overrides"], dict):
            deck.update(jigsaw_cfg["overrides"])
        # If user put known JIGSAW keys directly in the jigsaw section
        for key, val in jigsaw_cfg.items():
            if key not in ("standalone", "create_directory", "hours", "slurm", "overrides"):
                deck[key] = val

        # 9. Explicit keyword argument overrides
        deck.update(overrides)

        return deck

    def save(self, filename: str | Path) -> None:
        """Save this workflow configuration to a JSON file."""
        target = Path(filename).expanduser().resolve()
        target.parent.mkdir(parents=True, exist_ok=True)
        data = dict(self)
        if "config_file" in data:
            data["config_file"] = str(data["config_file"])
        with open(target, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=4)

    def to_dict(self) -> dict[str, Any]:
        """Convert to a plain Python dictionary."""
        return dict(self)


# ---------------------------------------------------------------------------
# JigsawConfigManager Class
# ---------------------------------------------------------------------------

class JigsawConfigManager:
    """Configuration manager for JIGSAW that handles defaults, templates, and serialization."""

    @staticmethod
    def get_default_config() -> dict[str, Any]:
        """Returns a copy of the default JIGSAW configuration dictionary."""
        return copy.deepcopy(DEFAULT_JIGSAW_CONFIG)

    @staticmethod
    def create_template_config(
        output_filename: str | Path,
        custom_values: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        """Create a JIGSAW configuration file with default values, optionally customized.

        Args:
            output_filename: Path to save the configuration file.
            custom_values: Optional dictionary of values to override defaults.

        Returns:
            The created configuration dictionary.
        """
        config = JigsawConfigManager.get_default_config()

        if custom_values:
            for key, value in custom_values.items():
                if isinstance(value, Path):
                    value = str(value)
                config[key] = value

        target = Path(output_filename).expanduser().resolve()
        target.parent.mkdir(parents=True, exist_ok=True)

        with open(target, "w", encoding="utf-8") as f:
            json.dump(config, f, indent=4)

        return config

    @staticmethod
    def load_config(filename: str | Path) -> dict[str, Any]:
        """Load a JIGSAW configuration from a JSON file.

        Args:
            filename: Path to the configuration file.

        Returns:
            The loaded configuration dictionary.
        """
        with open(Path(filename).expanduser().resolve(), "r", encoding="utf-8") as f:
            return json.load(f)


# ---------------------------------------------------------------------------
# Procedural & Convenience Functions
# ---------------------------------------------------------------------------

def load_workflow_config(filename: str | Path) -> WorkflowConfig:
    """Load a workflow JSON file, apply defaults, and resolve relative paths.

    Args:
        filename: Path to the workflow JSON file.

    Returns:
        WorkflowConfig instance (acts as a dict for drop-in backward compatibility).
    """
    return WorkflowConfig.load(filename)


def create_jigsaw_template_configuration_file(
    sFilename_configuration_json: str | Path,
    **kwargs: Any,
) -> dict[str, Any]:
    """Generate a JIGSAW config template file using keyword overrides.

    Args:
        sFilename_configuration_json: Path to save the configuration file.
        **kwargs: Additional configuration parameters to override defaults.

    Returns:
        The created configuration dictionary.
    """
    return JigsawConfigManager.create_template_config(
        sFilename_configuration_json,
        custom_values=kwargs,
    )


def read_jigsaw_configuration_file(
    sFilename_configuration_in: str | Path,
    iFlag_standalone_in: int = 1,
    iFlag_create_directory_in: int | None = None,
    iCase_index_in: int | None = None,
    sModel_in: str = "jigsaw",
    sDate_in: str | None = None,
    sWorkspace_output_in: str | Path | None = None,
):
    """Read a JIGSAW configuration from a JSON file and create a jigsawcase object.

    Args:
        sFilename_configuration_in: Path to the configuration JSON file.
        iFlag_standalone_in: Flag for standalone mode (default: 1).
        iFlag_create_directory_in: Flag to create output directory.
        iCase_index_in: Case index number.
        sModel_in: Model name (default: "jigsaw").
        sDate_in: Date string for the case.
        sWorkspace_output_in: Output workspace directory path.

    Returns:
        jigsawcase object initialized with the configuration.
    """
    from mpas_land_mesh.classes.jigsawcase import jigsawcase

    path_in = Path(sFilename_configuration_in).expanduser().resolve()
    if not path_in.is_file():
        print(f"{path_in} does not exist")
        return None

    with open(path_in, "r", encoding="utf-8") as f:
        aConfig = json.load(f)

    if iCase_index_in is not None:
        iCase_index = iCase_index_in
    else:
        iCase_index = int(aConfig.get("iCase_index", 1))

    if iFlag_standalone_in is not None:
        iFlag_standalone = iFlag_standalone_in
    else:
        iFlag_standalone = int(aConfig.get("iFlag_standalone", 1))

    sModel = sModel_in if sModel_in is not None else aConfig.get("sModel", "jigsaw")
    sDate = sDate_in if sDate_in is not None else aConfig.get("sDate")

    if sWorkspace_output_in is not None:
        sWorkspace_output = str(sWorkspace_output_in)
    else:
        sWorkspace_output = aConfig.get("sWorkspace_output", "./output")

    try:
        Path(sWorkspace_output).mkdir(parents=True, exist_ok=True)
    except Exception as e:
        print(f"The specified output workspace cannot be created: {e}")

    aConfig["iCase_index"] = iCase_index
    aConfig["iFlag_standalone"] = iFlag_standalone
    aConfig["sDate"] = sDate
    aConfig["sModel"] = sModel
    aConfig["sWorkspace_output"] = sWorkspace_output
    aConfig["sFilename_model_configuration"] = str(path_in)

    return jigsawcase(aConfig, iFlag_create_directory_in=iFlag_create_directory_in)


def create_jigsaw_case(
    config: dict[str, Any] | WorkflowConfig,
    generated_files: Mapping[str, str | Path] | None = None,
    output_workspace: str | Path | None = None,
    iFlag_create_directory_in: int = 1,
    **kwargs: Any,
):
    """Create and return an initialized jigsawcase instance.

    Accepts either a high-level WorkflowConfig (compiled automatically via
    build_jigsaw_config) or a low-level JIGSAW config dictionary.

    Args:
        config: WorkflowConfig or JIGSAW deck dictionary.
        generated_files: Optional mapping of preprocessed file paths.
        output_workspace: Optional output directory override.
        iFlag_create_directory_in: Flag to create output directory (default: 1).
        **kwargs: Additional parameters passed to build_jigsaw_config or jigsawcase.

    Returns:
        Initialized jigsawcase instance.
    """
    from mpas_land_mesh.classes.jigsawcase import jigsawcase

    if isinstance(config, WorkflowConfig):
        jigsaw_deck = config.build_jigsaw_config(
            generated_files=generated_files,
            output_workspace=output_workspace,
            **kwargs,
        )
    elif isinstance(config, dict):
        jigsaw_deck = copy.deepcopy(config)
        if kwargs:
            jigsaw_deck.update(kwargs)
    else:
        raise TypeError(f"Expected WorkflowConfig or dict, got {type(config).__name__}")

    return jigsawcase(jigsaw_deck, iFlag_create_directory_in=iFlag_create_directory_in)

