import os
import json
from pathlib import Path
from mpas_land_mesh.classes.jigsawcase import jigsawcase

class JigsawConfigManager:
    """Configuration manager for JIGSAW that handles defaults and serialization"""

    @staticmethod
    def get_default_config() -> dict:
        """Returns a dictionary with all default JIGSAW configuration values"""
        return {
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
            # Core mesh sizing and quality parameters from JigsawConfigManager
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

    @staticmethod
    def create_template_config(output_filename: str, custom_values: dict = None) -> dict:
        """Create a JIGSAW configuration file with default values, optionally customized

        Args:
            output_filename: Path to save the configuration file
            custom_values: Dictionary of values to override defaults

        Returns:
            The created configuration dictionary
        """
        # Get defaults
        config = JigsawConfigManager.get_default_config()

        # Apply customizations
        if custom_values:
            for key, value in custom_values.items():
                if isinstance(value, Path):
                    value = str(value)
                config[key] = value

        # Ensure directory exists
        os.makedirs(os.path.dirname(os.path.abspath(output_filename)), exist_ok=True)

        # Write to file
        with open(output_filename, "w") as f:
            json.dump(config, f, indent=4)

        return config

    @staticmethod
    def load_config(filename: str) -> dict:
        """Load a JIGSAW configuration from a JSON file

        Args:
            filename: Path to the configuration file

        Returns:
            The loaded configuration dictionary
        """
        with open(filename, "r") as f:
            return json.load(f)

def create_jigsaw_template_configuration_file(sFilename_configuration_json: str, **kwargs) -> dict:
    """Generate JIGSAW config template file using parameter keywords

    Args:
        sFilename_configuration_json: Path to save the configuration file
        **kwargs: Additional configuration parameters to override defaults

    Returns:
        The created configuration dictionary
    """

    # Create the configuration
    config = JigsawConfigManager.create_template_config(
        sFilename_configuration_json, custom_values=kwargs
    )

    # Ensure the output directory exists
    os.makedirs(os.path.dirname(sFilename_configuration_json), exist_ok=True)

    # Write the configuration to the specified JSON file
    with open(sFilename_configuration_json, "w") as f:
        json.dump(config, f, indent=4)

    return config

def read_jigsaw_configuration_file(sFilename_configuration_in: str,
                                   iFlag_standalone_in: int = 1,
    iFlag_create_directory_in: int = None,
    iCase_index_in: int = None,
    sModel_in: str = "jigsaw",
    sDate_in: str = None,
    sWorkspace_output_in: str = None):
    """Read a JIGSAW configuration from a JSON file and create a jigsawcase object

    Args:
        sFilename_configuration_in: Path to the configuration JSON file
        iFlag_standalone_in: Flag for standalone mode (default: 1)
        iFlag_create_directory_in: Flag to create output directory
        iCase_index_in: Case index number
        sModel_in: Model name (default: "jigsaw")
        sDate_in: Date string for the case
        sWorkspace_output_in: Output workspace directory path

    Returns:
        jigsawcase object initialized with the configuration
    """

    # Ensure input filenames are strings
    if isinstance(sFilename_configuration_in, Path):
        sFilename_configuration_in = str(sFilename_configuration_in)
    if isinstance(sWorkspace_output_in, Path):
        sWorkspace_output_in = str(sWorkspace_output_in)

    if not os.path.isfile(sFilename_configuration_in):
        print(sFilename_configuration_in + " does not exist")
        return

    with open(sFilename_configuration_in, "r") as f:
        aConfig = json.load(f)

    if iCase_index_in is not None:
        iCase_index = iCase_index_in
    else:
        iCase_index = int(aConfig["iCase_index"])

    if iFlag_standalone_in is not None:
        iFlag_standalone = iFlag_standalone_in
    else:
        iFlag_standalone = int(aConfig["iFlag_standalone"])

    if sModel_in is not None:
        sModel = sModel_in
    else:
        sModel = aConfig["sModel"]
        pass

    if sDate_in is not None:
        sDate = sDate_in
    else:
        sDate = aConfig["sDate"]
        pass

    if sWorkspace_output_in is not None:
        sWorkspace_output = sWorkspace_output_in
    else:
        sWorkspace_output = aConfig["sWorkspace_output"]
        # try to create this output folder first using

    try:
        print(
            "Creating the specified output workspace (if it does not exist): \n",
            sWorkspace_output,
        )
        Path(sWorkspace_output).mkdir(parents=True, exist_ok=True)
        print("The specified output workspace is: \n", sWorkspace_output)
    except ValueError:
        print("The specified output workspace cannot be created!")
        exit

    aConfig["iCase_index"] = iCase_index
    aConfig["iFlag_standalone"] = iFlag_standalone

    aConfig["sDate"] = sDate
    aConfig["sModel"] = sModel


    aConfig["sWorkspace_output"] = sWorkspace_output

    aConfig["sFilename_model_configuration"] = sFilename_configuration_in

    oJigsaw = jigsawcase(
        aConfig, iFlag_create_directory_in=iFlag_create_directory_in
    )

    return oJigsaw

