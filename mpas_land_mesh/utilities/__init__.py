"""
Core utilities for ULRM workflow

Note on import structure
------------------------
``object.py`` and ``io.py`` import mesh/river classes (edge, vertex, flowline, link)
from ``pyearthmesh`` and ``pyearthriver``. Import them directly when needed:

    from mpas_land_mesh.utilities.object import find_vertex_on_edge, ...
    from mpas_land_mesh.utilities.io import export_vertex_to_geojson, ...
"""

# --- config (pyearthmesh unified) ---
from pyearthmesh.config.config_manager import (
    DEFAULT_CONFIG,
    compile_workflow,
    create_workflow,
    dump_workflow_config,
    load_workflow_config,
)

# --- change_json_key_value ---
from .change_json_key_value import (
    change_json_key_value,
    change_json_keys_values,
)

# --- constants ---
from .constants import (
    IDL_threshold,
    earth_radius,
)

# --- gcsbuffer ---
from .gcsbuffer import (
    create_buffer_with_meter_distance,
    create_file_buffer_degrees,
    create_geometry_buffer_degrees,
    create_wkt_buffer_degrees,
    create_wkt_buffer_distance,
    meters_to_degrees,
)

# --- geometry ---
from .geometry import (
    calculate_angle_between_point,
    calculate_angle_between_vectors_degrees,
    calculate_distance_based_on_longitude_latitude,
    calculate_distance_to_plane,
    calculate_intersect_on_great_circle,
    calculate_polygon_area,
    calculate_polygon_file_area,
    calculate_spherical_triangle_area,
    convert_360_to_180,
    convert_longitude_latitude_to_sphere_3d,
    convert_sphere_3d_to_longitude_latitude,
    find_great_circle_intersection_with_meridian,
    find_minimal_enclosing_polygon,
    haversine,
    project_point_onto_plane,
    spherical_polygon_area,
    split_international_date_line_polygon_coordinates,
)

# --- raster ---
from .raster import (
    convert_vector_to_global_raster,
    create_raster_buffer_zone,
    fix_raster_antimeridian_issue,
    gdal_read_geotiff_file,
)

# --- spatial_reference ---
from .spatial_reference import (
    reproject_coordinates,
    reproject_coordinates_batch,
)

# --- system ---
from .system import (
    get_extension_from_path,
    get_python_environment,
)

# --- vector ---
from .vector import (
    SUPPORTED_VECTOR_FORMATS,
    add_field_to_vector_file,
    check_parquet_support,
    gdal_vector_format_support,
    get_available_vector_formats,
    get_extension_from_vector_format,
    get_field_and_value,
    get_vector_driver_from_filename,
    get_vector_driver_from_format,
    get_vector_format_from_extension,
    get_vector_format_from_filename,
    has_parquet_support,
    merge_features,
    print_supported_vector_formats,
    remove_small_polygon,
    write_wkt_to_vector_file,
)

# object and io are imported directly from their modules when needed:
#   from mpas_land_mesh.utilities.object import find_vertex_on_edge, slerp, ...
#   from mpas_land_mesh.utilities.io import export_vertex_to_geojson, ...

__all__ = [
    # constants
    "earth_radius",
    "IDL_threshold",
    # system
    "get_extension_from_path",
    # vector – format/driver helpers
    "SUPPORTED_VECTOR_FORMATS",
    "gdal_vector_format_support",
    "get_available_vector_formats",
    "print_supported_vector_formats",
    "get_vector_format_from_extension",
    "get_vector_format_from_filename",
    "get_extension_from_vector_format",
    "get_vector_driver_from_format",
    "get_vector_driver_from_filename",
    "check_parquet_support",
    "has_parquet_support",
    # vector – feature operations
    "get_field_and_value",
    "add_field_to_vector_file",
    "merge_features",
    "write_wkt_to_vector_file",
    "remove_small_polygon",
    # raster
    "convert_vector_to_global_raster",
    "create_raster_buffer_zone",
    "fix_raster_antimeridian_issue",
    # gcsbuffer
    "meters_to_degrees",
    "create_geometry_buffer_degrees",
    "create_wkt_buffer_degrees",
    "create_wkt_buffer_distance",
    "create_file_buffer_degrees",
    "create_buffer_with_meter_distance",
    # geometry
    "calculate_angle_between_point",
    "calculate_angle_between_vectors_degrees",
    "calculate_distance_based_on_longitude_latitude",
    "calculate_distance_to_plane",
    "project_point_onto_plane",
    "calculate_intersect_on_great_circle",
    "find_great_circle_intersection_with_meridian",
    "convert_longitude_latitude_to_sphere_3d",
    "convert_sphere_3d_to_longitude_latitude",
    "find_minimal_enclosing_polygon",
    "calculate_spherical_triangle_area",
    "haversine",
    "calculate_polygon_area",
    "calculate_polygon_file_area",
    "spherical_polygon_area",
    # spatial_reference
    "reproject_coordinates",
    "reproject_coordinates_batch",
    # config (pyearthmesh unified)
    "DEFAULT_CONFIG",
    "load_workflow_config",
    "dump_workflow_config",
    "create_workflow",
    "compile_workflow",
    # change_json_key_value
    "change_json_key_value",
    "change_json_keys_values",
]
