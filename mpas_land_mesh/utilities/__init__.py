"""
Core utilities for ULRM workflow

Note on import structure
------------------------
``object.py`` and ``io.py`` import from ``mpas_land_mesh.classes`` (edge, vertex,
flowline, link), which in turn import back from ``mpas_land_mesh.utilities``.
To break that circular dependency those two modules are **not** imported here at
package-init time; import them directly when needed:

    from mpas_land_mesh.utilities.object import find_vertex_on_edge, ...
    from mpas_land_mesh.utilities.io import export_vertex_to_geojson, ...
"""

# --- constants ---
from .constants import (
    earth_radius,
    IDL_threshold,
)

# --- system ---
from .system import (
    get_extension_from_path,
    get_python_environment,
)

# --- vector ---
from .vector import (
    SUPPORTED_VECTOR_FORMATS,
    gdal_vector_format_support,
    get_available_vector_formats,
    print_supported_vector_formats,
    get_vector_format_from_extension,
    get_vector_format_from_filename,
    get_extension_from_vector_format,
    get_vector_driver_from_format,
    get_vector_driver_from_filename,
    check_parquet_support,
    has_parquet_support,
    get_field_and_value,
    add_field_to_vector_file,
    merge_features,
    write_wkt_to_vector_file,
    remove_small_polygon,
)

# --- raster ---
from .raster import (
    convert_vector_to_global_raster,
    create_raster_buffer_zone,
    fix_raster_antimeridian_issue,
    gdal_read_geotiff_file,
)

# --- gcsbuffer ---
from .gcsbuffer import (
    meters_to_degrees,
    create_geometry_buffer_degrees,
    create_wkt_buffer_degrees,
    create_wkt_buffer_distance,
    create_file_buffer_degrees,
    create_buffer_with_meter_distance,
)

# --- geometry ---
from .geometry import (
    calculate_angle_between_point,
    calculate_angle_between_vectors_degrees,
    calculate_distance_based_on_longitude_latitude,
    calculate_distance_to_plane,
    project_point_onto_plane,
    calculate_intersect_on_great_circle,
    find_great_circle_intersection_with_meridian,
    convert_longitude_latitude_to_sphere_3d,
    convert_sphere_3d_to_longitude_latitude,
    find_minimal_enclosing_polygon,
    calculate_spherical_triangle_area,
    haversine,
    calculate_polygon_area,
    calculate_polygon_file_area,
    spherical_polygon_area,
    convert_360_to_180,
    split_international_date_line_polygon_coordinates,
)

# --- spatial_reference ---
from .spatial_reference import (
    reproject_coordinates,
    reproject_coordinates_batch,
)

# --- config_manager ---
from .config_manager import (
    jigsawcase,
    JigsawConfigManager,
    create_jigsaw_template_configuration_file,
    read_jigsaw_configuration_file,
)

# --- change_json_key_value ---
from .change_json_key_value import (
    change_json_key_value,
    change_json_keys_values,
)

# object and io are intentionally excluded from top-level imports to avoid
# circular imports (they depend on mpas_land_mesh.classes which depends on
# mpas_land_mesh.utilities).  Import them directly from their modules:
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
    # config_manager
    "jigsawcase",
    "JigsawConfigManager",
    "create_jigsaw_template_configuration_file",
    "read_jigsaw_configuration_file",
    # change_json_key_value
    "change_json_key_value",
    "change_json_keys_values",
]
