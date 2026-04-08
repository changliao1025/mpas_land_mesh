"""
Simple geometry buffer utilities

Provides simplified buffer operations without geodesic complexity
"""

from osgeo import ogr, osr

from mpas_land_mesh.utilities.constants import earth_radius

def meters_to_degrees(distance_meters, latitude=0.0):
    """
    Convert meters to approximate degrees at given latitude

    This is an approximation:
    - At equator: 1 degree ≈ 111 km
    - Latitude matters for longitude but not for latitude

    Args:
        distance_meters (float): Distance in meters
        latitude (float): Latitude in degrees (for longitude conversion)

    Returns:
        float: Approximate distance in degrees
    """
    import math

    # 1 degree latitude = ~111,320 meters (constant)

    #use earth_radius instead of 111320 to be more accurate
    meters_per_degree = 2 * math.pi * earth_radius / 360.0

    degrees_lat = distance_meters / meters_per_degree

    # For longitude, adjust by latitude
    # 1 degree longitude = meters_per_degree * cos(latitude) meters
    if latitude != 0.0:
        degrees_lon = distance_meters / (meters_per_degree * math.cos(math.radians(latitude)))
        # Return average for buffer (which expands in all directions)
        return (degrees_lat + degrees_lon) / 2.0

    return degrees_lat


def create_geometry_buffer_degrees(geometry, distance_degrees):
    """
    Create a simple buffer around a geometry using planar approximation

    Args:
        geometry (ogr.Geometry): Input geometry
        distance_degrees (float): Buffer distance in degrees

    Returns:
        ogr.Geometry: Buffered geometry
    """
    return geometry.Buffer(distance_degrees)


def create_wkt_buffer_degrees(wkt, distance_degrees, epsg=4326):
    """
    Create a simple buffer from WKT string

    Args:
        wkt (str): Input geometry as WKT
        distance_degrees (float): Buffer distance in degrees
        epsg (int): EPSG code (default: 4326)

    Returns:
        str: Buffered geometry as WKT
    """
    geometry = ogr.CreateGeometryFromWkt(wkt)
    buffered = geometry.Buffer(distance_degrees)
    return buffered.ExportToWkt()


def create_wkt_buffer_distance(wkt, distance_meter, epsg=4326, reference_latitude=0.0):
    """
    Create a simple buffer from WKT string

    Args:
        wkt (str): Input geometry as WKT
        distance_meters (float): Buffer distance in meters
        epsg (int): EPSG code (default: 4326)

    Returns:
        str: Buffered geometry as WKT
    """
    geometry = ogr.CreateGeometryFromWkt(wkt)
    distance_degrees = meters_to_degrees(distance_meter, reference_latitude)
    buffered = geometry.Buffer(distance_degrees)
    return buffered.ExportToWkt()

def create_file_buffer_degrees(input_file, output_file, distance_degrees):
    """
    Create buffer around all features in a file using simple planar method

    This is simpler and faster than geodesic buffering, suitable for:
    - Small buffer distances (< 1 degree)
    - Mid-latitude regions (not polar)
    - Quick prototyping

    For accurate geodesic buffers on a sphere, use the full
    pyearth.toolbox.geometry.create_gcs_buffer_zone function.

    Args:
        input_file (str): Input vector file path
        output_file (str): Output vector file path
        distance_degrees (float): Buffer distance in decimal degrees

    Example:
        # Buffer by ~10km at equator (approximately 0.09 degrees)
        create_simple_buffer_from_file('rivers.geojson', 'rivers_buffered.geojson', 0.09)
    """
    import os

    # Open input
    pDataset_in = ogr.Open(input_file)
    if pDataset_in is None:
        raise ValueError(f"Could not open {input_file}")

    pLayer_in = pDataset_in.GetLayer(0)
    pSpatialRef = pLayer_in.GetSpatialRef()
    geom_type = pLayer_in.GetGeomType()

    # Determine output geometry type
    if geom_type in [ogr.wkbLineString, ogr.wkbMultiLineString]:
        output_geom_type = ogr.wkbPolygon
    elif geom_type in [ogr.wkbPoint, ogr.wkbMultiPoint]:
        output_geom_type = ogr.wkbPolygon
    else:
        output_geom_type = ogr.wkbPolygon

    # Create output
    pDriver = ogr.GetDriverByName('GeoJSON')
    if os.path.exists(output_file):
        pDriver.DeleteDataSource(output_file)

    pDataset_out = pDriver.CreateDataSource(output_file)
    pLayer_out = pDataset_out.CreateLayer('buffered', pSpatialRef, output_geom_type)

    # Copy field definitions
    pLayerDefn_in = pLayer_in.GetLayerDefn()
    for i in range(pLayerDefn_in.GetFieldCount()):
        pFieldDefn = pLayerDefn_in.GetFieldDefn(i)
        pLayer_out.CreateField(pFieldDefn)

    # Buffer each feature
    for pFeature_in in pLayer_in:
        pGeometry = pFeature_in.GetGeometryRef()
        if pGeometry is None:
            continue

        # Create buffer
        pGeometry_buffered = pGeometry.Buffer(distance_degrees)

        # Create output feature
        pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
        pFeature_out.SetGeometry(pGeometry_buffered)

        # Copy attributes
        for i in range(pLayerDefn_in.GetFieldCount()):
            pFeature_out.SetField(
                pLayerDefn_in.GetFieldDefn(i).GetName(),
                pFeature_in.GetField(i)
            )

        pLayer_out.CreateFeature(pFeature_out)
        pFeature_out = None

    pDataset_in = None
    pDataset_out = None

    print(f"Created buffer: {output_file}")
    print(f"  Buffer distance: {distance_degrees} degrees")


def create_buffer_with_meter_distance(input_file, output_file, distance_meters,
                                      reference_latitude=0.0):
    """
    Create buffer using distance in meters (converted to degrees)

    Args:
        input_file (str): Input vector file
        output_file (str): Output vector file
        distance_meters (float): Buffer distance in meters
        reference_latitude (float): Reference latitude for conversion (default: equator)

    Example:
        # Buffer by 10 km
        create_buffer_with_meter_distance('rivers.geojson', 'buffered.geojson', 10000)
    """
    distance_degrees = meters_to_degrees(distance_meters, reference_latitude)

    print(f"Converting {distance_meters} meters to ~{distance_degrees:.6f} degrees")
    print(f"  (at latitude {reference_latitude}°)")

    create_simple_buffer_from_file(input_file, output_file, distance_degrees)


