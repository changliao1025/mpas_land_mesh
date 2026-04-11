"""
Vector file operations and utilities

Minimal implementations of vector processing functions
extracted from pyearth.toolbox and hexwatershed_utility
"""

import os
import logging
import time
from typing import Union, Optional, Tuple
import numpy as np
from osgeo import ogr, osr, gdal

from functools import lru_cache

from mpas_land_mesh.utilities.system import get_extension_from_path
from mpas_land_mesh.utilities.geometry import (
    calculate_polygon_area,
)
from mpas_land_mesh.utilities.constants import IDL_threshold


gdal.UseExceptions()

SUPPORTED_VECTOR_FORMATS = {
    ".geojson": "GeoJSON",
    ".json": "GeoJSON",
    ".shp": "ESRI Shapefile",
    ".gpkg": "GPKG",
    ".kml": "KML",
    ".gml": "GML",
    ".sqlite": "SQLite",
    ".parquet": "Parquet",
    ".geoparquet": "Parquet",
    ".csv": "CSV",
    ".tab": "MapInfo File",
    ".mif": "MapInfo File",
    ".dxf": "DXF",
    ".gdb": "OpenFileGDB",
    ".fgb": "FlatGeobuf",
}


__all__ = [
    "SUPPORTED_VECTOR_FORMATS",
    "gdal_vector_format_support",
    "get_available_vector_formats",
    "print_supported_vector_formats",
    "get_vector_format_from_extension",
    "get_vector_driver_from_filename",
    "get_vector_driver_from_format",
    "get_vector_format_from_filename",
    "get_extension_from_vector_format",
    "check_parquet_support",
    "has_parquet_support",
]


@lru_cache(maxsize=1)
def get_available_vector_formats() -> dict[str, str]:
    """
    Get supported vector formats that are actually available in the current GDAL/OGR installation.

    Returns:
    dict[str, str]: A dictionary of available vector formats mapping file extensions to OGR driver names.
    """
    available_formats = {}

    for ext, format_name in SUPPORTED_VECTOR_FORMATS.items():
        # Special handling for Parquet/GeoParquet
        if ext in (".parquet", ".geoparquet"):
            if check_parquet_support() is not None:
                available_formats[ext] = "Parquet"  # Use the canonical name
        else:
            # Check if the driver is available
            driver = ogr.GetDriverByName(format_name)
            if driver is not None:
                available_formats[ext] = format_name

    return available_formats


def gdal_vector_format_support() -> dict[str, str]:
    """
    Return a dictionary of supported vector formats that are available in the current GDAL/OGR installation.

    This function filters the static SUPPORTED_VECTOR_FORMATS to only include formats
    where the corresponding OGR driver is actually available at runtime.

    Returns:
    dict[str, str]: A dictionary of available vector formats mapping file extensions to OGR driver names.
    """
    return get_available_vector_formats()


def print_supported_vector_formats() -> None:
    """
    Print all supported vector formats that are available in the current GDAL/OGR installation.
    """
    print("Supported vector formats (available in current GDAL/OGR installation):")
    available_formats = get_available_vector_formats()
    for ext, driver in available_formats.items():
        print(f"  {ext}: {driver}")

    # Also show unavailable formats from the supported list
    unavailable_formats = {}
    for ext, format_name in SUPPORTED_VECTOR_FORMATS.items():
        if ext not in available_formats:
            unavailable_formats[ext] = format_name

    if unavailable_formats:
        print("\nFormats in supported list but not available in current installation:")
        for ext, driver in unavailable_formats.items():
            print(f"  {ext}: {driver} (driver not found)")


def get_vector_format_from_extension(sExtension: str) -> str:
    """
    Determine the OGR format string from file extension.

    Only accepts extensions that are both in the supported list and have available drivers.

    Parameters:
    sExtension (str): Input file extension (e.g., '.shp').

    Returns:
    str: Format string for OGR driver.

    Raises:
    ValueError: If the file extension is not supported or driver is not available.
    """
    ext = sExtension.lower()

    available_formats = get_available_vector_formats()
    if ext not in available_formats:
        if ext in SUPPORTED_VECTOR_FORMATS:
            raise ValueError(
                f"File extension '{ext}' is supported but driver '{SUPPORTED_VECTOR_FORMATS[ext]}' is not available in current GDAL/OGR installation"
            )
        else:
            raise ValueError(f"Unsupported vector file extension: {ext}")

    return available_formats[ext]


def get_vector_format_from_filename(filename: str) -> str:
    """
    Determine the OGR format string from a filename.

    Only accepts extensions that are both in the supported list and have available drivers.

    Parameters:
    filename (str): Input filename with extension.

    Returns:
    str: Format string for OGR driver.

    Raises:
    ValueError: If the file extension is not supported or driver is not available.
    """
    sExtension = get_extension_from_path(filename).lower()
    if sExtension == "":
        raise ValueError(
            f"Filename '{filename}' does not have an extension to determine vector format."
        )
    return get_vector_format_from_extension(sExtension)


def get_extension_from_vector_format(format_name: str) -> str:
    """
    Convert a GDAL vector format name to its corresponding file extension.

    Parameters:
    format_name (str): The GDAL/OGR driver name (e.g., 'GeoJSON', 'ESRI Shapefile', 'GPKG').

    Returns:
    str: The file extension (e.g., '.geojson', '.shp', '.gpkg').

    Raises:
    ValueError: If the format name is not found in the supported formats or if multiple
                extensions map to the same format (in which case the first found is returned).

    Examples:
    >>> get_extension_from_vector_format('GeoJSON')
    '.geojson'
    >>> get_extension_from_vector_format('ESRI Shapefile')
    '.shp'
    >>> get_extension_from_vector_format('GPKG')
    '.gpkg'
    """
    # Create a reverse mapping from format names to extensions
    format_to_extension = {}
    for ext, fmt in SUPPORTED_VECTOR_FORMATS.items():
        if fmt not in format_to_extension:
            format_to_extension[fmt] = ext

    # Look up the extension
    if format_name not in format_to_extension:
        raise ValueError(
            f"Format '{format_name}' not found in supported vector formats. "
            f"Available formats: {', '.join(sorted(set(SUPPORTED_VECTOR_FORMATS.values())))}"
        )

    return format_to_extension[format_name]


def get_vector_driver_from_format(file_format: str) :
    """
    Get the OGR driver based on the provided vector file format.

    Parameters:
    file_format (str): The vector file format (e.g., 'GeoJSON', 'ESRI Shapefile', 'GPKG').

    Returns:
    ogr.Driver: The OGR driver object corresponding to the file format.

    Raises:
    ValueError: If the driver is not found.
    """
    driver = ogr.GetDriverByName(file_format)
    if driver is None:
        raise ValueError(f"Driver for format '{file_format}' not found.")
    return driver


def get_vector_driver_from_filename(filename: str) :
    """
    Get the OGR driver based on file extension.

    Parameters:
    filename (str): Path to the input file.

    Returns:
    ogr.Driver: The OGR driver object corresponding to the file format.

    Raises:
    ValueError: If the file extension is not supported or driver is not found.
    """
    sExtension = get_extension_from_path(filename).lower()
    if sExtension == "":
        raise ValueError(
            f"Filename '{filename}' does not have an extension to determine vector format."
        )
    format_name = get_vector_format_from_extension(sExtension)
    return get_vector_driver_from_format(format_name)


@lru_cache(maxsize=1)
def check_parquet_support() -> Optional[str]:
    """Check whether an OGR driver supporting Parquet/Arrow is available.

    Returns
    -------
    Optional[str]
        The first matching OGR driver name (for example 'Parquet' or 'Arrow')
        if available, otherwise ``None``.

    Notes
    -----
    - Matches driver names that contain 'parquet' or 'arrow' (case-insensitive).
    - Result is cached for the life of the process since available drivers do
      not change at runtime.
    """
    for i in range(ogr.GetDriverCount()):
        drv = ogr.GetDriver(i)
        if drv is None:
            continue
        name = drv.GetName()
        if not name:
            continue
        lname = name.lower()
        if "parquet" in lname or "arrow" in lname:
            return name
    return None


def has_parquet_support() -> bool:
    """Return True if a Parquet/Arrow OGR driver is available."""
    return check_parquet_support() is not None


def get_field_and_value(sFilename_in):
    """
    Extract field names and values from the first feature of a vector file

    Args:
        sFilename_in (str): Path to input vector file

    Returns:
        tuple: (field_names, field_values)
    """
    pDriver = ogr.GetDriverByName('GeoJSON')
    pDataset = pDriver.Open(sFilename_in, 0)

    if pDataset is None:
        raise ValueError(f"Could not open {sFilename_in}")

    pLayer = pDataset.GetLayer(0)
    pFeature = pLayer.GetNextFeature()

    if pFeature is None:
        return [], []

    aField = []
    aValue = []

    pLayerDefn = pLayer.GetLayerDefn()
    nField = pLayerDefn.GetFieldCount()

    for i in range(nField):
        pFieldDefn = pLayerDefn.GetFieldDefn(i)
        sFieldName = pFieldDefn.GetName()
        aField.append(sFieldName)
        aValue.append(pFeature.GetField(sFieldName))

    pDataset = None

    return aField, aValue


def add_field_to_vector_file(sFilename_in, aField, aValue):
    """
    Add fields and values to all features in a vector file

    Args:
        sFilename_in (str): Path to vector file
        aField (list): List of field names
        aValue (list): List of field values
    """
    pDataset = ogr.Open(sFilename_in, 1)  # 1 = write mode

    if pDataset is None:
        raise ValueError(f"Could not open {sFilename_in} for writing")

    pLayer = pDataset.GetLayer(0)
    pLayerDefn = pLayer.GetLayerDefn()

    # Add fields if they don't exist
    for i, sFieldName in enumerate(aField):
        iFieldIndex = pLayerDefn.GetFieldIndex(sFieldName)

        if iFieldIndex == -1:
            # Field doesn't exist, create it
            pFieldDefn = ogr.FieldDefn(sFieldName, ogr.OFTString)
            pLayer.CreateField(pFieldDefn)

    # Update all features with the values
    pLayer.ResetReading()
    for pFeature in pLayer:
        for i, sFieldName in enumerate(aField):
            pFeature.SetField(sFieldName, str(aValue[i]))
        pLayer.SetFeature(pFeature)

    pDataset = None


def merge_features(sFilename_in, sFilename_out, iFlag_force=False):
    """
    Merge all features in a vector file into a single feature

    Args:
        sFilename_in (str): Input vector file
        sFilename_out (str): Output vector file
        iFlag_force (bool): Overwrite existing output file
    """
    if os.path.exists(sFilename_out) and not iFlag_force:
        print(f"Output file {sFilename_out} already exists")
        return

    pDataset_in = ogr.Open(sFilename_in, 0)
    if pDataset_in is None:
        raise ValueError(f"Could not open {sFilename_in}")

    pLayer_in = pDataset_in.GetLayer(0)
    pSpatialRef = pLayer_in.GetSpatialRef()

    # Create output
    pDriver = ogr.GetDriverByName('GeoJSON')

    if os.path.exists(sFilename_out):
        pDriver.DeleteDataSource(sFilename_out)

    pDataset_out = pDriver.CreateDataSource(sFilename_out)
    pLayer_out = pDataset_out.CreateLayer('merged', pSpatialRef, ogr.wkbMultiPolygon)

    # Collect all geometries
    aGeometry = []
    for pFeature in pLayer_in:
        pGeometry = pFeature.GetGeometryRef()
        if pGeometry is not None:
            aGeometry.append(pGeometry.Clone())

    # Create union of all geometries
    if len(aGeometry) > 0:
        pGeometry_merged = aGeometry[0]
        for i in range(1, len(aGeometry)):
            pGeometry_merged = pGeometry_merged.Union(aGeometry[i])

        # Create output feature
        pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
        pFeature_out.SetGeometry(pGeometry_merged)
        pLayer_out.CreateFeature(pFeature_out)
        pFeature_out = None

    pDataset_in = None
    pDataset_out = None


def write_wkt_to_vector_file(wkt, sFilename_out, iEPSG_in=4326):
    """
    Write WKT geometry to a vector file

    Args:
        wkt (str): Well-Known Text geometry
        sFilename_out (str): Output file path
        iEPSG_in (int): EPSG code for spatial reference
    """
    pDriver = ogr.GetDriverByName('GeoJSON')

    if os.path.exists(sFilename_out):
        pDriver.DeleteDataSource(sFilename_out)

    pSpatialRef = osr.SpatialReference()
    pSpatialRef.ImportFromEPSG(iEPSG_in)

    pDataset = pDriver.CreateDataSource(sFilename_out)
    pLayer = pDataset.CreateLayer('boundary', pSpatialRef, ogr.wkbPolygon)

    pGeometry = ogr.CreateGeometryFromWkt(wkt)

    pFeature = ogr.Feature(pLayer.GetLayerDefn())
    pFeature.SetGeometry(pGeometry)
    pLayer.CreateFeature(pFeature)

    pFeature = None
    pDataset = None


def remove_small_polygon(
    sFilename_vector_in: str,
    sFilename_vector_out: str,
    dThreshold_in: Union[float, int], # in square meters
    iFlag_algorithm: int = 2,
    verbose: bool = True,
    progress_interval: int = 1000,
) -> None:
    """
    Remove small polygons from a vector file based on area threshold.

    This function filters polygon geometries from an input vector file, keeping only
    those with areas greater than the specified threshold. The area calculation uses
    geodesic methods for accurate results on large regions. Both POLYGON and
    MULTIPOLYGON geometries are supported, with holes (inner rings) preserved.

    Parameters
    ----------
    sFilename_vector_in : str
        Path to the input vector file containing polygon geometries.
    sFilename_vector_out : str
        Path where the filtered output vector file will be created.
    dThreshold_in : float or int
        Minimum area threshold in square meters. Polygons with areas
        less than or equal to this value will be removed.
    iFlag_algorithm : int, optional
        Algorithm flag for area calculation (default is 2 for geodesic).
        - 1: Planar area calculation
        - 2: Geodesic area calculation (recommended for large areas)
    verbose : bool, optional
        If True, print progress information and statistics (default is True).
    progress_interval : int, optional
        Number of features to process before printing progress updates
        (default is 1000).

    Returns
    -------
    None
        The function creates a new vector file at `sFilename_vector_out`
        containing only polygons that meet the area threshold.

    Raises
    ------
    FileNotFoundError
        If the input file does not exist.
    ValueError
        If the threshold value is invalid or if vector format is unsupported.
    RuntimeError
        If GDAL operations fail (file opening, creation, or processing).

    Notes
    -----
    - The function assumes input geometries are in WGS84 (EPSG:4326) coordinate system
      for accurate geodesic area calculations.
    - Output file will contain additional fields: 'id' (sequential integer) and
      'area' (calculated area in square kilometers).
    - All original attribute fields are preserved except for 'id' and 'area' which
      are overwritten.
    - Degenerate polygons (less than 3 vertices) are automatically skipped.
    - For MULTIPOLYGON geometries, each individual polygon part is evaluated
      separately against the threshold.

    Examples
    --------
    Remove polygons smaller than 1 square kilometer:

    >>> remove_small_polygon(
    ...     'input_polygons.shp',
    ...     'filtered_polygons.shp',
    ...     1.0
    ... )

    Process with custom progress reporting:

    >>> remove_small_polygon(
    ...     'world_countries.geojson',
    ...     'large_countries.geojson',
    ...     1000.0,
    ...     verbose=True,
    ...     progress_interval=500
    ... )
    """
    # Input validation
    if not isinstance(sFilename_vector_in, str) or not sFilename_vector_in.strip():
        raise ValueError("Input filename must be a non-empty string")

    if not isinstance(sFilename_vector_out, str) or not sFilename_vector_out.strip():
        raise ValueError("Output filename must be a non-empty string")

    if not os.path.exists(sFilename_vector_in):
        raise FileNotFoundError(f"Input file does not exist: {sFilename_vector_in}")

    try:
        dThreshold = float(dThreshold_in)
    except (TypeError, ValueError) as e:
        raise ValueError(
            f"Threshold must be a numeric value, got {type(dThreshold_in)}: {e}"
        )

    if dThreshold <= 0:
        raise ValueError(f"Threshold must be positive, got {dThreshold}")

    if iFlag_algorithm not in [1, 2]:
        raise ValueError(f"Algorithm flag must be 1 or 2, got {iFlag_algorithm}")

    # Remove existing output file if it exists
    if os.path.exists(sFilename_vector_out):
        if verbose:
            logging.info(f"Removing existing output file: {sFilename_vector_out}")
        try:
            os.remove(sFilename_vector_out)
        except OSError as e:
            raise RuntimeError(f"Could not remove existing output file: {e}")

    # Determine formats and drivers
    try:
        sFormat_in = get_vector_format_from_filename(sFilename_vector_in)
        sFormat_out = get_vector_format_from_filename(sFilename_vector_out)
        pDriver_out = get_vector_driver_from_format(sFormat_out)
    except Exception as e:
        raise RuntimeError(f"Failed to determine vector format: {e}")

    if pDriver_out is None:
        available_formats = print_supported_vector_formats()
        raise ValueError(
            f"Output format '{sFormat_out}' is not supported. "
            f"Available formats: {available_formats}"
        )

    if verbose:
        logging.info(f"Input format: {sFormat_in}")
        logging.info(f"Output format: {sFormat_out}")
        logging.info(f"Area threshold: {dThreshold} m²")
        logging.info(
            f"Area calculation algorithm: {'Geodesic' if iFlag_algorithm == 2 else 'Planar'}"
        )

    # Set up spatial reference (WGS84 for geodesic calculations)
    try:
        pSrs = osr.SpatialReference()
        pSrs.ImportFromEPSG(4326)
    except Exception as e:
        raise RuntimeError(f"Failed to create spatial reference: {e}")

    # Open input datasource
    try:
        pDataSource_in = ogr.Open(sFilename_vector_in, 0)  # Read-only
        if pDataSource_in is None:
            raise RuntimeError(f"Could not open input file: {sFilename_vector_in}")
    except Exception as e:
        raise RuntimeError(f"GDAL error opening input file: {e}")

    try:
        pLayer_in = pDataSource_in.GetLayer()
        if pLayer_in is None:
            raise RuntimeError("Could not access input layer")

        # Get feature count for progress tracking
        nTotal_features = pLayer_in.GetFeatureCount()
        if nTotal_features < 0:
            logging.warning(
                "Could not determine total feature count, progress reporting may be inaccurate"
            )

        if verbose:
            logging.info(f"Processing {nTotal_features} features...")

        # Create output datasource
        try:
            pDataSource_out = pDriver_out.CreateDataSource(sFilename_vector_out)
            if pDataSource_out is None:
                raise RuntimeError(
                    f"Could not create output file: {sFilename_vector_out}"
                )
        except Exception as e:
            raise RuntimeError(f"GDAL error creating output file: {e}")

        try:
            # Create output layer
            pLayer_out = pDataSource_out.CreateLayer(
                "filtered_polygons", pSrs, ogr.wkbPolygon
            )
            if pLayer_out is None:
                raise RuntimeError("Could not create output layer")

            # Add standard fields
            field_id = ogr.FieldDefn("id", ogr.OFTInteger)
            field_area = ogr.FieldDefn("area", ogr.OFTReal)
            pLayer_out.CreateField(field_id)
            pLayer_out.CreateField(field_area)

            # Copy existing fields from input (excluding id and area if they exist)
            pLayerDefn_in = pLayer_in.GetLayerDefn()
            nFieldCount = pLayerDefn_in.GetFieldCount()

            for i in range(nFieldCount):
                pFieldDefn = pLayerDefn_in.GetFieldDefn(i)
                field_name = pFieldDefn.GetName().lower()
                if field_name not in ["id", "area"]:
                    try:
                        pLayer_out.CreateField(pFieldDefn)
                    except Exception as e:
                        logging.warning(f"Could not create field '{field_name}': {e}")

            pLayerDefn_out = pLayer_out.GetLayerDefn()

            # Process features
            start_time = time.time()
            lID = 1
            nProcessed = 0
            nKept = 0
            nSkipped = 0

            pLayer_in.ResetReading()  # Ensure we're at the beginning

            for pFeature_in in pLayer_in:
                nProcessed += 1

                if (
                    verbose
                    and progress_interval > 0
                    and nProcessed % progress_interval == 0
                ):
                    elapsed = time.time() - start_time
                    rate = nProcessed / elapsed if elapsed > 0 else 0
                    logging.info(
                        f"Processed {nProcessed}/{nTotal_features} features "
                        f"({rate:.1f} features/sec), kept {nKept}"
                    )

                pGeometry = pFeature_in.GetGeometryRef()
                if pGeometry is None:
                    nSkipped += 1
                    continue

                geometry_type = pGeometry.GetGeometryName()

                # Process POLYGON geometries
                if geometry_type == "POLYGON":
                    kept = _process_single_polygon(
                        pGeometry,
                        dThreshold,
                        iFlag_algorithm,
                        pSrs,
                        pLayerDefn_out,
                        pFeature_in,
                        pLayerDefn_in,
                        lID,
                        pLayer_out,
                    )
                    if kept:
                        lID += 1
                        nKept += 1

                # Process MULTIPOLYGON geometries
                elif geometry_type == "MULTIPOLYGON":
                    polygons_kept = _process_multipolygon(
                        pGeometry,
                        dThreshold,
                        iFlag_algorithm,
                        pSrs,
                        pLayerDefn_out,
                        pFeature_in,
                        pLayerDefn_in,
                        lID,
                        pLayer_out,
                    )
                    lID += polygons_kept
                    nKept += polygons_kept

                else:
                    nSkipped += 1
                    if verbose and nSkipped <= 5:  # Log first few non-polygon features
                        logging.debug(f"Skipping non-polygon geometry: {geometry_type}")

            # Final statistics
            total_time = time.time() - start_time

            if verbose:
                logging.info("Processing complete!")
                logging.info(f"Total features processed: {nProcessed}")
                logging.info(f"Features kept (area > {dThreshold} km²): {nKept}")
                logging.info(f"Features removed: {nProcessed - nKept - nSkipped}")
                logging.info(f"Features skipped (non-polygon): {nSkipped}")
                logging.info(f"Processing time: {total_time:.2f} seconds")
                logging.info(
                    f"Average processing rate: {nProcessed/total_time:.1f} features/sec"
                )
                logging.info(f"Output saved to: {sFilename_vector_out}")

        finally:
            # Clean up output datasource
            pDataSource_out = None

    finally:
        # Clean up input datasource
        pDataSource_in = None


def _process_single_polygon(
    pGeometry: ogr.Geometry,
    dThreshold: float,
    iFlag_algorithm: int,
    pSrs: osr.SpatialReference,
    pLayerDefn_out: ogr.FeatureDefn,
    pFeature_in: ogr.Feature,
    pLayerDefn_in: ogr.FeatureDefn,
    lID: int,
    pLayer_out: ogr.Layer,
) -> bool:
    """
    Process a single POLYGON geometry.

    Returns True if the polygon was kept, False otherwise.
    """
    try:
        # Extract outer ring coordinates
        pOuterRing = pGeometry.GetGeometryRef(0)
        if pOuterRing is None or pOuterRing.GetPointCount() < 3:
            return False

        aCoords_outer = np.array(
            [
                [pOuterRing.GetPoint(i)[0], pOuterRing.GetPoint(i)[1]]
                for i in range(pOuterRing.GetPointCount())
            ]
        )

        # Check if polygon is close to the antimeridian line (180° or -180°)
        # If it is, keep it regardless of area since it might be split and actual area is larger
        lon_coords = aCoords_outer[:, 0]

        is_near_antimeridian = np.any( np.abs(np.abs(lon_coords[:-1]) - 180.0) < IDL_threshold )

        if is_near_antimeridian:
            # Keep polygon near antimeridian without area check
            pass
        else:
            # Calculate area
            dArea = calculate_polygon_area(
                aCoords_outer[:, 0], aCoords_outer[:, 1], iFlag_algorithm
            )
            dAreakm = dArea / 1e6  # Convert to square kilometers

            if dArea <= dThreshold:
                return False

        # Recalculate area for output field (needed for both cases)
        dArea = calculate_polygon_area(
            aCoords_outer[:, 0], aCoords_outer[:, 1], iFlag_algorithm
        )
        dAreakm = dArea / 1e6  # Convert to square kilometers

        # Create output polygon with all rings
        pGeometry_out = ogr.Geometry(ogr.wkbPolygon)

        # Add outer ring
        pRing_outer = ogr.Geometry(ogr.wkbLinearRing)
        for coord in aCoords_outer:
            pRing_outer.AddPoint(coord[0], coord[1])
        pRing_outer.CloseRings()
        pGeometry_out.AddGeometry(pRing_outer)

        # Add inner rings (holes)
        for iRing in range(1, pGeometry.GetGeometryCount()):
            pInnerRing_src = pGeometry.GetGeometryRef(iRing)
            if pInnerRing_src is None:
                continue

            pInnerRing_new = ogr.Geometry(ogr.wkbLinearRing)
            for iPoint in range(pInnerRing_src.GetPointCount()):
                x, y, z = pInnerRing_src.GetPoint(iPoint)
                pInnerRing_new.AddPoint(x, y)
            pInnerRing_new.CloseRings()
            pGeometry_out.AddGeometry(pInnerRing_new)

        pGeometry_out.AssignSpatialReference(pSrs)

        # Create and populate output feature
        pFeature_out = ogr.Feature(pLayerDefn_out)
        pFeature_out.SetGeometry(pGeometry_out)
        pFeature_out.SetField("id", lID)
        pFeature_out.SetField("area", dAreakm)

        # Copy other fields
        for i in range(pLayerDefn_in.GetFieldCount()):
            field_name = pLayerDefn_in.GetFieldDefn(i).GetName()
            if field_name.lower() not in ["id", "area"]:
                try:
                    pFeature_out.SetField(field_name, pFeature_in.GetField(field_name))
                except Exception:
                    # Skip fields that can't be copied
                    pass

        pLayer_out.CreateFeature(pFeature_out)
        pFeature_out = None

        return True

    except Exception as e:
        logging.warning(f"Error processing polygon: {e}")
        return False


def _process_multipolygon(
    pGeometry: ogr.Geometry,
    dThreshold: float,
    iFlag_algorithm: int,
    pSrs: osr.SpatialReference,
    pLayerDefn_out: ogr.FeatureDefn,
    pFeature_in: ogr.Feature,
    pLayerDefn_in: ogr.FeatureDefn,
    lID_start: int,
    pLayer_out: ogr.Layer,
) -> int:
    """
    Process a MULTIPOLYGON geometry, creating separate features for each polygon part.

    Returns the number of polygons that were kept.
    """
    polygons_kept = 0

    try:
        for iPoly in range(pGeometry.GetGeometryCount()):
            pPolygon = pGeometry.GetGeometryRef(iPoly)
            if pPolygon is None or pPolygon.GetGeometryName() != "POLYGON":
                continue

            kept = _process_single_polygon(
                pPolygon,
                dThreshold,
                iFlag_algorithm,
                pSrs,
                pLayerDefn_out,
                pFeature_in,
                pLayerDefn_in,
                lID_start + polygons_kept,
                pLayer_out,
            )

            if kept:
                polygons_kept += 1

    except Exception as e:
        logging.warning(f"Error processing multipolygon: {e}")

    return polygons_kept
