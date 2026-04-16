"""
Coastline processing utilities

Simplified implementations for coastline operations
"""

import os
import logging
from osgeo import gdal, ogr, osr
from rtree.index import Index as RTreeindex
import logging
from pathlib import Path
from typing import Optional, Literal

import cartopy.feature as cfeature

from mpas_land_mesh.utilities.vector import (
    get_vector_driver_from_filename,
    get_vector_format_from_filename,
    remove_small_polygon,
    merge_features
)
from mpas_land_mesh.utilities.raster import convert_vector_to_global_raster, create_raster_buffer_zone, fix_raster_antimeridian_issue
from mpas_land_mesh.utilities.geometry import calculate_distance_based_on_longitude_latitude
from mpas_land_mesh.utilities.gcsbuffer import create_wkt_buffer_distance

gdal.UseExceptions()

# Set up logger for this module
logger = logging.getLogger(__name__)


# https://scitools.org.uk/cartopy/docs/latest/reference/generated/cartopy.feature.NaturalEarthFeature.html


def create_land_ocean_vector_mask_naturalearth(
    sFilename_out: str,
    sResolution_coastal: Literal["10m", "50m", "110m"] = "10m",
    sWorkspace_out: Optional[str] = None,
    sFormat: Optional[str] = None,
    iFlag_exclude_antarctica: bool = True,
    dLatitude_threshold_antarctica: float = -60.0,
) -> int:
    """
    Create a land-ocean vector mask from Natural Earth data.

    Downloads and processes Natural Earth physical land boundaries to create a
    vector mask distinguishing land from ocean. The mask is saved as individual
    polygon features that can be used for spatial analysis, clipping, or masking
    operations in GIS and climate modeling applications.

    Parameters
    ----------
    sFilename_out : str
        Output file path for the land-ocean mask.

        File extension determines output format if sFormat not specified:
        - '.geojson' → GeoJSON
        - '.shp' → ESRI Shapefile
        - '.gpkg' → GeoPackage
        - '.kml' → KML

        Examples:
        - '/path/to/land_mask.geojson'
        - '/path/to/land_mask.shp'
        - 'output/ocean_mask.gpkg'

        File will be overwritten if it exists.

    sResolution_coastal : {'10m', '50m', '110m'}, optional
        Resolution of Natural Earth coastal data, default '10m'.

        - '10m': Highest detail (1:10 million scale)
          * Most accurate coastlines
          * ~13 MB download on first use
          * Recommended for regional studies

        - '50m': Medium detail (1:50 million scale)
          * Good balance of detail and file size
          * ~3 MB download on first use
          * Suitable for continental/global studies

        - '110m': Lowest detail (1:110 million scale)
          * Simplified coastlines
          * ~1 MB download on first use
          * Suitable for global overviews

        Data is cached by Cartopy after first download.

    sWorkspace_out : str, optional
        Directory path for saving individual land geometry parts, default None.

        If specified:
        - Directory will be created if it doesn't exist
        - Each polygon/multipolygon part saved as separate GeoJSON
        - Files named: 'land_geometry_0.geojson', 'land_geometry_1.geojson', etc.
        - Useful for debugging or analyzing individual land masses

        If None:
        - Individual parts not saved (only combined output file created)
        - Faster processing, less disk usage

        Examples:
        - '/path/to/workspace/individual_parts'
        - 'output/geometries'

    sFormat : str, optional
        Output vector format name, default None (auto-detect from extension).

        Supported formats (driver name):
        - 'GeoJSON': JSON-based, human-readable, widely supported
        - 'ESRI Shapefile': Traditional GIS format
        - 'GPKG': GeoPackage, modern SQLite-based format
        - 'KML': Google Earth format
        - Other GDAL/OGR supported formats

        If None, format is auto-detected from sFilename_out extension.

        Examples:
        - 'GeoJSON'
        - 'ESRI Shapefile'
        - 'GPKG'

    iFlag_exclude_antarctica : bool, optional
        Whether to exclude Antarctica and sub-Antarctic islands, default True.

        - True: Exclude land masses with southern extent < dLatitude_threshold_antarctica
          * Reduces output file size and complexity
          * Useful for non-polar studies
          * Default excludes Antarctica (typically < -60°S)

        - False: Include all land masses including Antarctica
          * Complete global land mask
          * Necessary for polar studies

    dLatitude_threshold_antarctica : float, optional
        Latitude threshold for Antarctica exclusion (degrees), default -60.0.

        Land masses with minimum latitude (southern extent) below this value
        are excluded if iFlag_exclude_antarctica is True.

        Range: -90.0 to 0.0 (southern hemisphere)

        Common values:
        - -60.0: Standard Antarctic Circle approximation
        - -66.5: Actual Antarctic Circle
        - -50.0: Include sub-Antarctic islands

        Only used if iFlag_exclude_antarctica is True.

    Returns
    -------
    int
        Number of individual polygon features created in the output file.

        - Multipolygons are decomposed into individual polygons
        - Each polygon is a separate feature with ID and type attributes
        - Count excludes Antarctica if iFlag_exclude_antarctica is True

        Typical values:
        - 10m resolution: ~1500-2000 features
        - 50m resolution: ~500-1000 features
        - 110m resolution: ~200-500 features

    Raises
    ------
    TypeError
        If sFilename_out is not a string.
    ValueError
        - If sResolution_coastal not in ['10m', '50m', '110m']
        - If sFilename_out is empty string
        - If output directory is not writable
        - If sFormat is not supported by GDAL/OGR
    RuntimeError
        - If Natural Earth data download fails
        - If output file cannot be created
        - If GDAL driver not available
    OSError
        If workspace directory cannot be created.

    Warns
    -----
    UserWarning
        - If sFormat is specified but conflicts with file extension
        - If Natural Earth data is being downloaded (first time only)
        - If no geometries remain after Antarctica filtering

    Notes
    -----
    **Data Source:**

    - Uses Natural Earth Physical Land 1:10m, 1:50m, or 1:110m dataset
    - Downloaded automatically via Cartopy on first use
    - Data cached in: ~/.local/share/cartopy/ (Linux/Mac) or
      %LOCALAPPDATA%\\cartopy (Windows)
    - No API key or registration required

    **Processing Details:**

    1. Download Natural Earth land geometries via Cartopy
    2. Iterate through each land geometry (polygon or multipolygon)
    3. Check geometry type (POLYGON or MULTIPOLYGON)
    4. For each geometry:
       - Extract bounding box envelope
       - Apply Antarctica filter if enabled
       - Decompose multipolygons into individual polygons
       - Add each polygon as separate feature with attributes
    5. Optionally save individual parts as GeoJSON files
    6. Return total feature count

    **Output Attributes:**

    Each feature in the output file has these fields:
    - 'id': Integer, sequential feature ID (0, 1, 2, ...)
    - 'part_type': String, either 'polygon' or 'multipolygon_part'
      * 'polygon': Original geometry was a single polygon
      * 'multipolygon_part': Extracted from a multipolygon geometry

    **File Format Notes:**

    - GeoJSON: Single file, coordinates in WGS84
    - Shapefile: Multiple files (.shp, .shx, .dbf, .prj, etc.)
    - GeoPackage: Single SQLite database file

    **Performance:**

    - First run: Slower (downloads Natural Earth data)
    - Subsequent runs: Fast (uses cached data)
    - 10m resolution: ~30-60 seconds for global mask
    - 50m resolution: ~10-20 seconds
    - 110m resolution: ~5-10 seconds
    - Saving individual parts adds ~20-50% processing time

    Examples
    --------
    **Example 1: Basic usage with default settings (10m, GeoJSON)**

    >>> from pyearth.toolbox.data.ocean.define_land_ocean_mask import \\
    ...     create_land_ocean_vector_mask_naturalearth
    >>>
    >>> # Create high-resolution land mask (auto-detect GeoJSON format)
    >>> n_features = create_land_ocean_vector_mask_naturalearth(
    ...     sFilename_out='output/land_mask.geojson')
    >>> print(f"Created mask with {n_features} land features")
    Created mask with 1847 land features

    **Example 2: Medium resolution shapefile**

    >>> # Create medium-resolution shapefile
    >>> n_features = create_land_ocean_vector_mask_naturalearth(
    ...     sFilename_out='output/land_mask.shp',
    ...     sResolution_coastal='50m')
    >>> print(f"Created shapefile with {n_features} features")
    Created shapefile with 743 features

    **Example 3: Low resolution with GeoPackage format**

    >>> # Fast global overview in GeoPackage format
    >>> n_features = create_land_ocean_vector_mask_naturalearth(
    ...     sFilename_out='output/land_mask.gpkg',
    ...     sResolution_coastal='110m')
    >>> print(f"Created GeoPackage with {n_features} features")
    Created GeoPackage with 312 features

    **Example 4: Include Antarctica**

    >>> # Create mask including Antarctica for polar studies
    >>> n_features = create_land_ocean_vector_mask_naturalearth(
    ...     sFilename_out='output/land_mask_polar.geojson',
    ...     iFlag_exclude_antarctica=False)
    >>> print(f"Global mask (including Antarctica): {n_features} features")
    Global mask (including Antarctica): 1923 features

    **Example 5: Custom Antarctica threshold**

    >>> # Include sub-Antarctic islands (exclude only < -50°S)
    >>> n_features = create_land_ocean_vector_mask_naturalearth(
    ...     sFilename_out='output/land_mask_subantarctic.geojson',
    ...     iFlag_exclude_antarctica=True,
    ...     dLatitude_threshold_antarctica=-50.0)
    >>> print(f"Mask with sub-Antarctic islands: {n_features} features")
    Mask with sub-Antarctic islands: 1889 features

    **Example 6: Save individual land mass geometries**

    >>> # Create mask and save each land mass separately
    >>> n_features = create_land_ocean_vector_mask_naturalearth(
    ...     sFilename_out='output/land_mask.geojson',
    ...     sWorkspace_out='output/individual_parts',
    ...     sResolution_coastal='50m')
    >>> print(f"Created {n_features} features, saved individually to workspace")
    Created 743 features, saved individually to workspace
    >>> # Individual files: land_geometry_0.geojson, land_geometry_1.geojson, ...

    **Example 7: Explicit format specification**

    >>> # Force KML format even with .kml extension
    >>> n_features = create_land_ocean_vector_mask_naturalearth(
    ...     sFilename_out='output/land_mask.kml',
    ...     sFormat='KML',
    ...     sResolution_coastal='110m')
    >>> print(f"Created KML file with {n_features} features")
    Created KML file with 312 features

    **Example 8: Regional study (exclude Antarctica, high detail)**

    >>> # High-detail mask for non-polar regions
    >>> n_features = create_land_ocean_vector_mask_naturalearth(
    ...     sFilename_out='output/land_mask_regional.gpkg',
    ...     sResolution_coastal='10m',
    ...     iFlag_exclude_antarctica=True)
    >>> print(f"Regional high-detail mask: {n_features} features")
    Regional high-detail mask: 1847 features

    **Example 9: Debugging with individual parts (low resolution)**

    >>> # Create low-res mask and inspect individual land masses
    >>> import os
    >>> n_features = create_land_ocean_vector_mask_naturalearth(
    ...     sFilename_out='output/debug_mask.geojson',
    ...     sWorkspace_out='output/debug_parts',
    ...     sResolution_coastal='110m')
    >>> # Check individual files created
    >>> part_files = os.listdir('output/debug_parts')
    >>> print(f"Created {len(part_files)} individual geometry files")
    Created 312 individual geometry files

    **Example 10: Complete workflow for climate model**

    >>> # Create land mask for climate simulation
    >>> # High resolution, exclude Antarctica, save as GeoPackage
    >>> output_file = 'climate_model/land_ocean_mask.gpkg'
    >>> n_features = create_land_ocean_vector_mask_naturalearth(
    ...     sFilename_out=output_file,
    ...     sResolution_coastal='10m',
    ...     iFlag_exclude_antarctica=True,
    ...     dLatitude_threshold_antarctica=-60.0)
    >>>
    >>> # Verify output
    >>> from osgeo import ogr
    >>> ds = ogr.Open(output_file)
    >>> layer = ds.GetLayer(0)
    >>> print(f"Layer name: {layer.GetName()}")
    >>> print(f"Feature count: {layer.GetFeatureCount()}")
    >>> print(f"Extent: {layer.GetExtent()}")
    Layer name: land_ocean_mask
    Feature count: 1847
    Extent: (-180.0, 180.0, -60.0, 83.6)

    See Also
    --------
    cartopy.feature.NaturalEarthFeature : Access Natural Earth datasets
    osgeo.ogr.DataSource : GDAL vector data container
    pyearth.gis.gdal.gdal_vector_format_support : Vector format detection and driver utilities

    References
    ----------
    .. [1] Natural Earth. Free vector and raster map data at 1:10m, 1:50m, and
           1:110 million scales. https://www.naturalearthdata.com/
    .. [2] Cartopy documentation. Natural Earth feature interface.
           https://scitools.org.uk/cartopy/docs/latest/
    .. [3] GDAL/OGR Vector Data Model. https://gdal.org/user/vector_data_model.html

    """
    # ========================================================================
    # Input validation
    # ========================================================================

    # Validate output filename
    if not isinstance(sFilename_out, str):
        raise TypeError(f"sFilename_out must be a string, got {type(sFilename_out)}")

    if not sFilename_out or sFilename_out.isspace():
        raise ValueError("sFilename_out cannot be empty")

    # Validate resolution
    valid_resolutions = ["10m", "50m", "110m"]
    if sResolution_coastal not in valid_resolutions:
        raise ValueError(
            f"sResolution_coastal must be one of {valid_resolutions}, "
            f"got '{sResolution_coastal}'"
        )

    # Validate Antarctica threshold
    if not isinstance(dLatitude_threshold_antarctica, (int, float)):
        raise TypeError(
            f"dLatitude_threshold_antarctica must be numeric, "
            f"got {type(dLatitude_threshold_antarctica)}"
        )

    if not (-90.0 <= dLatitude_threshold_antarctica <= 0.0):
        raise ValueError(
            f"dLatitude_threshold_antarctica must be in range [-90, 0], "
            f"got {dLatitude_threshold_antarctica}"
        )

    logger.info(
        f"Creating land-ocean mask: resolution={sResolution_coastal}, "
        f"output={sFilename_out}, exclude_antarctica={iFlag_exclude_antarctica}"
    )

    # ========================================================================
    # Auto-detect format from output filename if not specified
    # ========================================================================

    if sFormat is None:
        try:
            sFormat = get_vector_format_from_filename(sFilename_out)
            logger.info(f"Auto-detected format from extension: {sFormat}")
        except Exception as e:
            raise ValueError(
                f"Cannot auto-detect format from '{sFilename_out}'. "
                f"Please specify sFormat parameter. Error: {e}"
            ) from e
    else:
        logger.info(f"Using specified format: {sFormat}")

    # ========================================================================
    # Get the appropriate GDAL driver
    # ========================================================================

    # Try to get driver from extension first, then fall back to driver by name
    try:
        pDriver_out = get_vector_driver_from_filename(sFilename_out)
        if pDriver_out is None:
            # Fallback to getting driver by name
            pDriver_out = ogr.GetDriverByName(sFormat)
            if pDriver_out is None:
                raise RuntimeError(
                    f"GDAL driver '{sFormat}' not available. "
                    f"Check GDAL installation and supported formats."
                )
    except Exception as e:
        # Fallback to getting driver by name
        logger.debug(
            f"get_vector_driver_from_extension failed: {e}, using GetDriverByName"
        )
        pDriver_out = ogr.GetDriverByName(sFormat)
        if pDriver_out is None:
            raise RuntimeError(
                f"GDAL driver '{sFormat}' not available. "
                f"Check GDAL installation and supported formats."
            )

    # Keep GeoJSON driver for intermediate files
    pDriver_geojson = ogr.GetDriverByName("GeoJSON")
    if pDriver_geojson is None:
        raise RuntimeError(
            "GeoJSON driver not available. Required for intermediate files."
        )

    # ========================================================================
    # Download Natural Earth data and extract geometries
    # ========================================================================

    logger.info(
        f"Downloading Natural Earth land data (resolution: {sResolution_coastal}). "
        f"This may take time on first run..."
    )

    try:
        # Create a land feature from Natural Earth
        land_feature = cfeature.NaturalEarthFeature(
            "physical", "land", sResolution_coastal
        )
        land_geometries = list(land_feature.geometries())
    except Exception as e:
        raise RuntimeError(
            f"Failed to download Natural Earth data at {sResolution_coastal} resolution. "
            f"Check internet connection and Cartopy configuration. Error: {e}"
        ) from e

    n_total_geometries = len(land_geometries)
    logger.info(f"Downloaded {n_total_geometries} land geometries from Natural Earth")

    # ========================================================================
    # Create output dataset
    # ========================================================================

    # Remove existing output file if it exists
    if os.path.exists(sFilename_out):
        logger.debug(f"Removing existing output file: {sFilename_out}")

        # For shapefiles, we need to delete all associated files
        if sFormat == "ESRI Shapefile":
            base_name = os.path.splitext(sFilename_out)[0]
            # Delete all shapefile components
            for ext in [".shp", ".shx", ".dbf", ".prj", ".cpg", ".qpj", ".sbn", ".sbx"]:
                file_to_delete = base_name + ext
                if os.path.exists(file_to_delete):
                    os.remove(file_to_delete)
                    logger.debug(f"Removed {file_to_delete}")
        else:
            os.remove(sFilename_out)

    # Create output directory if needed
    output_dir = os.path.dirname(sFilename_out)
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
            logger.debug(f"Created output directory: {output_dir}")
        except OSError as e:
            raise OSError(f"Cannot create output directory '{output_dir}': {e}") from e

    # Create output dataset
    try:
        pDataset_out = pDriver_out.CreateDataSource(sFilename_out)
        if pDataset_out is None:
            raise RuntimeError(f"Failed to create output file: {sFilename_out}")
    except Exception as e:
        raise RuntimeError(
            f"Cannot create output dataset '{sFilename_out}': {e}"
        ) from e

    # Create output layer with WGS84 spatial reference
    pSpatialRef = osr.SpatialReference()
    pSpatialRef.ImportFromEPSG(4326)  # WGS84

    pLayer_out = pDataset_out.CreateLayer(
        "land_ocean_mask", srs=pSpatialRef, geom_type=ogr.wkbPolygon
    )

    if pLayer_out is None:
        pDataset_out.Destroy()
        raise RuntimeError("Failed to create output layer")

    # Add fields to the output layer
    pField_id = ogr.FieldDefn("id", ogr.OFTInteger)
    pLayer_out.CreateField(pField_id)

    pField_part_type = ogr.FieldDefn("part_type", ogr.OFTString)
    pField_part_type.SetWidth(20)  # Set field width
    pLayer_out.CreateField(pField_part_type)

    logger.debug("Created output layer with 'id' and 'part_type' fields")

    # ========================================================================
    # Create workspace directory for individual parts if specified
    # ========================================================================

    if sWorkspace_out is not None:
        try:
            Path(sWorkspace_out).mkdir(parents=True, exist_ok=True)
            logger.info(f"Created workspace directory: {sWorkspace_out}")
        except OSError as e:
            logger.warning(
                f"Cannot create workspace directory '{sWorkspace_out}': {e}. "
                f"Individual parts will not be saved."
            )
            sWorkspace_out = None

    # ========================================================================
    # Process geometries
    # ========================================================================

    iCount = 0
    n_excluded_antarctica = 0
    n_polygons = 0
    n_multipolygons = 0

    logger.info("Processing land geometries...")

    for idx, land_geometry in enumerate(land_geometries):
        if idx % 100 == 0 and idx > 0:
            logger.debug(f"Processed {idx}/{n_total_geometries} geometries...")

        # Convert Shapely geometry to OGR geometry
        try:
            pGeometry_mesh = ogr.CreateGeometryFromWkb(land_geometry.wkb)
        except Exception as e:
            logger.warning(f"Failed to convert geometry {idx}: {e}")
            continue

        if pGeometry_mesh is None:
            logger.warning(f"Geometry {idx} is None, skipping")
            continue

        # Get envelope and geometry type
        envelope = pGeometry_mesh.GetEnvelope()  # (minX, maxX, minY, maxY)
        sGeometry_name = pGeometry_mesh.GetGeometryName()

        # ====================================================================
        # Handle POLYGON geometries
        # ====================================================================

        if sGeometry_name == "POLYGON":
            n_polygons += 1

            # Check Antarctica filter
            if (
                iFlag_exclude_antarctica
                and envelope[2] < dLatitude_threshold_antarctica
            ):
                n_excluded_antarctica += 1
                logger.debug(
                    f"Excluding polygon {idx} (min_lat={envelope[2]:.2f} < "
                    f"{dLatitude_threshold_antarctica})"
                )
                continue

            # Add to output layer
            pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
            pFeature_out.SetGeometry(pGeometry_mesh)
            pFeature_out.SetField("id", iCount)
            pFeature_out.SetField("part_type", "polygon")

            if pLayer_out.CreateFeature(pFeature_out) != 0:
                logger.warning(f"Failed to create feature for polygon {idx}")

            pFeature_out.Destroy()

            # Save individual part if workspace specified
            if sWorkspace_out is not None:
                sFilename_part = os.path.join(
                    sWorkspace_out, f"land_geometry_{iCount}.geojson"
                )

                try:
                    if os.path.exists(sFilename_part):
                        os.remove(sFilename_part)

                    pDataset_part = pDriver_geojson.CreateDataSource(sFilename_part)
                    pLayer_part = pDataset_part.CreateLayer(
                        "part", geom_type=ogr.wkbPolygon
                    )
                    pFeature_part = ogr.Feature(pLayer_part.GetLayerDefn())
                    pFeature_part.SetGeometry(pGeometry_mesh)
                    pLayer_part.CreateFeature(pFeature_part)
                    pFeature_part.Destroy()
                    pDataset_part.Destroy()
                except Exception as e:
                    logger.warning(f"Failed to save individual part {iCount}: {e}")

            iCount += 1

        # ====================================================================
        # Handle MULTIPOLYGON geometries
        # ====================================================================

        elif sGeometry_name == "MULTIPOLYGON":
            n_multipolygons += 1
            nPart1 = pGeometry_mesh.GetGeometryCount()

            logger.debug(f"Processing multipolygon {idx} with {nPart1} parts")

            for iPart in range(nPart1):
                pGeometry_part = pGeometry_mesh.GetGeometryRef(iPart)

                if pGeometry_part is None:
                    logger.warning(
                        f"Part {iPart} of multipolygon {idx} is None, skipping"
                    )
                    continue

                # Get envelope for this part
                envelope_part = pGeometry_part.GetEnvelope()

                # Check Antarctica filter
                if (
                    iFlag_exclude_antarctica
                    and envelope_part[2] < dLatitude_threshold_antarctica
                ):
                    n_excluded_antarctica += 1
                    logger.debug(
                        f"Excluding multipolygon part {iPart}/{nPart1} "
                        f"(min_lat={envelope_part[2]:.2f} < {dLatitude_threshold_antarctica})"
                    )
                    continue

                # Add to output layer
                pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
                pFeature_out.SetGeometry(pGeometry_part)
                pFeature_out.SetField("id", iCount)
                pFeature_out.SetField("part_type", "multipolygon_part")

                if pLayer_out.CreateFeature(pFeature_out) != 0:
                    logger.warning(
                        f"Failed to create feature for multipolygon part {iPart}"
                    )

                pFeature_out.Destroy()

                # Save individual part if workspace specified
                if sWorkspace_out is not None:
                    sFilename_part = os.path.join(
                        sWorkspace_out, f"land_geometry_{iCount}.geojson"
                    )

                    try:
                        if os.path.exists(sFilename_part):
                            os.remove(sFilename_part)

                        pDataset_part = pDriver_geojson.CreateDataSource(sFilename_part)
                        pLayer_part = pDataset_part.CreateLayer(
                            "part", geom_type=ogr.wkbPolygon
                        )
                        pFeature_part = ogr.Feature(pLayer_part.GetLayerDefn())
                        pFeature_part.SetGeometry(pGeometry_part)
                        pLayer_part.CreateFeature(pFeature_part)
                        pFeature_part.Destroy()
                        pDataset_part.Destroy()
                    except Exception as e:
                        logger.warning(f"Failed to save individual part {iCount}: {e}")

                iCount += 1

        else:
            logger.warning(
                f"Unknown geometry type '{sGeometry_name}' for geometry {idx}, skipping"
            )

    # ========================================================================
    # Clean up and finalize
    # ========================================================================

    # Close output dataset
    pDataset_out.Destroy()

    # ========================================================================
    # Log summary statistics
    # ========================================================================

    logger.info(
        f"Land-ocean mask creation complete:\n"
        f"  Output file: {sFilename_out}\n"
        f"  Total features created: {iCount}\n"
        f"  Original geometries: {n_total_geometries}\n"
        f"  Polygons processed: {n_polygons}\n"
        f"  Multipolygons processed: {n_multipolygons}\n"
        f"  Antarctica excluded: {n_excluded_antarctica}\n"
        f"  Resolution: {sResolution_coastal}\n"
        f"  Format: {sFormat}"
    )

    if sWorkspace_out is not None:
        logger.info(f"  Individual parts saved to: {sWorkspace_out}")

    if iCount == 0:
        logger.warning("No features were created! Check Antarctica filter settings.")

    return iCount


def create_land_ocean_mask_from_naturalearth(sWorkspace_coastline_output: str,
                                                                             dResolution_x_in: float, dResolution_y_in: float,
                                                                             dThreshold_area_island: float,
                                                                             dResolution_coastline_buffer: float,
                                                                             iRaster_buffer_pixel: int = 1) -> tuple:

    sFilename_naturalearth = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_naturalearth.geojson')
    create_land_ocean_vector_mask_naturalearth(sFilename_naturalearth)

    sFilename_vector_coastline = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island.geojson')
    remove_small_polygon(sFilename_naturalearth, sFilename_vector_coastline, dThreshold_area_island )

    sFilename_tif_wo_island = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island.tif')
    convert_vector_to_global_raster(sFilename_vector_coastline,
                                    sFilename_tif_wo_island,
                                    dResolution_x_in,
                                    dResolution_y_in,
                                    iFlag_boundary_only_in = 0,
                                    dFill_value_in = 2)


    if iRaster_buffer_pixel > 0:
        sFilename_tif_wo_island_buffered = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island_buffered.tif')
        create_raster_buffer_zone(sFilename_tif_wo_island, sFilename_tif_wo_island_buffered, 1, int(iRaster_buffer_pixel) )
        sFilename_tif_wo_island_buffered_fixed = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island_buffered_fixed.tif')
        #fix dateline issue
        fix_raster_antimeridian_issue(sFilename_tif_wo_island_buffered, sFilename_tif_wo_island_buffered_fixed, 1,2, iRaster_buffer_pixel )
    else:
        #no buffer is applied?
        sFilename_tif_wo_island_buffered_fixed = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island_buffered_fixed.tif')
        fix_raster_antimeridian_issue(sFilename_tif_wo_island, sFilename_tif_wo_island_buffered_fixed, 1, 2, 0)

    return sFilename_tif_wo_island_buffered_fixed, sFilename_vector_coastline


def geometries_bbox_overlap(bbox1: tuple, bbox2: tuple, tolerance: float = 1e-10) -> bool:
    """
    Check if two bounding boxes overlap with optional tolerance

    :param bbox1: (minX, maxX, minY, maxY) for first geometry
    :param bbox2: (minX, maxX, minY, maxY) for second geometry
    :param tolerance: small buffer for near-touching geometries
    :return: True if bounding boxes overlap
    """
    return not (bbox1[1] + tolerance < bbox2[0] or  # bbox1.maxX < bbox2.minX
                bbox2[1] + tolerance < bbox1[0] or  # bbox2.maxX < bbox1.minX
                bbox1[3] + tolerance < bbox2[2] or  # bbox1.maxY < bbox2.minY
                bbox2[3] + tolerance < bbox1[2])    # bbox2.maxY < bbox1.minY


def fix_naturalearth_hydrosheds_incompatibility(aFilename_hydrosheds_flowline: list,
                                                 sFilename_vector_naturalearth: str,
                                                 sFilename_vector_naturalearth_updated: str) -> None:
    """
    Fix incompatibility between Natural Earth land polygons and HydroSHEDS flowlines.

    This function identifies flowlines that are outside the Natural Earth land polygon
    and creates buffer zones around them to extend the land polygon, ensuring consistency
    between the land mask and river network. Uses RTree for efficient spatial indexing.

    Args:
        aFilename_hydrosheds_flowline: List of HydroSHEDS flowline file paths
        sFilename_vector_naturalearth: Path to Natural Earth land polygon file
        sFilename_vector_naturalearth_updated: Path for the updated Natural Earth file

    Returns:
        None
    """

    nFlowline = len(aFilename_hydrosheds_flowline)

    # Open the Natural Earth land polygon dataset
    pDataset_naturalearth = ogr.Open(sFilename_vector_naturalearth)
    if pDataset_naturalearth is None:
        logger.error(f'Could not open Natural Earth file: {sFilename_vector_naturalearth}')
        return None

    pLayer_naturalearth = pDataset_naturalearth.GetLayer(0)
    if pLayer_naturalearth is None:
        logger.error('No layer found in Natural Earth file')
        return None

    # Get spatial reference from Natural Earth layer
    pSpatial_reference = pLayer_naturalearth.GetSpatialRef()

    # Build spatial index for Natural Earth land polygons
    logger.info('Building spatial index for Natural Earth land polygons...')
    land_geometries = []
    land_spatial_index = None

    land_spatial_index = RTreeindex()

    # Index all land polygons
    nLand_features = pLayer_naturalearth.GetFeatureCount()
    for idx, pFeature_land in enumerate(pLayer_naturalearth):
        pGeometry_land = pFeature_land.GetGeometryRef()
        if pGeometry_land is None:
            continue

        land_geometries.append(pGeometry_land.Clone())
        envelope = pGeometry_land.GetEnvelope()  # (minX, maxX, minY, maxY)
        # RTree expects (minX, minY, maxX, maxY)
        bbox = (envelope[0], envelope[2], envelope[1], envelope[3])
        land_spatial_index.insert(idx, bbox)


    logger.info(f'Indexed {len(land_geometries)} land polygon features')

    # Create a list to store buffer geometries for flowlines outside land
    aBuffer_geometries = []

    for iFlowline in range(nFlowline):
        sFilename_flowline = aFilename_hydrosheds_flowline[iFlowline]
        logger.info(f'Processing flowline file {iFlowline + 1}/{nFlowline}: {os.path.basename(sFilename_flowline)}')
        # Open the flowline dataset
        pDataset_flowline = ogr.Open(sFilename_flowline)
        if pDataset_flowline is None:
            logger.warning(f'Could not open flowline file: {sFilename_flowline}')
            continue

        pLayer_flowline = pDataset_flowline.GetLayer(0)
        if pLayer_flowline is None:
            logger.warning(f'No layer found in flowline file: {sFilename_flowline}')
            continue

        # Process each flowline feature
        nFlowline_features = pLayer_flowline.GetFeatureCount()
        nOutside_count = 0

        for pFeature_flowline in pLayer_flowline:
            pGeometry_flowline = pFeature_flowline.GetGeometryRef()
            if pGeometry_flowline is None:
                continue

            # Get flowline bounding box for spatial query
            flowline_envelope = pGeometry_flowline.GetEnvelope()  # (minX, maxX, minY, maxY)

            # Find candidate land polygons using spatial index
            candidate_indices = []


            # RTree expects (minX, minY, maxX, maxY)
            query_bbox = (flowline_envelope[0], flowline_envelope[2], flowline_envelope[1], flowline_envelope[3])
            candidate_indices = list(land_spatial_index.intersection(query_bbox))


            # Check if flowline intersects with any candidate land polygon
            bIntersects_land = False

            # Start with the original flowline
            whole_flowline = pGeometry_flowline.Clone()

            for candidate_idx in candidate_indices:
                if candidate_idx >= len(land_geometries):
                    continue

                pGeometry_land = land_geometries[candidate_idx]
                if pGeometry_land is None:
                    continue

                # Check if flowline intersects with land polygon
                if whole_flowline.Intersects(pGeometry_land):
                    bIntersects_land = True
                    # Get the part of flowline that is outside this land polygon
                    outside_part = whole_flowline.Difference(pGeometry_land)
                    if outside_part is not None and not outside_part.IsEmpty():
                        # Update remaining flowline to only the outside part
                        outside_flowline = outside_part.Clone()
                        # Check if the outside part has significant length
                        outside_flowline.FlattenTo2D()
                        sGeometry_type = outside_flowline.GetGeometryName()
                        if sGeometry_type == 'LINESTRING':
                            nOutside_count += 1
                            sWkt= outside_flowline.ExportToWkt()
                            #calcuate the length of the outside_part to determine the buffer distance
                            point_start = outside_flowline.GetPoint(0)
                            point_end = outside_flowline.GetPoint(outside_flowline.GetPointCount() - 1)
                            distance = calculate_distance_based_on_longitude_latitude(point_start[0], point_start[1], point_end[0], point_end[1])
                            # Use midpoint latitude for accurate buffer distance conversion
                            reference_latitude = (point_start[1] + point_end[1]) / 2.0
                            sWkt_buffer_polygon = create_wkt_buffer_distance(sWkt, distance/2.0, epsg=4326, reference_latitude=reference_latitude)
                            #create a geometry from the buffer wkt
                            pBuffer_geometry = ogr.CreateGeometryFromWkt(sWkt_buffer_polygon)
                            aBuffer_geometries.append(pBuffer_geometry)
                        else:
                            if sGeometry_type == 'MULTILINESTRING':
                                for i in range(outside_flowline.GetGeometryCount()):
                                    outside_part_i = outside_flowline.GetGeometryRef(i)
                                    if outside_part_i is None or outside_part_i.IsEmpty():
                                        continue
                                    nOutside_count += 1
                                    sWkt= outside_part_i.ExportToWkt()
                                    point_start = outside_part_i.GetPoint(0)
                                    point_end = outside_part_i.GetPoint(outside_part_i.GetPointCount() - 1)
                                    distance = calculate_distance_based_on_longitude_latitude(point_start[0], point_start[1], point_end[0], point_end[1])
                                    #calcuate the length of the outside_part to determine the buffer distance
                                    # Use midpoint latitude for accurate buffer distance conversion
                                    reference_latitude = (point_start[1] + point_end[1]) / 2.0
                                    sWkt_buffer_polygon = create_wkt_buffer_distance(sWkt, distance/2.0, epsg=4326, reference_latitude=reference_latitude)
                                    #create a geometry from the buffer wkt
                                    pBuffer_geometry = ogr.CreateGeometryFromWkt(sWkt_buffer_polygon)
                                    aBuffer_geometries.append(pBuffer_geometry)
                        #use this distance to create a circle geometry as buffer
                        pass

                # If bIntersects_land is True but no outside_flowline_parts,
                # it means the entire flowline is inside land polygons - no buffer needed

        logger.info(f'  Found {nOutside_count} flowline features (or parts) requiring buffers out of {nFlowline_features} total')

        # Clean up flowline dataset
        pDataset_flowline = None

    # Create the updated Natural Earth file

    #first create a tmp file to store the original natural earth features
    sExtension = os.path.splitext(sFilename_vector_naturalearth)[1]
    #get the base name without extension
    sBase_name = os.path.splitext(os.path.basename(sFilename_vector_naturalearth))[0]

    sFilename_vector_naturalearth_tmp = os.path.join(os.path.dirname(sFilename_vector_naturalearth), sBase_name + '_tmp'+ sExtension)


    sFormat = 'GeoJSON'  # Default format

    pDriver = ogr.GetDriverByName(sFormat)
    if pDriver is None:
        logger.error(f'Driver {sFormat} not available')
        return None

    # Remove existing output file if it exists
    if os.path.exists(sFilename_vector_naturalearth_tmp):
        pDriver.DeleteDataSource(sFilename_vector_naturalearth_tmp)

    # Create output dataset
    pDataset_output = pDriver.CreateDataSource(sFilename_vector_naturalearth_tmp)
    if pDataset_output is None:
        logger.error('Could not create output dataset')
        return None

    # Create output layer
    pLayer_output = pDataset_output.CreateLayer('naturalearth_updated', pSpatial_reference, geom_type=ogr.wkbPolygon)
    if pLayer_output is None:
        logger.error('Could not create output layer')
        return None

    # Copy field definitions from input layer
    pLayerDefn_input = pLayer_naturalearth.GetLayerDefn()
    for i in range(pLayerDefn_input.GetFieldCount()):
        pFieldDefn = pLayerDefn_input.GetFieldDefn(i)
        if pLayer_output.CreateField(pFieldDefn) != 0:
            logger.warning(f'Failed to create field: {pFieldDefn.GetName()}')

    # Copy original Natural Earth features to output
    pLayer_naturalearth.ResetReading()
    for pFeature in pLayer_naturalearth:
        pFeature_out = ogr.Feature(pLayer_output.GetLayerDefn())
        pFeature_out.SetFrom(pFeature)
        if pLayer_output.CreateFeature(pFeature_out) != 0:
            logger.warning('Failed to create feature in output layer')
        pFeature_out = None


    # Add buffer geometries as new features
    if aBuffer_geometries:
        logger.info(f'Adding {len(aBuffer_geometries)} buffer zones for flowlines outside land')
        for pBuffer_geom in aBuffer_geometries:
            if pBuffer_geom is None or pBuffer_geom.IsEmpty():
                continue

            pFeature_out = ogr.Feature(pLayer_output.GetLayerDefn())
            pFeature_out.SetGeometry(pBuffer_geom)
            if pLayer_output.CreateFeature(pFeature_out) != 0:
                logger.warning('Failed to create buffer feature in output layer')
            pFeature_out = None

    # Clean up
    pDataset_output = None
    pDataset_naturalearth = None

    merge_features(sFilename_vector_naturalearth_tmp, sFilename_vector_naturalearth_updated)

    logger.info(f'Successfully created updated Natural Earth file: {sFilename_vector_naturalearth_updated}')
    logger.info(f'Added {len(aBuffer_geometries)} buffer zones for rivers outside original land polygons')

    return

