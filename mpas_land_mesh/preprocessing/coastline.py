"""
Coastline processing utilities

Simplified implementations for coastline operations
"""

import logging
import os
from pathlib import Path
from typing import Literal, Optional

_GDAL_IMPORT_ERROR: Optional[Exception] = None
try:
    from osgeo import gdal, ogr, osr
except Exception as exc:  # pragma: no cover - environment-dependent
    gdal = None
    ogr = None
    osr = None
    _GDAL_IMPORT_ERROR = exc

try:
    from rtree.index import Index as RTreeindex
except Exception:  # pragma: no cover - environment-dependent
    RTreeindex = None

from pyearth.gis.gdal.gdal_vector_format_support import (
    get_vector_driver_from_filename,
    get_vector_format_from_filename,
)
from pyearth.gis.geometry import calculate_distance_based_on_longitude_latitude
from pyearth.toolbox.analysis.image.raster_process import (
    create_raster_buffer_zone,
    fix_raster_antimeridian_issue,
)
from pyearth.toolbox.conversion.convert_vector_to_global_raster import (
    convert_vector_to_global_raster,
)
from pyearth.toolbox.management.vector.merge_features import merge_features
from pyearth.toolbox.management.vector.remove_internal_polygon import (
    remove_internal_polygon,
)
from pyearth.toolbox.management.vector.remove_small_polygon import (
    remove_small_polygon,
)
from pyearthbuffer import gcsbuffer

if gdal is not None:
    gdal.UseExceptions()

# Set up logger for this module
logger = logging.getLogger(__name__)


from pyearth.toolbox.data.ocean.define_land_ocean_mask import (
    create_land_ocean_vector_mask_naturalearth,
)


def create_land_ocean_mask_from_naturalearth(
    sWorkspace_coastline_output: str,
    dResolution_x_in: float,
    dResolution_y_in: float,
    dThreshold_area_island: float,
    dResolution_coastline_buffer: float,
    iRaster_buffer_pixel: int = 1,
) -> tuple:

    sFilename_naturalearth = os.path.join(
        sWorkspace_coastline_output, "land_ocean_mask_naturalearth.geojson"
    )
    create_land_ocean_vector_mask_naturalearth(sFilename_naturalearth)

    sFilename_vector_coastline = os.path.join(
        sWorkspace_coastline_output, "land_ocean_mask_wo_island.geojson"
    )
    remove_small_polygon(
        sFilename_naturalearth, sFilename_vector_coastline, dThreshold_area_island
    )

    sFilename_vector_coastline_wo_internal_ring = os.path.join(
        sWorkspace_coastline_output,
        "land_ocean_mask_wo_island_removed_internal.geojson",
    )
    remove_internal_polygon(
        sFilename_vector_coastline, sFilename_vector_coastline_wo_internal_ring
    )

    sFilename_tif_wo_island = os.path.join(
        sWorkspace_coastline_output, "land_ocean_mask_wo_island.tif"
    )
    convert_vector_to_global_raster(
        sFilename_vector_coastline_wo_internal_ring,
        sFilename_tif_wo_island,
        dResolution_x_in,
        dResolution_y_in,
        iFlag_boundary_only_in=0,
        dFill_value_in=2,
    )

    if iRaster_buffer_pixel > 0:
        sFilename_tif_wo_island_buffered = os.path.join(
            sWorkspace_coastline_output, "land_ocean_mask_wo_island_buffered.tif"
        )
        create_raster_buffer_zone(
            sFilename_tif_wo_island,
            sFilename_tif_wo_island_buffered,
            1,
            int(iRaster_buffer_pixel),
        )
        sFilename_tif_wo_island_buffered_fixed = os.path.join(
            sWorkspace_coastline_output, "land_ocean_mask_wo_island_buffered_fixed.tif"
        )
        # fix dateline issue
        fix_raster_antimeridian_issue(
            sFilename_tif_wo_island_buffered,
            sFilename_tif_wo_island_buffered_fixed,
            1,
            2,
            iRaster_buffer_pixel,
        )
    else:
        # no buffer is applied?
        sFilename_tif_wo_island_buffered_fixed = os.path.join(
            sWorkspace_coastline_output, "land_ocean_mask_wo_island_buffered_fixed.tif"
        )
        fix_raster_antimeridian_issue(
            sFilename_tif_wo_island, sFilename_tif_wo_island_buffered_fixed, 1, 2, 0
        )

    return (
        sFilename_tif_wo_island_buffered_fixed,
        sFilename_vector_coastline_wo_internal_ring,
    )


create_land_ocean_mask = create_land_ocean_mask_from_naturalearth


def geometries_bbox_overlap(
    bbox1: tuple, bbox2: tuple, tolerance: float = 1e-10
) -> bool:
    """
    Check if two bounding boxes overlap with optional tolerance

    :param bbox1: (minX, maxX, minY, maxY) for first geometry
    :param bbox2: (minX, maxX, minY, maxY) for second geometry
    :param tolerance: small buffer for near-touching geometries
    :return: True if bounding boxes overlap
    """
    return not (
        bbox1[1] + tolerance < bbox2[0]  # bbox1.maxX < bbox2.minX
        or bbox2[1] + tolerance < bbox1[0]  # bbox2.maxX < bbox1.minX
        or bbox1[3] + tolerance < bbox2[2]  # bbox1.maxY < bbox2.minY
        or bbox2[3] + tolerance < bbox1[2]
    )  # bbox2.maxY < bbox1.minY


def fix_naturalearth_hydrosheds_incompatibility(
    aFilename_hydrosheds_flowline: list,
    sFilename_vector_naturalearth: str,
    sFilename_vector_naturalearth_updated: str,
) -> None:
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

    if os.path.exists(sFilename_vector_naturalearth_updated):
        os.remove(sFilename_vector_naturalearth_updated)
        logger.info(
            f"Removed existing updated Natural Earth file: {sFilename_vector_naturalearth_updated}"
        )

    nFlowline = len(aFilename_hydrosheds_flowline)

    # Open the Natural Earth land polygon dataset
    pDataset_naturalearth = ogr.Open(sFilename_vector_naturalearth)
    if pDataset_naturalearth is None:
        logger.error(
            f"Could not open Natural Earth file: {sFilename_vector_naturalearth}"
        )
        return None

    pLayer_naturalearth = pDataset_naturalearth.GetLayer(0)
    if pLayer_naturalearth is None:
        logger.error("No layer found in Natural Earth file")
        return None

    # Get spatial reference from Natural Earth layer
    pSpatial_reference = pLayer_naturalearth.GetSpatialRef()

    # Build spatial index for Natural Earth land polygons
    logger.info("Building spatial index for Natural Earth land polygons...")
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

    logger.info(f"Indexed {len(land_geometries)} land polygon features")

    # Create a list to store buffer geometries for flowlines outside land
    aBuffer_geometries = []

    # -----------------------------------------------------------------------
    # Debug: prepare a directory and GeoJSON dataset to collect all
    # outside_flowline segments so they can be inspected in a GIS viewer.
    # -----------------------------------------------------------------------
    sWorkspace_debug = os.path.join(
        os.path.dirname(sFilename_vector_naturalearth), "debug_outside_flowlines"
    )
    os.makedirs(sWorkspace_debug, exist_ok=True)
    sFilename_debug_flowlines = os.path.join(
        sWorkspace_debug, "outside_flowlines.geojson"
    )
    sFilename_debug_buffers = os.path.join(
        sWorkspace_debug, "outside_flowline_buffers.geojson"
    )

    pDriver_debug = ogr.GetDriverByName("GeoJSON")
    pSrs_debug = osr.SpatialReference()
    pSrs_debug.ImportFromEPSG(4326)

    # Remove stale debug files from previous runs
    for _f in (sFilename_debug_flowlines, sFilename_debug_buffers):
        if os.path.exists(_f):
            pDriver_debug.DeleteDataSource(_f)

    pDataset_debug_fl = pDriver_debug.CreateDataSource(sFilename_debug_flowlines)
    pLayer_debug_fl = pDataset_debug_fl.CreateLayer(
        "outside_flowlines", pSrs_debug, ogr.wkbLineString
    )
    # Add informational fields
    pLayer_debug_fl.CreateField(ogr.FieldDefn("flowline_file", ogr.OFTString))
    pLayer_debug_fl.CreateField(ogr.FieldDefn("segment_idx", ogr.OFTInteger))
    pLayer_debug_fl.CreateField(ogr.FieldDefn("geom_type", ogr.OFTString))
    pLayer_debug_fl.CreateField(ogr.FieldDefn("distance_m", ogr.OFTReal))

    pDataset_debug_buf = pDriver_debug.CreateDataSource(sFilename_debug_buffers)
    pLayer_debug_buf = pDataset_debug_buf.CreateLayer(
        "outside_flowline_buffers", pSrs_debug, ogr.wkbPolygon
    )
    pLayer_debug_buf.CreateField(ogr.FieldDefn("flowline_file", ogr.OFTString))
    pLayer_debug_buf.CreateField(ogr.FieldDefn("segment_idx", ogr.OFTInteger))
    pLayer_debug_buf.CreateField(ogr.FieldDefn("distance_m", ogr.OFTReal))

    logger.info(f"Debug outside-flowline files will be written to: {sWorkspace_debug}")

    for iFlowline in range(nFlowline):
        sFilename_flowline = aFilename_hydrosheds_flowline[iFlowline]
        logger.info(
            f"Processing flowline file {iFlowline + 1}/{nFlowline}: {os.path.basename(sFilename_flowline)}"
        )
        # Open the flowline dataset
        pDataset_flowline = ogr.Open(sFilename_flowline)
        if pDataset_flowline is None:
            logger.warning(f"Could not open flowline file: {sFilename_flowline}")
            continue

        pLayer_flowline = pDataset_flowline.GetLayer(0)
        if pLayer_flowline is None:
            logger.warning(f"No layer found in flowline file: {sFilename_flowline}")
            continue

        # Process each flowline feature
        nFlowline_features = pLayer_flowline.GetFeatureCount()
        nOutside_count = 0
        sFlowline_basename = os.path.basename(sFilename_flowline)

        for pFeature_flowline in pLayer_flowline:
            pGeometry_flowline = pFeature_flowline.GetGeometryRef()
            if pGeometry_flowline is None:
                continue

            # Get flowline bounding box for spatial query
            flowline_envelope = (
                pGeometry_flowline.GetEnvelope()
            )  # (minX, maxX, minY, maxY)

            # Find candidate land polygons using spatial index
            candidate_indices = []

            # RTree expects (minX, minY, maxX, maxY)
            query_bbox = (
                flowline_envelope[0],
                flowline_envelope[2],
                flowline_envelope[1],
                flowline_envelope[3],
            )
            candidate_indices = list(land_spatial_index.intersection(query_bbox))

            # Progressively subtract every intersecting land polygon so the
            # final remainder captures ALL outside segments, not just the last.
            bIntersects_land = False
            whole_flowline = pGeometry_flowline.Clone()
            whole_flowline.FlattenTo2D()

            for candidate_idx in candidate_indices:
                if candidate_idx >= len(land_geometries):
                    continue
                if whole_flowline is None or whole_flowline.IsEmpty():
                    break

                pGeometry_land = land_geometries[candidate_idx]
                if pGeometry_land is None:
                    continue

                if whole_flowline.Intersects(pGeometry_land):
                    bIntersects_land = True
                    # Subtract this land polygon and keep the outside remainder
                    # for the next iteration (progressive subtraction).
                    outside_part = whole_flowline.Difference(pGeometry_land)
                    if outside_part is not None and not outside_part.IsEmpty():
                        outside_part.FlattenTo2D()
                        whole_flowline = outside_part  # <-- update for next iteration
                    else:
                        whole_flowline = None
                        break

            # Process the final outside remainder after all land polygons are subtracted
            if (
                bIntersects_land
                and whole_flowline is not None
                and not whole_flowline.IsEmpty()
            ):
                outside_flowline = whole_flowline
                sGeometry_type = outside_flowline.GetGeometryName()

                if sGeometry_type == "LINESTRING":
                    nOutside_count += 1
                    point_start = outside_flowline.GetPoint(0)
                    point_end = outside_flowline.GetPoint(
                        outside_flowline.GetPointCount() - 1
                    )
                    distance = calculate_distance_based_on_longitude_latitude(
                        point_start[0], point_start[1], point_end[0], point_end[1]
                    )
                    if distance <= 0:
                        continue
                    pBuffer_geometry = gcsbuffer(
                        outside_flowline, distance / 2.0, output_format="ogr"
                    )
                    aBuffer_geometries.append(pBuffer_geometry)

                    # --- debug ---
                    _pFeat_fl = ogr.Feature(pLayer_debug_fl.GetLayerDefn())
                    _pFeat_fl.SetGeometry(outside_flowline.Clone())
                    _pFeat_fl.SetField("flowline_file", sFlowline_basename)
                    _pFeat_fl.SetField("segment_idx", nOutside_count)
                    _pFeat_fl.SetField("geom_type", sGeometry_type)
                    _pFeat_fl.SetField("distance_m", distance)
                    pLayer_debug_fl.CreateFeature(_pFeat_fl)
                    _pFeat_fl = None

                    _pFeat_buf = ogr.Feature(pLayer_debug_buf.GetLayerDefn())
                    _pFeat_buf.SetGeometry(pBuffer_geometry.Clone())
                    _pFeat_buf.SetField("flowline_file", sFlowline_basename)
                    _pFeat_buf.SetField("segment_idx", nOutside_count)
                    _pFeat_buf.SetField("distance_m", distance)
                    pLayer_debug_buf.CreateFeature(_pFeat_buf)
                    _pFeat_buf = None
                    # --- end debug ---

                elif sGeometry_type == "MULTILINESTRING":
                    for i in range(outside_flowline.GetGeometryCount()):
                        outside_part_i = outside_flowline.GetGeometryRef(i)
                        if outside_part_i is None or outside_part_i.IsEmpty():
                            continue
                        nOutside_count += 1
                        point_start = outside_part_i.GetPoint(0)
                        point_end = outside_part_i.GetPoint(
                            outside_part_i.GetPointCount() - 1
                        )
                        distance = calculate_distance_based_on_longitude_latitude(
                            point_start[0], point_start[1], point_end[0], point_end[1]
                        )
                        if distance <= 0:
                            continue
                        pBuffer_geometry = gcsbuffer(
                            outside_part_i, distance / 2.0, output_format="ogr"
                        )
                        aBuffer_geometries.append(pBuffer_geometry)

                        # --- debug ---
                        _pFeat_fl = ogr.Feature(pLayer_debug_fl.GetLayerDefn())
                        _pFeat_fl.SetGeometry(outside_part_i.Clone())
                        _pFeat_fl.SetField("flowline_file", sFlowline_basename)
                        _pFeat_fl.SetField("segment_idx", nOutside_count)
                        _pFeat_fl.SetField("geom_type", "MULTILINESTRING_part")
                        _pFeat_fl.SetField("distance_m", distance)
                        pLayer_debug_fl.CreateFeature(_pFeat_fl)
                        _pFeat_fl = None

                        _pFeat_buf = ogr.Feature(pLayer_debug_buf.GetLayerDefn())
                        _pFeat_buf.SetGeometry(pBuffer_geometry.Clone())
                        _pFeat_buf.SetField("flowline_file", sFlowline_basename)
                        _pFeat_buf.SetField("segment_idx", nOutside_count)
                        _pFeat_buf.SetField("distance_m", distance)
                        pLayer_debug_buf.CreateFeature(_pFeat_buf)
                        _pFeat_buf = None
                        # --- end debug ---

                else:
                    print(
                        f"Warning: Unexpected geometry type for outside flowline: "
                        f"{sGeometry_type}, skipping buffer creation"
                    )

            elif not bIntersects_land:
                # The flowline does not intersect any land polygon at all —
                # the entire flowline is outside land and must be buffered.
                outside_flowline = pGeometry_flowline.Clone()
                outside_flowline.FlattenTo2D()
                sGeometry_type = outside_flowline.GetGeometryName()

                if sGeometry_type == "LINESTRING":
                    nOutside_count += 1
                    point_start = outside_flowline.GetPoint(0)
                    point_end = outside_flowline.GetPoint(
                        outside_flowline.GetPointCount() - 1
                    )
                    distance = calculate_distance_based_on_longitude_latitude(
                        point_start[0], point_start[1], point_end[0], point_end[1]
                    )
                    if distance <= 0:
                        continue
                    pBuffer_geometry = gcsbuffer(
                        outside_flowline, distance / 2.0, output_format="ogr"
                    )
                    aBuffer_geometries.append(pBuffer_geometry)

                    # --- debug ---
                    _pFeat_fl = ogr.Feature(pLayer_debug_fl.GetLayerDefn())
                    _pFeat_fl.SetGeometry(outside_flowline.Clone())
                    _pFeat_fl.SetField("flowline_file", sFlowline_basename)
                    _pFeat_fl.SetField("segment_idx", nOutside_count)
                    _pFeat_fl.SetField("geom_type", "LINESTRING_no_land")
                    _pFeat_fl.SetField("distance_m", distance)
                    pLayer_debug_fl.CreateFeature(_pFeat_fl)
                    _pFeat_fl = None

                    _pFeat_buf = ogr.Feature(pLayer_debug_buf.GetLayerDefn())
                    _pFeat_buf.SetGeometry(pBuffer_geometry.Clone())
                    _pFeat_buf.SetField("flowline_file", sFlowline_basename)
                    _pFeat_buf.SetField("segment_idx", nOutside_count)
                    _pFeat_buf.SetField("distance_m", distance)
                    pLayer_debug_buf.CreateFeature(_pFeat_buf)
                    _pFeat_buf = None
                    # --- end debug ---

                elif sGeometry_type == "MULTILINESTRING":
                    for i in range(outside_flowline.GetGeometryCount()):
                        outside_part_i = outside_flowline.GetGeometryRef(i)
                        if outside_part_i is None or outside_part_i.IsEmpty():
                            continue
                        nOutside_count += 1
                        point_start = outside_part_i.GetPoint(0)
                        point_end = outside_part_i.GetPoint(
                            outside_part_i.GetPointCount() - 1
                        )
                        distance = calculate_distance_based_on_longitude_latitude(
                            point_start[0], point_start[1], point_end[0], point_end[1]
                        )
                        if distance <= 0:
                            continue
                        pBuffer_geometry = gcsbuffer(
                            outside_part_i, distance / 2.0, output_format="ogr"
                        )
                        aBuffer_geometries.append(pBuffer_geometry)

                        # --- debug ---
                        _pFeat_fl = ogr.Feature(pLayer_debug_fl.GetLayerDefn())
                        _pFeat_fl.SetGeometry(outside_part_i.Clone())
                        _pFeat_fl.SetField("flowline_file", sFlowline_basename)
                        _pFeat_fl.SetField("segment_idx", nOutside_count)
                        _pFeat_fl.SetField("geom_type", "MULTILINESTRING_no_land")
                        _pFeat_fl.SetField("distance_m", distance)
                        pLayer_debug_fl.CreateFeature(_pFeat_fl)
                        _pFeat_fl = None

                        _pFeat_buf = ogr.Feature(pLayer_debug_buf.GetLayerDefn())
                        _pFeat_buf.SetGeometry(pBuffer_geometry.Clone())
                        _pFeat_buf.SetField("flowline_file", sFlowline_basename)
                        _pFeat_buf.SetField("segment_idx", nOutside_count)
                        _pFeat_buf.SetField("distance_m", distance)
                        pLayer_debug_buf.CreateFeature(_pFeat_buf)
                        _pFeat_buf = None
                        # --- end debug ---

        logger.info(
            f"  Found {nOutside_count} flowline features (or parts) requiring buffers out of {nFlowline_features} total"
        )

        # Clean up flowline dataset
        pDataset_flowline = None

    # Flush and close debug datasets
    pDataset_debug_fl = None
    pDataset_debug_buf = None
    logger.info(
        f"Debug outside-flowline segments saved to: {sFilename_debug_flowlines}"
    )
    logger.info(f"Debug outside-flowline buffers  saved to: {sFilename_debug_buffers}")

    # Create the updated Natural Earth file

    # first create a tmp file to store the original natural earth features
    sExtension = os.path.splitext(sFilename_vector_naturalearth)[1]
    # get the base name without extension
    sBase_name = os.path.splitext(os.path.basename(sFilename_vector_naturalearth))[0]

    sFilename_vector_naturalearth_tmp = os.path.join(
        os.path.dirname(sFilename_vector_naturalearth), sBase_name + "_tmp" + sExtension
    )
    sFilename_vector_naturalearth_tmp2 = os.path.join(
        os.path.dirname(sFilename_vector_naturalearth),
        sBase_name + "_tmp2" + sExtension,
    )

    sFormat = "GeoJSON"  # Default format

    pDriver = ogr.GetDriverByName(sFormat)
    if pDriver is None:
        logger.error(f"Driver {sFormat} not available")
        return None

    # Remove existing output file if it exists
    if os.path.exists(sFilename_vector_naturalearth_tmp):
        pDriver.DeleteDataSource(sFilename_vector_naturalearth_tmp)

    if os.path.exists(sFilename_vector_naturalearth_tmp2):
        pDriver.DeleteDataSource(sFilename_vector_naturalearth_tmp2)

    # Create output dataset
    pDataset_output = pDriver.CreateDataSource(sFilename_vector_naturalearth_tmp)
    if pDataset_output is None:
        logger.error("Could not create output dataset")
        return None

    # Create output layer
    pLayer_output = pDataset_output.CreateLayer(
        "naturalearth_updated", pSpatial_reference, geom_type=ogr.wkbPolygon
    )
    if pLayer_output is None:
        logger.error("Could not create output layer")
        return None

    # Copy field definitions from input layer
    pLayerDefn_input = pLayer_naturalearth.GetLayerDefn()
    for i in range(pLayerDefn_input.GetFieldCount()):
        pFieldDefn = pLayerDefn_input.GetFieldDefn(i)
        if pLayer_output.CreateField(pFieldDefn) != 0:
            logger.warning(f"Failed to create field: {pFieldDefn.GetName()}")

    # Copy original Natural Earth features to output
    pLayer_naturalearth.ResetReading()
    for pFeature in pLayer_naturalearth:
        pFeature_out = ogr.Feature(pLayer_output.GetLayerDefn())
        pFeature_out.SetFrom(pFeature)
        if pLayer_output.CreateFeature(pFeature_out) != 0:
            logger.warning("Failed to create feature in output layer")
        pFeature_out = None

    # Add buffer geometries as new features
    if aBuffer_geometries:
        logger.info(
            f"Adding {len(aBuffer_geometries)} buffer zones for flowlines outside land"
        )
        for pBuffer_geom in aBuffer_geometries:
            if pBuffer_geom is None or pBuffer_geom.IsEmpty():
                continue

            pFeature_out = ogr.Feature(pLayer_output.GetLayerDefn())
            pFeature_out.SetGeometry(pBuffer_geom)
            if pLayer_output.CreateFeature(pFeature_out) != 0:
                logger.warning("Failed to create buffer feature in output layer")
            pFeature_out = None

    # Clean up
    pDataset_output = None
    pDataset_naturalearth = None

    merge_features(
        sFilename_vector_naturalearth_tmp, sFilename_vector_naturalearth_tmp2
    )

    remove_internal_polygon(
        sFilename_vector_naturalearth_tmp2, sFilename_vector_naturalearth_updated
    )

    logger.info(
        f"Successfully created updated Natural Earth file: {sFilename_vector_naturalearth_updated}"
    )
    logger.info(
        f"Added {len(aBuffer_geometries)} buffer zones for rivers outside original land polygons"
    )

    return
