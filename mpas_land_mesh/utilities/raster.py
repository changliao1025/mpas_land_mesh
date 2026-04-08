"""
Raster processing utilities

Minimal implementations for raster operations
"""

import os
import numpy as np
from osgeo import gdal, ogr, osr


def convert_vector_to_global_raster(sFilename_vector_in, sFilename_raster_out,
                                     dResolution_x_in, dResolution_y_in,
                                     iFlag_boundary_only_in=1, dFill_value_in=1):
    """
    Convert vector file to global raster

    Args:
        sFilename_vector_in (str): Input vector file
        sFilename_raster_out (str): Output raster file
        dResolution_x_in (float): X resolution in degrees
        dResolution_y_in (float): Y resolution in degrees
        iFlag_boundary_only_in (int): Rasterize boundary only (1) or filled polygon (0)
        dFill_value_in (float): Value to burn into raster
    """
    # Global extent
    dLongitude_min = -180.0
    dLongitude_max = 180.0
    dLatitude_min = -90.0
    dLatitude_max = 90.0

    ncolumn = int((dLongitude_max - dLongitude_min) / dResolution_x_in)
    nrow = int((dLatitude_max - dLatitude_min) / dResolution_y_in)

    # Create output raster
    pDriver = gdal.GetDriverByName('GTiff')
    pDataset_out = pDriver.Create(sFilename_raster_out, ncolumn, nrow, 1, gdal.GDT_Byte)

    # Set geotransform
    pDataset_out.SetGeoTransform([dLongitude_min, dResolution_x_in, 0,
                                   dLatitude_max, 0, -dResolution_y_in])

    # Set projection
    pSpatialRef = osr.SpatialReference()
    pSpatialRef.ImportFromEPSG(4326)
    pDataset_out.SetProjection(pSpatialRef.ExportToWkt())

    # Initialize with zeros
    pBand = pDataset_out.GetRasterBand(1)
    pBand.SetNoDataValue(0)
    pBand.Fill(0)

    # Open vector file
    pDataset_vector = ogr.Open(sFilename_vector_in)
    pLayer = pDataset_vector.GetLayer(0)

    # Rasterize
    if iFlag_boundary_only_in == 1:
        # Boundary only
        gdal.RasterizeLayer(pDataset_out, [1], pLayer,
                           burn_values=[dFill_value_in],
                           options=['ALL_TOUCHED=TRUE'])
    else:
        # Filled polygon
        gdal.RasterizeLayer(pDataset_out, [1], pLayer,
                           burn_values=[dFill_value_in])

    pBand.FlushCache()
    pDataset_out = None
    pDataset_vector = None


def create_raster_buffer_zone(sFilename_raster_in, sFilename_raster_out,
                               iBuffer_value_in, iFill_value_in):
    """
    Create buffer zone around non-zero pixels in raster

    Args:
        sFilename_raster_in (str): Input raster
        sFilename_raster_out (str): Output raster
        iBuffer_value_in (int): Value for buffer pixels
        iFill_value_in (int): Value to buffer around
    """
    pDataset_in = gdal.Open(sFilename_raster_in)
    pBand_in = pDataset_in.GetRasterBand(1)

    aData_in = pBand_in.ReadAsArray()
    nrow, ncolumn = aData_in.shape

    # Create output array
    aData_out = np.copy(aData_in)

    # Find pixels with fill value
    mask = (aData_in == iFill_value_in)

    # Create buffer by expanding mask
    from scipy import ndimage
    struct = ndimage.generate_binary_structure(2, 1)  # 4-connectivity
    buffered = ndimage.binary_dilation(mask, structure=struct, iterations=1)

    # Set buffer pixels (that weren't already filled)
    buffer_mask = buffered & (~mask)
    aData_out[buffer_mask] = iBuffer_value_in

    # Create output raster
    pDriver = gdal.GetDriverByName('GTiff')
    pDataset_out = pDriver.Create(sFilename_raster_out, ncolumn, nrow, 1, gdal.GDT_Byte)

    pDataset_out.SetGeoTransform(pDataset_in.GetGeoTransform())
    pDataset_out.SetProjection(pDataset_in.GetProjection())

    pBand_out = pDataset_out.GetRasterBand(1)
    pBand_out.WriteArray(aData_out)
    pBand_out.SetNoDataValue(0)

    pBand_out.FlushCache()
    pDataset_out = None
    pDataset_in = None


def fix_raster_antimeridian_issue(sFilename_raster_in, sFilename_raster_out,
                                   iBuffer_value_in, iFill_value_in,
                                   iRaster_buffer_pixel=2):
    """
    Fix antimeridian discontinuity issues in global rasters

    Args:
        sFilename_raster_in (str): Input raster
        sFilename_raster_out (str): Output raster
        iBuffer_value_in (int): Buffer value
        iFill_value_in (int): Fill value
        iRaster_buffer_pixel (int): Number of pixels to buffer at edges
    """
    pDataset_in = gdal.Open(sFilename_raster_in)
    pBand_in = pDataset_in.GetRasterBand(1)

    aData = pBand_in.ReadAsArray()
    nrow, ncolumn = aData.shape

    # Check left and right edges for discontinuities
    left_edge = aData[:, 0:iRaster_buffer_pixel]
    right_edge = aData[:, -iRaster_buffer_pixel:]

    # If there are fill values on one edge and buffer values on the other,
    # propagate the buffer across the antimeridian
    for i in range(nrow):
        left_has_fill = np.any(left_edge[i, :] == iFill_value_in)
        right_has_fill = np.any(right_edge[i, :] == iFill_value_in)

        if left_has_fill and not right_has_fill:
            # Extend buffer from right to left
            aData[i, 0:iRaster_buffer_pixel] = iBuffer_value_in
        elif right_has_fill and not left_has_fill:
            # Extend buffer from left to right
            aData[i, -iRaster_buffer_pixel:] = iBuffer_value_in

    # Create output
    pDriver = gdal.GetDriverByName('GTiff')
    pDataset_out = pDriver.Create(sFilename_raster_out, ncolumn, nrow, 1, gdal.GDT_Byte)

    pDataset_out.SetGeoTransform(pDataset_in.GetGeoTransform())
    pDataset_out.SetProjection(pDataset_in.GetProjection())

    pBand_out = pDataset_out.GetRasterBand(1)
    pBand_out.WriteArray(aData)
    pBand_out.SetNoDataValue(0)

    pBand_out.FlushCache()
    pDataset_out = None
    pDataset_in = None


def gdal_read_geotiff_file(sFilename_in):
    """
    Read a GeoTIFF file and return its data and metadata.

    Local replacement for pyearth.gis.gdal.read.raster.gdal_read_geotiff_file.

    Args:
        sFilename_in (str): Path to the GeoTIFF file.

    Returns:
        dict with keys:
            - 'dataOut'     : 2D numpy array of raster data
            - 'originX'     : X coordinate of the upper-left corner
            - 'originY'     : Y coordinate of the upper-left corner
            - 'pixelWidth'  : Pixel width (positive, degrees per pixel in X)
            - 'pixelHeight' : Pixel height (negative for north-up rasters)
            - 'missingValue': NoData value (or -9999 if not set)
            - 'nrow'        : Number of rows
            - 'ncolumn'     : Number of columns
            - 'projection'  : WKT projection string
    """
    pDataset = gdal.Open(sFilename_in, gdal.GA_ReadOnly)
    if pDataset is None:
        raise FileNotFoundError(f"Cannot open raster file: {sFilename_in}")

    pBand = pDataset.GetRasterBand(1)
    aData = pBand.ReadAsArray()

    gt = pDataset.GetGeoTransform()
    originX = gt[0]       # X coordinate of upper-left corner
    originY = gt[3]       # Y coordinate of upper-left corner
    pixelWidth = gt[1]    # Pixel width (positive)
    pixelHeight = gt[5]   # Pixel height (negative for north-up)

    nodata = pBand.GetNoDataValue()
    missingValue = nodata if nodata is not None else -9999

    nrow, ncolumn = aData.shape
    projection = pDataset.GetProjection()

    pDataset = None  # close

    return {
        'dataOut': aData,
        'originX': originX,
        'originY': originY,
        'pixelWidth': pixelWidth,
        'pixelHeight': pixelHeight,
        'missingValue': missingValue,
        'nrow': nrow,
        'ncolumn': ncolumn,
        'projection': projection,
    }
