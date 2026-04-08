import os
import json
from osgeo import ogr, osr

from mpas_land_mesh.classes.edge import pyedge
from mpas_land_mesh.classes.link import pycelllink


def export_vertex_to_geojson(
    aVertex_in,
    sFilename_json_in,
    iFlag_projected_in=None,
    pSpatial_reference_in=None,
    aAttribute_data=None,
):
    """
    Export a list of vertex objects to a GeoJSON file.

    Args:
        aVertex_in (list): List of vertex objects with dLongitude_degree /
                           dLatitude_degree (geographic) or dx / dy (projected)
                           attributes.
        sFilename_json_in (str): Output GeoJSON file path.
        iFlag_projected_in (optional): If not None, use projected coordinates
                                       (dx, dy). Defaults to geographic.
        pSpatial_reference_in (optional): OGR SpatialReference. Defaults to
                                          WGS84 (EPSG:4326).
        aAttribute_data (optional): List of connectivity values (one per vertex).
    """

    if os.path.exists(sFilename_json_in):
        os.remove(sFilename_json_in)

    iFlag_projected_in = 0 if iFlag_projected_in is None else 1

    if pSpatial_reference_in is None:
        pSpatial_reference_in = osr.SpatialReference()
        pSpatial_reference_in.ImportFromEPSG(4326)  # WGS84 lat/lon

    iFlag_attribute = 1 if aAttribute_data is not None else 0
    aAttribute = aAttribute_data if aAttribute_data is not None else []

    pDriver = ogr.GetDriverByName("GeoJSON")
    pDataset_json = pDriver.CreateDataSource(sFilename_json_in)
    pLayer_json = pDataset_json.CreateLayer(
        "vertex", pSpatial_reference_in, ogr.wkbPoint
    )
    pLayer_json.CreateField(
        ogr.FieldDefn("pointid", ogr.OFTInteger64)
    )
    if iFlag_attribute == 1:
        pLayer_json.CreateField(
            ogr.FieldDefn("connectivity", ogr.OFTInteger64)
        )

    pLayerDefn = pLayer_json.GetLayerDefn()
    pFeature_out = ogr.Feature(pLayerDefn)

    for lID, pVertex in enumerate(aVertex_in):
        pPoint = ogr.Geometry(ogr.wkbPoint)
        if iFlag_projected_in == 1:
            pPoint.AddPoint(pVertex.dx, pVertex.dy)
        else:
            pPoint.AddPoint(pVertex.dLongitude_degree, pVertex.dLatitude_degree)

        pGeometry_out = ogr.CreateGeometryFromWkb(pPoint.ExportToWkb())
        pFeature_out.SetGeometry(pGeometry_out)
        pFeature_out.SetField("pointid", lID + 1)

        if iFlag_attribute == 1:
            pFeature_out.SetField("connectivity", int(aAttribute[lID]))

        pLayer_json.CreateFeature(pFeature_out)

    pDataset_json.FlushCache()
    pDataset_json = pLayer_json = pFeature_out = None

    return

def export_flowline_to_geojson(
    aFlowline_in,
    sFilename_json_in,
    iFlag_projected_in=None,
    pSpatial_reference_in=None,
    aAttribute_field=None,
    aAttribute_data=None,
    aAttribute_dtype=None,
):
    """
    convert a flowlist object list to json format.
    This function should be used for stream flowline only.
    """

    if os.path.exists(sFilename_json_in):
        os.remove(sFilename_json_in)

    nFlowline = len(aFlowline_in)
    if iFlag_projected_in is None:
        iFlag_projected_in = 0
    else:
        iFlag_projected_in = 1

    if pSpatial_reference_in is None:
        pSpatial_reference_in = osr.SpatialReference()
        pSpatial_reference_in.ImportFromEPSG(4326)  # WGS84 lat/lon
    else:
        pass

    iFlag_attribute = int(all([aAttribute_field, aAttribute_data, aAttribute_dtype]))

    if iFlag_attribute:
        nAttribute1, nAttribute2, nAttribute3 = map(
            len, [aAttribute_field, aAttribute_data, aAttribute_dtype]
        )
        if not nAttribute1 == nAttribute2 == nAttribute3 and nFlowline != len(
            aAttribute_data[0]
        ):
            print("The attribute is not correct, please check!")
            return

    pDriver_json = ogr.GetDriverByName("GeoJSON")
    pDataset_json = pDriver_json.CreateDataSource(sFilename_json_in)
    pLayer_json = pDataset_json.CreateLayer(
        "flowline", pSpatial_reference_in, ogr.wkbLineString
    )
    # Add one attribute
    pLayer_json.CreateField(
        ogr.FieldDefn("lineid", ogr.OFTInteger64)
    )  # long type for high resolution

    # add the other fields
    # if iFlag_attribute ==1:
    #    for i in range(nAttribute1):
    #        sField = aAttribute_field[i]
    #        dtype = aAttribute_dtype[i]
    #        if dtype == 'int':
    #            pLayer_json.CreateField(ogr.FieldDefn(sField, ogr.OFTInteger64))
    #            pass
    #        else:
    #            pLayer_json.CreateField(ogr.FieldDefn(sField, ogr.OFTReal))
    #            pass

    dtype_to_ogr = {"int": ogr.OFTInteger64, "float": ogr.OFTReal}
    if iFlag_attribute:
        for field, dtype in zip(aAttribute_field, aAttribute_dtype):
            pLayer_json.CreateField(ogr.FieldDefn(field, dtype_to_ogr[dtype]))

    pLayerDefn = pLayer_json.GetLayerDefn()
    pFeature_out = ogr.Feature(pLayerDefn)

    # lID = 0
    flag_to_attr = {1: ("dx", "dy"), 0: ("dLongitude_degree", "dLatitude_degree")}
    for lID in range(nFlowline):
        pFlowline = aFlowline_in[lID]
        pLine = ogr.Geometry(ogr.wkbLineString)
        for vertex in pFlowline.aVertex:
            pLine.AddPoint(
                *(getattr(vertex, attr) for attr in flag_to_attr[iFlag_projected_in])
            )

        pLine.FlattenTo2D()
        pFeature_out.SetGeometry(pLine)
        pFeature_out.SetField("lineid", lID + 1)

        if iFlag_attribute == 1:
            for sField, dtype, dummy in zip(
                aAttribute_field, aAttribute_dtype, aAttribute_data
            ):
                pFeature_out.SetField(
                    sField, int(dummy[lID]) if dtype == "int" else float(dummy[lID])
                )

        # Add new pFeature_shapefile to output Layer
        pLayer_json.CreateFeature(pFeature_out)
        pass

    pDataset_json.FlushCache()
    pDataset_json = pLayer_json = pFeature_out = None

    return

def export_flowline_info_to_json(
    aCell, aCell_intersect_in, aFlowline_in, sFilename_json_out
):
    # export the flowline topology to json
    ncell = len(aCell_intersect_in)
    nflowline = len(aFlowline_in)
    aLink = list()
    for i in range(1, nflowline + 1):
        pFlowline = aFlowline_in[i - 1]
        nVertex = pFlowline.nVertex
        nEdge = pFlowline.nEdge
        for j in range(1, nEdge + 1):
            pEdge = pFlowline.aEdge[j - 1]
            pVertex_start = pEdge.pVertex_start
            pVertex_end = pEdge.pVertex_end
            for k in range(ncell):
                if aCell_intersect_in[k].pVertex_center == pVertex_start:
                    pMpas_start = aCell_intersect_in[k]
                    pass
                if aCell_intersect_in[k].pVertex_center == pVertex_end:
                    pMpas_end = aCell_intersect_in[k]
                    pass
            pEdge_link = pyedge(pVertex_start, pVertex_end)
            pLink = pycelllink(pMpas_start, pMpas_end, pEdge_link)
            aLink.append(pLink)

    with open(sFilename_json_out, "w", encoding="utf-8") as f:
        sJson = json.dumps([json.loads(ob.tojson()) for ob in aLink], indent=4)
        f.write(sJson)
        f.close()

    return


