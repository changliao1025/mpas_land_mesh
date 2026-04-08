import os
import math
import numpy as np
from osgeo import ogr, osr, gdal

from mpas_land_mesh.mesh.mpas.convert_attributes import convert_gcs_attributes_to_cell
from mpas_land_mesh.mesh.jigsaw.run_jigsaw import run_jigsaw
from mpas_land_mesh.utilities.geometry import (
    convert_360_to_180,
    split_international_date_line_polygon_coordinates,
)

gdal.UseExceptions()


def calculate_potentiometric(dBed_elevation_in, dIce_thickness_in):
    """
    Calculate potentiometric surface elevation for ice-covered areas.

    The potentiometric surface is the hydraulic head beneath an ice sheet,
    calculated as bed elevation plus the pressure head from the ice column.

    Args:
        dBed_elevation_in (float): Bed elevation in meters
        dIce_thickness_in (float): Ice thickness in meters

    Returns:
        float: Potentiometric surface elevation in meters
    """
    dDensity_ice = 917.0   # kg/m^3
    dDensity_water = 1000.0  # kg/m^3
    dRatio = dDensity_ice / dDensity_water
    dPotentiometric = dBed_elevation_in + dRatio * dIce_thickness_in
    return dPotentiometric


def create_mpas_mesh(
    sFilename_output_in,
    iFlag_global_in=None,
    iFlag_use_mesh_dem_in=None,
    iFlag_save_mesh_in=None,
    iFlag_run_jigsaw_in=None,
    iFlag_antarctic_in=None,
    iFlag_arctic_in=None,
    iFlag_fill_hole_in=None,
    pBoundary_in=None,
    sWorkspace_jigsaw_in=None,
    sFilename_mpas_mesh_netcdf_in=None,
    sFilename_jigsaw_mesh_netcdf_in=None,
    sFilename_land_ocean_mask_in=None,
    aConfig_jigsaw_in=None,
    aFilename_river_network_in=None,
    aFilename_watershed_boundary_in=None,
    aFilename_lake_boundary_in=None,
    aFilename_coastline_in=None,
    iFlag_read_mesh_in=None,
    iFlag_generate_mesh_in=None,
):
    """
    Create a MPAS mesh from a JIGSAW-generated NetCDF file.

    Args:
        sFilename_output_in (str): Output GeoJSON filename for the mesh cells.
        iFlag_global_in (int): 1 for global mesh, 0 for regional. Default 0.
        iFlag_use_mesh_dem_in (int): 1 to include DEM elevation fields. Default 0.
        iFlag_save_mesh_in (int): 1 to save GeoJSON output. Default 1.
        iFlag_run_jigsaw_in (int): 1 to run JIGSAW before reading mesh. Default 0.
        iFlag_antarctic_in (int): 1 to process Antarctic region only. Default 0.
        iFlag_arctic_in (int): 1 to process Arctic region only. Default 0.
        iFlag_fill_hole_in (int): Reserved for future use.
        pBoundary_in (str): WKT string of the domain boundary polygon.
        sWorkspace_jigsaw_in (str): Path to JIGSAW workspace directory.
        sFilename_mpas_mesh_netcdf_in (str): Path to existing MPAS NetCDF mesh file.
        sFilename_jigsaw_mesh_netcdf_in (str): Path to JIGSAW output NetCDF file.
        sFilename_land_ocean_mask_in (str): Path to land/ocean mask file.
        aConfig_jigsaw_in (dict): JIGSAW configuration dictionary.
        aFilename_river_network_in (list): List of river network filenames.
        aFilename_watershed_boundary_in (list): List of watershed boundary filenames.
        aFilename_lake_boundary_in (list): List of lake boundary filenames.
        aFilename_coastline_in (list): List of coastline filenames.
        iFlag_read_mesh_in (int): Reserved for future use.
        iFlag_generate_mesh_in (int): Reserved for future use.

    Returns:
        list: List of pympas cell objects in the domain.
    """
    import netCDF4 as nc

    if iFlag_global_in is None:
        iFlag_global = 0
    else:
        iFlag_global = iFlag_global_in

    if iFlag_use_mesh_dem_in is None:
        iFlag_use_mesh_dem = 0
    else:
        iFlag_use_mesh_dem = iFlag_use_mesh_dem_in

    if iFlag_save_mesh_in is None:
        iFlag_save_mesh = 1
    else:
        iFlag_save_mesh = iFlag_save_mesh_in

    if iFlag_antarctic_in is None:
        iFlag_antarctic = 0
    else:
        iFlag_antarctic = iFlag_antarctic_in

    if iFlag_arctic_in is None:
        iFlag_arctic = 0
    else:
        iFlag_arctic = iFlag_arctic_in

    if sFilename_mpas_mesh_netcdf_in is None:
        print("This mesh file will be generated!")

    if pBoundary_in is None:
        pBoundary = None
    else:
        # for the reason that a geometry object will crash if the associated dataset is
        # closed, we must pass wkt string
        # https://gdal.org/api/python_gotchas.html
        pBoundary = ogr.CreateGeometryFromWkt(pBoundary_in)

    if os.path.exists(sFilename_output_in):
        os.remove(sFilename_output_in)

    iFlag_elevation_profile = 0
    iFlag_bed_elevation = 0
    iFlag_ice_thickness = 0

    pDriver_geojson = ogr.GetDriverByName("GeoJSON")
    pSpatial_reference_gcs = osr.SpatialReference()
    pSpatial_reference_gcs.ImportFromEPSG(4326)  # WGS84 lat/lon

    # geojson output
    if iFlag_save_mesh == 1:
        pDataset = pDriver_geojson.CreateDataSource(sFilename_output_in)
        pLayer = pDataset.CreateLayer("cell", pSpatial_reference_gcs, ogr.wkbPolygon)
        pLayer.CreateField(ogr.FieldDefn("cellid", ogr.OFTInteger64))
        pLayer.CreateField(ogr.FieldDefn("longitude", ogr.OFTReal))
        pLayer.CreateField(ogr.FieldDefn("latitude", ogr.OFTReal))
        pArea_field = ogr.FieldDefn("area", ogr.OFTReal)
        pArea_field.SetWidth(20)
        pArea_field.SetPrecision(2)
        pLayer.CreateField(pArea_field)
        if iFlag_use_mesh_dem == 1:
            pLayer.CreateField(ogr.FieldDefn("elevation_mean", ogr.OFTReal))
            pLayer.CreateField(ogr.FieldDefn("elevation_profile0", ogr.OFTReal))

        pLayerDefn = pLayer.GetLayerDefn()
        pFeature = ogr.Feature(pLayerDefn)

    if iFlag_run_jigsaw_in == 1:
        projector = [0.0, 0.0]
        geom, gprj, mesh, mprj = run_jigsaw(
            sWorkspace_jigsaw_in,
            projector,
            aConfig_in=aConfig_jigsaw_in,
            aFilename_river_network_in=aFilename_river_network_in,
            aFilename_watershed_boundary_in=aFilename_watershed_boundary_in,
            aFilename_lake_boundary_in=aFilename_lake_boundary_in,
            aFilename_coastline_in=aFilename_coastline_in,
        )

        # write output for ESM
        from mpas_land_mesh.mesh.jigsaw.saveesm import saveesm

        sFilename_culled_mesh, sFilename_invert_mesh = saveesm(
            sWorkspace_jigsaw_in,
            geom,
            mesh,
            sFilename_jigsaw_mesh_netcdf_in=sFilename_jigsaw_mesh_netcdf_in,
            sFilename_land_ocean_mask_in=sFilename_land_ocean_mask_in,
        )
        print("The generated MPAS mesh is: ", sFilename_invert_mesh)
        sFilename_mpas_mesh_netcdf_in = sFilename_invert_mesh
    else:
        # check whether the jigsaw mesh file exists
        if sFilename_jigsaw_mesh_netcdf_in is not None:
            if os.path.exists(sFilename_jigsaw_mesh_netcdf_in):
                pass
            else:
                print("JIGSAW Mesh file does not exist!")

    if os.path.exists(sFilename_mpas_mesh_netcdf_in):
        pass
    else:
        print("Mesh file does not exist!")
        return

    pDatasets_in = nc.Dataset(sFilename_mpas_mesh_netcdf_in)

    # read netcdf variables
    for sKey, aValue in pDatasets_in.variables.items():
        if sKey == "latCell":
            latCell0 = aValue
        if sKey == "lonCell":
            lonCell0 = aValue
        if sKey == "edgesOnCell":
            edgesOnCell0 = aValue
        if sKey == "cellsOnCell":
            cellsOnCell0 = aValue
        if sKey == "cellsOnEdge":
            cellsOnEdge0 = aValue
        if sKey == "verticesOnCell":
            verticesOnCell0 = aValue
        if sKey == "verticesOnEdge":
            verticesOnEdge0 = aValue
        if sKey == "indexToCellID":
            indexToCellID0 = aValue
        if sKey == "indexToEdgeID":
            indexToEdgeID0 = aValue
        if sKey == "indexToVertexID":
            indexToVertexID0 = aValue
        if sKey == "lonVertex":
            lonVertex0 = aValue
        if sKey == "latVertex":
            latVertex0 = aValue
        if sKey == "areaCell":
            areaCell0 = aValue
        if sKey == "bed_elevation":
            bed_elevation0 = aValue
            iFlag_bed_elevation = 1
        if sKey == "ice_thickness":
            ice_thickness0 = aValue
            iFlag_ice_thickness = 1
        if sKey == "dcEdge":
            dcEdge0 = aValue
        if sKey == "bed_elevation_profile":
            iFlag_elevation_profile = 1
            bed_elevation_profile0 = aValue

    aLatitudeVertex = latVertex0[:] / math.pi * 180
    aLongitudeVertex = lonVertex0[:] / math.pi * 180
    aLatitudeCell = latCell0[:] / math.pi * 180
    aLongitudeCell = lonCell0[:] / math.pi * 180
    aCellsOnCell = cellsOnCell0[:]
    aEdgesOnCell = edgesOnCell0[:]
    aVertexOnCell = verticesOnCell0[:]
    aVertexOnEdge0 = verticesOnEdge0[:]
    aIndexToCellID = indexToCellID0[:]
    if iFlag_bed_elevation == 1:
        aBed_elevation = bed_elevation0[:]
    if iFlag_ice_thickness == 1:
        aIce_thickness = ice_thickness0[:]
    aCellArea = areaCell0[:]
    aDcEdge = dcEdge0[:]
    if iFlag_elevation_profile == 1:
        aBed_elevation_profile = bed_elevation_profile0[:]

    ncell = len(aIndexToCellID)
    aMpas = list()
    aMpas_dict = dict()
    lCellIndex = 0

    aLongitudeCell_180 = convert_360_to_180(aLongitudeCell)
    aLongitudeVertex_180 = convert_360_to_180(aLongitudeVertex)

    # -------------------------------------------------------------------------
    # Inner helper: add a single MPAS cell into the list
    # -------------------------------------------------------------------------
    def add_cell_into_list(
        aList, i, lCellID, dArea, dElevation_mean, dElevation_profile0, aCoords_gcs
    ):
        iFlag_success = 1
        dLongitude_center = float(aLongitudeCell_180[i])
        dLatitude_center = float(aLatitudeCell[i])
        if dLongitude_center > 180:
            print("Warning: longitude > 180")

        aCellOnCellIndex = np.array(aCellsOnCell[i, :])
        aEdgesOnCellIndex = np.array(aEdgesOnCell[i, :])
        aVertexOnCellIndex = np.array(aVertexOnCell[i, :])

        dummy0 = np.where(aVertexOnCellIndex > 0)
        aVertexIndex = aVertexOnCellIndex[dummy0]
        if len(aVertexIndex) != len(set(aVertexIndex)):
            print("Duplicates found in aVertexIndex")
            iFlag_success = 0
            return iFlag_success, aList

        dummy1 = np.where(aEdgesOnCellIndex > 0)
        aEdgeIndex = aEdgesOnCellIndex[dummy1]
        if len(aEdgeIndex) != len(set(aEdgeIndex)):
            print("Duplicates found in aEdgeIndex")
            iFlag_success = 0
            return iFlag_success, aList

        dummy2 = np.where(aCellOnCellIndex > 0)
        aNeighborIndex = (aCellOnCellIndex[dummy2]).astype(int)
        aVertexIndexOnEdge = np.array(aVertexOnEdge0[aEdgeIndex - 1, :]).astype(int)

        if len(aVertexIndex) == len(aEdgeIndex) and len(aEdgeIndex) == len(
            aVertexIndexOnEdge
        ):
            pmpas = convert_gcs_attributes_to_cell(
                4,
                dLongitude_center,
                dLatitude_center,
                aCoords_gcs,
                aVertexIndex,
                aEdgeIndex,
                aVertexIndexOnEdge,
            )
            if pmpas is None:
                print("Warning: pmpas is None")
                iFlag_success = 0
                return iFlag_success, aList

            pmpas.dArea = dArea
            pmpas.calculate_edge_length()
            pmpas.dLength_flowline = pmpas.dLength  # Default
            pmpas.lCellID = lCellID
            pmpas.dElevation_mean = dElevation_mean
            pmpas.dElevation_profile0 = dElevation_profile0

            # setup neighbor information
            pmpas.aNeighbor = aNeighborIndex
            pmpas.nNeighbor = len(aNeighborIndex)
            if pmpas.nNeighbor != pmpas.nVertex:
                # this cell is next to the ocean boundary
                pmpas.nNeighbor_land = pmpas.nNeighbor
                pmpas.nNeighbor_ocean = pmpas.nVertex - pmpas.nNeighbor
                pmpas.aNeighbor_land = aNeighborIndex
                pmpas.nNeighbor_land = len(aNeighborIndex)
            else:
                # this cell is not at the land-ocean mask coastal line
                pmpas.nNeighbor_land = pmpas.nNeighbor
                pmpas.nNeighbor_ocean = 0
                pmpas.aNeighbor_land = aNeighborIndex
                pmpas.nNeighbor_land = len(aNeighborIndex)

            aDistance = list()
            for j in range(pmpas.nNeighbor):
                lEdgeID = aEdgeIndex[j]
                lIndex = lEdgeID - 1
                dDistance = aDcEdge[lIndex]
                aDistance.append(dDistance)

            pmpas.aNeighbor_distance = aDistance
            aList.append(pmpas)
            return iFlag_success, aList
        else:
            print(
                "Warning: len(aVertexIndex) != len(aVertexIndexOnEdge)",
                "cellID:",
                lCellID,
            )
            iFlag_success = 0
            return iFlag_success, aList

    # -------------------------------------------------------------------------
    # Main loop: iterate over all cells
    # -------------------------------------------------------------------------
    if iFlag_antarctic == 1:
        # Antarctic region: latitude < -60
        for i in range(ncell):
            dLongitude_center = float(aLongitudeCell_180[i])
            dLatitude_center = float(aLatitudeCell[i])
            aVertexOnCellIndex = np.array(aVertexOnCell[i, :])
            dummy0 = np.where(aVertexOnCellIndex > 0)
            aVertexIndex = aVertexOnCellIndex[dummy0]
            aLonVertex = aLongitudeVertex_180[aVertexIndex - 1]
            aLatVertex = aLatitudeVertex[aVertexIndex - 1]
            nVertex = len(aLonVertex)
            if nVertex < 3:
                print("Warning: nVertex < 3")
                continue

            iFlag = False
            ring = ogr.Geometry(ogr.wkbLinearRing)
            aCoords_gcs = np.full((nVertex, 2), -9999.0, dtype=float)
            for j in range(nVertex):
                x1 = aLonVertex[j]
                y1 = aLatVertex[j]
                ring.AddPoint(x1, y1)
                aCoords_gcs[j, 0] = x1
                aCoords_gcs[j, 1] = y1

            ring.CloseRings()
            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            pPolygon.AddGeometry(ring)

            if dLatitude_center < -60:
                iFlag = True

            if iFlag:
                if not pPolygon.IsValid():
                    print("Warning: invalid polygon")
                    continue

                lCellID = int(aIndexToCellID[i])
                dElevation_mean = float(aBed_elevation[i]) if iFlag_bed_elevation == 1 else 0.0
                dElevation_profile0 = (
                    float(aBed_elevation_profile[i, 0]) if iFlag_elevation_profile == 1 else 0.0
                )
                dThickness_ice = float(aIce_thickness[i]) if iFlag_ice_thickness == 1 else 0.0
                dArea = float(aCellArea[i])

                if dThickness_ice > 0:
                    dElevation_mean = calculate_potentiometric(dElevation_mean, dThickness_ice)
                    dElevation_profile0 = calculate_potentiometric(
                        dElevation_profile0, dThickness_ice
                    )

                iFlag_success, aMpas = add_cell_into_list(
                    aMpas, i, lCellID, dArea, dElevation_mean, dElevation_profile0, aCoords_gcs
                )
                if iFlag_success == 1:
                    aMpas_dict[lCellID] = lCellIndex
                    lCellIndex += 1

                if iFlag_save_mesh == 1:
                    pPolygon.FlattenTo2D()
                    pFeature.SetGeometry(pPolygon)
                    pFeature.SetField("cellid", int(lCellID))
                    pFeature.SetField("longitude", dLongitude_center)
                    pFeature.SetField("latitude", dLatitude_center)
                    pFeature.SetField("area", dArea)
                    if iFlag_use_mesh_dem == 1:
                        pFeature.SetField("elevation_mean", dElevation_mean)
                        pFeature.SetField("elevation_profile0", dElevation_profile0)
                    pLayer.CreateFeature(pFeature)

    else:
        if iFlag_arctic == 1:
            # Arctic region: latitude > 55
            for i in range(ncell):
                dLongitude_center = float(aLongitudeCell_180[i])
                dLatitude_center = float(aLatitudeCell[i])
                aVertexOnCellIndex = np.array(aVertexOnCell[i, :])
                dummy0 = np.where(aVertexOnCellIndex > 0)
                aVertexIndex = aVertexOnCellIndex[dummy0]
                aLonVertex = aLongitudeVertex_180[aVertexIndex - 1]
                aLatVertex = aLatitudeVertex[aVertexIndex - 1]
                nVertex = len(aLonVertex)
                if nVertex < 3:
                    print("Warning: nVertex < 3")
                    continue

                iFlag = False
                ring = ogr.Geometry(ogr.wkbLinearRing)
                aCoords_gcs = np.full((nVertex, 2), -9999.0, dtype=float)
                for j in range(nVertex):
                    x1 = aLonVertex[j]
                    y1 = aLatVertex[j]
                    ring.AddPoint(x1, y1)
                    aCoords_gcs[j, 0] = x1
                    aCoords_gcs[j, 1] = y1

                ring.CloseRings()
                pPolygon = ogr.Geometry(ogr.wkbPolygon)
                pPolygon.AddGeometry(ring)

                if dLatitude_center > 55.0:
                    iFlag = True

                if iFlag:
                    if not pPolygon.IsValid():
                        print("Warning: invalid polygon")
                        continue

                    lCellID = int(aIndexToCellID[i])
                    dElevation_mean = (
                        float(aBed_elevation[i]) if iFlag_bed_elevation == 1 else 0.0
                    )
                    dElevation_profile0 = (
                        float(aBed_elevation_profile[i, 0])
                        if iFlag_elevation_profile == 1
                        else 0.0
                    )
                    dThickness_ice = (
                        float(aIce_thickness[i]) if iFlag_ice_thickness == 1 else 0.0
                    )
                    dArea = float(aCellArea[i])

                    if dThickness_ice > 0:
                        dElevation_mean = calculate_potentiometric(
                            dElevation_mean, dThickness_ice
                        )
                        dElevation_profile0 = calculate_potentiometric(
                            dElevation_profile0, dThickness_ice
                        )

                    iFlag_success, aMpas = add_cell_into_list(
                        aMpas,
                        i,
                        lCellID,
                        dArea,
                        dElevation_mean,
                        dElevation_profile0,
                        aCoords_gcs,
                    )
                    if iFlag_success == 1:
                        aMpas_dict[lCellID] = lCellIndex
                        lCellIndex += 1

                    if iFlag_save_mesh == 1:
                        pPolygon.FlattenTo2D()
                        pFeature.SetGeometry(pPolygon)
                        pFeature.SetField("cellid", int(lCellID))
                        pFeature.SetField("longitude", dLongitude_center)
                        pFeature.SetField("latitude", dLatitude_center)
                        pFeature.SetField("area", dArea)
                        if iFlag_use_mesh_dem == 1:
                            pFeature.SetField("elevation_mean", dElevation_mean)
                            pFeature.SetField("elevation_profile0", dElevation_profile0)
                        pLayer.CreateFeature(pFeature)

        else:
            # General case: global or regional
            for i in range(ncell):
                dLongitude_center = float(aLongitudeCell_180[i])
                dLatitude_center = float(aLatitudeCell[i])

                aVertexOnCellIndex = np.array(aVertexOnCell[i, :])
                dummy0 = np.where(aVertexOnCellIndex > 0)
                aVertexIndex = aVertexOnCellIndex[dummy0]
                aLonVertex = aLongitudeVertex_180[aVertexIndex - 1]
                aLatVertex = aLatitudeVertex[aVertexIndex - 1]
                nVertex = len(aLonVertex)
                if nVertex < 3:
                    print("Warning: nVertex < 3")
                    continue

                ring = ogr.Geometry(ogr.wkbLinearRing)
                aCoords_gcs = np.full((nVertex, 2), -9999.0, dtype=float)
                for j in range(nVertex):
                    x1 = aLonVertex[j]
                    y1 = aLatVertex[j]
                    ring.AddPoint(x1, y1)
                    aCoords_gcs[j, 0] = x1
                    aCoords_gcs[j, 1] = y1

                # close the ring
                x1 = aLonVertex[0]
                y1 = aLatVertex[0]
                ring.AddPoint(x1, y1)
                pPolygon = ogr.Geometry(ogr.wkbPolygon)
                pPolygon.AddGeometry(ring)

                dLon_min = np.min(aCoords_gcs[:, 0])
                dLon_max = np.max(aCoords_gcs[:, 0])
                iFlag = False

                if iFlag_global == 1:
                    if dLatitude_center >= -60:
                        if dLatitude_center > 85:  # exclude north pole
                            continue
                        else:
                            if np.abs(dLon_min - dLon_max) > 100:
                                # polygon crosses international date line
                                aCoords_gcs_split = (
                                    split_international_date_line_polygon_coordinates(
                                        aCoords_gcs
                                    )
                                )
                                pMultiPolygon = ogr.Geometry(ogr.wkbMultiPolygon)
                                for aCoords_gcs_part in aCoords_gcs_split:
                                    nVertex_part = aCoords_gcs_part.shape[0]
                                    ring_part = ogr.Geometry(ogr.wkbLinearRing)
                                    for j in range(nVertex_part):
                                        x1 = aCoords_gcs_part[j, 0]
                                        y1 = aCoords_gcs_part[j, 1]
                                        ring_part.AddPoint(x1, y1)
                                    ring_part.CloseRings()
                                    pPolygon_part = ogr.Geometry(ogr.wkbPolygon)
                                    pPolygon_part.AddGeometry(ring_part)
                                    pMultiPolygon.AddGeometry(pPolygon_part)

                                if not pMultiPolygon.IsValid():
                                    print("Warning: invalid multipolygon")
                                    continue
                                else:
                                    pPolygon = pMultiPolygon
                            else:
                                iFlag = True
                    else:
                        continue  # remove Antarctic from global mesh for now
                else:
                    if np.abs(dLon_min - dLon_max) > 100:
                        # polygon crosses international date line – skip for regional
                        pass
                    else:
                        if not pPolygon.IsValid():
                            print("Warning: invalid polygon")
                            continue
                        else:
                            if pPolygon.Within(pBoundary):
                                iFlag = True
                            elif pPolygon.Intersects(pBoundary):
                                iFlag = True

                if iFlag:
                    lCellID = int(aIndexToCellID[i])
                    dElevation_mean = (
                        float(aBed_elevation[i]) if iFlag_bed_elevation == 1 else 0.0
                    )
                    dElevation_profile0 = (
                        float(aBed_elevation_profile[i, 0])
                        if iFlag_elevation_profile == 1
                        else 0.0
                    )
                    dThickness_ice = (
                        float(aIce_thickness[i]) if iFlag_ice_thickness == 1 else 0.0
                    )

                    if dThickness_ice > 0:
                        dElevation_mean = calculate_potentiometric(
                            dElevation_mean, dThickness_ice
                        )
                        if iFlag_elevation_profile == 1:
                            dElevation_profile0 = calculate_potentiometric(
                                dElevation_profile0, dThickness_ice
                            )

                    dArea = float(aCellArea[i])

                    iFlag_success, aMpas = add_cell_into_list(
                        aMpas,
                        i,
                        lCellID,
                        dArea,
                        dElevation_mean,
                        dElevation_profile0,
                        aCoords_gcs,
                    )
                    if iFlag_success == 1:
                        aMpas_dict[lCellID] = lCellIndex
                        lCellIndex += 1

                    if iFlag_save_mesh == 1:
                        pPolygon.FlattenTo2D()
                        dLongitude_center = float(aLongitudeCell_180[i])
                        dLatitude_center = float(aLatitudeCell[i])
                        pFeature.SetGeometry(pPolygon)
                        pFeature.SetField("cellid", int(lCellID))
                        pFeature.SetField("longitude", dLongitude_center)
                        pFeature.SetField("latitude", dLatitude_center)
                        pFeature.SetField("area", dArea)
                        if iFlag_use_mesh_dem == 1:
                            pFeature.SetField("elevation_mean", dElevation_mean)
                            pFeature.SetField("elevation_profile0", dElevation_profile0)
                        pLayer.CreateFeature(pFeature)

    # close the dataset
    if iFlag_save_mesh == 1:
        pFeature = None
        pLayer = None
        pDataset = None

    # For regional meshes, filter out cells whose neighbors are not in the domain
    if iFlag_global == 1:
        aMpas_out = aMpas
    else:
        aMpas_out = list()
        for pCell in aMpas:
            aNeighbor = pCell.aNeighbor
            aNeighbor_land_update = list()
            for lNeighbor in aNeighbor:
                if lNeighbor in aMpas_dict:
                    aNeighbor_land_update.append(lNeighbor)
            pCell.aNeighbor = aNeighbor_land_update
            pCell.nNeighbor = len(aNeighbor_land_update)
            pCell.aNeighbor_land = aNeighbor_land_update
            pCell.nNeighbor_land = len(aNeighbor_land_update)
            pCell.nNeighbor_ocean = pCell.nVertex - pCell.nNeighbor_land
            aMpas_out.append(pCell)

    return aMpas_out
