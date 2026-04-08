"""
Convert GCS (geographic coordinate system) attributes to MPAS mesh cell objects.

This module provides the convert_gcs_attributes_to_cell function for the MPAS
mesh type (iMesh_type == 4), adapted from pyflowline without external dependencies.
"""

import numpy as np

from mpas_land_mesh.classes.vertex import pyvertex
from mpas_land_mesh.classes.edge import pyedge
from mpas_land_mesh.classes.mpas import pympas


def convert_gcs_attributes_to_cell(
    iMesh_type_in,
    dLongitude_center_in,
    dLatitude_center_in,
    aCoordinates_gcs_in,
    aVertexID_in,
    aEdgeID_in,
    aVertexIndexOnEdge_in,
):
    """
    Convert GCS attributes to an MPAS mesh cell object.

    Only supports iMesh_type_in == 4 (MPAS). Returns None for other types.

    Args:
        iMesh_type_in (int): Mesh type identifier. Must be 4 for MPAS.
        dLongitude_center_in (float): Center longitude in degrees.
        dLatitude_center_in (float): Center latitude in degrees.
        aCoordinates_gcs_in (np.ndarray): Nx2 array of (lon, lat) vertex coordinates.
        aVertexID_in (np.ndarray): Array of vertex IDs.
        aEdgeID_in (np.ndarray): Array of edge IDs.
        aVertexIndexOnEdge_in (np.ndarray): Nx2 array of vertex indices for each edge.

    Returns:
        pympas or None: The constructed MPAS cell, or None on failure.
    """
    if iMesh_type_in != 4:
        print(f"convert_gcs_attributes_to_cell: unsupported mesh type {iMesh_type_in}")
        return None

    npoint = len(aVertexID_in)
    aVertex = []
    aEdge = []

    # Build vertex list
    for i in range(npoint):
        lon = float(aCoordinates_gcs_in[i][0])
        lat = float(aCoordinates_gcs_in[i][1])
        pVertex_dict = {
            "dLongitude_degree": lon,
            "dLatitude_degree": lat,
        }
        pVertex = pyvertex(pVertex_dict)
        pVertex.lVertexID = int(aVertexID_in[i])
        aVertex.append(pVertex)

    # Build edge list
    for j in range(npoint):
        aVertexID_dummy = aVertexIndexOnEdge_in[j, :]

        # Start vertex
        dummy_index = np.where(aVertexID_in == aVertexID_dummy[0])
        if len(dummy_index[0]) == 0:
            print("Vertex ID not found for edge start")
            return None
        pVertex1_dict = {
            "dLongitude_degree": float(aCoordinates_gcs_in[dummy_index, 0]),
            "dLatitude_degree": float(aCoordinates_gcs_in[dummy_index, 1]),
        }
        pVertex_start = pyvertex(pVertex1_dict)
        pVertex_start.lVertexID = int(aVertexID_dummy[0])

        # End vertex
        dummy_index = np.where(aVertexID_in == aVertexID_dummy[1])
        if len(dummy_index[0]) == 0:
            print("Vertex ID not found for edge end")
            return None
        pVertex2_dict = {
            "dLongitude_degree": float(aCoordinates_gcs_in[dummy_index, 0]),
            "dLatitude_degree": float(aCoordinates_gcs_in[dummy_index, 1]),
        }
        pVertex_end = pyvertex(pVertex2_dict)
        pVertex_end.lVertexID = int(aVertexID_dummy[1])

        pEdge = pyedge(pVertex_start, pVertex_end)
        pEdge.lEdgeID = int(aEdgeID_in[j])
        aEdge.append(pEdge)

    pMpas = pympas(dLongitude_center_in, dLatitude_center_in, aEdge, aVertex)
    return pMpas
