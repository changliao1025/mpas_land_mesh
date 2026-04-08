"""
Standalone pympas class for representing MPAS mesh cells.

This module defines the pympas class which provides MPAS-specific mesh cell
functionality for flowline network modeling without depending on pyflowline
or pyearth packages.

MPAS (Model for Prediction Across Scales) uses unstructured meshes with
variable polygon shapes (3-10 edges per cell).
"""

import json
from json import JSONEncoder
from typing import List, Optional, Any, Tuple
import numpy as np

from mpas_land_mesh.classes.vertex import pyvertex
from mpas_land_mesh.classes.edge import pyedge
from mpas_land_mesh.utilities.geometry import calculate_polygon_area


class MpasClassEncoder(JSONEncoder):
    """
    Custom JSON encoder for pympas objects.

    Handles numpy data types, pyvertex objects, and other complex types,
    converting them to native Python types for JSON serialization.
    """

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, (np.float32, np.float64)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, list):
            return obj
        if isinstance(obj, pyvertex):
            return json.loads(obj.tojson())
        if isinstance(obj, pyedge):
            return getattr(obj, "lEdgeID", str(obj))
        if isinstance(obj, pympas):
            return obj.lCellID
        return JSONEncoder.default(self, obj)


class pympas(object):
    """
    Standalone MPAS mesh cell class for flowline network representation.

    Represents a polygonal MPAS cell with variable number of vertices and
    edges (3-10). Does not inherit from pyflowline or pyearth classes.

    Attributes:
        lCellID (int): Cell identifier
        nFlowline (int): Number of flowlines in this cell
        nVertex (int): Number of vertices/edges
        nEdge (int): Number of edges
        dLength (float): Effective cell resolution (sqrt of area)
        dArea (float): Cell area in square meters
        dLongitude_center_degree (float): Center longitude in degrees
        dLatitude_center_degree (float): Center latitude in degrees
        dElevation_mean (float): Mean bed elevation
        dElevation_profile0 (float): Elevation profile value
        dLength_flowline (float): Total flowline length in cell
        nNeighbor (int): Total number of neighbors
        nNeighbor_land (int): Number of land neighbors
        nNeighbor_ocean (int): Number of ocean neighbors
        aNeighbor (list): Global IDs of all neighbors
        aNeighbor_land (list): Global IDs of land neighbors
        aNeighbor_distance (list): Distances to neighbors
        pBound (tuple): Bounding box (min_lon, min_lat, max_lon, max_lat)
        iFlag_watershed_boundary_burned (int): Watershed boundary flag
        aEdge (list): List of pyedge objects
        aVertex (list): List of pyvertex objects
    """

    def __init__(self, dLon: float, dLat: float, aEdge: List, aVertex: List) -> None:
        """
        Initialize an MPAS mesh cell object.

        Args:
            dLon (float): The longitude of the MPAS cell center in degrees
            dLat (float): The latitude of the MPAS cell center in degrees
            aEdge (List): A list of 3-10 edges that define the MPAS cell boundary
            aVertex (List): A list of 3-10 vertices that define the MPAS cell boundary

        Raises:
            ValueError: If the cell doesn't have 3-10 edges and vertices
            TypeError: If input parameters are not of the expected types
        """
        # Input validation
        if not isinstance(dLon, (int, float, np.number)):
            raise TypeError(f"dLon must be a number, got {type(dLon)}")
        if not isinstance(dLat, (int, float, np.number)):
            raise TypeError(f"dLat must be a number, got {type(dLat)}")
        if not isinstance(aEdge, list):
            raise TypeError(f"aEdge must be a list, got {type(aEdge)}")
        if not isinstance(aVertex, list):
            raise TypeError(f"aVertex must be a list, got {type(aVertex)}")

        # Validate coordinate ranges
        if not (-180 <= dLon <= 180):
            raise ValueError(f"Longitude must be between -180 and 180, got {dLon}")
        if not (-90 <= dLat <= 90):
            raise ValueError(f"Latitude must be between -90 and 90, got {dLat}")

        # Validate MPAS geometry requirements (3-10 edges for flexibility)
        nEdge = len(aEdge)
        nVertex = len(aVertex)

        if not (3 <= nEdge <= 10):
            raise ValueError(f"MPAS cell must have 3-10 edges, got {nEdge}")
        if not (3 <= nVertex <= 10):
            raise ValueError(f"MPAS cell must have 3-10 vertices, got {nVertex}")
        if nEdge != nVertex:
            raise ValueError(
                f"Number of edges ({nEdge}) must equal number of vertices ({nVertex})"
            )

        # Validate edges and vertices are not None
        for i, edge in enumerate(aEdge):
            if edge is None:
                raise ValueError(f"Edge at index {i} cannot be None")
        for i, vertex in enumerate(aVertex):
            if vertex is None:
                raise ValueError(f"Vertex at index {i} cannot be None")

        # Store geometry
        self.aEdge: List = aEdge
        self.aVertex: List = aVertex

        # Basic cell identification
        self.lCellID: int = -1

        # Flowline and geometry counts
        self.nFlowline: int = 0
        self.nVertex: int = nVertex
        self.nEdge: int = nEdge

        # Geometric properties
        self.dLength: float = 0.0
        self.dArea: float = 0.0
        self.dLongitude_center_degree: float = float(dLon)
        self.dLatitude_center_degree: float = float(dLat)

        # Elevation properties
        self.dElevation_mean: float = -9999.0
        self.dElevation_profile0: float = 0.0

        # Flowline properties
        self.dLength_flowline: float = 0.0

        # Status flags
        self.iFlag_intersected: int = -1
        self.iFlag_coast: int = 0

        # Downstream connectivity (after burning)
        self.lCellID_downstream_burned: int = -1
        self.iStream_order_burned: int = -1
        self.iStream_segment_burned: int = -1

        # Neighbor information
        self.nNeighbor: int = -1
        self.nNeighbor_land: int = -1
        self.nNeighbor_ocean: int = -1
        self.nNeighbor_land_virtual: int = -1

        # Neighbor lists
        self.aNeighbor_land_virtual: Optional[List] = None
        self.aNeighbor: Optional[List] = None
        self.aNeighbor_land: Optional[List] = None
        self.aNeighbor_ocean: Optional[List] = None
        self.aNeighbor_distance: Optional[List] = None

        # Spatial indexing
        self.pBound: Optional[Any] = None

        # MPAS-specific attributes
        self.iFlag_watershed_boundary_burned: int = 0

        # Calculate initial properties
        self.calculate_cell_area()
        if self.dArea > 0:
            self.calculate_edge_length()

    def __repr__(self) -> str:
        return (
            f"pympas(ID={self.lCellID}, "
            f"Center=({self.dLongitude_center_degree:.6f}, {self.dLatitude_center_degree:.6f}), "
            f"Edges={self.nEdge}, Area={self.dArea:.2f}m², Resolution={self.dLength:.2f}m, "
            f"Flowlines={self.nFlowline}, Neighbors={self.nNeighbor})"
        )

    def __str__(self) -> str:
        return (
            f"pympas(ID={self.lCellID}, "
            f"Center=({self.dLongitude_center_degree:.6f}, {self.dLatitude_center_degree:.6f}), "
            f"Edges={self.nEdge}, Area={self.dArea:.2f}m²)"
        )

    def __hash__(self) -> int:
        return hash(self.lCellID)

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, pympas):
            return False
        return self.lCellID == other.lCellID

    def calculate_cell_bound(self) -> Tuple[float, float, float, float]:
        """
        Calculate the bounding box of the MPAS cell.

        Handles the International Date Line crossing case.

        Returns:
            Tuple[float, float, float, float]: (min_lon, min_lat, max_lon, max_lat)
        """
        if not self.aVertex or len(self.aVertex) == 0:
            raise ValueError("Cannot calculate bounds: no vertices defined")

        dLat_min = 90.0
        dLat_max = -90.0
        dLon_min = 180.0
        dLon_max = -180.0

        for vertex in self.aVertex:
            if hasattr(vertex, "dLongitude_degree") and hasattr(vertex, "dLatitude_degree"):
                dLon_max = max(dLon_max, vertex.dLongitude_degree)
                dLon_min = min(dLon_min, vertex.dLongitude_degree)
                dLat_max = max(dLat_max, vertex.dLatitude_degree)
                dLat_min = min(dLat_min, vertex.dLatitude_degree)
            else:
                raise ValueError(
                    "Vertex must have dLongitude_degree and dLatitude_degree attributes"
                )

        # Handle International Date Line crossing
        if dLon_max - dLon_min > 180:
            tmp = dLon_max
            dLon_max = dLon_min + 360
            dLon_min = tmp

        self.pBound = (dLon_min, dLat_min, dLon_max, dLat_max)
        return self.pBound

    def has_this_edge(self, pEdge_in) -> int:
        """
        Check whether the MPAS cell contains a specific edge.

        Args:
            pEdge_in: The edge to be checked

        Returns:
            int: 1 if found, 0 if not found
        """
        if pEdge_in is None:
            raise ValueError("Edge to check cannot be None")

        iFlag_found = 0
        for pEdge in self.aEdge:
            if hasattr(pEdge, "is_overlap") and pEdge.is_overlap(pEdge_in):
                iFlag_found = 1
                break

        return iFlag_found

    def which_edge_cross_this_vertex(self, pVertex_in) -> Tuple[int, Optional[Any]]:
        """
        Find which edge overlaps with a vertex.

        Args:
            pVertex_in: The vertex to be checked

        Returns:
            Tuple[int, Optional[Any]]: (1 if found with edge object, 0 if not found with None)
        """
        if pVertex_in is None:
            raise ValueError("Vertex to check cannot be None")

        iFlag_found = 0
        pEdge_out = None

        for pEdge in self.aEdge:
            if hasattr(pEdge, "check_vertex_on_edge"):
                iFlag, dummy, diff = pEdge.check_vertex_on_edge(pVertex_in)
                if iFlag == 1:
                    iFlag_found = 1
                    pEdge_out = pEdge
                    break

        return iFlag_found, pEdge_out

    def calculate_cell_area(self) -> float:
        """
        Calculate the area of the MPAS cell using vertex coordinates.

        Returns:
            float: The area in square meters
        """
        if not self.aVertex or len(self.aVertex) == 0:
            raise ValueError("Cannot calculate area: no vertices defined")

        lons = []
        lats = []

        for vertex in self.aVertex:
            if hasattr(vertex, "dLongitude_degree") and hasattr(vertex, "dLatitude_degree"):
                lons.append(vertex.dLongitude_degree)
                lats.append(vertex.dLatitude_degree)
            else:
                raise ValueError(
                    "Vertex must have dLongitude_degree and dLatitude_degree attributes"
                )

        if not (3 <= len(lons) <= 10):
            raise ValueError(f"MPAS cell must have 3-10 vertices, got {len(lons)}")

        self.dArea = calculate_polygon_area(lons, lats)
        return self.dArea

    def calculate_edge_length(self) -> float:
        """
        Calculate the effective resolution/length of the MPAS cell.

        For MPAS cells, the effective length is the square root of the area.

        Returns:
            float: The effective cell resolution in meters
        """
        if self.dArea is None:
            raise ValueError("Area is not defined. Calculate area first.")
        if self.dArea < 0:
            raise ValueError(f"Area cannot be negative: {self.dArea}")
        if self.dArea == 0:
            self.calculate_cell_area()

        self.dLength = np.sqrt(self.dArea)
        return self.dLength

    def share_edge(self, other: "pympas") -> int:
        """
        Check whether this MPAS cell shares an edge with another MPAS cell.

        Args:
            other (pympas): The other MPAS cell to check

        Returns:
            int: 1 if cells share an edge, 0 if not
        """
        if other is None:
            raise ValueError("Other cell cannot be None")
        if not isinstance(other, pympas):
            raise TypeError(f"Other cell must be a pympas instance, got {type(other)}")

        iFlag_share = 0
        for pEdge in self.aEdge:
            for pEdge2 in other.aEdge:
                if hasattr(pEdge, "is_overlap") and hasattr(pEdge2, "is_overlap"):
                    if pEdge.is_overlap(pEdge2) == 1:
                        iFlag_share = 1
                        break
            if iFlag_share == 1:
                break

        return iFlag_share

    def set_cell_id(self, lCellID: int) -> None:
        """Set the cell ID."""
        self.lCellID = int(lCellID)

    def set_watershed_boundary_flag(self, iFlag: int) -> None:
        """
        Set the watershed boundary flag for the MPAS cell.

        Args:
            iFlag (int): Flag value (0 or 1)
        """
        if not isinstance(iFlag, (int, np.integer)):
            raise TypeError(f"Flag must be an integer, got {type(iFlag)}")
        if iFlag not in [0, 1]:
            raise ValueError(f"Flag must be 0 or 1, got {iFlag}")
        self.iFlag_watershed_boundary_burned = int(iFlag)

    def is_watershed_boundary(self) -> bool:
        """Check if the MPAS cell is on a watershed boundary."""
        return self.iFlag_watershed_boundary_burned == 1

    def crosses_international_dateline(self) -> bool:
        """Check if the MPAS cell crosses the International Date Line."""
        if self.pBound is None:
            self.calculate_cell_bound()
        dLon_min, _, dLon_max, _ = self.pBound
        return dLon_max > 180 or (dLon_max - dLon_min > 180)

    def get_corner_coordinates(self) -> List[Tuple[float, float]]:
        """
        Get the coordinates of all corners of the MPAS cell.

        Returns:
            List[Tuple[float, float]]: List of (longitude, latitude) tuples
        """
        corners = []
        for vertex in self.aVertex:
            if hasattr(vertex, "dLongitude_degree") and hasattr(vertex, "dLatitude_degree"):
                corners.append((vertex.dLongitude_degree, vertex.dLatitude_degree))
            else:
                raise ValueError(
                    "Vertex must have dLongitude_degree and dLatitude_degree attributes"
                )
        return corners

    def is_valid(self) -> bool:
        """
        Check if the MPAS cell has valid attributes.

        Returns:
            bool: True if cell has valid MPAS attributes
        """
        has_valid_edge_count = 3 <= len(self.aEdge) <= 10
        has_valid_vertex_count = 3 <= len(self.aVertex) <= 10
        has_matching_counts = len(self.aEdge) == len(self.aVertex)
        has_valid_lon = -180 <= self.dLongitude_center_degree <= 180
        has_valid_lat = -90 <= self.dLatitude_center_degree <= 90

        vertices_valid = all(
            hasattr(v, "dLongitude_degree") and hasattr(v, "dLatitude_degree")
            for v in self.aVertex
        )

        return (
            has_valid_edge_count
            and has_valid_vertex_count
            and has_matching_counts
            and has_valid_lon
            and has_valid_lat
            and vertices_valid
        )

    def copy(self) -> "pympas":
        """Create a deep copy of the MPAS cell."""
        new_mpas = pympas(
            self.dLongitude_center_degree,
            self.dLatitude_center_degree,
            self.aEdge.copy(),
            self.aVertex.copy(),
        )
        for attr_name, attr_value in self.__dict__.items():
            if hasattr(new_mpas, attr_name):
                if isinstance(attr_value, list) and attr_value is not None:
                    setattr(new_mpas, attr_name, attr_value.copy())
                else:
                    setattr(new_mpas, attr_name, attr_value)
        return new_mpas

    def tojson(self) -> str:
        """
        Convert MPAS cell object to a JSON string.

        Returns:
            str: JSON string representation of the MPAS cell
        """
        aSkip = [
            "aLine",
            "aPoint",
            "aEdge",
            "aFlowline",
            "pPoint_center",
            "dLongitude_radian",
            "dLatitude_radian",
            "wkt",
            "pVertex_start",
            "pVertex_end",
            "pBound",
        ]

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)

        sJson = json.dumps(
            obj, sort_keys=True, indent=4, ensure_ascii=True, cls=MpasClassEncoder
        )
        return sJson
