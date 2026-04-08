"""
Unified PyEdge class for representing edges/lines in flowline networks.

This module provides a unified edge class that combines line geometry operations
with flowline-specific edge attributes and functionality.
"""

import os
import json
from json import JSONEncoder
import numpy as np
from typing import List, Tuple, Optional, Dict, Any
from osgeo import ogr

# Import local classes
from mpas_land_mesh.classes.vertex import pyvertex

# Import local utility functions
from mpas_land_mesh.utilities.geometry import (
    calculate_intersect_on_great_circle,
)
from mpas_land_mesh.utilities.object import split_line_by_length
from mpas_land_mesh.utilities.geometry import (
    find_minimal_enclosing_polygon,
)
from mpas_land_mesh.utilities.geometry import (
    calculate_angle_between_point,
)
from mpas_land_mesh.utilities.geometry import (
    calculate_distance_to_plane,
)


# Constants for geometric calculations
ANGLE_THRESHOLD_COLLINEAR = (
    178.0  # Degrees - threshold for considering points collinear
)
ANGLE_THRESHOLD_PERPENDICULAR = 90.0  # Degrees - threshold for perpendicular check
DISTANCE_TOLERANCE = 1.0  # Meters - tolerance for point-on-line check
DISTANCE_PLANE_INITIAL = 9999.0  # Meters - initial large distance value


class EdgeClassEncoder(JSONEncoder):
    """
    Custom JSON encoder for pyedge objects.

    Handles numpy data types and pyvertex objects, converting them to
    native Python types for JSON serialization.
    """

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.float32):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, list):
            return obj
        if isinstance(obj, pyvertex):
            return json.loads(obj.tojson())
        return super().default(obj)


class pyedge(object):
    """
    Unified edge/line class for flowline network representation.

    This class combines line segment geometry (distance calculations, splitting, etc.)
    with edge-specific attributes used in flowline topology and stream network analysis.
    An edge represents a directed connection between two vertices, forming the basic
    building block of flowlines in the network.

    Attributes:
        pVertex_start (pyvertex): Starting vertex of the edge.
        pVertex_end (pyvertex): Ending vertex of the edge.
        dLength (float): The geodesic length of the line in meters.
        pBound (Tuple[float, float, float, float]): Bounding box (lon_min, lat_min, lon_max, lat_max).
        lEdgeID (int): Unique identifier for the edge (default: -1)
        lEdgeIndex (int): Array index for the edge in computational operations (default: -1)
        lIndex_upstream (int): Index of the upstream edge in the network (default: -1)
        lIndex_downstream (int): Index of the downstream edge in the network (default: -1)
        lLineID (int): Optional identifier for the line (default: -1, alias for lEdgeID).
        lLineIndex (int): Optional index for the line (default: -1, alias for lEdgeIndex).

    Args:
        pVertex_start_in (pyvertex or pyvertex): The starting vertex
        pVertex_end_in (pyvertex or pyvertex): The ending vertex

    Raises:
        TypeError: If vertex parameters are not pyvertex or pyvertex objects
        ValueError: If vertices are identical (zero-length edge)

    Example:
        >>> v1 = pyvertex({'dLongitude_degree': -77.0, 'dLatitude_degree': 38.0})
        >>> v2 = pyvertex({'dLongitude_degree': -77.1, 'dLatitude_degree': 38.1})
        >>> edge = pyedge(v1, v2)
        >>> edge.set_edge_id(100)
        >>> print(edge)
        pyedge(ID=100, Length=15588.45m)
    """

    def __init__(self, pVertex_start_in, pVertex_end_in) -> None:
        """
        Initialize a pyedge object.

        Converts vertex inputs to pyvertex objects if needed and initializes
        the geometry. Sets default values for edge-specific attributes.

        Args:
            pVertex_start_in (pyvertex or pyvertex): The starting vertex
            pVertex_end_in (pyvertex or pyvertex): The ending vertex

        Raises:
            TypeError: If inputs are not pyvertex or pyvertex objects
            ValueError: If start and end vertices are identical
        """
        # Validate inputs
        if pVertex_start_in is None or pVertex_end_in is None:
            raise ValueError("Start and end vertices cannot be None")

        # Convert pVertex_start_in and pVertex_end_in to pyvertex if they are not
        if hasattr(pVertex_start_in, "__dict__") and hasattr(
            pVertex_end_in, "__dict__"
        ):
            pVertex_start_in = pyvertex(pVertex_start_in.__dict__)
            pVertex_end_in = pyvertex(pVertex_end_in.__dict__)
        else:
            raise TypeError("Start and end vertices must be pyvertex or pyvertex objects")

        # Check if vertices are identical
        if pVertex_start_in == pVertex_end_in:
            raise ValueError(
                "Start and end vertices cannot be identical (zero-length edge)"
            )

        # Store the vertex objects BEFORE calling methods that depend on them
        self.pVertex_start = pVertex_start_in
        self.pVertex_end = pVertex_end_in

        self.dLength: float = self.calculate_length()
        self.pBound: Tuple[float, float, float, float] = self.calculate_line_bound()

        # Edge-specific attributes with defaults
        # lEdgeID: Unique identifier for the edge in the network
        self.lEdgeID: int = -1

        # lEdgeIndex: Used for array indexing in computational operations
        self.lEdgeIndex: int = -1

        # lIndex_upstream: Index of the upstream edge connected to this edge
        self.lIndex_upstream: int = -1

        # lIndex_downstream: Index of the downstream edge connected to this edge
        self.lIndex_downstream: int = -1

    @classmethod
    def create(cls, pVertex_start_in: pyvertex, pVertex_end_in: pyvertex) -> "pyedge":
        """
        Factory method to create a pyedge object.

        Args:
            pVertex_start_in (pyvertex): The starting vertex.
            pVertex_end_in (pyvertex): The ending vertex.

        Returns:
            pyedge: A pyedge object.
        """
        return cls(pVertex_start_in, pVertex_end_in)

    def calculate_line_bound(self) -> Tuple[float, float, float, float]:
        """
        Calculate the bounding box of the line.

        Returns:
            Tuple[float, float, float, float]: (lon_min, lat_min, lon_max, lat_max) in degrees.
        """
        dLon_max = max(
            self.pVertex_start.dLongitude_degree, self.pVertex_end.dLongitude_degree
        )
        dLon_min = min(
            self.pVertex_start.dLongitude_degree, self.pVertex_end.dLongitude_degree
        )
        dLat_max = max(
            self.pVertex_start.dLatitude_degree, self.pVertex_end.dLatitude_degree
        )
        dLat_min = min(
            self.pVertex_start.dLatitude_degree, self.pVertex_end.dLatitude_degree
        )
        self.pBound = (dLon_min, dLat_min, dLon_max, dLat_max)
        return self.pBound

    def calculate_length(self) -> float:
        """
        Calculate the length of the line.

        Returns:
            float: The geodesic distance between start and end points in meters.
        """
        self.dLength = self.pVertex_start.calculate_distance(self.pVertex_end)
        return self.dLength

    def check_shared_point(self, other: "pyedge") -> bool:
        """
        Check whether two edges share a common point.

        Args:
            other (pyedge): Another edge to compare with.

        Returns:
            bool: True if edges share at least one point, False otherwise.
        """
        return self.pVertex_start in (
            other.pVertex_start,
            other.pVertex_end,
        ) or self.pVertex_end in (other.pVertex_start, other.pVertex_end)

    def check_upstream(self, other: "pyedge") -> bool:
        """
        Check if another edge connects upstream to this edge.

        Args:
            other (pyedge): The potential upstream edge.

        Returns:
            bool: True if other.pVertex_end equals this.pVertex_start.
        """
        return self.pVertex_start == other.pVertex_end

    def check_downstream(self, other: "pyedge") -> bool:
        """
        Check if another edge connects downstream from this edge.

        Args:
            other (pyedge): The potential downstream edge.

        Returns:
            bool: True if this.pVertex_end equals other.pVertex_start.
        """
        return self.pVertex_end == other.pVertex_start

    def split_by_length(self, dLength_in: float) -> List["pyedge"]:
        """
        Split the edge into segments of specified maximum length.

        Args:
            dLength_in (float): Maximum length for each segment in meters.

        Returns:
            List[pyedge]: List of edge segments. Returns [self] if edge is shorter than threshold.

        Raises:
            ValueError: If dLength_in is not positive.
        """
        if dLength_in <= 0:
            raise ValueError("Length threshold must be positive.")

        if self.dLength <= dLength_in:
            return [self]
        else:
            return split_line_by_length(self, dLength_in)

    def reverse(self) -> "pyedge":
        """
        Create a reversed copy of the edge.

        Returns a new edge with start and end vertices swapped,
        and upstream/downstream indices swapped.

        Returns:
            pyedge: A new edge object with reversed direction
        """
        reversed_edge = pyedge(self.pVertex_end.copy(), self.pVertex_start.copy())
        reversed_edge.lEdgeID = self.lEdgeID
        reversed_edge.lEdgeIndex = self.lEdgeIndex
        reversed_edge.lIndex_upstream = self.lIndex_downstream
        reversed_edge.lIndex_downstream = self.lIndex_upstream
        reversed_edge.lLineID = self.lLineID
        reversed_edge.lLineIndex = self.lLineIndex
        return reversed_edge

    def is_overlap(self, pEdge_in: "pyedge") -> bool:
        """
        Check if two edges represent the same segment (possibly reversed).

        Args:
            pEdge_in (pyedge): Another edge to compare with.

        Returns:
            bool: True if edges have same endpoints (in any order), False otherwise.
        """
        return (
            self.pVertex_start == pEdge_in.pVertex_start
            and self.pVertex_end == pEdge_in.pVertex_end
        ) or (
            self.pVertex_start == pEdge_in.pVertex_end
            and self.pVertex_end == pEdge_in.pVertex_start
        )

    def check_vertex_on_edge(self, pVertex_in: pyvertex) -> Tuple[bool, float, float]:
        """
        Check if a vertex lies on this edge.

        Args:
            pVertex_in (pyvertex): The vertex to check.

        Returns:
            Tuple[bool, float, float]:
                - bool: True if vertex is on the edge
                - float: Distance along edge from start (-1 if not on edge)
                - float: Perpendicular distance to line plane
        """
        dDistance = -1.0
        dDistance_plane = DISTANCE_PLANE_INITIAL
        on_edge = False

        if pVertex_in != self.pVertex_start and pVertex_in != self.pVertex_end:
            d1 = self.pVertex_start.calculate_distance(pVertex_in)
            d2 = self.pVertex_end.calculate_distance(pVertex_in)
            d3 = d1 + d2 - self.dLength
            angle3deg = calculate_angle_between_point(
                self.pVertex_start.dLongitude_degree,
                self.pVertex_start.dLatitude_degree,
                pVertex_in.dLongitude_degree,
                pVertex_in.dLatitude_degree,
                self.pVertex_end.dLongitude_degree,
                self.pVertex_end.dLatitude_degree,
            )

            dDistance_plane = calculate_distance_to_plane(
                self.pVertex_start.dLongitude_degree,
                self.pVertex_start.dLatitude_degree,
                pVertex_in.dLongitude_degree,
                pVertex_in.dLatitude_degree,
                self.pVertex_end.dLongitude_degree,
                self.pVertex_end.dLatitude_degree,
            )

            if angle3deg > ANGLE_THRESHOLD_COLLINEAR and d3 < DISTANCE_TOLERANCE:
                on_edge = True
                dDistance = d1

        return on_edge, dDistance, dDistance_plane

    def calculate_distance_to_point(self, pVertex_in: pyvertex) -> Tuple[float, pyvertex]:
        """
        Calculate the minimum distance from a point to this line segment.

        Args:
            pPoint_in (pyvertex): The point to measure distance from.

        Returns:
            Tuple[float, pyvertex]:
                - float: Minimum distance in meters
                - pyvertex: The closest point on the line to pPoint_in
        """
        d1 = self.pVertex_start.calculate_distance(pVertex_in)
        d2 = self.pVertex_end.calculate_distance(pVertex_in)

        angle3deg = calculate_angle_between_point(
            self.pVertex_start.dLongitude_degree,
            self.pVertex_start.dLatitude_degree,
            pVertex_in.dLongitude_degree,
            pVertex_in.dLatitude_degree,
            self.pVertex_end.dLongitude_degree,
            self.pVertex_end.dLatitude_degree,
        )

        if angle3deg > ANGLE_THRESHOLD_PERPENDICULAR:
            dDistance_plane = calculate_distance_to_plane(
                self.pVertex_start.dLongitude_degree,
                self.pVertex_start.dLatitude_degree,
                pVertex_in.dLongitude_degree,
                pVertex_in.dLatitude_degree,
                self.pVertex_end.dLongitude_degree,
                self.pVertex_end.dLatitude_degree,
            )

            if dDistance_plane < d1 and dDistance_plane < d2:
                dLongitude_intersect, dLatitude_intersect = (
                    calculate_intersect_on_great_circle(
                        self.pVertex_start.dLongitude_degree,
                        self.pVertex_start.dLatitude_degree,
                        pVertex_in.dLongitude_degree,
                        pVertex_in.dLatitude_degree,
                        self.pVertex_end.dLongitude_degree,
                        self.pVertex_end.dLatitude_degree,
                    )
                )

                pVertex_out = pyvertex(
                    {
                        "dLongitude_degree": dLongitude_intersect,
                        "dLatitude_degree": dLatitude_intersect,
                    }
                )
                dDistance_min = pVertex_out.calculate_distance(pVertex_in)

                if dDistance_min < d1 and dDistance_min < d2:
                    return dDistance_min, pVertex_out

        if d1 < d2:
            return d1, self.pVertex_start
        else:
            return d2, self.pVertex_end

    def calculate_angle_with(self, other: "pyedge") -> float:
        """
        Calculate the angle between this edge and another edge.

        Uses the angle calculation function from the geometry module to compute
        the angle formed by the two line segments.

        Args:
            other (pyedge): Another edge to calculate angle with.

        Returns:
            float: Angle in degrees between the two edges.

        Raises:
            TypeError: If other is not a pyedge object.

        Example:
            >>> edge1 = pyedge(p1, p2)
            >>> edge2 = pyedge(p2, p3)
            >>> angle = edge1.calculate_angle_with(edge2)
        """
        if not isinstance(other, pyedge):
            raise TypeError(f"Expected pyedge, got {type(other)}")

        return calculate_angle_between_point(
            self.pVertex_start.dLongitude_degree,
            self.pVertex_start.dLatitude_degree,
            self.pVertex_end.dLongitude_degree,
            self.pVertex_end.dLatitude_degree,
            other.pVertex_start.dLongitude_degree,
            other.pVertex_start.dLatitude_degree,
            other.pVertex_end.dLongitude_degree,
            other.pVertex_end.dLatitude_degree,
        )

    def calculate_buffer_zone_polygon(
        self,
        dRadius: float,
        nPoint: int = 36,
        sFilename_out: Optional[str] = None,
        sFolder_out: Optional[str] = None,
    ) -> Tuple[str, List[pyvertex], List[pyvertex], List[pyvertex]]:
        """
        Calculate a buffer zone polygon around this line.

        Args:
            dRadius (float): Buffer radius in meters.
            nPoint (int): Number of points to use for circular approximation (default: 36).
            sFilename_out (Optional[str]): Output file path for the buffer polygon (unused, kept for compatibility).
            sFolder_out (Optional[str]): Output folder for intermediate files (unused, kept for compatibility).

        Returns:
            Tuple containing:
                - str: WKT representation of the buffer polygon
                - List[pyvertex]: Vertices of the buffer polygon
                - List[pyvertex]: Center points used
                - List[pyvertex]: All circle points generated

        Raises:
            ValueError: If dRadius is not positive or nPoint is less than 3.
        """
        if dRadius <= 0:
            raise ValueError("Buffer radius must be positive.")
        if nPoint < 3:
            raise ValueError("Number of points must be at least 3.")

        if self.dLength < dRadius * 1.5:
            aEdge = [self]
        else:
            aEdge = self.split_by_length(dRadius)

        aPoint_out = []
        aPoint_center = []
        aPoint_circle = []

        for i, pEdge in enumerate(aEdge):
            pVertex_start = pEdge.pVertex_start
            pVertex_end = pEdge.pVertex_end

            aPoint_center.extend([pVertex_start, pVertex_end])

            _, aPoint_start_buffer = pVertex_start.calculate_buffer_zone_circle(
                dRadius, nPoint
            )
            aPoint_circle.extend(aPoint_start_buffer)

            _, aPoint_end_buffer = pVertex_end.calculate_buffer_zone_circle(
                dRadius, nPoint
            )
            aPoint_circle.extend(aPoint_end_buffer)

        aLongitude_degree = [v.dLongitude_degree for v in aPoint_circle]
        aLatitude_degree = [v.dLatitude_degree for v in aPoint_circle]

        pPolygon_out = find_minimal_enclosing_polygon(
            aLongitude_degree, aLatitude_degree
        )

        aPoint_out = [
            pyvertex({"dLongitude_degree": p[0], "dLatitude_degree": p[1]})
            for p in pPolygon_out
        ]

        # Generate WKT from polygon points
        pGeometry = ogr.Geometry(ogr.wkbPolygon)
        pRing = ogr.Geometry(ogr.wkbLinearRing)
        for p in aPoint_out:
            pRing.AddPoint(p.dLongitude_degree, p.dLatitude_degree)
        pRing.CloseRings()
        pGeometry.AddGeometry(pRing)
        pGeometry.FlattenTo2D()
        sWkt_buffer_polygon = pGeometry.ExportToWkt()

        return sWkt_buffer_polygon, aPoint_out, aPoint_center, aPoint_circle

    def __eq__(self, other: object) -> bool:
        """
        Check if two edges are equivalent.

        Args:
            other: Object to compare with.

        Returns:
            bool: True if edges have identical start and end points in the same order.
        """
        if not isinstance(other, pyedge):
            return NotImplemented
        return (
            self.pVertex_start == other.pVertex_start
            and self.pVertex_end == other.pVertex_end
        )

    def __ne__(self, other: object) -> bool:
        """
        Check if two edges are not equivalent.

        Args:
            other: Object to compare with.

        Returns:
            bool: True if edges are not equivalent.
        """
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def __hash__(self) -> int:
        """
        Generate hash for the edge based on its endpoints.

        Returns:
            int: Hash value for the edge.
        """
        return hash((self.pVertex_start, self.pVertex_end))

    def __repr__(self) -> str:
        """
        Return a detailed string representation of the edge.

        Returns:
            str: Detailed representation including edge ID, index, and length
        """
        return (
            f"pyedge(ID={self.lEdgeID}, Index={self.lEdgeIndex}, "
            f"Length={self.dLength:.2f}m, "
            f"Upstream={self.lIndex_upstream}, Downstream={self.lIndex_downstream})"
        )

    def __str__(self) -> str:
        """
        Return a concise string representation of the edge.

        Returns:
            str: Concise representation with ID and length
        """
        return f"pyedge(ID={self.lEdgeID}, Length={self.dLength:.2f}m)"

    def tojson(self) -> str:
        """
        Convert edge object to a JSON string.

        Serializes all edge attributes including vertex information.
        Uses the custom EdgeClassEncoder to handle numpy data types
        and pyvertex objects.

        Returns:
            str: JSON string representation of the edge

        Example:
            >>> edge = pyedge(v1, v2)
            >>> json_str = edge.tojson()
        """
        aSkip = [
            "dLongitude_radian",
            "dLatitude_radian",
            "wkt",
            "pVertex_start",
            "pVertex_end",
        ]

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)

        sJson = json.dumps(
            obj, sort_keys=True, indent=4, ensure_ascii=True, cls=EdgeClassEncoder
        )
        return sJson

    def set_edge_id(self, lEdgeID: int) -> None:
        """
        Set the edge ID.

        Args:
            lEdgeID (int): New edge ID

        Raises:
            TypeError: If lEdgeID is not an integer
        """
        if not isinstance(lEdgeID, (int, np.integer)):
            raise TypeError(f"Edge ID must be an integer, got {type(lEdgeID)}")
        self.lEdgeID = int(lEdgeID)
        self.lLineID = int(lEdgeID)  # Keep in sync

    def set_edge_index(self, lEdgeIndex: int) -> None:
        """
        Set the edge array index.

        Args:
            lEdgeIndex (int): New edge index

        Raises:
            TypeError: If lEdgeIndex is not an integer
        """
        if not isinstance(lEdgeIndex, (int, np.integer)):
            raise TypeError(f"Edge index must be an integer, got {type(lEdgeIndex)}")
        self.lEdgeIndex = int(lEdgeIndex)
        self.lLineIndex = int(lEdgeIndex)  # Keep in sync

    def set_upstream_index(self, lIndex_upstream: int) -> None:
        """
        Set the upstream edge index.

        Args:
            lIndex_upstream (int): New upstream edge index

        Raises:
            TypeError: If lIndex_upstream is not an integer
        """
        if not isinstance(lIndex_upstream, (int, np.integer)):
            raise TypeError(
                f"Upstream index must be an integer, got {type(lIndex_upstream)}"
            )
        self.lIndex_upstream = int(lIndex_upstream)

    def set_downstream_index(self, lIndex_downstream: int) -> None:
        """
        Set the downstream edge index.

        Args:
            lIndex_downstream (int): New downstream edge index

        Raises:
            TypeError: If lIndex_downstream is not an integer
        """
        if not isinstance(lIndex_downstream, (int, np.integer)):
            raise TypeError(
                f"Downstream index must be an integer, got {type(lIndex_downstream)}"
            )
        self.lIndex_downstream = int(lIndex_downstream)

    def is_valid(self) -> bool:
        """
        Check if the edge has valid attributes.

        An edge is considered valid if:
        - It has non-zero length
        - Both vertices are valid
        - It has a valid ID or index

        Returns:
            bool: True if edge has valid attributes
        """
        has_valid_length = self.dLength > 0
        has_valid_vertices = (
            hasattr(self, "pVertex_start")
            and hasattr(self, "pVertex_end")
            and self.pVertex_start.is_valid()
            and self.pVertex_end.is_valid()
        )
        has_valid_id = self.lEdgeID > 0 or self.lEdgeIndex >= 0
        return has_valid_length and has_valid_vertices and has_valid_id

    def copy(self) -> "pyedge":
        """
        Create a deep copy of the edge.

        Returns:
            pyedge: A new edge object with the same attributes
        """
        new_edge = pyedge(self.pVertex_start.copy(), self.pVertex_end.copy())
        new_edge.lEdgeID = self.lEdgeID
        new_edge.lEdgeIndex = self.lEdgeIndex
        new_edge.lIndex_upstream = self.lIndex_upstream
        new_edge.lIndex_downstream = self.lIndex_downstream
        new_edge.lLineID = self.lLineID
        new_edge.lLineIndex = self.lLineIndex
        return new_edge

    def get_midpoint(self) -> pyvertex:
        """
        Calculate and return the midpoint of the edge.

        Returns:
            pyvertex: A new vertex at the midpoint of the edge
        """
        mid_lon = (
            self.pVertex_start.dLongitude_degree + self.pVertex_end.dLongitude_degree
        ) / 2.0
        mid_lat = (
            self.pVertex_start.dLatitude_degree + self.pVertex_end.dLatitude_degree
        ) / 2.0

        param = {"dLongitude_degree": mid_lon, "dLatitude_degree": mid_lat}
        return pyvertex(param)

    def is_connected_to(self, other: "pyedge") -> bool:
        """
        Check if this edge is connected to another edge.

        Two edges are connected if they share a common vertex.

        Args:
            other (pyedge): Another edge to check connection with

        Returns:
            bool: True if edges share a common vertex
        """
        if not isinstance(other, pyedge):
            return False

        return (
            self.pVertex_end == other.pVertex_start
            or self.pVertex_end == other.pVertex_end
            or self.pVertex_start == other.pVertex_start
            or self.pVertex_start == other.pVertex_end
        )

    def is_upstream_of(self, other: "pyedge") -> bool:
        """
        Check if this edge is directly upstream of another edge.

        Args:
            other (pyedge): Another edge to check

        Returns:
            bool: True if this edge's end vertex matches other's start vertex
        """
        if not isinstance(other, pyedge):
            return False
        return self.pVertex_end == other.pVertex_start

    def is_downstream_of(self, other: "pyedge") -> bool:
        """
        Check if this edge is directly downstream of another edge.

        Args:
            other (pyedge): Another edge to check

        Returns:
            bool: True if this edge's start vertex matches other's end vertex
        """
        if not isinstance(other, pyedge):
            return False
        return self.pVertex_start == other.pVertex_end


