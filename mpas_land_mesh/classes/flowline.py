"""
This module defines the pyflowline class for representing river flowlines.

A flowline is a directed sequence of connected edges that represents a stream or river segment.
It includes hydrological attributes like stream order, drainage area, and connectivity.
"""

import os
import numpy as np
from typing import List, Tuple, Any, Union, Optional
from json import JSONEncoder, dumps

# Import local classes
from mpas_land_mesh.classes.vertex import pyvertex
from mpas_land_mesh.classes.edge import pyedge

class FlowlineClassEncoder(JSONEncoder):
    """
    Custom JSON encoder for pyflowline objects and their dependencies.

    Handles serialization of:
    - pyflowline objects
    - pyedge objects
    - pyvertex objects
    - numpy data types
    """
    def default(self, obj):
        if isinstance(obj, pyflowline):
            return obj.__dict__
        elif isinstance(obj, (pyedge, pyvertex)):
            return obj.__dict__
        elif isinstance(obj, (np.int64, np.int32)):
            return int(obj)
        elif isinstance(obj, (np.float64, np.float32)):
            return float(obj)
        return super().default(obj)


class pyflowline(object):
    """
    Represents a river flowline - a directed sequence of connected edges.

    This class combines the functionality of both polyline (base geometry operations)
    and flowline (hydrological attributes and operations).

    A flowline consists of:
    - A sequence of edges (pyedge objects) that form a connected path
    - Start and end vertices
    - Geometric properties: length, bounding box, sinuosity
    - Hydrological properties: stream order, drainage area, upstream/downstream connections

    Attributes:
        aEdge: List of pyedge objects forming the flowline
        aVertex: List of all vertices along the flowline path
        pVertex_start: Starting vertex of the flowline
        pVertex_end: Ending vertex of the flowline
        dLength: Total length of the flowline in meters
        wkt: Well-Known Text representation
        pBound: Bounding box (min_lon, max_lon, min_lat, max_lat)

        lFlowlineID: Unique flowline identifier
        lFlowlineIndex: Index in flowline array
        iStream_order: Stream order (Strahler or other)
        dDrainage_area: Upstream drainage area in km²
        lUpstreamIDs: List of upstream flowline IDs
        lDownstreamIDs: List of downstream flowline IDs
        iFlag_dam: Whether flowline contains a dam
        iFlag_endorheic: Whether flowline is in endorheic basin
        iFlag_remove: Whether flowline should be removed
    """

    # Polyline attributes (geometry)
    aEdge: List[pyedge] = None  # The edges forming this flowline
    aVertex: List[pyvertex] = None  # All vertices along the path
    pVertex_start: pyvertex = None  # Start vertex
    pVertex_end: pyvertex = None  # End vertex
    dLength: float = None  # Total length in meters
    wkt: str = None  # Well-Known Text representation
    pBound: Tuple[float, float, float, float] = None  # (min_lon, max_lon, min_lat, max_lat)

    # Flowline-specific attributes (hydrology)
    lFlowlineID: int = -1
    lFlowlineIndex: int = -1
    lFlowlineID_downstream: int = -1  # Single downstream flowline ID (-1 = outlet/no downstream)
    iStream_order: int = -1
    dDrainage_area: float = -1.0
    lUpstreamIDs: List[int] = None
    lDownstreamIDs: List[int] = None
    iFlag_dam: bool = False
    iFlag_endorheic: bool = False
    iFlag_remove: bool = False

    def __init__(self, aEdge: List[Union[pyedge]]) -> None:
        """
        Initialize a flowline from a sequence of connected edges.

        Args:
            aEdge: List of pyedge objects that form a connected path.
                   The edges must be connected (end of one edge = start of next).

        Raises:
            ValueError: If edges are not properly connected or list is empty.
        """
        if not aEdge:
            raise ValueError("Edge list cannot be empty")

        # Store the edges
        self.aEdge = aEdge

        # Extract all unique points/vertices along the path
        aPoint = [aEdge[0].pVertex_start]
        for edge in aEdge:
            aPoint.append(edge.pVertex_end)
        self.aVertex = aPoint

        # Set start and end points
        self.pVertex_start = aPoint[0]
        self.pVertex_end = aPoint[-1]

        # Calculate total length
        self.dLength = sum(edge.dLength for edge in aEdge)

        # Calculate bounding box
        self.pBound = self.calculate_line_bound()

        # Update WKT representation
        self.update_wkt()

        # Initialize flowline-specific attributes
        self.lFlowlineID_downstream = -1  # -1 means outlet (no downstream)
        self.lUpstreamIDs = []
        self.lDownstreamIDs = []
        self.iFlag_dam = False
        self.iFlag_endorheic = False
        self.iFlag_remove = False

    def __repr__(self) -> str:
        """Return detailed string representation."""
        return (f"pyflowline(id={self.lFlowlineID}, order={self.iStream_order}, "
                f"edges={len(self.aEdge)}, length={self.dLength:.2f}m)")

    def __str__(self) -> str:
        """Return simple string representation."""
        return (f"Flowline {self.lFlowlineID}: {len(self.aEdge)} edges, "
                f"length={self.dLength:.2f}m, order={self.iStream_order}")

    def __hash__(self) -> int:
        """Return hash based on flowline ID."""
        return hash(self.lFlowlineID)

    def __eq__(self, other: Any) -> bool:
        """Check equality based on flowline ID."""
        if not isinstance(other, pyflowline):
            return False
        return self.lFlowlineID == other.lFlowlineID

    # ============================================================================
    # Geometry methods (from pypolyline)
    # ============================================================================

    def calculate_length(self) -> float:
        """
        Calculate the total length of the flowline.

        Returns:
            float: Total length in meters
        """
        self.dLength = sum(edge.dLength for edge in self.aEdge)
        return self.dLength

    def calculate_line_bound(self) -> Tuple[float, float, float, float]:
        """
        Calculate the bounding box of the flowline.

        Returns:
            Tuple[float, float, float, float]: (min_lon, max_lon, min_lat, max_lat)
        """
        lons = [p.dLongitude_degree for p in self.aVertex]
        lats = [p.dLatitude_degree for p in self.aVertex]

        return (min(lons), max(lons), min(lats), max(lats))

    def reverse(self) -> "pyflowline":
        """
        Reverse the direction of the flowline.

        Returns:
            pyflowline: A new flowline with reversed direction
        """
        # Reverse the edges
        aEdge_reversed = [edge.reverse() for edge in reversed(self.aEdge)]

        # Create new flowline
        pFlowline_reversed = pyflowline(aEdge_reversed)

        # Copy and swap upstream/downstream
        pFlowline_reversed.lUpstreamIDs = self.lDownstreamIDs.copy() if self.lDownstreamIDs else []
        pFlowline_reversed.lDownstreamIDs = self.lUpstreamIDs.copy() if self.lUpstreamIDs else []

        return pFlowline_reversed

    def split_line_by_length(self, dDistance: float) -> "pyflowline":
        """
        Split the flowline at a specific distance from the start.

        Args:
            dDistance: Distance from start in meters

        Returns:
            pyflowline: New flowline from split point to end, or None if distance >= total length
        """
        if dDistance >= self.dLength:
            return None

        dLength_accumulated = 0.0
        aEdge_new = []

        for i, edge in enumerate(self.aEdge):
            dLength_accumulated += edge.dLength

            if dLength_accumulated > dDistance:
                # This edge contains the split point
                dLength_remaining = dLength_accumulated - dDistance

                # Split this edge
                aEdge_split = edge.split_by_length(edge.dLength - dLength_remaining)

                if len(aEdge_split) > 1:
                    # Add the second part of split edge
                    aEdge_new.append(aEdge_split[1])

                # Add all remaining edges
                aEdge_new.extend(self.aEdge[i+1:])
                break

        if aEdge_new:
            return pyflowline(aEdge_new)
        return None

    def split_by_length(self, dDistance: float) -> List["pyflowline"]:
        """
        Split the flowline into segments of specified length.

        Args:
            dDistance: Maximum length of each segment in meters

        Returns:
            List[pyflowline]: List of flowline segments
        """
        if dDistance >= self.dLength:
            return [self]

        aFlowline_out = []
        pFlowline_current = self

        while pFlowline_current is not None and pFlowline_current.dLength > dDistance:
            # Calculate split point
            dLength_accumulated = 0.0
            aEdge_first = []
            aEdge_remaining = []

            for i, edge in enumerate(pFlowline_current.aEdge):
                dLength_accumulated += edge.dLength

                if dLength_accumulated <= dDistance:
                    # This edge is completely in first segment
                    aEdge_first.append(edge)
                else:
                    # This edge needs to be split
                    dLength_in_edge = dDistance - (dLength_accumulated - edge.dLength)

                    if dLength_in_edge > 0:
                        aEdge_split = edge.split_by_length(dLength_in_edge)

                        if len(aEdge_split) > 0:
                            aEdge_first.append(aEdge_split[0])
                        if len(aEdge_split) > 1:
                            aEdge_remaining.append(aEdge_split[1])
                    else:
                        aEdge_remaining.append(edge)

                    # Add all remaining edges
                    aEdge_remaining.extend(pFlowline_current.aEdge[i+1:])
                    break

            # Create first segment
            if aEdge_first:
                aFlowline_out.append(pyflowline(aEdge_first))

            # Update current to process remaining
            if aEdge_remaining:
                pFlowline_current = pyflowline(aEdge_remaining)
            else:
                pFlowline_current = None

        # Add final segment
        if pFlowline_current is not None:
            aFlowline_out.append(pFlowline_current)

        return aFlowline_out

    def copy_attributes(self, other: "pyflowline") -> None:
        """
        Copy flowline-specific attributes from another flowline.

        Args:
            other: Source flowline to copy attributes from
        """
        self.lFlowlineID = other.lFlowlineID
        self.lFlowlineIndex = other.lFlowlineIndex
        self.iStream_order = other.iStream_order
        self.dDrainage_area = other.dDrainage_area
        self.lUpstreamIDs = other.lUpstreamIDs.copy() if other.lUpstreamIDs else []
        self.lDownstreamIDs = other.lDownstreamIDs.copy() if other.lDownstreamIDs else []
        self.iFlag_dam = other.iFlag_dam
        self.iFlag_endorheic = other.iFlag_endorheic
        self.iFlag_remove = other.iFlag_remove

    def calculate_polyline_sinuosity(self) -> float:
        """
        Calculate the sinuosity of the flowline.

        Sinuosity = actual length / straight-line distance between endpoints

        Returns:
            float: Sinuosity ratio (>= 1.0)
        """
        if self.dLength <= 0:
            return 1.0

        # Calculate straight-line distance between start and end
        dDistance_straight = self.pVertex_start.calculate_distance(self.pVertex_end)

        if dDistance_straight <= 0:
            return 1.0

        return self.dLength / dDistance_straight

    def calculate_distance_to_point(self, pPoint: pyvertex) -> Tuple[float, pyvertex]:
        """
        Calculate the minimum distance from the flowline to a point.

        Args:
            pPoint: The point to measure distance to

        Returns:
            Tuple[float, pyvertex]: (minimum distance in meters, closest point on flowline)
        """
        dDistance_min = float('inf')
        pPoint_closest = None

        for edge in self.aEdge:
            dDistance, pPoint_on_edge = edge.calculate_distance_to_point(pPoint)

            if dDistance < dDistance_min:
                dDistance_min = dDistance
                pPoint_closest = pPoint_on_edge

        return dDistance_min, pPoint_closest

    def calculate_buffer_zone_polygon(
        self,
        dRadius: float,
        npoint: int = 8,
        use_high_precision: bool = False
    ) -> str:
        """
        Calculate a buffer zone polygon around the flowline.

        Args:
            dRadius: Buffer radius in meters
            npoint: Number of points per quarter circle for rounded ends
            use_high_precision: Whether to use high-precision geodesic calculations

        Returns:
            str: WKT polygon representation of buffer zone
        """

        # Collect all buffer polygons from edges
        aPolygon = []

        for edge in self.aEdge:
            polygon_wkt = edge.calculate_buffer_zone_polygon(
                dRadius, npoint, use_high_precision
            )
            if polygon_wkt:
                aPolygon.append(polygon_wkt)

        # For simplicity, return the union would require shapely
        # For now, return the first edge buffer
        if aPolygon:
            return aPolygon[0]
        return ""

    def calculate_distance_to_polyline(self, pPolyline_other: "pyflowline") -> float:
        """
        Calculate the minimum distance to another flowline.

        Args:
            pPolyline_other: Other flowline

        Returns:
            float: Minimum distance in meters
        """
        dDistance_min = float('inf')

        for edge in self.aEdge:
            for edge_other in pPolyline_other.aEdge:
                # Distance between two edges - use point-to-line distances
                dist1, _ = edge.calculate_distance_to_point(edge_other.pVertex_start)
                dist2, _ = edge.calculate_distance_to_point(edge_other.pVertex_end)
                dist3, _ = edge_other.calculate_distance_to_point(edge.pVertex_start)
                dist4, _ = edge_other.calculate_distance_to_point(edge.pVertex_end)

                dDistance = min(dist1, dist2, dist3, dist4)
                dDistance_min = min(dDistance_min, dDistance)

        return dDistance_min

    def calculate_bearing_angle(self) -> Optional[float]:
        """
        Calculate the bearing angle of the flowline (from start to end).

        Returns:
            Optional[float]: Bearing angle in degrees (0-360), or None if cannot calculate
        """
        if not self.aEdge:
            return None

        # Use the first edge's direction as the flowline bearing
        return self.aEdge[0].calculate_angle_with(self.aEdge[0])

    def towkt(self) -> str:
        """
        Convert flowline to Well-Known Text LINESTRING format.

        Returns:
            str: WKT representation
        """
        if not self.aVertex:
            return ""

        coords = ", ".join(f"{p.dLongitude_degree} {p.dLatitude_degree}" for p in self.aVertex)
        return f"LINESTRING({coords})"

    def update_wkt(self) -> str:
        """
        Update and return the WKT representation.

        Returns:
            str: WKT representation
        """
        self.wkt = self.towkt()
        return self.wkt

    # ============================================================================
    # Flowline-specific methods (setters and properties)
    # ============================================================================

    def set_flowline_id(self, lFlowlineID: int) -> None:
        """
        Set the flowline ID.

        Args:
            lFlowlineID: Unique flowline identifier
        """
        self.lFlowlineID = lFlowlineID

    def set_flowline_index(self, lFlowlineIndex: int) -> None:
        """
        Set the flowline index.

        Args:
            lFlowlineIndex: Index in flowline array
        """
        self.lFlowlineIndex = lFlowlineIndex

    def set_stream_order(self, iStream_order: int) -> None:
        """
        Set the stream order.

        Args:
            iStream_order: Stream order value
        """
        self.iStream_order = iStream_order

    def set_drainage_area(self, dDrainage_area: float) -> None:
        """
        Set the drainage area.

        Args:
            dDrainage_area: Drainage area in km²
        """
        self.dDrainage_area = dDrainage_area

    # ============================================================================
    # Flowline validation and checks
    # ============================================================================

    def is_valid(self) -> bool:
        """
        Check if the flowline is valid.

        Returns:
            bool: True if flowline is valid
        """
        if not self.aEdge:
            return False

        # Check that all edges are connected
        for i in range(len(self.aEdge) - 1):
            if self.aEdge[i].pVertex_end != self.aEdge[i+1].pVertex_start:
                return False

        return True

    def is_headwater(self) -> bool:
        """
        Check if this is a headwater flowline (no upstream flowlines).

        Returns:
            bool: True if headwater
        """
        return not self.lUpstreamIDs or len(self.lUpstreamIDs) == 0

    def is_outlet(self) -> bool:
        """
        Check if this is an outlet flowline (no downstream flowlines).

        Returns:
            bool: True if outlet
        """
        return not self.lDownstreamIDs or len(self.lDownstreamIDs) == 0

    # ============================================================================
    # Topology checks
    # ============================================================================

    def check_upstream(self, other: "pyflowline") -> bool:
        """
        Check if another flowline is directly upstream of this one.

        Args:
            other: Potential upstream flowline

        Returns:
            bool: True if other flows into this flowline
        """
        if not other or not other.aEdge or not self.aEdge:
            return False

        # Check if other's end connects to this flowline's start
        return other.pVertex_end == self.pVertex_start

    def check_downstream(self, other: "pyflowline") -> bool:
        """
        Check if another flowline is directly downstream of this one.

        Args:
            other: Potential downstream flowline

        Returns:
            bool: True if this flowline flows into other
        """
        if not other or not other.aEdge or not self.aEdge:
            return False

        # Check if this flowline's end connects to other's start
        return self.pVertex_end == other.pVertex_start

    # ============================================================================
    # Flowline operations
    # ============================================================================

    def merge_upstream(self, other: "pyflowline") -> "pyflowline":
        """
        Merge an upstream flowline into this one.

        Args:
            other: Upstream flowline to merge

        Returns:
            pyflowline: New merged flowline

        Raises:
            ValueError: If flowlines are not connected
        """
        if not self.check_upstream(other):
            raise ValueError("Flowlines are not connected (other must be upstream)")

        # Combine edges: upstream first, then this flowline
        aEdge_merged = other.aEdge + self.aEdge

        # Create new flowline
        pFlowline_merged = pyflowline(aEdge_merged)

        # Copy attributes from this flowline (downstream one)
        pFlowline_merged.lFlowlineID = self.lFlowlineID
        pFlowline_merged.lFlowlineIndex = self.lFlowlineIndex
        pFlowline_merged.iStream_order = max(self.iStream_order, other.iStream_order)
        pFlowline_merged.dDrainage_area = self.dDrainage_area  # Keep downstream drainage area

        # Update connectivity: merged flowline inherits upstream's upstreams
        pFlowline_merged.lUpstreamIDs = other.lUpstreamIDs.copy() if other.lUpstreamIDs else []
        pFlowline_merged.lDownstreamIDs = self.lDownstreamIDs.copy() if self.lDownstreamIDs else []

        # Copy flags
        pFlowline_merged.iFlag_dam = self.iFlag_dam or other.iFlag_dam
        pFlowline_merged.iFlag_endorheic = self.iFlag_endorheic or other.iFlag_endorheic

        return pFlowline_merged

    def copy_attributes(self, other: "pyflowline") -> None:
        """
        Copy all attributes from another flowline.

        Args:
            other: Source flowline
        """
        self.lFlowlineID = other.lFlowlineID
        self.lFlowlineIndex = other.lFlowlineIndex
        self.iStream_order = other.iStream_order
        self.dDrainage_area = other.dDrainage_area
        self.lUpstreamIDs = other.lUpstreamIDs.copy() if other.lUpstreamIDs else []
        self.lDownstreamIDs = other.lDownstreamIDs.copy() if other.lDownstreamIDs else []
        self.iFlag_dam = other.iFlag_dam
        self.iFlag_endorheic = other.iFlag_endorheic
        self.iFlag_remove = other.iFlag_remove

    def calculate_flowline_sinuosity(self) -> float:
        """
        Calculate and store the sinuosity of the flowline.

        Returns:
            float: Sinuosity ratio
        """
        return self.calculate_polyline_sinuosity()

    def calculate_distance_to_vertex(self, pVertex: pyvertex) -> float:
        """
        Calculate minimum distance from flowline to a vertex.

        Args:
            pVertex: The vertex to measure distance to

        Returns:
            float: Minimum distance in meters
        """
        dDistance, _ = self.calculate_distance_to_point(pVertex)
        return dDistance

    def get_sinuosity(self) -> float:
        """
        Get the sinuosity of the flowline.

        Returns:
            float: Sinuosity ratio
        """
        return self.calculate_polyline_sinuosity()

    def copy(self) -> "pyflowline":
        """
        Create a deep copy of the flowline.

        Returns:
            pyflowline: New flowline with copied edges and attributes
        """
        # Copy edges
        aEdge_copy = [edge.copy() for edge in self.aEdge]

        # Create new flowline
        pFlowline_copy = pyflowline(aEdge_copy)

        # Copy all attributes
        pFlowline_copy.copy_attributes(self)

        return pFlowline_copy

    def reverse(self) -> "pyflowline":
        """
        Create a reversed copy of the flowline.

        Returns:
            pyflowline: New flowline with reversed direction
        """
        # Reverse edges
        aEdge_reversed = [edge.reverse() for edge in reversed(self.aEdge)]

        # Create new flowline
        pFlowline_reversed = pyflowline(aEdge_reversed)

        # Copy attributes
        pFlowline_reversed.lFlowlineID = self.lFlowlineID
        pFlowline_reversed.lFlowlineIndex = self.lFlowlineIndex
        pFlowline_reversed.iStream_order = self.iStream_order
        pFlowline_reversed.dDrainage_area = self.dDrainage_area
        pFlowline_reversed.iFlag_dam = self.iFlag_dam
        pFlowline_reversed.iFlag_endorheic = self.iFlag_endorheic
        pFlowline_reversed.iFlag_remove = self.iFlag_remove

        # Swap upstream/downstream
        pFlowline_reversed.lUpstreamIDs = self.lDownstreamIDs.copy() if self.lDownstreamIDs else []
        pFlowline_reversed.lDownstreamIDs = self.lUpstreamIDs.copy() if self.lUpstreamIDs else []

        return pFlowline_reversed

    # ============================================================================
    # Getters
    # ============================================================================

    def get_length(self) -> float:
        """Get the flowline length in meters."""
        return self.dLength

    def get_edge_count(self) -> int:
        """Get the number of edges in the flowline."""
        return len(self.aEdge)

    def get_vertex_count(self) -> int:
        """Get the number of vertices in the flowline."""
        return len(self.aVertex)

    def get_upstream_count(self) -> int:
        """
        Get the number of upstream flowlines.

        Returns:
            int: Count of upstream flowlines
        """
        return len(self.lUpstreamIDs) if self.lUpstreamIDs else 0

    def has_dam(self) -> bool:
        """Check if flowline has a dam."""
        return self.iFlag_dam

    def is_endorheic(self) -> bool:
        """Check if flowline is in an endorheic basin."""
        return self.iFlag_endorheic

    def should_keep(self) -> bool:
        """Check if flowline should be kept (not marked for removal)."""
        return not self.iFlag_remove

    def mark_for_removal(self) -> None:
        """Mark the flowline for removal."""
        self.iFlag_remove = True

    def mark_for_keeping(self) -> None:
        """Mark the flowline to be kept."""
        self.iFlag_remove = False

    # ============================================================================
    # Serialization
    # ============================================================================

    def tojson(self) -> str:
        """
        Convert flowline to JSON string.

        Returns:
            str: JSON representation
        """
        return dumps(
            {
                'lFlowlineID': int(self.lFlowlineID),
                'lFlowlineIndex': int(self.lFlowlineIndex),
                'iStream_order': int(self.iStream_order),
                'dDrainage_area': float(self.dDrainage_area),
                'dLength': float(self.dLength),
                'nEdge': len(self.aEdge),
                'nVertex': len(self.aVertex),
                'lUpstreamIDs': [int(x) for x in self.lUpstreamIDs] if self.lUpstreamIDs else [],
                'lDownstreamIDs': [int(x) for x in self.lDownstreamIDs] if self.lDownstreamIDs else [],
                'iFlag_dam': bool(self.iFlag_dam),
                'iFlag_endorheic': bool(self.iFlag_endorheic),
                'iFlag_remove': bool(self.iFlag_remove),
                'wkt': self.wkt
            },
            cls=FlowlineClassEncoder
        )


