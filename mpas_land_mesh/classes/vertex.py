"""
PyVertex class for representing vertices/points in flowline networks.

This module provides a unified vertex class that combines point geometry
with flowline-specific vertex attributes and functionality.
"""

import json
from json import JSONEncoder
import numpy as np
from typing import List, Tuple, Dict, Any, Optional
from osgeo import ogr

# Import local distance calculation function
from mpas_land_mesh.utilities.geometry import (
    calculate_distance_based_on_longitude_latitude,
)

iPrecision_default = 8  # ~1 mm in latitude degrees (1e-8 deg ~= 1.11 mm)


class VertexClassEncoder(JSONEncoder):
    """
    Custom JSON encoder for pyvertex objects.

    Handles numpy data types and converts them to native Python types
    for JSON serialization.
    """

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, pyvertex):
            return json.loads(obj.tojson())
        return JSONEncoder.default(self, obj)


class pyvertex(object):
    """
    Unified vertex/point class for flowline network representation.

    This class combines point geometry (coordinates, distance calculations, etc.)
    with vertex-specific attributes used in flowline topology and stream network
    analysis. A vertex represents a point in the flowline network that may be
    shared between multiple flowlines or serve as a confluence point.

    Attributes:
        dLongitude_degree (float): Longitude in degrees
        dLatitude_degree (float): Latitude in degrees
        dLongitude_radian (float): Longitude in radians
        dLatitude_radian (float): Latitude in radians
        dX_meter (float): X coordinate in meters (ECEF)
        dY_meter (float): Y coordinate in meters (ECEF)
        dZ_meter (float): Z coordinate in meters (ECEF)
        dElevation (float): Elevation value (default: 0.0)
        wkt (str): WKT representation of the point
        lVertexIndex (int): Array index for the vertex (default: -1)
        lVertexID (int): Unique identifier for the vertex (default: 1)
        lFlowlineID (int): Associated flowline ID, primarily used for
                          intersection operations (default: -1)

    Args:
        aParameter (dict): Dictionary containing vertex parameters.
                          Required keys: 'dLongitude_degree', 'dLatitude_degree'
                          Optional keys: 'lVertexIndex', 'lVertexID', 'lFlowlineID',
                                       'x', 'y', 'z', 'dElevation'

    Example:
        >>> vertex_params = {
        ...     'dLongitude_degree': -77.0,
        ...     'dLatitude_degree': 38.0,
        ...     'lVertexID': 5
        ... }
        >>> vertex = pyvertex(vertex_params)
        >>> print(vertex)
        pyvertex(ID=5, Lon=-77.0, Lat=38.0)
    """

    def __init__(self, aParameter: Dict[str, Any]) -> None:
        """
        Initialize a vertex object.

        Args:
            aParameter (dict): A dictionary containing vertex parameters.
                               Expected keys: 'dLongitude_degree', 'dLatitude_degree'.
                               Optional keys: 'x', 'y', 'z', 'dElevation',
                                            'lVertexIndex', 'lVertexID', 'lFlowlineID'.

        Raises:
            ValueError: If required parameters are missing
            TypeError: If aParameter is not a dictionary
        """
        if not isinstance(aParameter, dict):
            raise TypeError(f"aParameter must be a dictionary, got {type(aParameter)}")

        # Initialize coordinate attributes
        self.dX_meter = float(aParameter.get("x", -9999.0))
        self.dY_meter = float(aParameter.get("y", -9999.0))
        self.dZ_meter = float(aParameter.get("z", -9999.0))
        self.dElevation = float(aParameter.get("dElevation", 0.0))

        # dLongitude and dLatitude are always required
        if (
            "dLongitude_degree" not in aParameter
            or "dLatitude_degree" not in aParameter
        ):
            raise ValueError(
                "Initialization of pyvertex failed: 'dLongitude_degree' and 'dLatitude_degree' are required."
            )

        self.dLongitude_degree = float(aParameter["dLongitude_degree"])
        self.dLatitude_degree = float(aParameter["dLatitude_degree"])

        self.dLongitude_radian = np.radians(self.dLongitude_degree)
        self.dLatitude_radian = np.radians(self.dLatitude_degree)

        # Recalculate x, y, z based on dLongitude and dLatitude if default values are used
        if (
            self.dX_meter == -9999.0
            and self.dY_meter == -9999.0
            and self.dZ_meter == -9999.0
        ):
            self.dX_meter, self.dY_meter, self.dZ_meter = self.calculate_xyz()

        # Vertex-specific attributes with defaults
        # lVertexIndex: Used for array indexing in computational operations
        self.lVertexIndex = int(aParameter.get("lVertexIndex", -1))

        # lVertexID: Unique identifier for the vertex in the network
        self.lVertexID = int(aParameter.get("lVertexID", 1))

        # lFlowlineID: Associated flowline ID (used during intersection operations)
        self.lFlowlineID = int(aParameter.get("lFlowlineID", -1))

        self.wkt = self.towkt()

        return

    def __repr__(self) -> str:
        """
        Return a detailed string representation of the vertex.

        Returns:
            str: Detailed representation including vertex ID and coordinates
        """
        return (
            f"pyvertex(ID={self.lVertexID}, Index={self.lVertexIndex}, "
            f"Lon={self.dLongitude_degree:.6f}, Lat={self.dLatitude_degree:.6f}, "
            f"FlowlineID={self.lFlowlineID})"
        )

    def __str__(self) -> str:
        """
        Return a concise string representation of the vertex.

        Returns:
            str: Concise representation with ID and coordinates
        """
        return f"pyvertex(ID={self.lVertexID}, Lon={self.dLongitude_degree:.6f}, Lat={self.dLatitude_degree:.6f})"

    def toNvector(self, use_high_precision=False):
        """
        Convert geographic coordinates to n-vector representation.

        Uses high precision (float128) for trigonometric calculations when enabled
        to maintain accuracy during spherical interpolation operations.

        Note: replicated in LatLon_NvectorEllipsoidal

        Args:
            use_high_precision (bool): Use float128 for calculations (default: False)

        Returns:
            pynvector: An n-vector object with unit length
        """
        from mpas_land_mesh.classes.nvector import pynvector

        # Use high precision if requested
        if use_high_precision:
            dtype = np.float128
            a = dtype(self.dLatitude_radian)
            b = dtype(self.dLongitude_radian)
        else:
            a = self.dLatitude_radian
            b = self.dLongitude_radian

        # Calculate n-vector components
        c = np.sin(a)
        e = np.cos(a)
        d = np.sin(b)
        f = np.cos(b)

        # Right-handed vector: x -> 0°E,0°N; y -> 90°E,0°N, z -> 90°N
        x = e * f
        y = e * d
        z = c

        point = dict()
        point["x"] = x
        point["y"] = y
        point["z"] = z
        pNvector = pynvector(point, use_high_precision=use_high_precision)
        return pNvector

    def __hash__(self) -> int:
        """
        Return hash value for the vertex.

        Uses rounded coordinates to allow vertices to be used in sets and
        as dictionary keys.

        Returns:
            int: Hash value based on coordinates
        """
        return hash(
            (
                round(self.dLongitude_degree, iPrecision_default),
                round(self.dLatitude_degree, iPrecision_default),
            )
        )

    def __eq__(self, other: Any) -> bool:
        """
        Check whether two vertices are equivalent.
        Uses numpy.isclose for robust floating-point comparison.

        Vertices are considered equal if their coordinates match within
        a precision threshold, regardless of their IDs or indices.

        Args:
            other: Another object to compare with

        Returns:
            bool: True if vertices have the same coordinates
        """
        if not isinstance(other, pyvertex):
            return NotImplemented

        dThreshold_in = 10.0 ** (-1 * iPrecision_default)
        return np.isclose(
            self.dLongitude_degree, other.dLongitude_degree, atol=dThreshold_in, rtol=0
        ) and np.isclose(
            self.dLatitude_degree, other.dLatitude_degree, atol=dThreshold_in, rtol=0
        )

    def __ne__(self, other: Any) -> bool:
        """
        Check whether two vertices are not equivalent

        Args:
            other (pyvertex): The other vertex

        Returns:
            bool: True if not equivalent, False if equivalent
        """
        return not self.__eq__(other)

    def calculate_distance(self, other: "pyvertex") -> float:
        """
        Calculate the distance between two vertices

        Args:
            other (pyvertex): The other vertex

        Returns:
            float: The great circle distance in meters
        """
        dDistance = 0.0
        lon1 = self.dLongitude_degree
        lat1 = self.dLatitude_degree
        lon2 = other.dLongitude_degree
        lat2 = other.dLatitude_degree
        dDistance = calculate_distance_based_on_longitude_latitude(
            lon1, lat1, lon2, lat2
        )
        return dDistance

    def calculate_buffer_zone_point(self, dRadius: float, dBearing: float = 90) -> str:
        """
        Calculate a single buffer point at a given radius and bearing.

        Args:
            dRadius (float): Buffer radius in meters
            dBearing (float): Bearing angle in degrees (default: 90)

        Returns:
            str: WKT representation of the buffer point
        """
        # Create a geodesic object
        from geographiclib.geodesic import Geodesic

        geod = Geodesic.WGS84  # the default is WGS84
        # Calculate the geodesic buffer
        pPoint_buffer = geod.Direct(
            self.dLatitude_degree, self.dLongitude_degree, dBearing, dRadius
        )
        # Extract the latitude and longitude of the buffer point
        # create a point object using the buffer point
        point0 = dict()
        point0["dLongitude_degree"] = pPoint_buffer["lon2"]
        point0["dLatitude_degree"] = pPoint_buffer["lat2"]
        pPoint_out = pyvertex(point0)

        # convert the point to a wkt string
        sWkt_buffer_point = pPoint_out.towkt()

        return sWkt_buffer_point

    def calculate_buffer_zone_circle(
        self,
        dRadius: float,
        nPoint: int = 360,
        iFlag_support_antimeridian_in: int = 0,
        sFilename_out: Optional[str] = None,
    ) -> Tuple[str, List["pyvertex"]]:
        """
        Calculate a buffer zone circle around a point using geodesic distances.

        Args:
            dRadius (float): Buffer radius in meters.
            nPoint (int): Number of points to approximate the circle (default: 360).
            iFlag_support_antimeridian_in (int): Flag to support antimeridian crossing (default: 0).
            sFilename_out (Optional[str]): Output file path for the buffer polygon.

        Returns:
            Tuple containing:
                - str: WKT representation of the buffer polygon
                - List[pyvertex]: Points around the circle boundary
        """
        if nPoint < 3:
            raise ValueError("nPoint must be at least 3 to form a polygon.")

        # Create a geodesic object
        from geographiclib.geodesic import Geodesic

        geod = Geodesic.WGS84  # the default is WGS84
        aPoint = []

        # Calculate the geodesic buffer
        for dBearing in np.linspace(0.0, 360.0, num=nPoint, endpoint=False):
            pPoint_buffer = geod.Direct(
                self.dLatitude_degree, self.dLongitude_degree, dBearing, dRadius
            )
            point0 = dict()
            point0["dLongitude_degree"] = pPoint_buffer["lon2"]
            point0["dLatitude_degree"] = pPoint_buffer["lat2"]
            pPoint_out = pyvertex(point0)
            aPoint.append(pPoint_out)

        if (
            np.max([p.dLongitude_degree for p in aPoint])
            - np.min([p.dLongitude_degree for p in aPoint])
            > 180
        ):
            iFlag_cross_idl = 1
        else:
            iFlag_cross_idl = 0

        if iFlag_support_antimeridian_in == 1:
            # split int two parts if cross the antimeridian
            aPoint_part1 = [p for p in aPoint if p.dLongitude_degree >= 0]
            aPoint_part2 = [p for p in aPoint if p.dLongitude_degree < 0]
            # create a multi polygon wkt string
            pGeometry = ogr.Geometry(ogr.wkbMultiPolygon)
            if len(aPoint_part1) >= 3:
                pPolygon1 = ogr.Geometry(ogr.wkbPolygon)
                pRing1 = ogr.Geometry(ogr.wkbLinearRing)
                for p in aPoint_part1:
                    pRing1.AddPoint(p.dLongitude_degree, p.dLatitude_degree)
                pRing1.CloseRings()
                pPolygon1.AddGeometry(pRing1)
                pGeometry.AddGeometry(pPolygon1)
            if len(aPoint_part2) >= 3:
                pPolygon2 = ogr.Geometry(ogr.wkbPolygon)
                pRing2 = ogr.Geometry(ogr.wkbLinearRing)
                for p in aPoint_part2:
                    pRing2.AddPoint(p.dLongitude_degree, p.dLatitude_degree)
                pRing2.CloseRings()
                pPolygon2.AddGeometry(pRing2)
                pGeometry.AddGeometry(pPolygon2)
            pGeometry.FlattenTo2D()  # Ensure the geometry is 2D for WKT export
            sWkt_buffer_polygon = pGeometry.ExportToWkt()
        else:
            if iFlag_cross_idl == 1:
                if self.dLongitude_degree >= 0:
                    # drop points with longitude < 0
                    aPoint = [p for p in aPoint if p.dLongitude_degree >= 0]
                else:
                    # drop points with longitude > 0
                    aPoint = [p for p in aPoint if p.dLongitude_degree <= 0]
                pass
            else:
                pass


            # Use OGR for robust polygon WKT generation.
            pGeometry = ogr.Geometry(ogr.wkbPolygon)
            pRing = ogr.Geometry(ogr.wkbLinearRing)
            for p in aPoint:
                pRing.AddPoint(p.dLongitude_degree, p.dLatitude_degree)
            pRing.CloseRings()
            pGeometry.AddGeometry(pRing)
            pGeometry.FlattenTo2D()  # Ensure the geometry is 2D for WKT export
            sWkt_buffer_polygon = pGeometry.ExportToWkt()

        return sWkt_buffer_polygon, aPoint

    def calculate_xyz(self) -> Tuple[float, float, float]:
        """
        Calculate the x, y, z based on dLongitude and dLatitude

        Returns:
            tuple: The x, y, z in meters (ECEF coordinates)
        """
        dX_meter = 0.0
        dY_meter = 0.0
        dZ_meter = 0.0
        dRadius = 6371000.0  # earth radius in meter
        dX_meter = (
            dRadius * np.cos(self.dLatitude_radian) * np.cos(self.dLongitude_radian)
        )
        dY_meter = (
            dRadius * np.cos(self.dLatitude_radian) * np.sin(self.dLongitude_radian)
        )
        dZ_meter = dRadius * np.sin(self.dLatitude_radian)
        return dX_meter, dY_meter, dZ_meter

    def tojson(self) -> str:
        """
        Convert vertex object to a JSON string.

        Serializes all vertex attributes except internal computed values
        like radians. Uses the custom VertexClassEncoder to handle
        numpy data types.

        Returns:
            str: JSON string representation of the vertex

        Example:
            >>> vertex = pyvertex({'dLongitude_degree': -77.0, 'dLatitude_degree': 38.0})
            >>> json_str = vertex.tojson()
        """
        aSkip = ["dLongitude_radian", "dLatitude_radian"]

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)

        sJson = json.dumps(
            obj, sort_keys=True, indent=4, ensure_ascii=True, cls=VertexClassEncoder
        )
        return sJson

    def towkt(self) -> str:
        """
        Convert a vertex object to a WKT string

        Returns:
            str: A WKT string
        """
        sWKT = "POINT ("
        sWKT += str(self.dLongitude_degree) + " "
        sWKT += str(self.dLatitude_degree) + ")"
        self.wkt = sWKT
        return sWKT

    def set_vertex_id(self, lVertexID: int) -> None:
        """
        Set the vertex ID.

        Args:
            lVertexID (int): New vertex ID

        Raises:
            TypeError: If lVertexID is not an integer
        """
        if not isinstance(lVertexID, (int, np.integer)):
            raise TypeError(f"Vertex ID must be an integer, got {type(lVertexID)}")
        self.lVertexID = int(lVertexID)

    def set_vertex_index(self, lVertexIndex: int) -> None:
        """
        Set the vertex array index.

        Args:
            lVertexIndex (int): New vertex index

        Raises:
            TypeError: If lVertexIndex is not an integer
        """
        if not isinstance(lVertexIndex, (int, np.integer)):
            raise TypeError(
                f"Vertex index must be an integer, got {type(lVertexIndex)}"
            )
        self.lVertexIndex = int(lVertexIndex)

    def set_flowline_id(self, lFlowlineID: int) -> None:
        """
        Set the associated flowline ID.

        Args:
            lFlowlineID (int): New flowline ID

        Raises:
            TypeError: If lFlowlineID is not an integer
        """
        if not isinstance(lFlowlineID, (int, np.integer)):
            raise TypeError(f"Flowline ID must be an integer, got {type(lFlowlineID)}")
        self.lFlowlineID = int(lFlowlineID)

    def is_valid(self) -> bool:
        """
        Check if the vertex has valid attributes.

        Returns:
            bool: True if vertex has valid coordinates and at least one valid ID
        """
        has_valid_coords = (
            -180 <= self.dLongitude_degree <= 180
            and -90 <= self.dLatitude_degree <= 90
        )
        has_valid_id = self.lVertexID > 0 or self.lVertexIndex >= 0
        return has_valid_coords and has_valid_id

    def copy(self) -> "pyvertex":
        """
        Create a deep copy of the vertex.

        Returns:
            pyvertex: A new vertex object with the same attributes
        """
        param = {
            "dLongitude_degree": self.dLongitude_degree,
            "dLatitude_degree": self.dLatitude_degree,
            "x": self.dX_meter,
            "y": self.dY_meter,
            "z": self.dZ_meter,
            "dElevation": self.dElevation,
            "lVertexIndex": self.lVertexIndex,
            "lVertexID": self.lVertexID,
            "lFlowlineID": self.lFlowlineID,
        }
        return pyvertex(param)

