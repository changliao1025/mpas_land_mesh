
import math
from typing import Union, Tuple, List, Optional
import numpy as np
from osgeo import osr, gdal, ogr

from mpas_land_mesh.utilities.geometry  import (
        calculate_distance_based_on_longitude_latitude
    )

from mpas_land_mesh.utilities.constants import earth_radius, IDL_TOLERANCE, IDL_OFFSET

def calculate_angle_between_point(
    dLongitude1_in: float,
    dLatitude1_in: float,
    dLongitude2_in: float,
    dLatitude2_in: float,
    dLongitude3_in: float,
    dLatitude3_in: float,
    iFlag_radian: bool = False,
) -> float:
    """

    Calculates the angle between three points on a sphere.

    Args:
        dLongitude1_in (float): Longitude of the first point.
        dLatitude1_in (float): Latitude of the first point.
        dLongitude2_in (float): Longitude of the second point (the middle one).
        dLatitude2_in (float): Latitude of the second point (the middle one).
        dLongitude3_in (float): Longitude of the third point.
        dLatitude3_in (float): Latitude of the third point.
        iFlag_radian (bool, optional): If True, input coordinates are in radians. Defaults to False (degrees).

    Returns:
        float: The angle in degrees between the vectors from the middle point to the other two.
    """

    # Determine whether the downstream conversion routine should treat inputs as radians.
    iFlag_convert_radian = True if iFlag_radian else None

    # The points in 3D space
    a3 = convert_longitude_latitude_to_sphere_3d(
        dLongitude1_in, dLatitude1_in, iFlag_convert_radian
    )
    b3 = convert_longitude_latitude_to_sphere_3d(
        dLongitude2_in, dLatitude2_in, iFlag_convert_radian
    )
    c3 = convert_longitude_latitude_to_sphere_3d(
        dLongitude3_in, dLatitude3_in, iFlag_convert_radian
    )

    # Vectors in 3D space
    a3vec = a3 - b3
    c3vec = c3 - b3

    angle3deg = calculate_angle_between_vectors_degrees(a3vec, c3vec)
    return angle3deg


def calculate_angle_between_vectors_degrees(
    vector1: np.ndarray, vector2: np.ndarray
) -> float:
    """
    Return the angle between two vectors in any dimension space, in degrees.

    Parameters:
    vector1 (np.ndarray): The first vector.
    vector2 (np.ndarray): The second vector.

    Returns:
    float: The angle between the two vectors in degrees.
    """
    dot_product = np.dot(vector1, vector2)
    norm_product = np.linalg.norm(vector1) * np.linalg.norm(vector2)
    cosine_similarity = np.clip(dot_product / norm_product, -1.0, 1.0)
    angle_in_radians = np.arccos(cosine_similarity)
    angle_in_degrees = np.degrees(angle_in_radians)

    return angle_in_degrees


def calculate_distance_based_on_longitude_latitude(
    aLongitude_from: Union[float, list, np.ndarray],
    aLatitude_from: Union[float, list, np.ndarray],
    aLongitude_to: Union[float, list, np.ndarray],
    aLatitude_to: Union[float, list, np.ndarray],
    bUnits_are_radians: bool = False,
    dRadius_in: float = None,
) -> Union[float, np.ndarray]:
    """
    Calculate the great circle distances between arrays of points
    on the earth using the haversine formula.
    https://en.wikipedia.org/wiki/Great-circle_distance

    This function is a vectorized version that can operate on arrays of coordinates efficiently.

    Args:
        aLongitude_from (Union[float, list, np.ndarray]): The longitudes of the start points
        aLatitude_from (Union[float, list, np.ndarray]): The latitudes of the start points
        aLongitude_to (Union[float, list, np.ndarray]): The longitudes of the end points
        aLatitude_to (Union[float, list, np.ndarray]): The latitudes of the end points
        bUnits_are_radians (bool, optional): If True, input is in radians. Defaults to False (degrees).
        dRadius_in (float, optional): Custom earth radius in meters. If None, default earth_radius is used.

    Returns:
        Union[float, np.ndarray]: The great circle distances in meters (or in radians if bUnits_are_radians is True)
    """
    # Convert inputs to numpy arrays to ensure vectorized operations
    aLongitude_from = np.asarray(aLongitude_from)
    aLatitude_from = np.asarray(aLatitude_from)
    aLongitude_to = np.asarray(aLongitude_to)
    aLatitude_to = np.asarray(aLatitude_to)

    if not bUnits_are_radians:
        dLongitude_radian_from = np.radians(aLongitude_from)
        dLatitude_radian_from = np.radians(aLatitude_from)
        dLongitude_radian_to = np.radians(aLongitude_to)
        dLatitude_radian_to = np.radians(aLatitude_to)
    else:  # already in radian
        dLongitude_radian_from = aLongitude_from
        dLatitude_radian_from = aLatitude_from
        dLongitude_radian_to = aLongitude_to
        dLatitude_radian_to = aLatitude_to

    # haversine formula
    dlon = dLongitude_radian_to - dLongitude_radian_from
    dlat = dLatitude_radian_to - dLatitude_radian_from

    a = (
        np.sin(dlat / 2) ** 2
        + np.cos(dLatitude_radian_from)
        * np.cos(dLatitude_radian_to)
        * np.sin(dlon / 2) ** 2
    )
    c = 2 * np.arcsin(np.sqrt(a))

    if bUnits_are_radians:
        return c
    else:
        radius = dRadius_in if dRadius_in is not None else earth_radius
        d = c * radius
        return d


def calculate_distance_to_plane(
    dLongitude1_in: float,
    dLatitude1_in: float,
    dLongitude2_in: float,
    dLatitude2_in: float,
    dLongitude3_in: float,
    dLatitude3_in: float,
    iFlag_radian: bool = False,
) -> float:
    """
    Calculate the distance from a point to a plane defined by three points in 3D space.

    This function computes the perpendicular distance from a geographic point (point 2)
    to a plane that passes through the Earth's center and three defining points
    (points 1, 2, and 3). The calculation uses 3D Cartesian coordinates on a unit sphere.

    Args:
        dLongitude1_in: Longitude of first plane-defining point (degrees or radians)
        dLatitude1_in: Latitude of first plane-defining point (degrees or radians)
        dLongitude2_in: Longitude of query point (degrees or radians)
        dLatitude2_in: Latitude of query point (degrees or radians)
        dLongitude3_in: Longitude of third plane-defining point (degrees or radians)
        dLatitude3_in: Latitude of third plane-defining point (degrees or radians)
        iFlag_radian: If True, input coordinates are in radians; if False, in degrees (default: False)

    Returns:
        float: Perpendicular distance from point 2 to the plane defined by points 1, 2, and 3
               in unit sphere coordinates (dimensionless, typically < 2)

    Note:
        - Point 2 is the query point whose distance to the plane is calculated
        - Points 1 and 3 (along with point 2) define the plane
        - The plane passes through the Earth's center (origin)
        - Returns 0.0 if the three points are collinear (cannot define a unique plane)
        - Distance is in unit sphere coordinates, not meters
        - For great circle distance checking, see related functions

    Example:
        >>> # Check if a point lies on a great circle defined by two other points
        >>> distance = calculate_distance_to_plane(
        ...     -120.0, 37.0,  # Point 1: defines plane
        ...     -122.0, 38.0,  # Point 2: query point
        ...     -118.0, 36.0   # Point 3: defines plane
        ... )
        >>> if distance < 0.001:
        ...     print("Point 2 is approximately on the great circle through points 1 and 3")

    Algorithm:
        1. Convert all three geographic coordinates to 3D Cartesian (x,y,z) on unit sphere
        2. Create two vectors in the plane:
           - v1: from point 1 to point 2
           - v2: from point 1 to point 3
        3. Calculate plane normal vector using cross product: n = v1 × v2
        4. Check if normal is zero (points are collinear) → return 0.0
        5. Compute plane equation: Ax + By + Cz + D = 0 where:
           - (A, B, C) = normal vector
           - D = -(A*x1 + B*y1 + C*z1)
        6. Calculate perpendicular distance using point-to-plane formula:
           distance = |Ax2 + By2 + Cz2 + D| / sqrt(A² + B² + C²)

    References:
        - Point-to-plane distance: https://mathworld.wolfram.com/Point-PlaneDistance.html
        - Great circle checking: https://stackoverflow.com/questions/8204998/
    """
    # Convert to radians if input is in degrees
    if iFlag_radian:
        lon1_rad, lat1_rad = dLongitude1_in, dLatitude1_in
        lon2_rad, lat2_rad = dLongitude2_in, dLatitude2_in
        lon3_rad, lat3_rad = dLongitude3_in, dLatitude3_in
    else:
        lon1_rad, lat1_rad = np.radians([dLongitude1_in, dLatitude1_in])
        lon2_rad, lat2_rad = np.radians([dLongitude2_in, dLatitude2_in])
        lon3_rad, lat3_rad = np.radians([dLongitude3_in, dLatitude3_in])

    # Convert geographic coordinates to 3D Cartesian coordinates on unit sphere
    x1, y1, z1 = convert_longitude_latitude_to_sphere_3d(
        lon1_rad, lat1_rad
    )  # Plane point 1
    x2, y2, z2 = convert_longitude_latitude_to_sphere_3d(
        lon2_rad, lat2_rad
    )  # Query point
    x3, y3, z3 = convert_longitude_latitude_to_sphere_3d(
        lon3_rad, lat3_rad
    )  # Plane point 2

    # Calculate two vectors in the plane
    # v1: from point 1 to point 2
    v1_x, v1_y, v1_z = x2 - x1, y2 - y1, z2 - z1
    # v2: from point 1 to point 3
    v2_x, v2_y, v2_z = x3 - x1, y3 - y1, z3 - z1

    # Compute the normal vector using cross product: n = v1 × v2
    normal_x = v1_y * v2_z - v1_z * v2_y
    normal_y = v1_z * v2_x - v1_x * v2_z
    normal_z = v1_x * v2_y - v1_y * v2_x

    # Check if the normal vector is zero (points are collinear)
    # If collinear, the three points don't define a unique plane
    if abs(normal_x) < 1e-10 and abs(normal_y) < 1e-10 and abs(normal_z) < 1e-10:
        return 0.0

    # Calculate the plane equation coefficients: Ax + By + Cz + D = 0
    A, B, C = normal_x, normal_y, normal_z
    D = -(A * x1 + B * y1 + C * z1)

    # Calculate the perpendicular distance from point 2 to the plane
    # Formula: distance = |Ax + By + Cz + D| / sqrt(A² + B² + C²)
    numerator = abs(A * x2 + B * y2 + C * z2 + D)
    denominator = math.sqrt(A**2 + B**2 + C**2)
    distance = numerator / denominator

    return distance


def project_point_onto_plane(p: np.ndarray, normal: np.ndarray) -> np.ndarray:
    """Project a 3D point onto a plane defined by its normal vector.

    Given a point and a plane (defined by a normal vector passing through
    the origin), this function calculates the orthogonal projection of the
    point onto the plane.

    Parameters
    ----------
    p : np.ndarray
        3D point to be projected, shape (3,).
        Example: [x, y, z]
    normal : np.ndarray
        Normal vector of the plane, shape (3,).
        Should be a unit vector (magnitude = 1) for correct results.
        The plane passes through the origin with this normal direction.

    Returns
    -------
    np.ndarray
        Projected point on the plane, shape (3,).
        This is the point on the plane closest to the input point.

    Raises
    ------
    ValueError
        If normal vector has zero magnitude.

    Notes
    -----
    - The plane is defined by all points r such that: normal · r = 0
      (plane through origin perpendicular to normal)
    - Projection formula: p_proj = p - (p · n) * n, where n is unit normal
    - If normal is not unit length, it will be normalized automatically
    - Distance from point to plane is: |p · n| (when n is unit vector)

    Examples
    --------
    >>> # Project point onto xy-plane (normal = z-axis)
    >>> p = np.array([1, 2, 3])
    >>> normal = np.array([0, 0, 1])
    >>> proj = project_point_onto_plane(p, normal)
    >>> np.allclose(proj, [1, 2, 0])
    True

    See Also
    --------
    calculate_distance_to_plane : Calculate distance from point to plane
    """
    # Validate and normalize the normal vector
    normal_magnitude = np.linalg.norm(normal)

    if normal_magnitude < 1e-10:
        raise ValueError(
            "Normal vector has zero or near-zero magnitude. "
            f"Got magnitude: {normal_magnitude}. Cannot define a plane."
        )

    # Normalize to unit vector if not already
    if not np.isclose(normal_magnitude, 1.0):
        normal_unit = normal / normal_magnitude
    else:
        normal_unit = normal

    # Calculate the projection
    # Distance from point to plane along normal direction
    distance_to_plane = np.dot(p, normal_unit)

    # Subtract the normal component to get projection
    projection = p - distance_to_plane * normal_unit

    return projection


def calculate_intersect_on_great_circle(
    dLongitude1_in: float,
    dLatitude1_in: float,
    dLongitude2_in: float,
    dLatitude2_in: float,
    dLongitude3_in: float,
    dLatitude3_in: float,
    iFlag_radian: bool = False,
) -> tuple:
    """Calculate the closest point on a great circle arc to a query point.

    This function finds the point on the great circle arc between point1 and
    point3 that is closest to point2 (the query point). This is useful for
    finding the nearest location on a line to a given point in mesh operations.

    Parameters
    ----------
    dLongitude1_in : float
        Longitude of first point on the great circle arc.
        In degrees by default, or radians if iFlag_radian=True.
    dLatitude1_in : float
        Latitude of first point on the great circle arc.
        In degrees by default, or radians if iFlag_radian=True.
    dLongitude2_in : float
        Longitude of query point (point to find nearest location to).
        In degrees by default, or radians if iFlag_radian=True.
    dLatitude2_in : float
        Latitude of query point (point to find nearest location to).
        In degrees by default, or radians if iFlag_radian=True.
    dLongitude3_in : float
        Longitude of second point on the great circle arc.
        In degrees by default, or radians if iFlag_radian=True.
    dLatitude3_in : float
        Latitude of second point on the great circle arc.
        In degrees by default, or radians if iFlag_radian=True.
    iFlag_radian : bool, optional
        If True, input coordinates are in radians. If False (default),
        input coordinates are in degrees.

    Returns
    -------
    tuple
        (longitude, latitude) of the closest point on the great circle arc.
        Always returned in degrees.

    Notes
    -----
        - The returned point is the orthogonal projection onto the great circle,
            which may lie outside the arc segment between point1 and point3
        - All calculations are performed on a unit sphere
        - The great circle is defined by the plane through point1, point3,
            and Earth's center

    Examples
    --------
    >>> # Find closest point on equator arc to North Pole
    >>> lon, lat = calculate_intersect_on_great_circle(
    ...     0, 0,      # Point 1: (0°E, 0°N) on equator
    ...     0, 90,     # Query: North Pole
    ...     90, 0      # Point 3: (90°E, 0°N) on equator
    ... )
    >>> # Result should be close to (0°, 0°) - nearest point on equator

    Algorithm
    ---------
    1. Convert all points to 3D Cartesian coordinates on unit sphere
    2. Calculate great circle plane normal: n = normalize(cross(point1, point3))
    3. Project query point onto plane: projected = query - dot(query, n) * n
    4. Normalize projected point to sphere surface
    5. Convert back to longitude/latitude

    See Also
    --------
    calculate_distance_to_line : Calculate perpendicular distance to line segment
    calculate_distance_to_plane : Calculate distance from point to plane
    """
    # Convert all points to 3D Cartesian coordinates
    point1_3d = convert_longitude_latitude_to_sphere_3d(
        dLongitude1_in, dLatitude1_in, iFlag_radian
    )
    query_point_3d = convert_longitude_latitude_to_sphere_3d(
        dLongitude2_in, dLatitude2_in, iFlag_radian
    )
    point3_3d = convert_longitude_latitude_to_sphere_3d(
        dLongitude3_in, dLatitude3_in, iFlag_radian
    )

    # Calculate the normal vector to the great circle plane
    # The plane passes through point1, point3, and Earth's center
    plane_normal = np.cross(point1_3d, point3_3d)
    plane_normal /= np.linalg.norm(plane_normal)

    # Project the query point onto the great circle plane
    distance_to_plane = np.dot(query_point_3d, plane_normal)
    projected_point = query_point_3d - distance_to_plane * plane_normal

    # Normalize to sphere surface (the projection may not be on unit sphere)
    projected_point /= np.linalg.norm(projected_point)

    # Convert back to longitude/latitude (always in degrees)
    longitude_intersect, latitude_intersect = convert_sphere_3d_to_longitude_latitude(
        *projected_point
    )

    return float(longitude_intersect), float(latitude_intersect)


def find_great_circle_intersection_with_meridian(
    lon1: float, lat1: float, lon2: float, lat2: float, target_lon: float
) -> tuple:
    """Find the latitude on a great circle arc at a specified longitude.

    Given two points defining a great circle, find the latitude where
    the great circle crosses a specified longitude meridian. This is commonly
    used for polygon splitting operations at the International Date Line (±180°).

    Parameters
    ----------
    lon1 : float
        Longitude of the first point in degrees.
    lat1 : float
        Latitude of the first point in degrees.
    lon2 : float
        Longitude of the second point in degrees.
    lat2 : float
        Latitude of the second point in degrees.
    target_lon : float
        The target longitude in degrees where intersection is sought.

    Returns
    -------
    tuple
        (target_lon, target_lat) where target_lat is the latitude at the
        target longitude.

    Raises
    ------
    ValueError
        If the great circle does not cross the target longitude meridian,
        or if the two input points define the same location.

    Notes
    -----
    - A great circle may cross a meridian at 0, 1, or 2 points
    - This function returns one intersection point
    - The target longitude is returned unchanged in the result tuple
    - All angles are in degrees

    Examples
    --------
    >>> # Great circle from (0°E, 0°N) to (180°E, 0°N) crosses 90°E at equator
    >>> lon, lat = calculate_great_circle_latitude_at_longitude(0, 0, 180, 0, 90)
    >>> np.allclose([lon, lat], [90, 0])
    True

    >>> # Great circle crossing International Date Line
    >>> lon, lat = calculate_great_circle_latitude_at_longitude(170, 45, -170, 50, 180)
    >>> # Returns latitude where circle crosses ±180°

    Algorithm
    ---------
    1. Convert both points to 3D Cartesian coordinates
    2. Calculate great circle plane normal: n = cross(p1, p2)
    3. Solve plane equation for latitude at target longitude:
       n_x * cos(lat) * cos(lon) + n_y * cos(lat) * sin(lon) + n_z * sin(lat) = 0
    4. Simplify to: tan(lat) = -(n_x*cos(lon) + n_y*sin(lon)) / n_z

    See Also
    --------
    calculate_intersect_on_great_circle : Find closest point on great circle
    split_polygon_cross_idl : Uses this function for polygon splitting

    References
    ----------
    .. [1] Williams, E. "Aviation Formulary V1.47"
           http://www.edwilliams.org/avform147.htm
    """
    # Convert to radians
    lon1_rad = np.radians(lon1)
    lat1_rad = np.radians(lat1)
    lon2_rad = np.radians(lon2)
    lat2_rad = np.radians(lat2)
    target_lon_rad = np.radians(target_lon)

    # Convert both points to 3D Cartesian coordinates
    p1 = convert_longitude_latitude_to_sphere_3d(lon1, lat1, iFlag_radian=False)
    p2 = convert_longitude_latitude_to_sphere_3d(lon2, lat2, iFlag_radian=False)

    # Calculate the great circle plane normal
    normal = np.cross(p1, p2)
    normal_magnitude = np.linalg.norm(normal)

    # Check if points are the same or antipodal
    if normal_magnitude < 1e-10:
        raise ValueError(
            f"Points ({lon1}, {lat1}) and ({lon2}, {lat2}) are too close "
            "or antipodal to define a unique great circle."
        )

    normal = normal / normal_magnitude

    # The great circle plane equation is: normal · r = 0
    # For a point at longitude=target_lon: r = [cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)]
    # Substituting into plane equation and solving for lat:
    # n_x * cos(lat) * cos(target_lon) + n_y * cos(lat) * sin(target_lon) + n_z * sin(lat) = 0
    # cos(lat) * (n_x * cos(target_lon) + n_y * sin(target_lon)) + n_z * sin(lat) = 0

    nx, ny, nz = normal
    cos_target = np.cos(target_lon_rad)
    sin_target = np.sin(target_lon_rad)

    # Coefficient of cos(lat)
    A = nx * cos_target + ny * sin_target
    # Coefficient of sin(lat)
    B = nz

    # We have: A * cos(lat) + B * sin(lat) = 0
    # This can be written as: tan(lat) = -A / B (if B != 0)
    #
    # If B ~= 0, the great-circle plane contains the Earth's rotation axis,
    # so the great circle passes through both poles. In that case the target
    # meridian intersection is a pole, and tiny non-zero A values are usually
    # just floating-point noise.
    tolerance = 1e-8

    if abs(B) < tolerance:
        mean_lat = 0.5 * (lat1 + lat2)
        fallback_lat = 90.0 if mean_lat >= 0.0 else -90.0

        if abs(A) < tolerance:
            # The target meridian lies nearly in the great-circle plane.
            # Return a stable representative point on the relevant pole side.
            return (target_lon, fallback_lat)

        # Pole-passing great circle: every meridian intersects it at the poles.
        return (target_lon, fallback_lat)

    # Calculate latitude
    lat_rad = np.arctan2(-A, B)
    lat_deg = np.degrees(lat_rad)

    return (target_lon, float(lat_deg))

def convert_longitude_latitude_to_sphere_3d(
    dLongitude_in: float, dLatitude_in: float, iFlag_radian: bool = False
) -> np.ndarray:
    """Convert longitude/latitude to 3D Cartesian coordinates on a unit sphere.

    Transforms geographic coordinates (longitude, latitude) to 3D Cartesian
    coordinates (x, y, z) on a sphere with radius 1. This is useful for
    spherical geometry calculations, great circle distances, and 3D visualizations.

    Parameters
    ----------
    dLongitude_in : float
        Longitude coordinate. In degrees by default, or radians if iFlag_radian=True.
        Valid range: [-180, 180] degrees or [-π, π] radians.
    dLatitude_in : float
        Latitude coordinate. In degrees by default, or radians if iFlag_radian=True.
        Valid range: [-90, 90] degrees or [-π/2, π/2] radians.
    iFlag_radian : bool, optional
        If True, input coordinates are in radians. If False (default),
        input coordinates are in degrees.

    Returns
    -------
    np.ndarray
        3D Cartesian coordinates [x, y, z] on a unit sphere (radius = 1).
        - x = cos(lat) * cos(lon)
        - y = cos(lat) * sin(lon)
        - z = sin(lat)

    Raises
    ------
    ValueError
        If latitude is outside valid range [-90, 90] degrees or [-π/2, π/2] radians.
        If longitude is outside valid range [-180, 180] degrees or [-π, π] radians.
    TypeError
        If input coordinates cannot be converted to float.

    Notes
    -----
    - The returned vector has magnitude (norm) equal to 1
    - This uses the standard geographic coordinate system:
      * Longitude: 0° at Prime Meridian, positive East, negative West
      * Latitude: 0° at Equator, positive North, negative South
    - The 3D coordinate system:
      * x-axis: points to (0°N, 0°E) - Equator at Prime Meridian
      * y-axis: points to (0°N, 90°E) - Equator at 90° East
      * z-axis: points to North Pole (90°N)

    Examples
    --------
    >>> # North Pole
    >>> coords = convert_longitude_latitude_to_sphere_3d(0, 90)
    >>> np.allclose(coords, [0, 0, 1])
    True

    >>> # Equator at Prime Meridian
    >>> coords = convert_longitude_latitude_to_sphere_3d(0, 0)
    >>> np.allclose(coords, [1, 0, 0])
    True

    >>> # Equator at 90° East
    >>> coords = convert_longitude_latitude_to_sphere_3d(90, 0)
    >>> np.allclose(coords, [0, 1, 0])
    True

    >>> # Using radians
    >>> coords = convert_longitude_latitude_to_sphere_3d(np.pi/2, 0, iFlag_radian=True)
    >>> np.allclose(coords, [0, 1, 0])
    True

    >>> # Verify unit sphere (magnitude = 1)
    >>> coords = convert_longitude_latitude_to_sphere_3d(45, 45)
    >>> np.allclose(np.linalg.norm(coords), 1.0)
    True
    """
    # Validate and convert inputs to float
    try:
        longitude = float(dLongitude_in)
        latitude = float(dLatitude_in)
    except (TypeError, ValueError) as e:
        raise TypeError(
            f"Coordinates must be numeric values. "
            f"Got longitude={dLongitude_in}, latitude={dLatitude_in}. Error: {e}"
        )

    # Convert to radians if needed and validate ranges
    if not iFlag_radian:
        # Input in degrees - validate and convert
        if not -90 <= latitude <= 90:
            raise ValueError(
                f"Latitude must be in range [-90, 90] degrees. Got {latitude}°"
            )
        if not -180 <= longitude <= 180:
            raise ValueError(
                f"Longitude must be in range [-180, 180] degrees. Got {longitude}°"
            )
        longitude_rad, latitude_rad = np.radians([longitude, latitude])
    else:
        # Input in radians - validate
        if not -np.pi / 2 <= latitude <= np.pi / 2:
            raise ValueError(
                f"Latitude must be in range [-π/2, π/2] radians. Got {latitude} rad"
            )
        if not -np.pi <= longitude <= np.pi:
            raise ValueError(
                f"Longitude must be in range [-π, π] radians. Got {longitude} rad"
            )
        longitude_rad = longitude
        latitude_rad = latitude

    # Convert to 3D Cartesian coordinates on unit sphere
    # x = cos(lat) * cos(lon)
    # y = cos(lat) * sin(lon)
    # z = sin(lat)
    cos_lat = np.cos(latitude_rad)

    x = cos_lat * np.cos(longitude_rad)
    y = cos_lat * np.sin(longitude_rad)
    z = np.sin(latitude_rad)

    return np.array([x, y, z])


def convert_sphere_3d_to_longitude_latitude(
    x: Union[float, np.ndarray],
    y: Union[float, np.ndarray],
    z: Union[float, np.ndarray],
) -> Tuple[float, float]:
    """Convert 3D Cartesian coordinates to longitude/latitude coordinates.

    Transforms 3D Cartesian coordinates (x, y, z) to geographic coordinates
    (longitude, latitude) on a sphere. The input coordinates are automatically
    normalized to a unit sphere before conversion.

    This is the inverse operation of convert_longitude_latitude_to_sphere_3d.

    Parameters
    ----------
    x : float or np.ndarray
        X coordinate in 3D Cartesian space.
        - Points to (0°N, 0°E) - Equator at Prime Meridian
        - For unit sphere: x = cos(lat) * cos(lon)
    y : float or np.ndarray
        Y coordinate in 3D Cartesian space.
        - Points to (0°N, 90°E) - Equator at 90° East
        - For unit sphere: y = cos(lat) * sin(lon)
    z : float or np.ndarray
        Z coordinate in 3D Cartesian space.
        - Points to North Pole (90°N)
        - For unit sphere: z = sin(lat)

    Returns
    -------
    Tuple[float, float]
        Geographic coordinates (longitude, latitude) in decimal degrees:
        - longitude: Range [-180, 180] degrees (negative = West, positive = East)
        - latitude: Range [-90, 90] degrees (negative = South, positive = North)

    Raises
    ------
    ValueError
        If all coordinates are zero (degenerate case, no unique point).
        If inputs contain NaN or infinite values.
    TypeError
        If inputs cannot be converted to numeric values.

    Notes
    -----
    - Input coordinates are normalized before conversion (divided by their magnitude)
    - Works for any non-zero radius; automatically projects to unit sphere
    - Uses atan2 for longitude to handle all quadrants correctly
    - Uses asin for latitude (assumes normalized coordinates)
    - Z-coordinate is clamped to [-1, 1] to handle floating-point precision errors

    Mathematical Formulas
    ---------------------
    After normalization to unit sphere:
        longitude = atan2(y, x)
        latitude = asin(z)

    Examples
    --------
    >>> # North Pole
    >>> lon, lat = convert_sphere_3d_to_longitude_latitude(0, 0, 1)
    >>> np.allclose([lon, lat], [0, 90])
    True

    >>> # Equator at Prime Meridian
    >>> lon, lat = convert_sphere_3d_to_longitude_latitude(1, 0, 0)
    >>> np.allclose([lon, lat], [0, 0])
    True

    >>> # Equator at 90° East
    >>> lon, lat = convert_sphere_3d_to_longitude_latitude(0, 1, 0)
    >>> np.allclose([lon, lat], [90, 0])
    True

    >>> # Equator at 180° (or -180°)
    >>> lon, lat = convert_sphere_3d_to_longitude_latitude(-1, 0, 0)
    >>> np.allclose(abs(lon), 180) and np.allclose(lat, 0)
    True

    >>> # Point at 45°N, 45°E (on unit sphere)
    >>> x = np.cos(np.radians(45)) * np.cos(np.radians(45))
    >>> y = np.cos(np.radians(45)) * np.sin(np.radians(45))
    >>> z = np.sin(np.radians(45))
    >>> lon, lat = convert_sphere_3d_to_longitude_latitude(x, y, z)
    >>> np.allclose([lon, lat], [45, 45], atol=1e-10)
    True

    >>> # Works with any radius (auto-normalized)
    >>> lon, lat = convert_sphere_3d_to_longitude_latitude(2, 0, 0)
    >>> np.allclose([lon, lat], [0, 0])
    True

    >>> # Verify round-trip conversion
    >>> xyz = convert_longitude_latitude_to_sphere_3d(45, 30)
    >>> lon, lat = convert_sphere_3d_to_longitude_latitude(*xyz)
    >>> np.allclose([lon, lat], [45, 30])
    True

    See Also
    --------
    convert_longitude_latitude_to_sphere_3d : Inverse transformation (lon/lat to xyz)
    numpy.arctan2 : 2-argument arctangent for correct quadrant
    numpy.arcsin : Arcsine function

    References
    ----------
    .. [1] Spherical coordinate system:
           https://en.wikipedia.org/wiki/Spherical_coordinate_system
    .. [2] Geographic coordinate system:
           https://en.wikipedia.org/wiki/Geographic_coordinate_system
    """
    # Validate and convert inputs
    try:
        x_val = float(x) if np.isscalar(x) else np.asarray(x, dtype=float)
        y_val = float(y) if np.isscalar(y) else np.asarray(y, dtype=float)
        z_val = float(z) if np.isscalar(z) else np.asarray(z, dtype=float)
    except (TypeError, ValueError) as e:
        raise TypeError(
            f"Coordinates must be numeric values. "
            f"Got x={x}, y={y}, z={z}. Error: {e}"
        )

    # Check for NaN or infinite values
    if np.any(np.isnan([x_val, y_val, z_val])):
        raise ValueError(
            "Coordinates cannot contain NaN values. "
            f"Got x={x_val}, y={y_val}, z={z_val}"
        )
    if np.any(np.isinf([x_val, y_val, z_val])):
        raise ValueError(
            "Coordinates cannot contain infinite values. "
            f"Got x={x_val}, y={y_val}, z={z_val}"
        )

    # Calculate magnitude (norm)
    norm = np.sqrt(x_val**2 + y_val**2 + z_val**2)

    # Check for degenerate case (all zeros)
    if norm == 0 or np.isclose(norm, 0, atol=1e-15):
        raise ValueError(
            "Cannot convert origin (0, 0, 0) to longitude/latitude. "
            "All coordinates are zero (or extremely close to zero)."
        )

    # Normalize coordinates to unit sphere
    x_norm = x_val / norm
    y_norm = y_val / norm
    z_norm = z_val / norm

    # Convert to spherical coordinates
    # Longitude: use atan2 for correct quadrant handling
    lon_rad = np.arctan2(y_norm, x_norm)

    # Latitude: use asin (z is already normalized)
    # Clamp z_norm to [-1, 1] to handle floating-point errors
    z_clamped = np.clip(z_norm, -1.0, 1.0)
    lat_rad = np.arcsin(z_clamped)

    # Convert radians to degrees
    lon_deg = np.degrees(lon_rad)
    lat_deg = np.degrees(lat_rad)

    return float(lon_deg), float(lat_deg)


def find_minimal_enclosing_polygon(aLongitude_degree, aLatitude_degree):
    from scipy.spatial import ConvexHull
    # Convert aLatitude_degree/aLongitude_degree points to Cartesian coordinates
    aVertex_cartesian = list()
    for i in range(len(aLongitude_degree)):
        x = aLongitude_degree[i]
        y = aLatitude_degree[i]
        aVertex_cartesian.append((x, y))

    aVertex_cartesian = np.array(aVertex_cartesian)
    # Create a ConvexHull object from the Cartesian coordinates
    hull = ConvexHull(aVertex_cartesian)
    # Get the indices of the points on the convex hull
    indices = hull.vertices
    # Get the points on the convex hull in Cartesian coordinates
    aVertex_on_hull_cartesian = [aVertex_cartesian[i] for i in indices]
    # Convert the points on the convex hull back to aLatitude_degree/aLongitude_degree coordinates
    aVertex_on_hull_latlon = []
    for p in aVertex_on_hull_cartesian:
        dLongitude_degree = p[0]
        dLatitude_degree = p[1]
        aVertex_on_hull_latlon.append((dLongitude_degree, dLatitude_degree))

    return aVertex_on_hull_latlon


def calculate_spherical_triangle_area(
    aLongitude_in: Union[List[float], np.ndarray],
    aLatitude_in: Union[List[float], np.ndarray],
    iFlag_radian: bool = False,
    dRadius_in: Optional[float] = None,
) -> float:
    """Calculate the area of a spherical triangle using L'Huilier's theorem.

    Computes the area of a triangle on a spherical surface given the
    longitude/latitude coordinates of its three points. Uses L'Huilier's
    formula for calculating the spherical excess, which is numerically
    stable for triangles of all sizes.

    Parameters
    ----------
    aLongitude_in : list or np.ndarray
        Longitudes of the three triangle points.
        Must contain exactly 3 values.
        In degrees by default, or radians if iFlag_radian=True.
    aLatitude_in : list or np.ndarray
        Latitudes of the three triangle points.
        Must contain exactly 3 values.
        In degrees by default, or radians if iFlag_radian=True.
    iFlag_radian : bool, optional
        If True, input coordinates are in radians and output is spherical
        excess in radians. If False (default), input is in degrees and
        output is area in square meters (or square units of dRadius_in).
    dRadius_in : float, optional
        Sphere radius in meters. If None (default), uses Earth's mean
        radius from global variables (approximately 6371229 m).
        Only used when iFlag_radian=False.

    Returns
    -------
    float
        If iFlag_radian=True: Spherical excess in radians.
        If iFlag_radian=False: Triangle area in square meters
        (or square units of dRadius_in).

    Raises
    ------
    ValueError
        If aLongitude_in or aLatitude_in don't contain exactly 3 values.
        If the three points are collinear (degenerate triangle).

    Notes
    -----
    - The algorithm uses L'Huilier's theorem for numerical stability
    - All calculations are performed on a perfect sphere
    - The three points must not be collinear
    - Triangle points should be ordered (clockwise or counterclockwise)

    Algorithm
    ---------
     1. Calculate great circle distances between each pair of points:
         a = distance(point0, point1)
         b = distance(point1, point2)
         c = distance(point2, point0)

    2. Compute semi-perimeter:
       s = (a + b + c) / 2

    3. Apply L'Huilier's formula for spherical excess:
       tan(E/4) = sqrt(tan(s/2) * tan((s-a)/2) * tan((s-b)/2) * tan((s-c)/2))
       E = 4 * arctan(tan(E/4))

    4. Calculate area:
       Area = E * R²

    Examples
    --------
    >>> # Equilateral triangle on Earth's surface
    >>> # Vertices at (0°, 0°), (1°, 0°), (0.5°, 0.866°)
    >>> lons = [0, 1, 0.5]
    >>> lats = [0, 0, 0.866]
    >>> area = calculate_spherical_triangle_area(lons, lats)
    >>> # Returns area in square meters

    >>> # Using radians and getting spherical excess
    >>> import numpy as np
    >>> lons_rad = np.radians([0, 10, 5])
    >>> lats_rad = np.radians([0, 0, 8.66])
    >>> excess = calculate_spherical_triangle_area(
    ...     lons_rad, lats_rad, iFlag_radian=True
    ... )
    >>> # Returns spherical excess in radians

    >>> # Custom sphere radius
    >>> area = calculate_spherical_triangle_area(
    ...     [0, 1, 0.5], [0, 0, 1], dRadius_in=1000.0
    ... )
    >>> # Returns area in square units of radius 1000

    See Also
    --------
    calculate_polygon_area : Calculate area of multi-sided spherical polygon
    calculate_distance_based_on_longitude_latitude : Great circle distance

    References
    ----------
    .. [1] L'Huilier, S.-A.-J. "Mémoire sur la polyèdrométrie", 1812.
    .. [2] Girard's theorem: https://mathworld.wolfram.com/GirardsSphericalExcessFormula.html
    .. [3] L'Huilier's theorem: https://mathworld.wolfram.com/LHuiliersTheorem.html
    """
    # Validate inputs
    aLongitude_in = np.asarray(aLongitude_in)
    aLatitude_in = np.asarray(aLatitude_in)

    if len(aLongitude_in) != 3:
        raise ValueError(
            f"Triangle requires exactly 3 vertices. "
            f"Got {len(aLongitude_in)} longitude values."
        )

    if len(aLatitude_in) != 3:
        raise ValueError(
            f"Triangle requires exactly 3 vertices. "
            f"Got {len(aLatitude_in)} latitude values."
        )

    # Convert to radians if needed
    if iFlag_radian:
        aLongitude_radian = aLongitude_in
        aLatitude_radian = aLatitude_in
    else:
        aLongitude_radian = np.radians(aLongitude_in)
        aLatitude_radian = np.radians(aLatitude_in)

    # Calculate the three great circle arc lengths (edges of the triangle)
    # Edge a: from point 0 to point 1
    a = calculate_distance_based_on_longitude_latitude(
        aLongitude_radian[0],
        aLatitude_radian[0],
        aLongitude_radian[1],
        aLatitude_radian[1],
        iFlag_radian=True,
    )

    # Edge b: from point 1 to point 2
    b = calculate_distance_based_on_longitude_latitude(
        aLongitude_radian[1],
        aLatitude_radian[1],
        aLongitude_radian[2],
        aLatitude_radian[2],
        iFlag_radian=True,
    )

    # Edge c: from point 2 to point 0
    c = calculate_distance_based_on_longitude_latitude(
        aLongitude_radian[2],
        aLatitude_radian[2],
        aLongitude_radian[0],
        aLatitude_radian[0],
        iFlag_radian=True,
    )

    # Check for degenerate triangle (collinear points)
    if np.isclose(a + b, c) or np.isclose(b + c, a) or np.isclose(a + c, b):
        raise ValueError(
            "Triangle points are collinear (degenerate triangle). "
            f"Edge lengths: a={a}, b={b}, c={c}"
        )

    # Calculate semi-perimeter
    s = 0.5 * (a + b + c)

    # L'Huilier's formula for spherical excess
    # tan(E/4) = sqrt(tan(s/2) * tan((s-a)/2) * tan((s-b)/2) * tan((s-c)/2))
    tan_quarter_excess = np.sqrt(
        np.tan(0.5 * s)
        * np.tan(0.5 * (s - a))
        * np.tan(0.5 * (s - b))
        * np.tan(0.5 * (s - c))
    )

    # Spherical excess in radians
    spherical_excess = 4.0 * np.arctan(tan_quarter_excess)

    # Return based on flag
    if iFlag_radian:
        # Return spherical excess in radians
        return float(spherical_excess)
    else:
        # Calculate area in square meters (or square units of radius)
        if dRadius_in is not None:
            area = spherical_excess * dRadius_in**2
        else:
            area = spherical_excess * earth_radius**2

        return float(area)


def haversine(x: float) -> float:
    """Calculate the haversine function: hav(x) = (1 - cos(x)) / 2.

    The haversine function is commonly used in spherical trigonometry
    calculations, particularly for great circle distances and areas.

    Parameters
    ----------
    x : float
        Angle in radians.

    Returns
    -------
    float
        The haversine of x, in range [0, 1].
        - hav(0) = 0
        - hav(π) = 1
        - hav(π/2) = 0.5

    Notes
    -----
    This is equivalent to: sin²(x/2)

    The haversine function was historically important for navigation
    before the advent of electronic calculators, as it avoids issues
    with numerical precision when calculating small angles.

    Examples
    --------
    >>> haversine(0)
    0.0
    >>> haversine(np.pi)
    1.0
    >>> np.isclose(haversine(np.pi/2), 0.5)
    True

    See Also
    --------
    spherical_polygon_area : Uses haversine in Karney's algorithm
    """
    return (1.0 - math.cos(x)) / 2.0


def calculate_polygon_area(
    aLongitude_in: Union[list, np.ndarray],
    aLatitude_in: Union[list, np.ndarray],
    iFlag_algorithm: int = 2,
    iFlag_radian: bool = False,
    dRadius_in: Optional[float] = None,
    dLine_threshold: Optional[float] = None,
) -> float:
    """Calculate area of a spherical polygon on Earth's surface.

    Computes the area of a polygon on a spherical Earth using one of three
    available algorithms. The polygon is automatically closed if the first
    and last points differ.

    Parameters
    ----------
    aLongitude_in : list or np.ndarray
        Longitude coordinates of polygon vertices.
        In degrees by default (or radians if iFlag_radian=True).
    aLatitude_in : list or np.ndarray
        Latitude coordinates of polygon vertices.
        In degrees by default (or radians if iFlag_radian=True).
        Must have same length as aLongitude_in.
    iFlag_algorithm : int, optional
        Algorithm selection (default: 2, recommended):
        - 0: Green's Theorem line integral (fast, good for regular polygons)
        - 1: L'Huilier's theorem triangulation (good for convex polygons)
        - 2: Karney's method (most accurate, handles edge cases)
    iFlag_radian : bool, optional
        If True, input coordinates are in radians. If False (default),
        coordinates are in degrees.
    dRadius_in : float, optional
        Sphere radius in meters. If None (default), uses Earth's mean
        radius from global variables (approximately 6371229 m).
    dLine_threshold : float, optional
        Threshold for detecting degenerate line-like polygons.
        If max_edge_length / sum_other_edges > (1 - threshold),
        the polygon is considered degenerate and area = 0.
        If None, no check is performed.

    Returns
    -------
    float
        Polygon area in square meters (or square units of dRadius_in).
        Returns 0.0 for degenerate polygons when dLine_threshold is set.

    Raises
    ------
    ValueError
        If less than 3 points are provided (not a valid polygon).
    ValueError
        If aLongitude_in and aLatitude_in have different lengths.

    Notes
    -----
    - The polygon is automatically closed if not already closed
    - All algorithms assume a perfect sphere, not the WGS84 ellipsoid
    - For more accurate ellipsoidal calculations, consider using nvector API
    - Algorithm 2 (Karney) is recommended for general use
    - The area is always positive (uses absolute value internally)

    Examples
    --------
    >>> # Square on equator (approx 111km × 111km)
    >>> lon = [0, 1, 1, 0]
    >>> lat = [0, 0, 1, 1]
    >>> area = calculate_polygon_area(lon, lat, iFlag_algorithm=2)
    >>> # Result: approximately 12,364 km²

    >>> # Triangle with custom radius
    >>> lon = [0, 1, 0.5]
    >>> lat = [0, 0, 1]
    >>> area = calculate_polygon_area(lon, lat, dRadius_in=1000.0)

    >>> # Detect degenerate polygon
    >>> lon = [0, 100, 0.01]  # Nearly a line
    >>> lat = [0, 0, 0]
    >>> area = calculate_polygon_area(lon, lat, dLine_threshold=0.01)
    >>> # Returns 0.0 if polygon is too narrow

    Algorithm Details
    -----------------
    Algorithm 0 (Green's Theorem):
        Based on line integral formulation. Computes colatitudes and
        azimuths, then integrates (1-cos(colat)) * daz. Fast but may
        have precision issues for very small or irregular polygons.

    Algorithm 1 (L'Huilier):
        Triangulates polygon from first vertex, computes spherical
        excess for each triangle using L'Huilier's formula, then sums.
        Good accuracy for convex polygons.

    Algorithm 2 (Karney):
        Uses Karney's spherical polygon area formula from JPL.
        Most robust and accurate, especially for edge cases like
        polygons crossing poles or anti-meridian.

    See Also
    --------
    calculate_polygon_file_area : Calculate area from GeoJSON file
    calculate_spherical_triangle_area : Area of single spherical triangle
    spherical_polygon_area : Karney's algorithm implementation

    References
    ----------
    .. [1] Karney, C.F.F. "Algorithms for geodesics", 2013.
           https://trs.jpl.nasa.gov/handle/2014/41271
    .. [2] https://mathworld.wolfram.com/SphericalPolygon.html
    """
    # Convert to numpy arrays for easier manipulation
    aLongitude_in = np.asarray(aLongitude_in)
    aLatitude_in = np.asarray(aLatitude_in)

    # Validate inputs
    npoint = len(aLongitude_in)
    if npoint < 3:
        raise ValueError(
            f"A polygon requires at least 3 points. Got {npoint} point(s)."
        )

    if len(aLatitude_in) != npoint:
        raise ValueError(
            f"Longitude and latitude arrays must have same length. "
            f"Got {npoint} longitudes and {len(aLatitude_in)} latitudes."
        )

    if len(aLatitude_in) != npoint:
        raise ValueError(
            f"Longitude and latitude arrays must have same length. "
            f"Got {npoint} longitudes and {len(aLatitude_in)} latitudes."
        )

    # Close polygon if not already closed
    if aLatitude_in[0] != aLatitude_in[-1] or aLongitude_in[0] != aLongitude_in[-1]:
        aLatitude_in = np.append(aLatitude_in, aLatitude_in[0])
        aLongitude_in = np.append(aLongitude_in, aLongitude_in[0])
        npoint = len(aLongitude_in)

    # Convert to radians if needed
    if iFlag_radian:
        aLongitude_radian_in = aLongitude_in
        aLatitude_radian_in = aLatitude_in
    else:
        aLongitude_radian_in = np.radians(aLongitude_in)
        aLatitude_radian_in = np.radians(aLatitude_in)

    # Optional: Check for degenerate line-like polygons
    if dLine_threshold is not None:
        aLength = np.zeros(npoint - 1)
        for i in range(npoint - 1):
            dLength = calculate_distance_based_on_longitude_latitude(
                aLongitude_in[i],
                aLatitude_in[i],
                aLongitude_in[i + 1],
                aLatitude_in[i + 1],
            )
            aLength[i] = dLength

        dLength_max = np.max(aLength)
        dlength_rest = np.sum(aLength) - dLength_max

        # Check if polygon is too narrow (close to a line)
        if dlength_rest > 0 and (dLength_max / dlength_rest) > (1 - dLine_threshold):
            # Degenerate polygon - essentially a line
            return 0.0

    # Algorithm selection
    if iFlag_algorithm == 0:
        # Algorithm 0: Green's Theorem Line Integral
        # Get colatitude (surface distance as angle from origin)
        a = (
            np.sin(aLatitude_radian_in / 2) ** 2
            + np.cos(aLatitude_radian_in) * np.sin(aLongitude_radian_in / 2) ** 2
        )
        colat = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

        # Azimuth of each point from arbitrary origin
        az = np.arctan2(
            np.cos(aLatitude_radian_in) * np.sin(aLongitude_radian_in),
            np.sin(aLatitude_radian_in),
        ) % (2 * np.pi)

        # Calculate azimuth differences
        daz = np.diff(az)
        daz = (daz + np.pi) % (2 * np.pi) - np.pi

        # Average surface distance for each step
        deltas = np.diff(colat) / 2
        colat = colat[0:-1] + deltas

        # Integral: (1 - cos(colat)) * daz
        integrands = (1 - np.cos(colat)) * daz

        # Sum and take absolute value
        area = abs(np.sum(integrands))

        # Could be area of inside or outside - take minimum
        area = min(area, 1 - area)

    elif iFlag_algorithm == 1:
        # Algorithm 1: L'Huilier's Theorem (Spherical Triangulation)
        dLongitude_root = aLongitude_radian_in[0]
        dLatitude_root = aLatitude_radian_in[0]
        dLongitude_b = aLongitude_radian_in[1]
        dLatitude_b = aLatitude_radian_in[1]

        nTriangle = npoint - 2
        aArea = np.zeros(nTriangle)

        for i in np.arange(1, nTriangle + 1, 1):
            # Define triangle vertices
            dLongitude_a = dLongitude_b
            dLatitude_a = dLatitude_b
            dLongitude_b = aLongitude_radian_in[i + 1]
            dLatitude_b = aLatitude_radian_in[i + 1]

            # Calculate spherical triangle area
            aLongitude_temp = [dLongitude_root, dLongitude_a, dLongitude_b]
            aLatitude_temp = [dLatitude_root, dLatitude_a, dLatitude_b]
            dArea_triangle = calculate_spherical_triangle_area(
                aLongitude_temp, aLatitude_temp, iFlag_radian=True
            )
            aArea[i - 1] = dArea_triangle

        area = np.sum(aArea)

    elif iFlag_algorithm == 2:
        # Algorithm 2: Karney's Method (JPL/NASA - Most Accurate)
        if dRadius_in is not None:
            dRadius = dRadius_in
        else:
            dRadius = earth_radius

        dArea_m = spherical_polygon_area(
            aLatitude_radian_in, aLongitude_radian_in, dRadius
        )
        return float(dArea_m)

    else:
        raise ValueError(
            f"Invalid algorithm selection: {iFlag_algorithm}. "
            "Valid options are 0 (Green), 1 (L'Huilier), or 2 (Karney)."
        )

    # Convert from fraction of sphere to square meters
    if iFlag_radian:
        # Return as fraction of sphere
        return float(area)
    else:
        # Convert to square meters
        if dRadius_in is not None:
            dArea_m = area * dRadius_in**2
        else:
            dArea_m = area * earth_radius**2

        return float(dArea_m)


def calculate_polygon_file_area(sFilename_polygon_in: str) -> float:
    """Calculate total area of all polygons in a GeoJSON file.

    Reads a GeoJSON file and computes the sum of areas for all polygon
    features. Supports both POLYGON and MULTIPOLYGON geometry types.

    Parameters
    ----------
    sFilename_polygon_in : str
        Path to the GeoJSON file containing polygon geometries.

    Returns
    -------
    float
        Total area of all polygons in square meters.
        Returns 0.0 if no valid polygons are found.

    Raises
    ------
    FileNotFoundError
        If the specified file cannot be opened.

    Notes
    -----
    - Only POLYGON and MULTIPOLYGON geometries are processed
    - Empty geometries are skipped
    - Uses algorithm 2 (Karney) for all area calculations
    - Requires GDAL/OGR to be installed

    Examples
    --------
    >>> # Calculate total area from a GeoJSON file
    >>> area = calculate_polygon_file_area('watersheds.geojson')
    >>> print(f"Total area: {area / 1e6:.2f} km²")

    See Also
    --------
    calculate_polygon_area : Calculate area of single polygon
    get_geometry_coordinates : Extract coordinates from OGR geometry
    """

    pDriver = ogr.GetDriverByName("GeoJSON")
    pDataSource = pDriver.Open(sFilename_polygon_in, 0)

    if pDataSource is None:
        raise FileNotFoundError(
            f"Could not open file: {sFilename_polygon_in}. "
            "Please check the file path and ensure it's a valid GeoJSON file."
        )

    pLayer = pDataSource.GetLayer()
    dArea = 0.0

    pFeature = pLayer.GetNextFeature()
    while pFeature:
        pGeometry = pFeature.GetGeometryRef()

        if pGeometry is None:
            pFeature = pLayer.GetNextFeature()
            continue

        if pGeometry.IsEmpty():
            pFeature = pLayer.GetNextFeature()
            continue

        sGeometryType = pGeometry.GetGeometryName()

        if sGeometryType == "POLYGON":
            aCoords_gcs = get_geometry_coordinates(pGeometry)
            dArea += calculate_polygon_area(aCoords_gcs[:, 0], aCoords_gcs[:, 1])
        elif sGeometryType == "MULTIPOLYGON":
            for i in range(pGeometry.GetGeometryCount()):
                pGeometry_temp = pGeometry.GetGeometryRef(i)
                aCoords_gcs = get_geometry_coordinates(pGeometry_temp)
                dArea += calculate_polygon_area(aCoords_gcs[:, 0], aCoords_gcs[:, 1])
        # Skip unsupported geometry types

        pFeature = pLayer.GetNextFeature()

    pDataSource = None
    return dArea


def spherical_polygon_area(
    lat: Union[list, np.ndarray], lon: Union[list, np.ndarray], r: float
) -> float:
    """Calculate area of spherical polygon using Karney's method.

    This is the most accurate algorithm for spherical polygon area calculation,
    based on research from JPL/NASA. It handles edge cases well including
    polygons that cross the anti-meridian or contain poles.

    The algorithm computes the spherical excess using a robust formulation
    based on L'Huilier's theorem applied to each edge's spherical triangle.

    Parameters
    ----------
    lat : list or np.ndarray
        Latitudes of all vertices in radians.
        Must be in range [-π/2, π/2].
    lon : list or np.ndarray
        Longitudes of all vertices in radians.
        Must be in range [-π, π] or [0, 2π].
    r : float
        Spherical radius in meters (typically Earth's mean radius ≈ 6371229 m).

    Returns
    -------
    float
        Area of the spherical polygon in square units of r.
        Always returns a positive value.

    Notes
    -----
    - Input coordinates must be in radians
    - Polygon should be closed (first point = last point)
    - The algorithm is numerically stable even for very small or large polygons
    - Handles degenerate cases (e.g., consecutive identical vertices)

    Algorithm
    ---------
    For each edge of the polygon:
    1. Compute the spherical triangle formed by:
       - The two edge vertices
       - The North pole (or origin point)
    2. Calculate the spherical excess using L'Huilier's formula:
       E = 4 * arctan(sqrt(tan(s/2) * tan((s-a)/2) * tan((s-b)/2) * tan((s-c)/2)))
       where s = (a+b+c)/2 is the semi-perimeter
    3. Sum the signed excesses (sign depends on edge direction)
    4. Multiply by r² to get area

    Examples
    --------
    >>> # Spherical cap at North pole (small circle at 45°N)
    >>> lats = np.radians([45, 45, 45, 45])
    >>> lons = np.radians([0, 90, 180, 270])
    >>> r = 6371229.0  # Earth's mean radius
    >>> area = spherical_polygon_area(lats, lons, r)

    See Also
    --------
    calculate_polygon_area : High-level interface with algorithm selection
    haversine : Haversine function used in this algorithm

    References
    ----------
    .. [1] Karney, C.F.F. "Algorithms for geodesics", Journal of Geodesy, 2013.
           https://trs.jpl.nasa.gov/handle/2014/41271
    .. [2] Beyer, W.H. "CRC Standard Mathematical Tables", 28th ed.
    """
    lam1 = lam2 = beta1 = beta2 = cosB1 = cosB2 = 0.0
    hav = 0.0
    sum_excess = 0.0

    for j in range(len(lat)):
        k = (j + 1) % len(lat)

        if j == 0:
            lam1 = lon[j]
            beta1 = lat[j]
            lam2 = lon[j + 1]
            beta2 = lat[j + 1]
            cosB1 = math.cos(beta1)
            cosB2 = math.cos(beta2)
        else:
            # Reuse previous endpoint as new starting point
            lam1 = lam2
            beta1 = beta2
            lam2 = lon[k]
            beta2 = lat[k]
            cosB1 = cosB2
            cosB2 = math.cos(beta2)

        # Skip if edge has zero length (same longitude)
        if lam1 != lam2:
            # Calculate haversine distance between the two points
            hav = haversine(beta2 - beta1) + cosB1 * cosB2 * haversine(lam2 - lam1)
            a = 2 * math.asin(math.sqrt(hav))

            # Compute complementary latitudes (colatitudes)
            b = math.pi / 2 - beta2
            c = math.pi / 2 - beta1

            # Semi-perimeter of spherical triangle
            s = 0.5 * (a + b + c)

            # L'Huilier's formula for spherical excess
            # tan(E/4) = sqrt(tan(s/2) * tan((s-a)/2) * tan((s-b)/2) * tan((s-c)/2))
            t = (
                math.tan(s / 2)
                * math.tan((s - a) / 2)
                * math.tan((s - b) / 2)
                * math.tan((s - c) / 2)
            )

            excess = abs(4 * math.atan(math.sqrt(abs(t))))

            # Sign of excess depends on edge direction
            if lam2 < lam1:
                excess = -excess

            sum_excess += excess

    # Area = |spherical excess| * r²
    return abs(sum_excess) * r * r


def convert_360_to_180(aLongitude_in):
    """
    Convert longitude values from [0, 360] range to [-180, 180] range.

    Args:
        aLongitude_in: Scalar or array-like of longitude values in [0, 360]

    Returns:
        numpy array or scalar with longitudes in [-180, 180]
    """
    aLon = np.array(aLongitude_in, dtype=float)
    aLon = np.where(aLon > 180.0, aLon - 360.0, aLon)
    return aLon


def split_international_date_line_polygon_coordinates(aCoords_gcs):
    """
    Split a polygon that crosses the International Date Line into two parts.

    When a polygon has vertices spanning more than 180 degrees of longitude
    (i.e., it crosses the antimeridian at ±180°), this function splits it
    into two sub-polygons on either side of the date line.

    Args:
        aCoords_gcs (numpy.ndarray): Nx2 array of (longitude, latitude) coordinates
                                     in [-180, 180] range.

    Returns:
        list of numpy.ndarray: List of coordinate arrays for each sub-polygon.
    """
    aCoords = np.array(aCoords_gcs, dtype=float)
    lons = aCoords[:, 0]
    lats = aCoords[:, 1]
    n = len(lons)

    # Separate into west (lon < 0) and east (lon >= 0) parts
    west_coords = []
    east_coords = []

    for i in range(n):
        lon = lons[i]
        lat = lats[i]
        if lon < 0:
            west_coords.append([lon, lat])
            # Also add a mirrored point on the east side at the date line
        else:
            east_coords.append([lon, lat])

    # Build two polygons: one on the west side (lon < 0 → shift to > 180)
    # and one on the east side (lon >= 0)
    # Strategy: split at the antimeridian by adjusting coordinates
    west_part = []
    east_part = []

    for i in range(n):
        lon = lons[i]
        lat = lats[i]
        if lon < 0:
            # This vertex is on the west side; for the west polygon keep as-is
            # For the east polygon, shift to positive (lon + 360)
            west_part.append([lon, lat])
        else:
            # This vertex is on the east side
            east_part.append([lon, lat])

    # Close the polygons if they have enough points
    result = []
    if len(west_part) >= 3:
        west_arr = np.array(west_part, dtype=float)
        result.append(west_arr)
    if len(east_part) >= 3:
        east_arr = np.array(east_part, dtype=float)
        result.append(east_arr)

    # Fallback: if splitting didn't produce valid sub-polygons, return original
    if len(result) == 0:
        result = [aCoords]

    return result



def calculate_signed_area_shoelace(coords: np.ndarray) -> float:
    x, y = coords[:, 0], coords[:, 1]

    # Vectorized shoelace formula - more efficient than loops
    # Handle the wrap-around (last point to first point) implicitly
    x_rolled = np.roll(x, -1)  # x[i+1] for all i, with wrap-around
    y_rolled = np.roll(y, -1)  # y[i+1] for all i, with wrap-around

    signed_area = 0.5 * np.sum(x * y_rolled - x_rolled * y)
    return signed_area

def remove_duplicate_closure(coords: NDArray[np.floating]) -> NDArray[np.floating]:
    """Remove duplicated closing vertex when first and last are equal."""
    if len(coords) > 1 and np.allclose(coords[0], coords[-1], atol=1.0e-12):
        return coords[:-1]
    return coords

def calculate_signed_area_spherical_polar(
    coords: np.ndarray, pole: str = "north"
) -> float:
    """Calculate signed spherical area for polygons enclosing a pole.

    The polygon is projected with Lambert azimuthal equal-area (LAEA) centered
    on the requested pole, then planar signed area is computed with the
    shoelace formula. With unit sphere radius, output area is in steradians.

    Parameters
    ----------
    coords : np.ndarray
        Polygon coordinates as (lon, lat) degrees.
    pole : str, default="north"
        Polar center for projection, either "north" or "south".

    Returns
    -------
    float
        Signed spherical area (steradians on unit sphere).
    """
    if not isinstance(coords, np.ndarray) or coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError("coords must be a 2D numpy array with shape (n, 2)")

    if len(coords) < 3:
        return 0.0

    arr = remove_duplicate_closure(coords)
    if len(arr) < 3:
        return 0.0

    lon_rad = np.deg2rad(arr[:, 0])
    lat_rad = np.deg2rad(arr[:, 1])

    cos_lat = np.cos(lat_rad)
    sin_lat = np.sin(lat_rad)

    pole_lc = pole.lower()
    if pole_lc == "north":
        # LAEA centered at +90 deg latitude
        denom = 1.0 + sin_lat
        # Avoid divide-by-zero at antipode (not expected for pole-enclosing cells)
        denom = np.maximum(denom, 1.0e-15)
        k = np.sqrt(2.0 / denom)
        x = k * cos_lat * np.sin(lon_rad)
        y = -k * cos_lat * np.cos(lon_rad)
    elif pole_lc == "south":
        # LAEA centered at -90 deg latitude
        denom = 1.0 - sin_lat
        denom = np.maximum(denom, 1.0e-15)
        k = np.sqrt(2.0 / denom)
        x = k * cos_lat * np.sin(lon_rad)
        y = k * cos_lat * np.cos(lon_rad)
    else:
        raise ValueError("pole must be either 'north' or 'south'")

    projected = np.column_stack((x, y))
    return calculate_signed_area_shoelace(projected)


def _point_on_segment_2d(
    point: NDArray[np.floating],
    seg_start: NDArray[np.floating],
    seg_end: NDArray[np.floating],
    tol: float = 1.0e-12,
) -> bool:
    """Return True if a 2D point lies on a line segment within tolerance."""
    px, py = point
    x1, y1 = seg_start
    x2, y2 = seg_end

    dx = x2 - x1
    dy = y2 - y1

    # Degenerate segment
    if abs(dx) < tol and abs(dy) < tol:
        return np.hypot(px - x1, py - y1) <= tol

    # Cross-product distance to the infinite line
    cross = (px - x1) * dy - (py - y1) * dx
    if abs(cross) > tol:
        return False

    # Dot-product bounds check for segment extents
    dot = (px - x1) * dx + (py - y1) * dy
    if dot < -tol:
        return False

    seg_len_sq = dx * dx + dy * dy
    if dot - seg_len_sq > tol:
        return False

    return True

def _point_in_polygon_2d(
    point: NDArray[np.floating],
    polygon: NDArray[np.floating],
    include_boundary: bool = False,
    tol: float = 1.0e-12,
) -> bool:
    """2D point-in-polygon using ray casting."""
    poly = remove_duplicate_closure(polygon)
    if len(poly) < 3:
        return False

    x, y = point

    # Boundary check first
    for i in range(len(poly)):
        p1 = poly[i]
        p2 = poly[(i + 1) % len(poly)]
        if _point_on_segment_2d(point, p1, p2, tol=tol):
            return include_boundary

    inside = False
    for i in range(len(poly)):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % len(poly)]

        intersects = (y1 > y) != (y2 > y)
        if intersects:
            x_intersect = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
            if x_intersect > x:
                inside = not inside

    return inside

def polygon_includes_pole(
    coords: NDArray[np.floating],
    pole: str = "south",
    include_boundary: bool = False,
    tol: float = 1.0e-12,
) -> bool:
    """Check whether a polygon includes the requested pole in its interior.

    The test projects lon/lat vertices into a local polar plane where the pole
    maps to the origin, then runs a standard 2D point-in-polygon query.
    """
    if coords is None:
        return False

    arr = np.asarray(coords, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 2 or len(arr) < 3:
        return False

    arr = remove_duplicate_closure(arr)
    if len(arr) < 3:
        return False

    lons = arr[:, 0]
    lats = arr[:, 1]

    lon_rad = np.deg2rad(lons)
    pole_lc = pole.lower()
    if pole_lc == "south":
        radial = np.maximum(0.0, 90.0 + lats)
    elif pole_lc == "north":
        radial = np.maximum(0.0, 90.0 - lats)
    else:
        raise ValueError("pole must be either 'south' or 'north'")

    projected = np.column_stack((radial * np.cos(lon_rad), radial * np.sin(lon_rad)))
    origin = np.array([0.0, 0.0])

    return _point_in_polygon_2d(
        origin, projected, include_boundary=include_boundary, tol=tol
    )

def unwrap_longitudes(coords: np.ndarray) -> np.ndarray:
    """Unwrap longitude coordinates to handle International Date Line crossings.

    Adjusts longitude values to ensure they are all within 180 degrees of the
    first coordinate, preventing artificial jumps when calculating polygon areas.

    Parameters
    ----------
    coords : np.ndarray
        Array of shape (n, 2) representing polygon coordinates in (longitude, latitude) format.

    Returns
    -------
    np.ndarray
        Array with unwrapped longitude coordinates, same shape as input.

    Notes
    -----
    This function uses the first longitude point as a reference and adjusts all
    subsequent longitudes to be within ±180 degrees of this reference point by
    adding or subtracting 360 degrees as needed.

    Examples
    --------
    >>> coords = np.array([[-170, 10], [170, 20], [-160, 30]])
    >>> unwrapped = unwrap_longitudes(coords)
    >>> # Longitudes adjusted to avoid large jumps across IDL
    """
    if not isinstance(coords, np.ndarray) or coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError("coords must be a 2D numpy array with shape (n, 2)")

    coords_unwrapped = coords.copy()
    lons = coords_unwrapped[:, 0]
    ref_lon = lons[0]

    # Vectorized approach for better performance
    diff = lons - ref_lon
    # Adjust longitudes that are more than 180 degrees away
    lons[diff > 180] -= 360
    lons[diff < -180] += 360

    return coords_unwrapped



def check_cross_international_date_line_polygon(
    coords: np.ndarray,
) -> Tuple[bool, np.ndarray]:
    """Check if polygon coordinates cross the International Date Line.

    This function distinguishes between actual IDL edge crossings and vertices
    that happen to lie exactly on the IDL meridian (±180°). Vertices on the IDL
    without edge crossings are adjusted slightly to avoid numerical issues.

    Parameters
    ----------
    coords : np.ndarray
        Array of polygon coordinates in (longitude, latitude) format.

    Returns
    -------
    Tuple[bool, np.ndarray]
        A tuple containing:
        - bool: True if the polygon has actual edges crossing the IDL, False otherwise
        - np.ndarray or None:
          - If False: Adjusted coordinates with IDL vertices moved slightly off the meridian
          - If True: None (original coordinates should be used for splitting)

    Notes
    -----
    **IDL Detection Logic:**

    1. **Vertex Detection**: Identifies vertices exactly on ±180° meridian
    2. **Edge Crossing Detection**: Finds edges that actually cross the IDL
       (excludes edges involving IDL vertices)
    3. **Coordinate Adjustment**: If no edge crossings but IDL vertices exist,
       moves these vertices slightly (±1e-6°) based on neighboring vertices
    4. **Hemisphere Assignment**: IDL vertices are moved to the hemisphere
       containing the majority of their neighbors

    **Return Behavior:**

    - Returns (False, adjusted_coords) for polygons with IDL vertices but no crossings
    - Returns (True, None) for polygons with actual IDL edge crossings
    - Returns (False, coords) for polygons that don't involve the IDL

    This approach ensures that single IDL vertices don't cause false crossing
    detection while preserving the ability to split truly crossing polygons.
    """
    try:
        if (
            not isinstance(coords, np.ndarray)
            or coords.ndim != 2
            or coords.shape[1] != 2
        ):
            raise ValueError("coords must be a 2D numpy array with shape (n, 2)")

        if len(coords) < 3:
            return False, coords.copy()

        # Validate coordinate ranges
        lons, lats = coords[:, 0], coords[:, 1]
        if not (
            np.all((-180 <= lons) & (lons <= 180))
            and np.all((-90 <= lats) & (lats <= 90))
        ):
            return False, coords.copy()

    except (ValueError, IndexError, TypeError):
        return False, coords.copy()  # Invalid coordinates can't cross IDL

    # Make a copy to avoid modifying the original
    coords_updated = coords.copy()
    lons = coords_updated[:, 0]

    # Check if any vertices are exactly on the IDL
    # Exclude the last point to avoid duplication in closed polygons (first == last)
    # Use module-level IDL_TOLERANCE and IDL_OFFSET constants
    idl_vertices = np.zeros_like(lons, dtype=bool)
    idl_vertices[:-1] = np.abs(np.abs(lons[:-1]) - 180.0) < IDL_TOLERANCE

    # Track hemisphere support using only non-IDL vertices.
    non_idl_lons = lons[~idl_vertices]
    has_eastern = np.any((non_idl_lons > 0) & (non_idl_lons < 180.0))
    has_western = np.any((non_idl_lons < 0) & (non_idl_lons > -180.0))
    spans_both_hemispheres = has_eastern and has_western

    idl_touch_positive = np.any(np.abs(lons[idl_vertices] - 180.0) < IDL_TOLERANCE)
    idl_touch_negative = np.any(np.abs(lons[idl_vertices] + 180.0) < IDL_TOLERANCE)
    touches_both_idl_sides = idl_touch_positive and idl_touch_negative

    # A polygon enclosing either pole necessarily crosses the IDL.
    if polygon_includes_pole(coords_updated, pole="north") or polygon_includes_pole(
        coords_updated, pole="south"
    ):
        return True, None

    # If there are vertices on IDL, check if there are actual edge crossings
    if np.any(idl_vertices):
        lons_next = np.roll(lons, -1)

        # Find eastward crossings: positive to negative longitude, excluding IDL vertices
        eastward_crossings = (
            (lons > 0)
            & (lons < 180.0)
            & (lons_next < 0)
            & ~idl_vertices
            & ~np.roll(idl_vertices, -1)
        )

        # Find westward crossings: negative to positive longitude, excluding IDL vertices
        westward_crossings = (
            (lons < 0) & (lons_next > 0) & ~idl_vertices & ~np.roll(idl_vertices, -1)
        )

        # If no actual edge crossings but polygon spans both hemispheres, it crosses IDL
        if touches_both_idl_sides or spans_both_hemispheres:
            # This is a true IDL crossing with vertices on the meridian
            return True, None
        elif not (np.any(eastward_crossings) or np.any(westward_crossings)):
            # Adjust vertices that are exactly on the IDL by slightly moving them
            idl_indices = np.where(idl_vertices)[0]
            for idx in idl_indices:
                # Determine which hemisphere to move to based on neighboring vertices
                prev_idx = (idx - 1) % len(coords_updated)
                next_idx = (idx + 1) % len(coords_updated)

                prev_lon = coords_updated[prev_idx, 0]
                next_lon = coords_updated[next_idx, 0]

                # Count non-IDL neighbors to determine preferred hemisphere
                neighbor_lons = []
                if abs(abs(prev_lon) - 180.0) > IDL_TOLERANCE:  # Previous vertex not on IDL
                    neighbor_lons.append(prev_lon)
                if abs(abs(next_lon) - 180.0) > IDL_TOLERANCE:  # Next vertex not on IDL
                    neighbor_lons.append(next_lon)

                if neighbor_lons:
                    # Move to hemisphere containing majority of neighbors
                    positive_neighbors = sum(1 for lon in neighbor_lons if lon > 0)
                    if positive_neighbors >= len(neighbor_lons) / 2:
                        # Move to eastern hemisphere
                        coords_updated[idx, 0] = 180.0 - IDL_OFFSET
                    else:
                        # Move to western hemisphere
                        coords_updated[idx, 0] = -180.0 + IDL_OFFSET
                else:
                    # If neighbors are also on IDL, check overall polygon distribution
                    positive_lons = np.sum(lons > 0)
                    negative_lons = np.sum(lons < 0)

                    if positive_lons >= negative_lons:
                        coords_updated[idx, 0] = 180.0 - IDL_OFFSET
                    else:
                        coords_updated[idx, 0] = -180.0 + IDL_OFFSET

            # After adjustment, no crossing
            return False, coords_updated
        else:
            # There are actual edge crossings
            return True, None
    else:
        # Original logic for cases without IDL vertices
        # Check for large jumps between consecutive points (> 180 degrees)
        lon_diffs = np.abs(np.diff(lons))
        max_jump = np.max(lon_diffs)

        # Also check the wrap-around from last to first point
        wrap_jump = abs(lons[-1] - lons[0])

        crossing = max_jump > 180 or wrap_jump > 180
        return crossing, None


def check_counter_clockwise(coords: np.ndarray) -> bool:
    """Check if polygon coordinates are in counter-clockwise order.

    Parameters
    ----------
    coords : np.ndarray
        Array of shape (n, 2) representing polygon coordinates.

    Returns
    -------
    bool
        True if vertices are in counter-clockwise order, False otherwise.
    """

    if not isinstance(coords, np.ndarray) or coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError("coords must be a 2D numpy array with shape (n, 2)")

    if len(coords) < 3:
        return True  # Degenerate case

    if polygon_includes_pole(coords, pole="north"):
        signed_area = calculate_signed_area_spherical_polar(coords, pole="north")
    elif polygon_includes_pole(coords, pole="south"):
        signed_area = calculate_signed_area_spherical_polar(coords, pole="south")
    else:
        # Check if polygon crosses the International Date Line
        iFlag_cross, _ = check_cross_international_date_line_polygon(coords)
        if iFlag_cross:
            coords_unwrapped = unwrap_longitudes(coords)
            # Calculate signed area using optimized shoelace formula
            signed_area = calculate_signed_area_shoelace(coords_unwrapped)
        else:
            # Standard case: calculate signed area directly
            signed_area = calculate_signed_area_shoelace(coords)
    return signed_area > 0

def get_geometry_coordinates(
    geometry: ogr.Geometry,
    enforce_ccw: bool = False,
) -> Union[np.ndarray, List[np.ndarray]]:
    """Extract coordinates from an OGR geometry object.

    Dispatches to appropriate handler based on geometry type. For polygon
    geometries, coordinates can optionally be enforced to be in
    counter-clockwise (CCW) order.

    Parameters
    ----------
    geometry : ogr.Geometry
        OGR geometry object from which to extract coordinates.
    enforce_ccw : bool, default=False
        If True, polygon and multipolygon exterior rings are returned in
        counter-clockwise order. If False, coordinates are returned in the
        original geometry order.

    Returns
    -------
    np.ndarray or List[np.ndarray]
        For POINT, LINESTRING, POLYGON, LINEARRING: numpy array of shape (n, 2)
        For MULTIPOLYGON: list of numpy arrays, one per polygon part

    Raises
    ------
    ValueError
        If geometry is None, invalid, or of unsupported type.
    RuntimeError
        If coordinate extraction fails.

    Notes
    -----
    Supported geometry types:
    - POINT: Returns single point as array of shape (1, 2)
    - LINESTRING: Returns array of line vertices
    - LINEARRING: Returns array of ring vertices
    - POLYGON: Returns exterior ring coordinates in CCW order
    - MULTIPOLYGON: Returns list of exterior rings, each in CCW order

    Currently unsupported (will raise ValueError):
    - MULTIPOINT, MULTILINESTRING, GEOMETRYCOLLECTION

    Z-coordinates (3D) are ignored; only X and Y are extracted.

    Examples
    --------
    >>> from osgeo import ogr
    >>> # Create a point geometry
    >>> point = ogr.Geometry(ogr.wkbPoint)
    >>> point.AddPoint(10.0, 20.0)
    >>> coords = get_geometry_coordinates(point)
    >>> coords.shape
    (1, 2)

    See Also
    --------
    get_polygon_exterior_coords : Extract polygon exterior ring
    get_multipolygon_exterior_coords : Extract multipolygon exterior rings
    """
    # Validate geometry
    if geometry is None:
        raise ValueError("Geometry object cannot be None")

    try:
        geometry_type = geometry.GetGeometryName()
    except Exception as e:
        raise ValueError(f"Invalid geometry object: {e}")

    # Dispatch to appropriate handler
    if geometry_type == "POINT":
        return get_point_coords(geometry)
    elif geometry_type == "LINESTRING":
        return get_linestring_coords(geometry)
    elif geometry_type == "POLYGON":
        return get_polygon_exterior_coords(geometry, enforce_ccw=enforce_ccw)
    elif geometry_type == "LINEARRING":
        return get_linearring_coords(geometry)
    elif geometry_type == "MULTIPOLYGON":
        return get_multipolygon_exterior_coords(geometry, enforce_ccw=enforce_ccw)
    else:
        raise ValueError(
            f"Unsupported geometry type: '{geometry_type}'. "
            "Supported types: POINT, LINESTRING, POLYGON, LINEARRING, MULTIPOLYGON"
        )


def get_polygon_exterior_coords(
    polygon_geometry: ogr.Geometry, enforce_ccw: bool = False
) -> np.ndarray:
    """Extract exterior ring coordinates from a polygon.

    Parameters
    ----------
    polygon_geometry : ogr.Geometry
        OGR Polygon geometry object.
    enforce_ccw : bool, default=False
        If True, reverse the extracted exterior ring when needed so the
        returned coordinates are counter-clockwise.

    Returns
    -------
    np.ndarray
        Array of shape (n, 2) containing exterior ring coordinates.

    Raises
    ------
    ValueError
        If geometry is not a valid polygon or has no exterior ring.
    RuntimeError
        If coordinate extraction fails.

    Notes
    -----
    - Only extracts the exterior ring; holes/interior rings are ignored
    - Optionally reverses coordinate order if clockwise
    - Returns (longitude, latitude) or (x, y) coordinate pairs

    Examples
    --------
    >>> from osgeo import ogr
    >>> # Create a square polygon
    >>> ring = ogr.Geometry(ogr.wkbLinearRing)
    >>> ring.AddPoint(0, 0)
    >>> ring.AddPoint(1, 0)
    >>> ring.AddPoint(1, 1)
    >>> ring.AddPoint(0, 1)
    >>> ring.AddPoint(0, 0)
    >>> polygon = ogr.Geometry(ogr.wkbPolygon)
    >>> polygon.AddGeometry(ring)
    >>> coords = get_polygon_exterior_coords(polygon)
    >>> coords.shape[1]
    2
    """
    if polygon_geometry is None:
        raise ValueError("Polygon geometry cannot be None")

    try:
        # Get the exterior ring (index 0)
        ring = polygon_geometry.GetGeometryRef(0)
        if ring is None:
            raise ValueError("Polygon has no exterior ring")

        n_points = ring.GetPointCount()
        if n_points == 0:
            raise ValueError("Polygon exterior ring has no points")

        # Extract coordinates
        exterior_coords = []
        for i in range(n_points):
            point = ring.GetPoint(i)
            exterior_coords.append((point[0], point[1]))

    except Exception as e:
        raise RuntimeError(f"Failed to extract polygon coordinates: {e}")

    # Convert to numpy array
    coords_array = np.array(exterior_coords)

    # Ensure counter-clockwise order when requested
    if enforce_ccw and not check_counter_clockwise(coords_array):
        coords_array = coords_array[::-1]

    return coords_array


def get_multipolygon_exterior_coords(
    multipolygon_geometry: ogr.Geometry,
    enforce_ccw: bool = False,
) -> List[np.ndarray]:
    """Extract exterior ring coordinates from all parts of a multipolygon.

    Each polygon part's coordinates can optionally be enforced to be in
    counter-clockwise order.

    Parameters
    ----------
    multipolygon_geometry : ogr.Geometry
        OGR MultiPolygon geometry object.
    enforce_ccw : bool, default=False
        If True, reverse each exterior ring when needed so returned
        coordinates are counter-clockwise.

    Returns
    -------
    List[np.ndarray]
        List of numpy arrays, one per polygon part. Each array has shape (n, 2)
        and contains coordinates in counter-clockwise order.

    Raises
    ------
    ValueError
        If geometry is not a valid multipolygon or has no parts.
    RuntimeError
        If coordinate extraction fails for any part.

    Notes
    -----
    - Only extracts exterior rings; holes/interior rings are ignored
    - Each polygon part can be independently enforced to CCW order
    - Empty parts are skipped with a warning

    Examples
    --------
    >>> from osgeo import ogr
    >>> # Create a multipolygon with 2 parts
    >>> multi = ogr.Geometry(ogr.wkbMultiPolygon)
    >>> # Add first polygon
    >>> ring1 = ogr.Geometry(ogr.wkbLinearRing)
    >>> ring1.AddPoint(0, 0)
    >>> ring1.AddPoint(1, 0)
    >>> ring1.AddPoint(1, 1)
    >>> ring1.AddPoint(0, 0)
    >>> poly1 = ogr.Geometry(ogr.wkbPolygon)
    >>> poly1.AddGeometry(ring1)
    >>> multi.AddGeometry(poly1)
    >>> coords_list = get_multipolygon_exterior_coords(multi)
    >>> len(coords_list)
    1
    """
    if multipolygon_geometry is None:
        raise ValueError("MultiPolygon geometry cannot be None")

    # Validate geometry type
    try:
        geometry_type = multipolygon_geometry.GetGeometryName()
        if geometry_type != "MULTIPOLYGON":
            raise ValueError(f"Expected MULTIPOLYGON geometry, got '{geometry_type}'")
    except Exception as e:
        raise ValueError(f"Invalid geometry object: {e}")

    try:
        n_parts = multipolygon_geometry.GetGeometryCount()
        if n_parts == 0:
            raise ValueError("MultiPolygon has no parts")

        exterior_coords_list = []

        for i in range(n_parts):
            # Get polygon part
            polygon = multipolygon_geometry.GetGeometryRef(i)
            if polygon is None:
                continue  # Skip invalid parts

            # Get exterior ring
            ring = polygon.GetGeometryRef(0)
            if ring is None:
                continue  # Skip polygons without rings

            n_points = ring.GetPointCount()
            if n_points == 0:
                continue  # Skip empty rings

            # Extract coordinates for this polygon part
            part_coords = []
            for j in range(n_points):
                point = ring.GetPoint(j)
                part_coords.append((point[0], point[1]))

            # Convert to numpy array
            coords_array = np.array(part_coords)

            # Ensure counter-clockwise order when requested
            if enforce_ccw and not check_counter_clockwise(coords_array):
                coords_array = coords_array[::-1]

            exterior_coords_list.append(coords_array)

    except Exception as e:
        raise RuntimeError(f"Failed to extract multipolygon coordinates: {e}")

    if len(exterior_coords_list) == 0:
        raise ValueError("No valid polygon parts found in multipolygon")

    return exterior_coords_list


def get_linestring_coords(linestring_geometry: ogr.Geometry) -> np.ndarray:
    """Extract coordinates from a linestring geometry.

    Parameters
    ----------
    linestring_geometry : ogr.Geometry
        OGR LineString geometry object.

    Returns
    -------
    np.ndarray
        Array of shape (n, 2) containing linestring vertex coordinates.

    Raises
    ------
    ValueError
        If geometry is not a valid linestring or has no points.
    RuntimeError
        If coordinate extraction fails.

    Examples
    --------
    >>> from osgeo import ogr
    >>> line = ogr.Geometry(ogr.wkbLineString)
    >>> line.AddPoint(0, 0)
    >>> line.AddPoint(1, 1)
    >>> coords = get_linestring_coords(line)
    >>> coords.shape
    (2, 2)
    """
    if linestring_geometry is None:
        raise ValueError("LineString geometry cannot be None")

    try:
        n_points = linestring_geometry.GetPointCount()
        if n_points == 0:
            raise ValueError("LineString has no points")

        # Extract coordinates using list comprehension
        coords = [
            (linestring_geometry.GetPoint(i)[0], linestring_geometry.GetPoint(i)[1])
            for i in range(n_points)
        ]

    except Exception as e:
        raise RuntimeError(f"Failed to extract linestring coordinates: {e}")

    return np.array(coords)


def get_point_coords(point_geometry: ogr.Geometry) -> np.ndarray:
    """Extract coordinates from a point geometry.

    Parameters
    ----------
    point_geometry : ogr.Geometry
        OGR Point geometry object.

    Returns
    -------
    np.ndarray
        Array of shape (1, 2) containing the point coordinates.

    Raises
    ------
    ValueError
        If geometry is not a valid point.
    RuntimeError
        If coordinate extraction fails.

    Examples
    --------
    >>> from osgeo import ogr
    >>> point = ogr.Geometry(ogr.wkbPoint)
    >>> point.AddPoint(10.5, 20.3)
    >>> coords = get_point_coords(point)
    >>> coords.shape
    (1, 2)
    >>> coords[0]
    array([10.5, 20.3])
    """
    if point_geometry is None:
        raise ValueError("Point geometry cannot be None")

    try:
        point = point_geometry.GetPoint()
        if point is None or len(point) < 2:
            raise ValueError("Invalid point geometry")

    except Exception as e:
        raise RuntimeError(f"Failed to extract point coordinates: {e}")

    return np.array([(point[0], point[1])])


def get_linearring_coords(linearring_geometry: ogr.Geometry) -> np.ndarray:
    """Extract coordinates from a linear ring geometry.

    Parameters
    ----------
    linearring_geometry : ogr.Geometry
        OGR LinearRing geometry object.

    Returns
    -------
    np.ndarray
        Array of shape (n, 2) containing ring vertex coordinates.

    Raises
    ------
    ValueError
        If geometry is not a valid linear ring or has no points.
    RuntimeError
        If coordinate extraction fails.

    Notes
    -----
    - A linear ring is a closed linestring (first point = last point)
    - This function does not enforce CCW order (use for non-polygon contexts)

    Examples
    --------
    >>> from osgeo import ogr
    >>> ring = ogr.Geometry(ogr.wkbLinearRing)
    >>> ring.AddPoint(0, 0)
    >>> ring.AddPoint(1, 0)
    >>> ring.AddPoint(1, 1)
    >>> ring.AddPoint(0, 0)
    >>> coords = get_linearring_coords(ring)
    >>> coords.shape
    (4, 2)
    """
    if linearring_geometry is None:
        raise ValueError("LinearRing geometry cannot be None")

    try:
        n_points = linearring_geometry.GetPointCount()
        if n_points == 0:
            raise ValueError("LinearRing has no points")

        # Extract coordinates using list comprehension
        coords = [
            (linearring_geometry.GetPoint(i)[0], linearring_geometry.GetPoint(i)[1])
            for i in range(n_points)
        ]

    except Exception as e:
        raise RuntimeError(f"Failed to extract linear ring coordinates: {e}")

    return np.array(coords)
