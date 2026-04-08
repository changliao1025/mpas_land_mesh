"""
Utility function to find vertices that lie on a given edge.

Copied from pyflowline/algorithms/auxiliary/find_index_in_list.py and adapted
to remove external pyflowline dependencies. Uses rtree spatial indexing when
available, with a pure-Python fallback.
"""


import numpy as np
import logging
from rtree.index import Index as RTreeindex
from mpas_land_mesh.utilities.spatial_reference import reproject_coordinates

# NOTE: mpas_land_mesh.classes imports (pyvertex, pyedge, pyflowline) are done
# lazily inside each function to avoid circular imports:
#   object.py  ←  edge.py  ←  object.py  (via flowline.py → edge.py)


# Configure logger for this module
logger = logging.getLogger(__name__)

def find_vertex_on_edge(aVertex_in, pEdge_in):
    """
    Find vertices from aVertex_in that lie on pEdge_in.

    Args:
        aVertex_in (list): List of vertex objects with dLongitude_degree and
                           dLatitude_degree attributes.
        pEdge_in: Edge object with pVertex_start, pVertex_end, and
                  check_vertex_on_edge() method.

    Returns:
        tuple: (iFlag_exist, npoint, aIndex_order)
            - iFlag_exist (int): 1 if any vertex found on edge, 0 otherwise
            - npoint (int): Number of vertices found on edge
            - aIndex_order (list): Ordered indices of vertices on edge
                                   (ordered by distance from edge start)
    """
    iFlag_exist = 0
    aIndex = []
    aIndex_order = []
    aDistance = []
    nVertex = len(aVertex_in)
    npoint = 0

    if nVertex == 0:
        return iFlag_exist, npoint, aIndex_order

    pVertex_start = pEdge_in.pVertex_start
    pVertex_end = pEdge_in.pVertex_end
    x1 = pVertex_start.dLongitude_degree
    y1 = pVertex_start.dLatitude_degree
    x2 = pVertex_end.dLongitude_degree
    y2 = pVertex_end.dLatitude_degree


    # Use R-tree spatial index for efficient lookup
    index_vertex = RTreeindex()
    for i in range(nVertex):
        x = aVertex_in[i].dLongitude_degree
        y = aVertex_in[i].dLatitude_degree
        pBound = (x - 1e-5, y - 1e-5, x + 1e-5, y + 1e-5)
        index_vertex.insert(i, pBound)
    left = np.min([x1, x2])
    right = np.max([x1, x2])
    bottom = np.min([y1, y2])
    top = np.max([y1, y2])
    pBound = (left, bottom, right, top)
    aIntersect = list(index_vertex.intersection(pBound))
    for k in aIntersect:
        pVertex = aVertex_in[k]
        iFlag_overlap, dDistance, diff = pEdge_in.check_vertex_on_edge(pVertex)
        if iFlag_overlap == 1:
            iFlag_exist = 1
            aDistance.append(dDistance)
            aIndex.append(k)
            npoint += 1


    # Re-order by distance from edge start
    if iFlag_exist == 1:
        x = np.array(aDistance)
        b = np.argsort(x)
        c = np.array(aIndex)
        d = c[b]
        aIndex_order = list(d)

    return iFlag_exist, npoint, aIndex_order


def slerp(n1, n2, t):
    """
    Spherical linear interpolation between two n-vectors.

    Args:
        n1: First n-vector
        n2: Second n-vector
        t: Interpolation parameter (0 to 1)

    Returns:
        Interpolated n-vector
    """
    from mpas_land_mesh.classes.nvector import pynvector

    use_high_precision = getattr(n1, "use_high_precision", False)
    dtype = np.float128 if use_high_precision else np.float64
    eps = dtype(1e-18 if use_high_precision else 1e-12)

    dot_product = dtype(np.clip(n1.dot(n2), -1.0, 1.0))

    # Near-parallel vectors: lerp is stable and avoids division by tiny sin(theta).
    if dot_product > dtype(1.0) - eps:
        result = n1 * (dtype(1.0) - dtype(t)) + n2 * dtype(t)
        return result.normalize()

    # Near-antiparallel vectors: avoid lerp because midpoint can collapse to zero.
    if dot_product < -dtype(1.0) + eps:
        abs_x = abs(n1.dX)
        abs_y = abs(n1.dY)
        abs_z = abs(n1.dZ)

        # Build a deterministic orthogonal direction to n1.
        if abs_x <= abs_y and abs_x <= abs_z:
            ortho = pynvector(
                {"x": dtype(0.0), "y": -n1.dZ, "z": n1.dY},
                use_high_precision=use_high_precision,
            )
        elif abs_y <= abs_x and abs_y <= abs_z:
            ortho = pynvector(
                {"x": -n1.dZ, "y": dtype(0.0), "z": n1.dX},
                use_high_precision=use_high_precision,
            )
        else:
            ortho = pynvector(
                {"x": -n1.dY, "y": n1.dX, "z": dtype(0.0)},
                use_high_precision=use_high_precision,
            )

        ortho = ortho.normalize()
        theta = dtype(np.pi) * dtype(t)
        return (n1 * np.cos(theta) + ortho * np.sin(theta)).normalize()

    # Standard slerp
    theta = np.arccos(dot_product)
    sin_theta = np.sin(theta)

    if abs(sin_theta) <= eps:
        result = n1 * (dtype(1.0) - dtype(t)) + n2 * dtype(t)
        return result.normalize()

    factor1 = np.sin((1 - t) * theta) / sin_theta
    factor2 = np.sin(t * theta) / sin_theta

    return (n1 * factor1 + n2 * factor2).normalize()


def split_polyline_by_length(aFlowline_in, dDistance):

    aPolyline_out = list()
    nPolyline = len(aFlowline_in)
    for i in range(nPolyline):
        pPolyline = aFlowline_in[i]
        pPolyline_out = pPolyline.split_by_length(dDistance)
        aPolyline_out.append(pPolyline_out)

    return aPolyline_out


def split_line_by_length(pLine_in, dLength_in, tolerance=1e-6, use_high_precision=True):
    """
    Split a line into smaller segments with maximum length constraint.

    Uses optimal segmentation with spherical interpolation. Supports high-precision
    mode (float128) to reduce numerical errors in the coordinate conversion chain.

    Args:
        pLine_in: Input line to split
        dLength_in: Maximum length for each segment (meters)
        tolerance: Relative tolerance for length comparison (default: 1e-6)
        use_high_precision: Use float128 for conversions (default: True)

    Returns:
        List of line segments

    Raises:
        ValueError: If length threshold is not positive
        AttributeError: If line doesn't have required attributes

    Note:
        When use_high_precision=True, the conversion chain uses numpy.float128
        to maintain precision during spherical interpolation. This significantly
        reduces degenerate segment warnings but has a small performance cost (~10-30%).

        The adaptive segmentation limits the number of segments based on numerical
        precision to prevent creating segments that cannot be distinguished.
    """
    from mpas_land_mesh.classes.edge import pyedge

    # Input validation
    if dLength_in <= 0:
        raise ValueError("Length threshold must be positive")

    if not hasattr(pLine_in, "dLength"):
        raise AttributeError("Line must have dLength attribute")

    dLength_total = pLine_in.dLength

    # Early return for lines that are already short enough
    if dLength_total <= dLength_in * (1 + tolerance):
        return [pLine_in]

    # Calculate optimal number of segments
    nSegments = int(np.ceil(dLength_total / dLength_in))

    # Adaptive segmentation: limit based on numerical precision
    # Estimate the angular span on the sphere
    dRadius_earth = 6371000.0  # meters (mean Earth radius)
    dAngle_total = dLength_total / dRadius_earth  # approximate angle in radians

    # Minimum distinguishable angle based on precision
    # float128: ~33 decimal digits, can distinguish ~1e-30 radians
    # float64: ~15 decimal digits, can distinguish ~1e-15 radians
    # Use conservative estimates to ensure points remain distinguishable
    dAngle_min = 1e-14 if use_high_precision else 1e-12
    nSegments_max = max(1, int(dAngle_total / dAngle_min))

    if nSegments > nSegments_max:
        logger.info(
            f"Reducing segments from {nSegments} to {nSegments_max} "
            f"due to precision limits (line length={dLength_total:.2f}m, "
            f"angle={np.degrees(dAngle_total):.6e}°)"
        )
        nSegments = nSegments_max

    pPoint_start = pLine_in.pPoint_start
    pPoint_end = pLine_in.pPoint_end

    # Use high precision for n-vector conversion if requested
    n1 = pPoint_start.toNvector(use_high_precision=use_high_precision)
    n2 = pPoint_end.toNvector(use_high_precision=use_high_precision)

    aLine_out = []

    # Build segments from the last accepted endpoint so skipped degenerate
    # interpolation intervals are merged into the next valid segment.
    dThreshold_segment_m = 1e-3  # 1 mm
    pPoint_start_seg = pPoint_start

    for i in range(nSegments):
        t1 = i / nSegments
        t2 = (i + 1) / nSegments

        if i == nSegments - 1:
            pPoint_end_seg = pPoint_end
        else:
            mid2 = slerp(n1, n2, t2)
            pPoint_end_seg = mid2.toLatLon()

        # Check for degenerate segments (nearly same start and end points)
        dSeg_candidate = pPoint_start_seg.calculate_distance(pPoint_end_seg)
        if dSeg_candidate <= dThreshold_segment_m:
            logger.debug(
                f"Skipping degenerate interval {i}/{nSegments}: "
                f"coords=({pPoint_start_seg.dLongitude_degree:.15f}, "
                f"{pPoint_start_seg.dLatitude_degree:.15f}), "
                f"distance={dSeg_candidate:.6e}m, "
                f"t1={t1:.6f}, t2={t2:.6f}, "
                f"precision={'float128' if use_high_precision else 'float64'}"
            )
            print(
                f"Skipping degenerate interval {i}/{nSegments}: "
                f"coords=({pPoint_start_seg.dLongitude_degree:.15f}, "
                f"{pPoint_start_seg.dLatitude_degree:.15f}), "
                f"distance={dSeg_candidate:.6e}m, "
                f"t1={t1:.6f}, t2={t2:.6f}, "
                f"precision={'float128' if use_high_precision else 'float64'}"
            )
            continue

        # Create line segment with error handling
        try:
            pLine = pyedge(pPoint_start_seg, pPoint_end_seg)
            aLine_out.append(pLine)
            pPoint_start_seg = pPoint_end_seg
        except ValueError as e:
            logger.debug(f"Skipping interval {i}/{nSegments}: {str(e)}")
            print(f"Skipping interval {i}/{nSegments}: {str(e)}")
            continue

    return aLine_out


def convert_gcs_coordinates_to_flowline(aCoordinates_in):
    """convert coordinates to flowline, but we cannot setup index and id yet

    Args:
        aCoordinates_in (_type_): _description_

    Returns:
        _type_: _description_
    """
    from mpas_land_mesh.classes.vertex import pyvertex
    from mpas_land_mesh.classes.edge import pyedge
    from mpas_land_mesh.classes.flowline import pyflowline

    npoint = len(aCoordinates_in)
    aVertex = [
        pyvertex({"dLongitude_degree": x, "dLatitude_degree": y})
        for x, y in aCoordinates_in
    ]
    aEdge = [
        pyedge(aVertex[j], aVertex[j + 1])
        for j in range(npoint - 1)
        if aVertex[j] != aVertex[j + 1]
    ]
    if len(aEdge) == 0:
        print("No edge is created")
        return None

    return pyflowline(aEdge)


def convert_pcs_coordinates_to_flowline(aCoordinates_in, pProjection_in):
    from mpas_land_mesh.classes.vertex import pyvertex
    from mpas_land_mesh.classes.edge import pyedge
    from mpas_land_mesh.classes.flowline import pyflowline

    npoint = len(aCoordinates_in)

    aVertex = list()
    for i in range(npoint):
        x = aCoordinates_in[i][0]
        y = aCoordinates_in[i][1]
        dummy = dict()
        dummy["x"] = x
        dummy["y"] = y
        lon, lat = reproject_coordinates(x, y, pProjection_in)
        dummy["lon"] = lon
        dummy["lat"] = lat
        pVertex = pyvertex(dummy)
        aVertex.append(pVertex)

    aEdge = list()
    for j in range(npoint - 1):
        pEdge = pyedge(aVertex[j], aVertex[j + 1])
        aEdge.append(pEdge)

    pFlowline = pyflowline(aEdge)

    return pFlowline
