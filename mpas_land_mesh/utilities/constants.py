"""Physical and geometric constants."""

from pyearth.gis.idl_handler import IDL_OFFSET, IDL_TOLERANCE

earth_radius = 6371000  # in meters
IDL_threshold = 1e-3  # Degrees from antimeridian

# Constants with explanations
KM2_TO_M2 = 1.0e6
ISLAND_AREA_MULTIPLIER = 10  # At least 10 grid cells
DRAINAGE_AREA_MULTIPLIER = 100  # At least 100 grid cells

__all__ = [
    "DRAINAGE_AREA_MULTIPLIER",
    "IDL_OFFSET",
    "IDL_TOLERANCE",
    "IDL_threshold",
    "ISLAND_AREA_MULTIPLIER",
    "KM2_TO_M2",
    "earth_radius",
]
