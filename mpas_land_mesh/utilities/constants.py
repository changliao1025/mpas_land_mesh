earth_radius = 6371000  # in meters

IDL_threshold = 1e-3  # Degrees from antimeridian

# Constants with explanations
KM2_TO_M2 = 1.0E6
ISLAND_AREA_MULTIPLIER = 10  # At least 10 grid cells
DRAINAGE_AREA_MULTIPLIER = 100  # At least 100 grid cells
# Module-level constants for IDL detection and adjustment
IDL_TOLERANCE = 1e-6  # Tolerance for detecting points on the IDL (±180°)
IDL_OFFSET = 1e-7     # Offset to move points off IDL (< IDL_TOLERANCE so they remain detectable)