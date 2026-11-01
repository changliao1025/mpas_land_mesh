"""
Preprocessing utilities for river networks and coastlines.
"""

from hexwatershed_utility.preprocess.feature.river import (
    crop_river_network_by_boundary,
    get_outlet_location,
    simplify_river_network,
    tag_river_outlet,
)
from hexwatershed_utility.preprocess.feature.river.simplify_river_network import (
    convert_geometry_flowline,
    precompute_flowline_geometries,
    precompute_flowline_geometries_by_segment,
)

# Alias for backward compatibility
simplify_hydrorivers_network = simplify_river_network

__all__ = [
    # river network (delegated to hexwatershed_utility)
    "convert_geometry_flowline",
    "crop_river_network_by_boundary",
    "get_outlet_location",
    "precompute_flowline_geometries",
    "precompute_flowline_geometries_by_segment",
    "simplify_hydrorivers_network",
    "simplify_river_network",
    "tag_river_outlet",
]

try:
    from .coastline import (
        create_land_ocean_mask,
        create_land_ocean_mask_from_naturalearth,
        fix_naturalearth_hydrosheds_incompatibility,
        geometries_bbox_overlap,
    )

    __all__ += [
        "create_land_ocean_mask",
        "create_land_ocean_mask_from_naturalearth",
        "geometries_bbox_overlap",
        "fix_naturalearth_hydrosheds_incompatibility",
    ]
except ImportError:
    pass
