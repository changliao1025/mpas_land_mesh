"""
Preprocessing utilities for river networks and coastlines.

Note on import structure
------------------------
``coastlines.py`` imports ``pyearth`` and ``pyearthbuffer`` at module level.
Those packages may not be installed in every environment, so the coastline
symbols are imported inside a ``try/except`` block.  If the import fails the
coastline functions are simply not available at the package level; import them
directly when needed::

    from mpas_land_mesh.preprocessing.coastlines import (
        create_land_ocean_mask_from_naturalearth,
        fix_naturalearth_hydrosheds_incompatibility,
        ...
    )
"""

from .river_network import (
    convert_geometry_flowline,
    precompute_flowline_geometries,
    precompute_flowline_geometries_by_segment,
    simplify_hydrorivers_network,
    get_outlet_location,
    tag_river_outlet,
)

__all__ = [
    # river_network
    'convert_geometry_flowline',
    'precompute_flowline_geometries',
    'precompute_flowline_geometries_by_segment',
    'simplify_hydrorivers_network',
    'get_outlet_location',
    'tag_river_outlet',
]

try:
    from .coastline import (
        create_land_ocean_mask,
        create_land_ocean_mask_from_naturalearth,
        geometries_bbox_overlap,
        fix_naturalearth_hydrosheds_incompatibility,
    )
    __all__ += [
        'create_land_ocean_mask',
        'create_land_ocean_mask_from_naturalearth',
        'geometries_bbox_overlap',
        'fix_naturalearth_hydrosheds_incompatibility',
    ]
except ImportError:
    pass
