"""
Confluence class for representing river confluence points.

Copied from pyflowline/classes/confluence.py and adapted to remove
external pyflowline and pyearth dependencies.
"""

import json
from json import JSONEncoder
import numpy as np

from mpas_land_mesh.classes.vertex import pyvertex
from mpas_land_mesh.utilities.geometry import calculate_angle_between_point


class ConfluenceClassEncoder(JSONEncoder):
    """Confluence Class Encoder

    Args:
        JSONEncoder (_type_): _description_
    """

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.float32):
            return float(obj)
        if isinstance(obj, list):
            pass
        if isinstance(obj, pyvertex):
            return json.loads(obj.tojson())
        # Handle pyflowline objects by their ID attribute
        if hasattr(obj, "lFlowlineID"):
            return obj.lFlowlineID

        return JSONEncoder.default(self, obj)


class pyconfluence:
    """
    The pyconfluence class
    Returns:
        object: A confluence object
    """

    def __init__(self, pVertex_center, aFlowline_upstream_in, aFlowline_downstream_in):
        """
        Initialize a pyconfluence object

        Args:
            pVertex_center (pyvertex): The center vertex
            aFlowline_upstream_in (list [pyflowline]): A list of upstream flowlines
            aFlowline_downstream_in (list [pyflowline]): A list of downstream flowlines
        """
        # Initialize instance attributes with default values
        self.lIndex = -1
        self.lConfluenceID = -1
        self.dAngle_upstream = 0.0
        self.nUpstream = len(aFlowline_upstream_in)
        self.nDownstream = len(aFlowline_downstream_in)

        try:
            self.pVertex_confluence = pVertex_center
            self.aFlowline_upstream = aFlowline_upstream_in
            self.aFlowline_downstream = (
                aFlowline_downstream_in  # to support braided flowlines in the future
            )
        except:
            print("Initialization of confluence failed!")

        return

    def calculate_branching_angle(self):
        """
        Calculate the confluence branching angle
        (https://www.pnas.org/doi/10.1073/pnas.1215218109)

        Returns:
            float: The branching angle in degree
        """
        # normally there are 2 edges meet at confluence
        if len(self.aFlowline_upstream) == 2:
            pFlowline1 = self.aFlowline_upstream[0]
            pFlowline2 = self.aFlowline_upstream[1]
            nedge1 = pFlowline1.nEdge
            nedge2 = pFlowline2.nEdge
            x1 = pFlowline1.aEdge[nedge1 - 1].pVertex_start.dLongitude_degree
            y1 = pFlowline1.aEdge[nedge1 - 1].pVertex_start.dLatitude_degree
            x2 = self.pVertex_confluence.dLongitude_degree
            y2 = self.pVertex_confluence.dLatitude_degree
            x3 = pFlowline2.aEdge[nedge2 - 1].pVertex_start.dLongitude_degree
            y3 = pFlowline2.aEdge[nedge2 - 1].pVertex_start.dLatitude_degree
            self.dAngle_upstream = calculate_angle_between_point(x1, y1, x2, y2, x3, y3)
        else:
            print("multiple upstream")
            print(len(self.aFlowline_upstream))
            self.dAngle_upstream = 0.0

        return self.dAngle_upstream

    def tojson(self):
        """
        Convert a pyconfluence object to json

        Returns:
            json str: A json string
        """
        aSkip = ["aFlowline_upstream"]
        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)
        sJson = json.dumps(
            obj, sort_keys=True, indent=4, ensure_ascii=True, cls=ConfluenceClassEncoder
        )
        return sJson
