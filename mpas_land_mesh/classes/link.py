"""
PyCellLink class for representing connections between mesh cells in flowline networks.

Copied from pyflowline/classes/link.py and adapted to remove
external pyflowline dependencies.
"""

import json
from json import JSONEncoder
from typing import Optional, Any, Tuple
import numpy as np


class LinkClassEncoder(JSONEncoder):
    """
    Custom JSON encoder for pycelllink objects.

    Handles numpy data types, mesh cell objects, and other complex types,
    converting them to native Python types for JSON serialization.
    """

    def default(self, obj):
        if isinstance(obj, (np.integer, np.int32, np.int64)):
            return int(obj)
        if isinstance(obj, (np.float32, np.float64)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, list):
            return obj
        # Handle edge-like objects
        if hasattr(obj, "lEdgeID"):
            return obj.lEdgeID
        # Handle vertex-like objects
        if hasattr(obj, "lVertexID"):
            return obj.lVertexID
        # Handle flowline-like objects
        if hasattr(obj, "lFlowlineID"):
            return obj.lFlowlineID
        # Handle meshcell-like objects
        if hasattr(obj, "lCellID"):
            return obj.lCellID
        if isinstance(obj, pycelllink):
            return obj.lLinkID
        return JSONEncoder.default(self, obj)


class pycelllink:
    """
    Cell link class for representing connections between mesh cells.

    Represents the topological connection between two mesh cells through a
    shared edge. Cell links are fundamental for defining mesh connectivity,
    neighbor relationships, and flow routing in mesh-based models.

    Attributes:
        lLinkIndex (int): Sequential index of the link in a collection
        lLinkID (int): Unique identifier for the link
        pCell_start: Starting/source mesh cell
        pCell_end: Ending/target mesh cell
        pEdge_link: Edge object that connects the two cells
        dLink_length (float): Length of the connecting edge
        dLink_width (float): Width of the connection (if applicable)
        iLink_type (int): Type of link (internal, boundary, etc.)
        bLink_active (bool): Whether the link is active/valid
    """

    def __init__(self, pCell_start_in, pCell_end_in, pEdge_link_in) -> None:
        """
        Initialize a cell link object.

        Args:
            pCell_start_in: The starting/source mesh cell object
            pCell_end_in: The ending/target mesh cell object
            pEdge_link_in: The edge object that connects the two cells
        """
        # Input validation
        if pCell_start_in is None:
            raise ValueError("Starting cell cannot be None")
        if pCell_end_in is None:
            raise ValueError("Ending cell cannot be None")
        if pEdge_link_in is None:
            raise ValueError("Linking edge cannot be None")

        # Check that cells are different (a cell shouldn't link to itself)
        if (
            hasattr(pCell_start_in, "lCellID")
            and hasattr(pCell_end_in, "lCellID")
            and pCell_start_in.lCellID == pCell_end_in.lCellID
            and pCell_start_in.lCellID != -1
        ):
            raise ValueError("Cannot create link: start and end cells are the same")

        # Initialize link attributes
        self.lLinkIndex: int = 0
        self.lLinkID: int = 0

        # Core link components
        self.pCell_start = pCell_start_in
        self.pCell_end = pCell_end_in
        self.pEdge_link = pEdge_link_in

        # Additional link properties
        self.dLink_length: float = -1.0  # Length of the connecting edge
        self.dLink_width: float = -1.0   # Width of the connection
        self.iLink_type: int = 0         # Type of link (0=internal, 1=boundary, etc.)
        self.bLink_active: bool = True   # Whether the link is active

        # Calculate link properties if possible
        self._calculate_link_properties()

    def __repr__(self) -> str:
        start_id = (
            getattr(self.pCell_start, "lCellID", "Unknown")
            if self.pCell_start
            else "None"
        )
        end_id = (
            getattr(self.pCell_end, "lCellID", "Unknown") if self.pCell_end else "None"
        )
        return (
            f"pycelllink(ID={self.lLinkID}, Index={self.lLinkIndex}, "
            f"Start={start_id}, End={end_id}, Length={self.dLink_length:.2f}m, "
            f"Type={self.iLink_type}, Active={self.bLink_active})"
        )

    def __str__(self) -> str:
        start_id = (
            getattr(self.pCell_start, "lCellID", "Unknown")
            if self.pCell_start
            else "None"
        )
        end_id = (
            getattr(self.pCell_end, "lCellID", "Unknown") if self.pCell_end else "None"
        )
        return (
            f"pycelllink(ID={self.lLinkID}, Start={start_id}, End={end_id}, "
            f"Length={self.dLink_length:.2f}m)"
        )

    def __hash__(self) -> int:
        start_hash = hash(id(self.pCell_start)) if self.pCell_start else 0
        end_hash = hash(id(self.pCell_end)) if self.pCell_end else 0
        edge_hash = hash(id(self.pEdge_link)) if self.pEdge_link else 0
        return hash((start_hash, end_hash, edge_hash))

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, pycelllink):
            return NotImplemented
        same_forward = (
            self.pCell_start == other.pCell_start
            and self.pCell_end == other.pCell_end
            and self.pEdge_link == other.pEdge_link
        )
        same_reverse = (
            self.pCell_start == other.pCell_end
            and self.pCell_end == other.pCell_start
            and self.pEdge_link == other.pEdge_link
        )
        return same_forward or same_reverse

    def _calculate_link_properties(self) -> None:
        """Calculate link properties such as length and width."""
        try:
            if self.pEdge_link and hasattr(self.pEdge_link, "calculate_length"):
                self.dLink_length = self.pEdge_link.calculate_length()
            elif self.pEdge_link and hasattr(self.pEdge_link, "dLength"):
                self.dLink_length = self.pEdge_link.dLength

            if (
                self.pCell_start
                and self.pCell_end
                and hasattr(self.pCell_start, "iFlag_coast")
                and hasattr(self.pCell_end, "iFlag_coast")
            ):
                if self.pCell_start.iFlag_coast == 1 or self.pCell_end.iFlag_coast == 1:
                    self.iLink_type = 1
                else:
                    self.iLink_type = 0
        except Exception:
            pass

    def set_link_id(self, lLinkID: int) -> None:
        if not isinstance(lLinkID, (int, np.integer)):
            raise TypeError(f"Link ID must be an integer, got {type(lLinkID)}")
        self.lLinkID = int(lLinkID)

    def set_link_index(self, lLinkIndex: int) -> None:
        if not isinstance(lLinkIndex, (int, np.integer)):
            raise TypeError(f"Link index must be an integer, got {type(lLinkIndex)}")
        self.lLinkIndex = int(lLinkIndex)

    def set_link_type(self, iLink_type: int) -> None:
        if not isinstance(iLink_type, (int, np.integer)):
            raise TypeError(f"Link type must be an integer, got {type(iLink_type)}")
        if iLink_type < 0:
            raise ValueError(f"Link type must be non-negative, got {iLink_type}")
        self.iLink_type = int(iLink_type)

    def set_active(self, bActive: bool) -> None:
        if not isinstance(bActive, bool):
            raise TypeError(f"Active status must be a boolean, got {type(bActive)}")
        self.bLink_active = bActive

    def reverse_direction(self) -> None:
        """Reverse the direction of the link."""
        temp_cell = self.pCell_start
        self.pCell_start = self.pCell_end
        self.pCell_end = temp_cell

    def get_connected_cells(self) -> Tuple[Any, Any]:
        return (self.pCell_start, self.pCell_end)

    def get_other_cell(self, cell) -> Optional[Any]:
        if cell is None:
            raise ValueError("Reference cell cannot be None")
        if self.pCell_start == cell:
            return self.pCell_end
        elif self.pCell_end == cell:
            return self.pCell_start
        else:
            return None

    def contains_cell(self, cell) -> bool:
        if cell is None:
            raise ValueError("Cell to check cannot be None")
        return self.pCell_start == cell or self.pCell_end == cell

    def get_link_properties(self) -> dict:
        return {
            "link_id": self.lLinkID,
            "link_index": self.lLinkIndex,
            "link_type": self.iLink_type,
            "active": self.bLink_active,
            "length": self.dLink_length,
            "width": self.dLink_width,
            "start_cell_id": (
                getattr(self.pCell_start, "lCellID", None) if self.pCell_start else None
            ),
            "end_cell_id": (
                getattr(self.pCell_end, "lCellID", None) if self.pCell_end else None
            ),
            "edge_id": (
                getattr(self.pEdge_link, "lEdgeID", None) if self.pEdge_link else None
            ),
        }

    def is_boundary_link(self) -> bool:
        return self.iLink_type == 1

    def is_valid(self) -> bool:
        has_valid_cells = self.pCell_start is not None and self.pCell_end is not None
        has_valid_edge = self.pEdge_link is not None
        has_reasonable_length = self.dLink_length >= 0
        return (
            has_valid_cells
            and has_valid_edge
            and has_reasonable_length
            and self.bLink_active
        )

    def copy(self) -> "pycelllink":
        new_link = pycelllink(self.pCell_start, self.pCell_end, self.pEdge_link)
        new_link.lLinkIndex = self.lLinkIndex
        new_link.lLinkID = self.lLinkID
        new_link.dLink_length = self.dLink_length
        new_link.dLink_width = self.dLink_width
        new_link.iLink_type = self.iLink_type
        new_link.bLink_active = self.bLink_active
        return new_link

    def tojson(self) -> str:
        """Convert cell link object to a JSON string."""
        obj = self.__dict__.copy()
        obj_clean = {}
        for key, value in obj.items():
            if key in ["pCell_start", "pCell_end", "pEdge_link"]:
                if hasattr(value, "lCellID"):
                    obj_clean[key + "_id"] = value.lCellID
                elif hasattr(value, "lEdgeID"):
                    obj_clean[key + "_id"] = value.lEdgeID
                else:
                    obj_clean[key + "_id"] = str(value) if value else None
            else:
                obj_clean[key] = value

        sJson = json.dumps(
            obj_clean, sort_keys=True, indent=4, ensure_ascii=True, cls=LinkClassEncoder
        )
        return sJson
