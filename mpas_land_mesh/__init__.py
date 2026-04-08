"""
Unified Land-River Mesh (ULRM) Workflow Module

A lightweight, standalone module for running land-river mesh workflows
without requiring full dependency installations.

Version: 1.0.0
"""

__version__ = "1.0.0"
__author__ = "ULRM Team"

from . import utilities
from . import preprocessing

__all__ = ['utilities', 'preprocessing']
