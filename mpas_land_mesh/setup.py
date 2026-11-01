"""
Setup script for ULRM Workflow module
"""

from pathlib import Path

from setuptools import find_packages, setup

here = Path(__file__).resolve().parent
readme_path = here / "README.md"
if not readme_path.exists():
    readme_path = here.parent / "README.md"

if readme_path.exists():
    with open(readme_path, encoding="utf-8") as fh:
        long_description = fh.read()
else:
    long_description = "Lightweight module for unified land-river mesh workflows"

setup(
    name="ulrm-workflow",
    version="1.0.0",
    author="ULRM Team",
    description="Lightweight module for unified land-river mesh workflows",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: GIS",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
        "pyyaml>=6.0",
        "pyearth>=0.2.1",
        "pyearthmesh>=0.1.1",
        "pyearthriver>=0.1.0",
        "pyearthbuffer>=0.1.0",
        "geographiclib>=2.0",
    ],
    extras_require={
        "gdal": ["gdal>=3.0.0"],
        "dev": ["pytest>=6.0", "black", "flake8"],
    },
)
