# Installation Guide

## Prerequisites

### System Requirements
- Python 3.8 or higher
- GDAL library (version 3.0 or higher)

### Quick Start

#### Option 1: Using Conda (Recommended)

```bash
# Create a new environment
conda create -n ulrm python=3.9
conda activate ulrm

# Install dependencies
conda install -c conda-forge gdal numpy scipy

# Install ULRM workflow module
cd mpas_land_mesh
pip install -e .
```

#### Option 2: Using pip

```bash
# Install Python dependencies
pip install numpy scipy

# Install GDAL (system-dependent, see below)

# Install ULRM workflow module
cd mpas_land_mesh
pip install -e .
```

### GDAL Installation by Platform

#### Linux (Ubuntu/Debian)
```bash
sudo apt-get update
sudo apt-get install gdal-bin libgdal-dev
pip install gdal==$(gdal-config --version)
```

#### macOS
```bash
brew install gdal
pip install gdal==$(gdal-config --version)
```

#### Windows
```bash
# Use conda (easiest)
conda install -c conda-forge gdal

# Or download OSGeo4W installer
# https://trac.osgeo.org/osgeo4w/
```

## Verification

Test your installation:

```python
from osgeo import gdal, ogr
print(f"GDAL version: {gdal.__version__}")

import mpas_land_mesh
print(f"ULRM Workflow version: {mpas_land_mesh.__version__}")
```

## Running the Example

```bash
cd mpas_land_mesh/examples
python run_minimal_workflow.py
```

## Troubleshooting

### GDAL Import Error
If you get `ImportError: No module named 'osgeo'`:
- Ensure GDAL is installed system-wide
- Try using conda instead of pip
- Check Python version matches GDAL version

### Version Mismatch
If you get version compatibility errors:
- Reinstall with matching versions: `pip install gdal==$(gdal-config --version)`
- Use conda to manage versions automatically

### Missing scipy.ndimage
The raster buffer functions require scipy:
```bash
pip install scipy
```

## Next Steps

See [`USAGE.md`](USAGE.md) for detailed usage instructions.
