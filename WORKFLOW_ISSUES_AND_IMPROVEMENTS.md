# Workflow Issues and Improvements Analysis

## Overview
This document identifies potential issues and improvement opportunities in the minimal workflow (`examples/run_minimal_workflow.py`) and related modules.

**Note:** Completed improvements have been moved to [`COMPLETED_IMPROVEMENTS.md`](COMPLETED_IMPROVEMENTS.md).

---

## Critical Issues

### 1. **Hardcoded User-Specific Paths**
**Location:** [`examples/run_minimal_workflow.py:55-56`](examples/run_minimal_workflow.py:55)

```python
sWorkspace_input = '/qfs/people/liao313/workspace/python/unified_land_river_mesh/data/global/input'
sWorkspace_output = '/compyfs/liao313/04model/pyhexwatershed/global/'
```

**Issue:** Hardcoded absolute paths specific to user `liao313` will fail for other users.

**Impact:** High - Workflow will not run on different systems or for different users.

**Recommendation:**
- Use environment variables or configuration files
- Provide relative paths or make paths configurable via command-line arguments
- Add path validation with helpful error messages

```python
# Suggested improvement:
import os
sWorkspace_input = os.environ.get('MPAS_INPUT_DIR', './data/input')
sWorkspace_output = os.environ.get('MPAS_OUTPUT_DIR', './output')
```

---

### 2. **Hardcoded HPC-Specific Paths in Job Creation**
**Location:** [`mpas_land_mesh/classes/jigsawcase.py:189-193`](mpas_land_mesh/classes/jigsawcase.py:189)

```python
sLine = (
    "#!/qfs/people/liao313/.conda/envs/"
    + sConda_env_name
    + "/bin/"
    + "python3"
    + "\n"
)
```

**Issue:** Hardcoded user path in shebang line makes HPC job scripts non-portable.

**Impact:** High - Job submission will fail for other users.

**Recommendation:**
- Use `#!/usr/bin/env python3` for portability
- Or dynamically detect the Python interpreter path

---

### 3. **Missing Error Handling for File Operations**
**Location:** Multiple locations in [`examples/run_minimal_workflow.py`](examples/run_minimal_workflow.py:99)

```python
sFilename_flowline_hydrosheds_in = '/compyfs/liao313/00raw/hydrology/hydrosheds/hydroriver/...'
# No check if file exists before use
```

**Issue:** No validation that input files exist before processing.

**Impact:** Medium - Cryptic errors when files are missing.

**Recommendation:**
```python
if not os.path.exists(sFilename_flowline_hydrosheds_in):
    raise FileNotFoundError(f"Required input file not found: {sFilename_flowline_hydrosheds_in}")
```

---

### 4. **Inconsistent Flag Logic**
**Location:** [`examples/run_minimal_workflow.py:48-49`](examples/run_minimal_workflow.py:48)

```python
iFlag_simplify_hydrosheds_river_network = 0
iFlag_process_coastline = 0
```

**Issue:** Both flags are set to 0 (disabled), but the workflow continues as if data exists. The `else` blocks (lines 122-124, 149-151) have only `pass` statements with no validation.

**Impact:** Medium - Silent failures if expected output files don't exist.

**Recommendation:**
```python
else:
    # Validate that required files exist
    if not os.path.exists(sFilename_flowline_hydrosheds_out):
        raise FileNotFoundError(f"Simplified river network not found. Set iFlag_simplify_hydrosheds_river_network=1 or provide existing file.")
```

---

### 5. **Potential Variable Reference Before Assignment**
**Location:** [`mpas_land_mesh/preprocessing/coastlines.py:52`](mpas_land_mesh/preprocessing/coastlines.py:52)

```python
if iRaster_buffer_pixel > 0:
    sFilename_tif_wo_island_buffered = ...
    # ... processing ...
else:
    sFilename_tif_wo_island_buffered_fixed = ...
    fix_raster_antimeridian_issue(sFilename_tif_wo_island_buffered, ...)  # ERROR!
```

**Issue:** In the `else` block, `sFilename_tif_wo_island_buffered` is used but was never defined.

**Impact:** High - Runtime error when `iRaster_buffer_pixel <= 0`.

**Recommendation:**
```python
else:
    sFilename_tif_wo_island_buffered_fixed = sFilename_tif_wo_island
    # No buffering needed
```

---

## Medium Priority Issues

### 6. **Unclear Variable Naming Convention**
**Location:** Throughout codebase

**Issue:** Inconsistent naming with Hungarian notation (`s` for string, `i` for int, `d` for double, `a` for array).

**Impact:** Low - Reduces code readability for developers unfamiliar with this convention.

**Recommendation:**
- Document the naming convention in a style guide
- Or migrate to more Pythonic naming (e.g., `filename_flowline` instead of `sFilename_flowline`)

---

### 7. **Redundant Directory Creation**
**Location:** [`examples/run_minimal_workflow.py:60-66`](examples/run_minimal_workflow.py:60)

```python
if os.path.exists(sWorkspace_river_network_output) is False:
    os.makedirs(sWorkspace_river_network_output)
```

**Issue:** Using `is False` is not Pythonic; `os.makedirs` already has `exist_ok` parameter.

**Impact:** Low - Code verbosity.

**Recommendation:**
```python
os.makedirs(sWorkspace_river_network_output, exist_ok=True)
```

---

### 8. **Incomplete Return Value Handling**
**Location:** [`mpas_land_mesh/preprocessing/coastlines.py:54`](mpas_land_mesh/preprocessing/coastlines.py:54)

```python
return sFilename_tif_wo_island_buffered_fixed, sFilename_wo_island
```

**Issue:** Function returns two values, but in the workflow (line 131), only one is captured:

```python
sFilename_tif_wo_island, sFilename_vector_coastline = create_land_ocean_mask_from_naturalearth(...)
```

**Impact:** Low - Works correctly, but variable names are misleading.

**Recommendation:**
- Rename variables to match what's actually returned
- Or update function to return values in expected order

---

### 9. **Debug Code Left in Production**
**Location:** [`mpas_land_mesh/preprocessing/river_networks.py:615-617`](mpas_land_mesh/preprocessing/river_networks.py:615)

```python
if iSegment_upstream == 624:
    print('debug')
```

**Issue:** Debug statements with hardcoded segment IDs left in code (appears multiple times).

**Impact:** Low - Clutters output, indicates incomplete cleanup.

**Recommendation:**
- Remove debug statements or convert to proper logging
- Use logging module with debug level

---

### 10. **Recursion Limit Modification**
**Location:** [`mpas_land_mesh/preprocessing/river_networks.py:7-8`](mpas_land_mesh/preprocessing/river_networks.py:7)

```python
import sys
sys.setrecursionlimit(100000)
```

**Issue:** Globally modifying recursion limit can mask infinite recursion bugs and cause stack overflow.

**Impact:** Medium - Potential for hard-to-debug crashes.

**Recommendation:**
- Refactor recursive algorithms to use iteration where possible
- Document why high recursion limit is needed
- Add safeguards in recursive functions

---

### 11. **Missing Configuration Validation**
**Location:** [`examples/run_minimal_workflow.py:163-207`](examples/run_minimal_workflow.py:163)

**Issue:** Many configuration values are set via `change_json_key_value` but no validation that required keys exist or values are valid.

**Impact:** Medium - Silent failures or cryptic errors during mesh generation.

**Recommendation:**
```python
# Validate configuration after all changes
required_keys = ['dResolution_ocean', 'dResolution_land', 'sFilename_river_network_vector']
for key in required_keys:
    if key not in config:
        raise ValueError(f"Required configuration key missing: {key}")
```

---

### 12. **Inconsistent String Formatting**
**Location:** Throughout codebase

**Issue:** Mix of `.format()`, f-strings, and `+` concatenation.

```python
sBasin_id = '{:04d}'.format(i)  # Line 140
sLine = "#!/bin/bash\n"  # Line 266
sLine = "#SBATCH --job-name=" + self.sCase + "\n"  # Line 270
```

**Impact:** Low - Inconsistent style.

**Recommendation:**
- Standardize on f-strings (Python 3.6+) for consistency and readability

---

## Low Priority / Enhancement Suggestions

### 13. **No Progress Indicators for Long Operations**
**Location:** [`mpas_land_mesh/preprocessing/river_networks.py:735`](mpas_land_mesh/preprocessing/river_networks.py:735)

**Issue:** Long-running loops have minimal progress feedback.

**Impact:** Low - User doesn't know if process is hung or progressing.

**Recommendation:**
```python
from tqdm import tqdm
for i in tqdm(range(nFlowline_outlet), desc="Processing basins"):
    # ... processing ...
```

---

### 14. **Commented-Out Code**
**Location:** [`mpas_land_mesh/preprocessing/river_networks.py:186`](mpas_land_mesh/preprocessing/river_networks.py:186)

```python
#change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_RRS18to6_ocean", "true")
```

**Issue:** Commented-out code should be removed or explained.

**Impact:** Low - Code clutter.

**Recommendation:**
- Remove if not needed
- Or add comment explaining why it's disabled

---

### 15. **Missing Docstrings for Some Functions**
**Location:** Various helper functions

**Issue:** Some functions lack docstrings (e.g., `convert_geometry_flowline`, `tag_upstream`).

**Impact:** Low - Reduces code maintainability.

**Recommendation:**
- Add comprehensive docstrings following NumPy or Google style

---

### 16. **No Unit Tests**
**Location:** Project structure

**Issue:** No visible test directory or test files.

**Impact:** Medium - Difficult to verify correctness and prevent regressions.

**Recommendation:**
- Create `tests/` directory
- Add unit tests for core functions
- Add integration tests for workflow

---

## Configuration Management Issues

### 17. **Configuration Scattered Across Multiple Files**
**Location:** [`examples/run_minimal_workflow.py`](examples/run_minimal_workflow.py:27)

**Issue:** Configuration parameters are hardcoded in the workflow script rather than in a separate config file.

**Impact:** Medium - Hard to manage different scenarios, requires code changes for different runs.

**Recommendation:**
```python
# Create config.yaml
resolution:
  ocean: 100
  land: 100
  river_network: 10
  coastline: 10

paths:
  input: ./data/input
  output: ./output

# Load in script
import yaml
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

---

### 18. **Boolean Flags as Integers**
**Location:** Throughout codebase

```python
iFlag_simplify_hydrosheds_river_network = 0  # Should be False
iFlag_process_coastline = 0  # Should be False
```

**Issue:** Using 0/1 instead of True/False for boolean flags.

**Impact:** Low - Less Pythonic, slightly less readable.

**Recommendation:**
```python
flag_simplify_hydrosheds_river_network = False
flag_process_coastline = False
```

---

## Performance Issues

### 19. **Inefficient Dictionary Lookups**
**Location:** [`mpas_land_mesh/preprocessing/river_networks.py:762-764`](mpas_land_mesh/preprocessing/river_networks.py:762)

**Issue:** Creating dictionaries for O(1) lookup is good, but could be done earlier to avoid repeated linear searches.

**Impact:** Low - Already optimized in most places.

**Recommendation:**
- Ensure all lookups use dictionaries consistently
- Profile code to identify remaining bottlenecks

---

### 20. **Repeated File I/O in Loop**
**Location:** [`mpas_land_mesh/utilities/change_json_key_value.py:29-50`](mpas_land_mesh/utilities/change_json_key_value.py:29)

**Issue:** Each call to `change_json_key_value` reads and writes the entire JSON file.

**Impact:** Medium - Inefficient when making many changes (lines 165-200 in workflow).

**Recommendation:**
```python
# Use change_json_keys_values (plural) for batch updates
key_value_pairs = {
    "iFlag_geom": "true",
    "iFlag_geom_river_network": "true",
    "dResolution_ocean": dResolution_ocean,
    # ... all changes at once
}
change_json_keys_values(sFilename_jigsaw_configuration_copy, key_value_pairs)
```

---

## Documentation Issues

### 21. **Incomplete README/Documentation**
**Location:** Project root

**Issue:** No clear getting-started guide or example of how to run the workflow with custom data.

**Impact:** Medium - High barrier to entry for new users.

**Recommendation:**
- Add comprehensive README with:
  - Installation instructions
  - Quick start guide
  - Example with sample data
  - Configuration options
  - Troubleshooting section

---

### 22. **Missing Dependency Documentation**
**Location:** [`mpas_land_mesh/requirements.txt`](mpas_land_mesh/requirements.txt:1)

**Issue:** External dependencies like `pyearth`, `pyflowline`, `pyearthbuffer` are used but may not be clearly documented.

**Impact:** Medium - Users may not know what to install.

**Recommendation:**
- Ensure all dependencies are in requirements.txt
- Add version constraints
- Document any system-level dependencies (GDAL, etc.)

---

## Summary

This document tracks **22 pending issues** for improvement. For details on the **4 completed improvements**, see [`COMPLETED_IMPROVEMENTS.md`](COMPLETED_IMPROVEMENTS.md).

### Priority Breakdown

**Critical Issues (5):**
- Issues #1-5: Path portability, error handling, flag logic, variable references

**Medium Priority (7):**
- Issues #6-12: Code quality, naming conventions, configuration validation

**Low Priority / Enhancements (5):**
- Issues #13-16, plus configuration and performance issues

**Configuration Management (2):**
- Issues #17-18: Config file support, boolean flag conventions

**Performance (2):**
- Issues #19-20: Dictionary lookups, file I/O optimization

**Documentation (2):**
- Issues #21-22: README improvements, dependency documentation

---

## Conclusion

The workflow is functional but has several portability, maintainability, and robustness issues. The most critical issues involve hardcoded paths and missing error handling. Addressing the critical and medium priority issues would significantly improve the code quality and user experience.
