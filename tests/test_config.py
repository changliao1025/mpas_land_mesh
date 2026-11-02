import json
import os
import platform
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock

# Gracefully mock optional native C/Fortran libraries if running in minimal Python environment
for mod in ["numpy", "numpy.typing", "osgeo", "osgeo.gdal", "osgeo.ogr", "osgeo.osr"]:
    if mod not in sys.modules:
        try:
            __import__(mod)
        except ImportError:
            m = MagicMock()
            m.__path__ = []
            sys.modules[mod] = m

from pyearth.gis.idl_handler import IdlHandler, IdlStrategy
from pyearthmesh.config.config_manager import (
    DEFAULT_CONFIG,
    compile_workflow,
    create_workflow,
    dump_workflow_config,
    load_workflow_config,
)
from pyearthmesh.config.workflows.mpas_workflow import MPASWorkflow

from mpas_land_mesh.config import load_ulrm_config
from mpas_land_mesh.utilities import DEFAULT_CONFIG as UTILS_DEFAULT_CONFIG
from mpas_land_mesh.utilities import compile_workflow as utils_compile_workflow
from mpas_land_mesh.utilities import create_workflow as utils_create_workflow
from mpas_land_mesh.utilities import dump_workflow_config as utils_dump_workflow_config
from mpas_land_mesh.utilities import load_workflow_config as utils_load_workflow_config
from mpas_land_mesh.workflow import ULRMWorkflow


class TestWorkflowConfig(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory()
        self.base_dir = Path(self.test_dir.name)

    def tearDown(self):
        self.test_dir.cleanup()

    def test_reexported_utilities(self):
        """Verify mpas_land_mesh.utilities re-exports pyearthmesh config functions."""
        self.assertIs(utils_load_workflow_config, load_workflow_config)
        self.assertIs(utils_dump_workflow_config, dump_workflow_config)
        self.assertIs(utils_create_workflow, create_workflow)
        self.assertIs(utils_compile_workflow, compile_workflow)
        self.assertEqual(UTILS_DEFAULT_CONFIG, DEFAULT_CONFIG)

    def test_idl_handler_and_functions(self):
        """Verify IDL crossing detection and splitting using pyearth IdlHandler."""
        import numpy as np

        handler = IdlHandler(strategy=IdlStrategy.CHECK_ONLY)
        self.assertIsNotNone(handler)

        # Polygon crossing the IDL
        crossing_coords = np.array(
            [
                [170.0, -10.0],
                [170.0, 10.0],
                [-170.0, 10.0],
                [-170.0, -10.0],
                [170.0, -10.0],
            ]
        )
        result = handler.detect(crossing_coords)
        self.assertTrue(result.crosses_idl)

        split_handler = IdlHandler(strategy=IdlStrategy.SPLIT)
        split_result = split_handler.detect(crossing_coords)
        self.assertTrue(split_result.crosses_idl)
        self.assertIsNotNone(split_result.sub_polygons)
        sub_polys = [sp for sp in split_result.sub_polygons if len(sp) >= 3]
        self.assertEqual(len(sub_polys), 2)
        for poly in sub_polys:
            self.assertGreaterEqual(len(poly), 3)

        # Normal non-crossing polygon
        normal_coords = np.array(
            [
                [10.0, 10.0],
                [20.0, 10.0],
                [20.0, 20.0],
                [10.0, 20.0],
                [10.0, 10.0],
            ]
        )
        normal_result = handler.detect(normal_coords)
        self.assertFalse(normal_result.crosses_idl)

    def test_workflow_config_defaults_json(self):
        """Verify load_workflow_config loads JSON, merges defaults, and normalizes dates."""
        cfg_path = self.base_dir / "workflow.json"
        raw_data = {
            "workflow": {
                "case_index": 5,
                "mesh_type": "mpas",
                "process_coastline": True,
            },
            "resolutions": {
                "ocean_km": 60,
                "land_km": 30,
            },
            "paths": {
                "output_workspace": "./my_output",
            },
        }
        with open(cfg_path, "w", encoding="utf-8") as f:
            json.dump(raw_data, f)

        config = load_workflow_config(cfg_path)

        # Check dictionary structure
        self.assertIsInstance(config, dict)
        self.assertEqual(config["workflow"]["case_index"], 5)
        self.assertTrue(config["workflow"]["process_coastline"])
        self.assertFalse(config["workflow"]["simplify_hydrosheds_river_network"])
        self.assertIsNotNone(config["workflow"]["date"])

        # Check resolutions with merged defaults
        self.assertEqual(config["resolutions"]["ocean_km"], 60)
        self.assertEqual(config["resolutions"]["land_km"], 30)
        self.assertEqual(config["resolutions"]["coastline_km"], 10)
        self.assertEqual(config["resolutions"]["river_network_km"], 10)

        # Check path resolution relative to config file directory
        expected_output = str((self.base_dir / "my_output").resolve())
        self.assertEqual(config["paths"]["output_workspace"], expected_output)

    def test_workflow_config_loading_yaml_with_comments(self):
        """Verify load_workflow_config parses YAML with inline and block comments."""
        yaml_path = self.base_dir / "workflow_config.yaml"
        yaml_content = """# ==============================================================================
# Tier 1 Specification: Commented MPAS Workflow
# ==============================================================================
workflow:
  date: "20261031"                                # Run date
  mesh_type: mpas                                 # Target MPAS mesh family
  case_index: 7                                   # Case index
  model: jigsaw                                   # Backend engine
  process_coastline: true                         # Enable coastline buffer
  run_on_hpc: false                               # Local execution

resolutions:
  ocean_km: 75.0                                  # Open ocean resolution
  land_km: 25.0                                   # Land resolution
  river_network_km: 5.0                           # River corridor refinement

paths:
  output_workspace: ./annotated_output            # Resolved relative to this file

jigsaw:
  iFlag_geom_coastline: true                      # Enforce coastline geometry
  iFlag_spac_coastline: true                      # Enforce coastline spacing
"""
        yaml_path.write_text(yaml_content, encoding="utf-8")

        config = load_workflow_config(yaml_path)
        self.assertIsInstance(config, dict)
        self.assertEqual(config["workflow"]["mesh_type"], "mpas")
        self.assertEqual(config["workflow"]["case_index"], 7)
        self.assertEqual(config["workflow"]["date"], "20261031")
        self.assertTrue(config["workflow"]["process_coastline"])
        self.assertEqual(config["resolutions"]["ocean_km"], 75.0)
        self.assertEqual(config["resolutions"]["land_km"], 25.0)
        self.assertEqual(config["resolutions"]["river_network_km"], 5.0)

        expected_output = str((self.base_dir / "annotated_output").resolve())
        self.assertEqual(config["paths"]["output_workspace"], expected_output)

        # Verify compilation from YAML
        compiled = compile_workflow(yaml_path)
        self.assertEqual(compiled["iCase_index"], 7)
        self.assertEqual(compiled["dResolution_ocean"], 75.0)
        self.assertEqual(compiled["dResolution_land"], 25.0)
        self.assertTrue(compiled["iFlag_geom_coastline"])
        self.assertTrue(compiled["iFlag_spac_coastline"])

    def test_platform_specific_paths(self):
        """Verify platform-specific paths are resolved according to current OS."""
        current_os = platform.system()
        other_os = "Windows" if current_os != "Windows" else "Linux"

        cfg_path = self.base_dir / "paths.json"
        raw_data = {
            "paths": {
                "river_network_vector": {
                    current_os: "./rivers_current.shp",
                    other_os: "/other/path/rivers.shp",
                    "default": "./default_rivers.shp",
                },
                "region_geometry": {
                    "default": "./default_region.geojson",
                },
            }
        }
        with open(cfg_path, "w", encoding="utf-8") as f:
            json.dump(raw_data, f)

        config = load_workflow_config(cfg_path)
        expected_river = str((self.base_dir / "rivers_current.shp").resolve())
        expected_region = str((self.base_dir / "default_region.geojson").resolve())

        self.assertEqual(config["paths"]["river_network_vector"], expected_river)
        self.assertEqual(config["paths"]["region_geometry"], expected_region)

    def test_unknown_setting_raises_error(self):
        """Verify unknown keys in configuration sections raise ValueError."""
        cfg_path = self.base_dir / "invalid.yaml"
        cfg_path.write_text("workflow:\n  unknown_key: 123\n", encoding="utf-8")

        with self.assertRaises(ValueError):
            load_workflow_config(cfg_path)

    def test_compile_workflow(self):
        """Verify compile_workflow translates frontend spec into native MPAS parameters."""
        cfg_path = self.base_dir / "workflow.yaml"
        yaml_content = """workflow:
  date: "20260901"
  case_index: 2
  mesh_type: mpas

resolutions:
  ocean_km: 30
  land_km: 15
  river_network_km: 5
  coastline_km: 10

paths:
  output_workspace: ./output

jigsaw:
  iFlag_geom_coastline: true
  iFlag_spac_coastline: true
"""
        cfg_path.write_text(yaml_content, encoding="utf-8")

        compiled = compile_workflow(cfg_path)

        self.assertIsInstance(compiled, dict)
        self.assertEqual(compiled["iCase_index"], 2)
        self.assertEqual(compiled["sDate"], "20260901")
        self.assertEqual(compiled["dResolution_ocean"], 30.0)
        self.assertEqual(compiled["dResolution_land"], 15.0)
        self.assertEqual(compiled["dResolution_river_network"], 5.0)
        self.assertEqual(compiled["dResolution_coastline"], 10.0)
        self.assertTrue(compiled["iFlag_geom_coastline"])
        self.assertTrue(compiled["iFlag_spac_coastline"])

    def test_create_workflow_and_hpc_job(self):
        """Verify create_workflow instantiates MPASWorkflow and delegates HPC script generation."""
        cfg_path = self.base_dir / "workflow.yaml"
        yaml_content = """workflow:
  date: "20260901"
  case_index: 1
  mesh_type: mpas

resolutions:
  coastline_km: 10

paths:
  output_workspace: ./output

jigsaw:
  standalone: true
  hours: 3
  slurm: slurm
"""
        cfg_path.write_text(yaml_content, encoding="utf-8")

        workflow = create_workflow(cfg_path)

        self.assertIsInstance(workflow, MPASWorkflow)
        self.assertEqual(workflow.case.iCase_index, 1)
        self.assertEqual(workflow.case.sDate, "20260901")

        # Test HPC job generation
        workflow.create_hpc_job()
        case_dir = Path(workflow.case.sWorkspace_output)
        self.assertTrue(case_dir.is_dir())
        self.assertTrue((case_dir / "run_mpas.py").is_file())
        self.assertTrue((case_dir / "submit.job").is_file())

    def test_directory_and_stem_auto_discovery(self):
        """Verify load_workflow_config auto-discovers workflow_config.yaml in a directory or stem."""
        case_dir = self.base_dir / "test_case"
        case_dir.mkdir()
        yaml_file = case_dir / "workflow_config.yaml"
        yaml_file.write_text(
            "workflow:\n  mesh_type: mpas\n  case_index: 9\n", encoding="utf-8"
        )

        # Pass directory directly
        cfg_from_dir = load_workflow_config(case_dir)
        self.assertEqual(cfg_from_dir["workflow"]["case_index"], 9)
        self.assertEqual(cfg_from_dir["config_file"], yaml_file.resolve())

        # Pass stem without suffix
        stem = case_dir / "workflow_config"
        cfg_from_stem = load_workflow_config(stem)
        self.assertEqual(cfg_from_stem["workflow"]["case_index"], 9)

    def test_dump_workflow_config_roundtrip(self):
        """Verify dump_workflow_config exports valid YAML and JSON files that round-trip cleanly."""
        spec_data = {
            "workflow": {
                "date": "20261031",
                "mesh_type": "mpas",
                "case_index": 4,
                "model": "jigsaw",
                "debug": True,
                "run_on_hpc": False,
            },
            "resolutions": {"ocean_km": 100, "land_km": 20},
            "paths": {"output_workspace": str(self.base_dir / "dumped_out")},
        }

        # Test YAML dump & reload
        yaml_file = self.base_dir / "dumped.yaml"
        dump_workflow_config(spec_data, yaml_file, fmt="yaml")
        self.assertTrue(yaml_file.is_file())
        loaded_yaml = load_workflow_config(yaml_file)
        self.assertEqual(loaded_yaml["workflow"]["case_index"], 4)
        self.assertEqual(loaded_yaml["resolutions"]["land_km"], 20)

        # Test JSON dump & reload
        json_file = self.base_dir / "dumped.json"
        dump_workflow_config(spec_data, json_file, fmt="json")
        self.assertTrue(json_file.is_file())
        loaded_json = load_workflow_config(json_file)
        self.assertEqual(loaded_json["workflow"]["case_index"], 4)
        self.assertEqual(loaded_json["resolutions"]["land_km"], 20)

    def test_yaml_and_json_parity_for_all_examples(self):
        """Verify that every example's workflow_config.yaml and workflow_config.json compile identically."""
        workspace_root = Path(__file__).resolve().parent.parent

        example_dirs = [
            workspace_root / "examples/minimal_workflow",
            workspace_root / "examples/ocn30coast10lnd10",
            workspace_root / "examples/qinghaihu",
            workspace_root / "examples/huanghe",
        ]

        for example_dir in example_dirs:
            yaml_path = example_dir / "workflow_config.yaml"
            json_path = example_dir / "workflow_config.json"

            self.assertTrue(yaml_path.exists(), f"Missing {yaml_path}")
            self.assertTrue(json_path.exists(), f"Missing {json_path}")

            cfg_yaml = load_ulrm_config(yaml_path)
            cfg_json = load_ulrm_config(json_path)

            wf_yaml = ULRMWorkflow(cfg_yaml)
            wf_json = ULRMWorkflow(cfg_json)

            mesh_ready_yaml = wf_yaml.compile_mesh_ready_config(save_yaml=False)
            mesh_ready_json = wf_json.compile_mesh_ready_config(save_yaml=False)

            compiled_yaml = compile_workflow(mesh_ready_yaml)
            compiled_json = compile_workflow(mesh_ready_json)

            # Both YAML and JSON should produce identical compiled keys and values
            # (ignoring sFilename_model_configuration which holds the specific file path)
            for k in compiled_yaml:
                if k == "sFilename_model_configuration":
                    continue
                self.assertIn(
                    k,
                    compiled_json,
                    f"Key {k} missing in JSON compiled output for {example_dir.name}",
                )
                self.assertEqual(
                    compiled_yaml[k],
                    compiled_json[k],
                    f"Value mismatch for key '{k}' in {example_dir.name}: "
                    f"{compiled_yaml[k]} (YAML) != {compiled_json[k]} (JSON)",
                )


if __name__ == "__main__":
    unittest.main()
