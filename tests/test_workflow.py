import json
import os
import platform
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

# Gracefully mock optional native C/Fortran libraries if running in minimal Python environment
for mod in ["numpy", "numpy.typing", "osgeo", "osgeo.gdal", "osgeo.ogr", "osgeo.osr"]:
    if mod not in sys.modules:
        try:
            __import__(mod)
        except ImportError:
            m = MagicMock()
            m.__path__ = []
            sys.modules[mod] = m

import yaml

from mpas_land_mesh.config import (
    DEFAULT_ULRM_CONFIG,
    _deep_merge,
    _resolve_paths,
    _select_platform_path,
    dump_ulrm_config,
    load_ulrm_config,
)
from mpas_land_mesh.workflow import MPASLandMeshWorkflow, ULRMWorkflow, run_workflow


class TestULRMConfig(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory()
        self.base_dir = Path(self.test_dir.name).resolve()

    def tearDown(self):
        self.test_dir.cleanup()

    def test_default_config_structure(self):
        """Verify default configuration has expected sections."""
        self.assertIn("workflow", DEFAULT_ULRM_CONFIG)
        self.assertIn("preprocessing", DEFAULT_ULRM_CONFIG)
        self.assertIn("resolutions", DEFAULT_ULRM_CONFIG)
        self.assertIn("paths", DEFAULT_ULRM_CONFIG)
        self.assertIn("jigsaw", DEFAULT_ULRM_CONFIG)

    def test_deep_merge(self):
        """Verify recursive dictionary merging."""
        defaults = {"a": 1, "b": {"x": 10, "y": 20}}
        overrides = {"b": {"y": 99, "z": 30}, "c": 3}
        merged = _deep_merge(defaults, overrides)
        self.assertEqual(merged["a"], 1)
        self.assertEqual(merged["b"]["x"], 10)
        self.assertEqual(merged["b"]["y"], 99)
        self.assertEqual(merged["b"]["z"], 30)
        self.assertEqual(merged["c"], 3)

    def test_select_platform_path(self):
        """Verify platform path selection."""
        self.assertEqual(_select_platform_path(""), "")
        self.assertEqual(_select_platform_path(None), "")
        self.assertEqual(_select_platform_path("/path/to/file"), "/path/to/file")

        current_os = platform.system()
        os_dict = {
            current_os: f"/custom/{current_os}/path",
            "default": "/default/path",
        }
        self.assertEqual(_select_platform_path(os_dict), f"/custom/{current_os}/path")

        fallback_dict = {"OtherOS": "/other", "default": "/default/path"}
        if current_os != "OtherOS":
            self.assertEqual(_select_platform_path(fallback_dict), "/default/path")

    def test_load_ulrm_config_dict(self):
        """Verify loading and resolving config from dictionary."""
        raw = {
            "workflow": {"case_index": 42},
            "preprocessing": {"simplify_rivers": True},
            "resolutions": {"ocean_km": 60.0, "land_km": 15.0},
        }
        cfg = load_ulrm_config(raw)
        self.assertEqual(cfg["workflow"]["case_index"], 42)
        self.assertTrue(cfg["preprocessing"]["simplify_rivers"])
        self.assertEqual(cfg["resolutions"]["ocean_km"], 60.0)
        # Default fields should be preserved
        self.assertIsNotNone(cfg["workflow"]["date"])
        self.assertIn("paths", cfg)

    def test_load_ulrm_config_yaml_file(self):
        """Verify loading from YAML file."""
        cfg_file = self.base_dir / "test_config.yaml"
        sample_data = {
            "workflow": {"case_index": 5},
            "preprocessing": {"process_coastline": True},
            "paths": {
                "output_workspace": "./out",
            },
        }
        with open(cfg_file, "w", encoding="utf-8") as f:
            yaml.safe_dump(sample_data, f)

        loaded = load_ulrm_config(cfg_file)
        self.assertEqual(loaded["workflow"]["case_index"], 5)
        self.assertTrue(loaded["preprocessing"]["process_coastline"])
        expected_out = str((self.base_dir / "out").resolve())
        self.assertEqual(loaded["paths"]["output_workspace"], expected_out)

    def test_load_ulrm_config_legacy_workflow_flags(self):
        """Verify legacy workflow flags map into preprocessing section."""
        legacy_data = {
            "workflow": {
                "simplify_hydrosheds_river_network": True,
                "process_coastline": True,
                "largest_outlets": 25,
            }
        }
        cfg = load_ulrm_config(legacy_data)
        self.assertTrue(cfg["preprocessing"]["simplify_rivers"])
        self.assertTrue(cfg["preprocessing"]["process_coastline"])
        self.assertEqual(cfg["preprocessing"]["largest_outlets"], 25)

    def test_dump_ulrm_config(self):
        """Verify serializing config to YAML and JSON."""
        data = {"workflow": {"case_index": 7}, "resolutions": {"land_km": 20.0}}
        yaml_path = self.base_dir / "saved.yaml"
        json_path = self.base_dir / "saved.json"

        dump_ulrm_config(data, yaml_path)
        self.assertTrue(yaml_path.is_file())
        with open(yaml_path, encoding="utf-8") as f:
            loaded_yaml = yaml.safe_load(f)
        self.assertEqual(loaded_yaml["workflow"]["case_index"], 7)

        dump_ulrm_config(data, json_path)
        self.assertTrue(json_path.is_file())
        with open(json_path, encoding="utf-8") as f:
            loaded_json = json.load(f)
        self.assertEqual(loaded_json["workflow"]["case_index"], 7)


class TestULRMWorkflow(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory()
        self.base_dir = Path(self.test_dir.name).resolve()
        self.out_dir = self.base_dir / "output"
        self.out_dir.mkdir(parents=True, exist_ok=True)

        self.sample_config = {
            "workflow": {
                "case_index": 1,
                "date": "20261102",
                "mesh_type": "mpas",
                "model": "jigsaw",
                "debug": True,
                "run_on_hpc": False,
            },
            "preprocessing": {
                "simplify_rivers": False,
                "process_coastline": False,
                "largest_outlets": 5,
            },
            "resolutions": {
                "ocean_km": 100.0,
                "land_km": 50.0,
                "river_network_km": 10.0,
                "coastline_km": 10.0,
            },
            "paths": {
                "input_workspace": str(self.base_dir / "input"),
                "output_workspace": str(self.out_dir),
                "river_network_vector": str(self.base_dir / "river.geojson"),
                "river_network_raster": str(self.base_dir / "river.tif"),
                "coastline_vector": str(self.base_dir / "coast.geojson"),
                "coastline_raster": str(self.base_dir / "coast.tif"),
            },
            "jigsaw": {
                "standalone": True,
                "create_directory": True,
                "hours": 5,
                "slurm": "slurm",
            },
        }

    def tearDown(self):
        self.test_dir.cleanup()

    def test_workflow_init(self):
        """Verify workflow initialization and alias."""
        wf = ULRMWorkflow(self.sample_config)
        self.assertEqual(Path(wf.output_workspace).resolve(), self.out_dir)
        self.assertIs(MPASLandMeshWorkflow, ULRMWorkflow)

    def test_preprocess_skip(self):
        """Verify preprocess passes through existing direct paths when flags are False."""
        wf = ULRMWorkflow(self.sample_config)
        paths = wf.preprocess()

        self.assertEqual(
            Path(paths["river_network_vector"]).resolve(),
            (self.base_dir / "river.geojson").resolve(),
        )
        self.assertEqual(
            Path(paths["coastline_raster"]).resolve(),
            (self.base_dir / "coast.tif").resolve(),
        )
        self.assertEqual(Path(paths["output_workspace"]).resolve(), self.out_dir)
        self.assertIn("nrow_space", paths)
        self.assertIn("ncolumn_space", paths)

    @patch("mpas_land_mesh.workflow.simplify_river_network")
    @patch("mpas_land_mesh.workflow.convert_vector_to_global_raster")
    @patch("mpas_land_mesh.workflow.create_land_ocean_mask_from_naturalearth")
    @patch("mpas_land_mesh.workflow.extend_land_mask_to_rivers")
    @patch("mpas_land_mesh.workflow.merge_features")
    def test_preprocess_with_execution(
        self,
        mock_merge,
        mock_extend,
        mock_coastline,
        mock_rasterize,
        mock_simplify,
    ):
        """Verify river and coastline pipelines execute when enabled."""
        raw_river_file = self.base_dir / "input" / "HydroRIVERS_v10.shp"
        raw_river_file.parent.mkdir(parents=True, exist_ok=True)
        raw_river_file.touch()

        mock_coastline.return_value = (
            str(self.out_dir / "land_mask.tif"),
            str(self.out_dir / "coast.geojson"),
        )

        cfg = dict(self.sample_config)
        cfg["preprocessing"] = {
            "simplify_rivers": True,
            "process_coastline": True,
            "largest_outlets": 2,
            "source": "naturalearth",
        }
        cfg["paths"]["hydrosheds_rivers"] = str(raw_river_file)

        wf = ULRMWorkflow(cfg)
        paths = wf.preprocess()

        mock_simplify.assert_called_once()
        mock_rasterize.assert_called_once()
        mock_coastline.assert_called_once()
        self.assertIsNotNone(paths)

    def test_compile_mesh_ready_config(self):
        """Verify compile_mesh_ready_config produces valid pyearthmesh config dict and file."""
        wf = ULRMWorkflow(self.sample_config)
        mesh_ready = wf.compile_mesh_ready_config(save_yaml=True)

        self.assertIn("workflow", mesh_ready)
        self.assertIn("resolutions", mesh_ready)
        self.assertIn("paths", mesh_ready)
        self.assertIn("jigsaw", mesh_ready)
        self.assertEqual(mesh_ready["workflow"]["mesh_type"], "mpas")
        self.assertEqual(mesh_ready["workflow"]["case_index"], 1)

        saved_yaml = self.out_dir / "mesh_ready_workflow.yaml"
        self.assertTrue(saved_yaml.is_file())

    @patch("mpas_land_mesh.workflow.create_workflow")
    def test_create_mesh_workflow(self, mock_create):
        """Verify create_mesh_workflow forwards compiled config to pyearthmesh."""
        mock_instance = MagicMock()
        mock_create.return_value = mock_instance

        wf = ULRMWorkflow(self.sample_config)
        result = wf.create_mesh_workflow()

        mock_create.assert_called_once()
        self.assertEqual(result, mock_instance)

    @patch("mpas_land_mesh.workflow.create_workflow")
    def test_run_local(self, mock_create):
        """Verify run executes setup and run for local execution."""
        mock_instance = MagicMock()
        mock_create.return_value = mock_instance

        wf = ULRMWorkflow(self.sample_config)
        wf.run()

        mock_instance.setup.assert_called_once()
        mock_instance.run.assert_called_once()

    @patch("mpas_land_mesh.workflow.create_workflow")
    def test_run_hpc(self, mock_create):
        """Verify run generates HPC SLURM scripts when run_on_hpc is True."""
        mock_instance = MagicMock()
        mock_instance.output_workspace = str(self.out_dir)
        mock_create.return_value = mock_instance

        cfg = dict(self.sample_config)
        cfg["workflow"] = dict(self.sample_config["workflow"])
        cfg["workflow"]["run_on_hpc"] = True

        wf = ULRMWorkflow(cfg)
        wf.run()

        mock_instance.create_hpc_job.assert_called_once_with(hours=5)
        mock_instance.setup.assert_not_called()
        mock_instance.run.assert_not_called()

    @patch("mpas_land_mesh.workflow.ULRMWorkflow.run")
    def test_run_workflow_helper(self, mock_run):
        """Verify top-level run_workflow function."""
        mock_run.return_value = "success"
        res = run_workflow(self.sample_config)
        mock_run.assert_called_once()
        self.assertEqual(res, "success")


if __name__ == "__main__":
    unittest.main()
