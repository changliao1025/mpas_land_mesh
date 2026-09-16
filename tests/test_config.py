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

from mpas_land_mesh.utilities.config import (
    DEFAULT_CONFIG,
    DEFAULT_JIGSAW_CONFIG,
    DEFAULT_WORKFLOW_CONFIG,
    JigsawConfigManager,
    WorkflowConfig,
    create_jigsaw_case,
    create_jigsaw_template_configuration_file,
    load_workflow_config,
    read_jigsaw_configuration_file,
)
from mpas_land_mesh.classes.jigsawcase import jigsawcase


class TestUnifiedConfig(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory()
        self.base_dir = Path(self.test_dir.name)

    def tearDown(self):
        self.test_dir.cleanup()

    def test_workflow_config_defaults(self):
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

        # Check it behaves as a dict
        self.assertIsInstance(config, dict)
        self.assertIsInstance(config, WorkflowConfig)
        self.assertEqual(config["workflow"]["case_index"], 5)
        self.assertTrue(config["workflow"]["process_coastline"])
        self.assertFalse(config["workflow"]["simplify_hydrosheds_river_network"])  # Default
        self.assertIsNotNone(config["workflow"]["date"])  # Auto-generated

        # Check resolutions with merged defaults
        self.assertEqual(config["resolutions"]["ocean_km"], 60)
        self.assertEqual(config["resolutions"]["land_km"], 30)
        self.assertEqual(config["resolutions"]["coastline_km"], 10)  # Default
        self.assertEqual(config["resolutions"]["river_network_km"], 10)  # Default

        # Check property access
        self.assertEqual(config.workflow["case_index"], 5)
        self.assertEqual(config.resolutions["ocean_km"], 60)
        self.assertEqual(config.config_file, cfg_path.resolve())

        # Check path resolution relative to config file directory
        expected_output = str((self.base_dir / "my_output").resolve())
        self.assertEqual(config.paths["output_workspace"], expected_output)

    def test_platform_specific_paths(self):
        current_os = platform.system()
        other_os = "Windows" if current_os != "Windows" else "Linux"

        cfg_path = self.base_dir / "paths.json"
        raw_data = {
            "paths": {
                "hydrosheds_rivers": {
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

        config = WorkflowConfig.load(cfg_path)
        expected_river = str((self.base_dir / "rivers_current.shp").resolve())
        expected_region = str((self.base_dir / "default_region.geojson").resolve())

        self.assertEqual(config.paths["hydrosheds_rivers"], expected_river)
        self.assertEqual(config.paths["region_geometry"], expected_region)

    def test_unknown_setting_raises_error(self):
        cfg_path = self.base_dir / "invalid.json"
        raw_data = {
            "workflow": {
                "unknown_key": 123,
            }
        }
        with open(cfg_path, "w", encoding="utf-8") as f:
            json.dump(raw_data, f)

        with self.assertRaises(ValueError):
            load_workflow_config(cfg_path)

    def test_build_jigsaw_config(self):
        cfg_path = self.base_dir / "workflow.json"
        raw_data = {
            "workflow": {
                "date": "20260901",
                "case_index": 2,
                "process_coastline": True,
                "simplify_hydrosheds_river_network": True,
            },
            "resolutions": {
                "ocean_km": 30,
                "land_km": 15,
                "river_network_km": 5,
                "coastline_km": 10,
            },
            "paths": {
                "output_workspace": str(self.base_dir / "output"),
            },
            "jigsaw": {
                "standalone": True,
                "overrides": {
                    "optm_iter": 48,
                    "mesh_rad2": 1.2,
                },
            },
        }
        with open(cfg_path, "w", encoding="utf-8") as f:
            json.dump(raw_data, f)

        config = WorkflowConfig.load(cfg_path)

        generated = {
            "river_network_vector": str(self.base_dir / "rivers.geojson"),
            "river_network_raster": str(self.base_dir / "rivers.tif"),
            "coastline_raster": str(self.base_dir / "coastline.tif"),
        }

        jigsaw_deck = config.build_jigsaw_config(generated_files=generated)

        # Verify grid dimension calculations
        dResolution_x_in = 30.0 / 3600.0 * 10.0
        expected_nrow = int(180.0 / dResolution_x_in)
        expected_ncolumn = int(360.0 / dResolution_x_in)

        self.assertEqual(jigsaw_deck["nrow_space"], expected_nrow)
        self.assertEqual(jigsaw_deck["ncolumn_space"], expected_ncolumn)
        self.assertEqual(jigsaw_deck["dResolution_ocean"], 30.0)
        self.assertEqual(jigsaw_deck["dResolution_land"], 15.0)
        self.assertEqual(jigsaw_deck["dResolution_river_network"], 5.0)
        self.assertEqual(jigsaw_deck["dResolution_coastline"], 10.0)

        # Verify feature flags
        self.assertTrue(jigsaw_deck["iFlag_geom"])
        self.assertTrue(jigsaw_deck["iFlag_spac"])
        self.assertTrue(jigsaw_deck["iFlag_spac_coastline"])
        self.assertTrue(jigsaw_deck["iFlag_geom_river_network"])
        self.assertTrue(jigsaw_deck["iFlag_spac_river_network"])

        # Verify generated file paths
        self.assertEqual(jigsaw_deck["sFilename_river_network_vector"], generated["river_network_vector"])
        self.assertEqual(jigsaw_deck["sFilename_river_network_raster"], generated["river_network_raster"])
        self.assertEqual(jigsaw_deck["sFilename_coastline_raster"], generated["coastline_raster"])

        # Verify metadata
        self.assertEqual(jigsaw_deck["iCase_index"], 2)
        self.assertEqual(jigsaw_deck["sDate"], "20260901")
        self.assertEqual(jigsaw_deck["sWorkspace_output"], str(self.base_dir / "output"))

        # Verify overrides
        self.assertEqual(jigsaw_deck["optm_iter"], 48)
        self.assertEqual(jigsaw_deck["mesh_rad2"], 1.2)

    def test_create_jigsaw_case_and_hpc_job(self):
        cfg_path = self.base_dir / "workflow.json"
        raw_data = {
            "workflow": {
                "date": "20260901",
                "case_index": 1,
            },
            "resolutions": {
                "coastline_km": 10,
            },
            "paths": {
                "output_workspace": str(self.base_dir / "output"),
            },
            "jigsaw": {
                "standalone": True,
            },
        }
        with open(cfg_path, "w", encoding="utf-8") as f:
            json.dump(raw_data, f)

        config = load_workflow_config(cfg_path)
        case = create_jigsaw_case(config, iFlag_create_directory_in=1)

        self.assertIsInstance(case, jigsawcase)
        self.assertEqual(case.iCase_index, 1)
        self.assertEqual(case.sDate, "20260901")

        # Test setup
        case.jigsaw_setup()
        self.assertTrue(Path(case.sWorkspace_output).is_dir())
        self.assertTrue((Path(case.sWorkspace_output) / "tmp").is_dir())
        self.assertTrue((Path(case.sWorkspace_output) / "out").is_dir())

        # Test HPC job generation creates jigsaw_configuration.json
        case._jigsaw_create_hpc_job(hours_in=2)
        self.assertTrue((Path(case.sWorkspace_output) / "jigsaw_configuration.json").is_file())
        self.assertTrue((Path(case.sWorkspace_output) / "run_jigsaw.py").is_file())
        self.assertTrue((Path(case.sWorkspace_output) / "submit.job").is_file())

    def test_jigsaw_config_manager(self):
        template_file = self.base_dir / "jigsaw_template.json"
        cfg = create_jigsaw_template_configuration_file(
            str(template_file),
            dResolution_ocean=50.0,
            iFlag_geom=True,
        )
        self.assertTrue(template_file.is_file())
        self.assertEqual(cfg["dResolution_ocean"], 50.0)
        self.assertTrue(cfg["iFlag_geom"])

        loaded = JigsawConfigManager.load_config(str(template_file))
        self.assertEqual(loaded["dResolution_ocean"], 50.0)

        # Test read_jigsaw_configuration_file creates jigsawcase
        case = read_jigsaw_configuration_file(
            str(template_file),
            iCase_index_in=3,
            sWorkspace_output_in=str(self.base_dir / "out_legacy"),
        )
        self.assertIsInstance(case, jigsawcase)
        self.assertEqual(case.iCase_index, 3)

    def test_backward_compatibility_imports(self):
        # Test importing from workflow_config
        from mpas_land_mesh.utilities.workflow_config import (
            DEFAULT_CONFIG as WC_DEFAULT,
            load_workflow_config as wc_load,
        )
        self.assertEqual(WC_DEFAULT, DEFAULT_WORKFLOW_CONFIG)

        # Test importing from config_manager
        from mpas_land_mesh.utilities.config_manager import (
            JigsawConfigManager as CM_JCM,
            create_jigsaw_template_configuration_file as cm_create,
            read_jigsaw_configuration_file as cm_read,
        )
        self.assertIs(CM_JCM, JigsawConfigManager)
        self.assertIs(cm_create, create_jigsaw_template_configuration_file)
        self.assertIs(cm_read, read_jigsaw_configuration_file)

    def test_load_existing_example_configs(self):
        workspace_root = Path(__file__).resolve().parent.parent

        examples = [
            workspace_root / "examples/minimal_workflow/minimal_workflow.json",
            workspace_root / "examples/ocn30coast10lnd10/ocn30coast10lnd10.json",
            workspace_root / "examples/qinghaihu/ocn100coast10lnd5riv4lak3.json",
        ]

        for example_path in examples:
            if example_path.exists():
                cfg = load_workflow_config(example_path)
                self.assertIsInstance(cfg, WorkflowConfig)
                deck = cfg.build_jigsaw_config()
                self.assertIsInstance(deck, dict)
                self.assertIn("ncolumn_space", deck)
                self.assertIn("nrow_space", deck)
                self.assertIn("dResolution_coastline", deck)


if __name__ == "__main__":
    unittest.main()
