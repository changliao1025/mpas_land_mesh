"""Unified Land-River Mesh (ULRM) Workflow Pipeline.

Orchestrates raw data preprocessing (river network simplification, coastline/ocean
mask generation, coastal flowline reconciliation) and compiles a mesh-ready
specification for execution by pyearthmesh.
"""

from __future__ import annotations

import glob
import logging
import os
from pathlib import Path
from typing import Any, Union

from hexwatershed_utility.preprocess.feature.coastline import (
    create_land_ocean_mask_from_hydrobasins,
    create_land_ocean_mask_from_naturalearth,
    create_land_ocean_mask_from_shoreline,
    extend_land_mask_to_rivers,
)
from hexwatershed_utility.preprocess.feature.river import simplify_river_network
from pyearth.toolbox.conversion.convert_vector_to_global_raster import (
    convert_vector_to_global_raster,
)
from pyearth.toolbox.management.vector.fields import (
    add_field_to_vector_file,
    get_field_and_value,
)
from pyearth.toolbox.management.vector.merge_features import merge_features
from pyearthmesh.config.config_manager import create_workflow

from mpas_land_mesh.config import dump_ulrm_config, load_ulrm_config
from mpas_land_mesh.utilities.constants import (
    DRAINAGE_AREA_MULTIPLIER,
    ISLAND_AREA_MULTIPLIER,
    KM2_TO_M2,
)

logger = logging.getLogger(__name__)


class ULRMWorkflow:
    """Unified Land-River Mesh (ULRM) end-to-end workflow runner."""

    def __init__(self, config_or_path: str | Path | dict[str, Any]) -> None:
        """Initialize workflow with a path to YAML/JSON or a config dict."""
        self.config: dict[str, Any] = load_ulrm_config(config_or_path)
        self.preprocessed_paths: dict[str, Any] = {}
        self._mesh_ready_config: dict[str, Any] | None = None

    @property
    def output_workspace(self) -> str:
        return self.config["paths"]["output_workspace"]

    def preprocess(self) -> dict[str, Any]:
        """Execute river network and coastline preprocessing pipelines."""
        resolutions = self.config["resolutions"]
        preprocessing = self.config["preprocessing"]
        paths = self.config["paths"]

        dResolution_ocean = float(resolutions["ocean_km"])
        dResolution_land = float(resolutions["land_km"])
        dResolution_river_network = float(
            resolutions.get("river_network_km") or dResolution_land
        )
        dResolution_coastline = float(
            resolutions.get("coastline_km") or dResolution_land
        )

        island_mult = float(
            preprocessing.get("island_area_multiplier", ISLAND_AREA_MULTIPLIER)
        )
        drainage_mult = float(
            preprocessing.get("drainage_area_multiplier", DRAINAGE_AREA_MULTIPLIER)
        )
        nOutlet_largest = int(preprocessing.get("largest_outlets", 10))
        iRaster_buffer_pixel = int(preprocessing.get("raster_buffer_pixel", 2))

        dThreshold_area_island = (
            dResolution_ocean * dResolution_ocean * island_mult * KM2_TO_M2
        )
        dResolution_coastline_buffer = dResolution_coastline * 1.0e3
        dDrainage_area_threshold = (
            dResolution_land * dResolution_land * drainage_mult * KM2_TO_M2
        )
        dDistance_tolerance = dResolution_river_network * 1.0e3

        dResolution_x_in = 30.0 / 3600.0 * dResolution_coastline
        dResolution_y_in = dResolution_x_in

        sDistance_tolerance = f"{dDistance_tolerance:.2E}"
        sDrainage_area_threshold = f"{dDrainage_area_threshold:.2E}"
        sCoastline_buffer = f"{dResolution_coastline_buffer:.1E}"
        sThreshold_area_island = f"{dThreshold_area_island:.1E}"

        sWorkspace_output = paths["output_workspace"]
        sWorkspace_river_network_output = os.path.join(
            sWorkspace_output,
            "river_network",
            f"{sDistance_tolerance}_{sDrainage_area_threshold}",
        )
        sWorkspace_coastline_output = os.path.join(
            sWorkspace_output,
            "coastline",
            f"{sCoastline_buffer}_{sThreshold_area_island}",
        )
        os.makedirs(sWorkspace_river_network_output, exist_ok=True)
        os.makedirs(sWorkspace_coastline_output, exist_ok=True)

        # -------------------------------------------------------------------
        # 1. River Network Preprocessing
        # -------------------------------------------------------------------
        sFilename_flowline_raw_in = (
            paths.get("river_network_raw")
            or paths.get("raw_rivers")
            or paths.get("hydrosheds_rivers")
            or ""
        )
        bSimplify_rivers = bool(preprocessing.get("simplify_rivers", False))

        if (
            not sFilename_flowline_raw_in
            and bSimplify_rivers
            and paths.get("river_network_vector")
        ):
            sFilename_flowline_raw_in = paths["river_network_vector"]

        if sFilename_flowline_raw_in:
            raw_stem = Path(sFilename_flowline_raw_in).stem
            if "simplified" in raw_stem.lower():
                base_name = (
                    f"{raw_stem}_{sDistance_tolerance}_{sDrainage_area_threshold}"
                )
            else:
                base_name = f"{raw_stem}_simplified_{sDistance_tolerance}_{sDrainage_area_threshold}"
        else:
            base_name = f"river_network_simplified_{sDistance_tolerance}_{sDrainage_area_threshold}"

        sFilename_flowline_simplified_out = os.path.join(
            sWorkspace_river_network_output, f"{base_name}.geojson"
        )
        sFilename_river_network_raster = os.path.join(
            sWorkspace_river_network_output, "river_network_raster.tif"
        )
        sFilename_flowline_out = sFilename_flowline_simplified_out

        if (
            bSimplify_rivers
            and sFilename_flowline_raw_in
            and os.path.exists(sFilename_flowline_raw_in)
        ):
            logger.info("Simplifying river network from: %s", sFilename_flowline_raw_in)

            # Optional column mappings for general hydrographic datasets (HydroRIVERS, NHDPlus, MERIT, etc.)
            river_cols = preprocessing.get("river_columns") or {}
            sCol_id = preprocessing.get("col_id") or river_cols.get("col_id")
            sCol_main_riv = preprocessing.get("col_main_riv") or river_cols.get(
                "col_main_riv"
            )
            sCol_stream_order = preprocessing.get("col_stream_order") or river_cols.get(
                "col_stream_order"
            )
            sCol_drainage_area = preprocessing.get(
                "col_drainage_area"
            ) or river_cols.get("col_drainage_area")
            sCol_next_down = preprocessing.get("col_next_down") or river_cols.get(
                "col_next_down"
            )
            sCol_endorheic = preprocessing.get("col_endorheic") or river_cols.get(
                "col_endorheic"
            )
            dDrainage_multiplier = float(
                preprocessing.get("drainage_area_multiplier_unit")
                or preprocessing.get("drainage_area_unit_multiplier")
                or 1.0e6
            )
            iFlag_endorheic = int(preprocessing.get("endorheic", 0))

            simplify_river_network(
                sFilename_flowline_hydroshed_in=sFilename_flowline_raw_in,
                sFilename_flowline_hydroshed_out=sFilename_flowline_simplified_out,
                dDistance_tolerance_in=dDistance_tolerance,
                dDrainage_area_threshold_in=dDrainage_area_threshold,
                iFlag_endorheic_in=iFlag_endorheic,
                nOutlet_largest=nOutlet_largest,
                sCol_id=sCol_id,
                sCol_main_riv=sCol_main_riv,
                sCol_stream_order=sCol_stream_order,
                sCol_drainage_area=sCol_drainage_area,
                sCol_next_down=sCol_next_down,
                sCol_endorheic=sCol_endorheic,
                dDrainage_area_multiplier=dDrainage_multiplier,
            )
            logger.info(
                "Rasterizing simplified river network to: %s",
                sFilename_river_network_raster,
            )
            convert_vector_to_global_raster(
                sFilename_flowline_simplified_out,
                sFilename_river_network_raster,
                dResolution_x_in,
                dResolution_y_in,
            )
            sFilename_flowline_out = sFilename_flowline_simplified_out
        elif os.path.exists(sFilename_flowline_simplified_out) and os.path.exists(
            sFilename_river_network_raster
        ):
            logger.info(
                "Reusing existing river network outputs in: %s",
                sWorkspace_river_network_output,
            )
            sFilename_flowline_out = sFilename_flowline_simplified_out
        elif paths.get("river_network_vector"):
            sFilename_flowline_out = paths["river_network_vector"]
            sFilename_river_network_raster = paths.get("river_network_raster", "")
            if (
                not sFilename_river_network_raster
                or not os.path.exists(sFilename_river_network_raster)
            ) and os.path.exists(sFilename_flowline_out):
                target_raster = os.path.join(
                    sWorkspace_river_network_output, "river_network_raster.tif"
                )
                logger.info(
                    "Rasterizing provided river network vector to: %s", target_raster
                )
                try:
                    convert_vector_to_global_raster(
                        sFilename_flowline_out,
                        target_raster,
                        dResolution_x_in,
                        dResolution_y_in,
                    )
                    sFilename_river_network_raster = target_raster
                except Exception as e:
                    logger.warning("Could not auto-rasterize river network: %s", e)
        elif sFilename_flowline_raw_in and os.path.exists(sFilename_flowline_raw_in):
            logger.info(
                "Using raw river network directly without simplification: %s",
                sFilename_flowline_raw_in,
            )
            sFilename_flowline_out = sFilename_flowline_raw_in
            target_raster = os.path.join(
                sWorkspace_river_network_output, "river_network_raster.tif"
            )
            if not os.path.exists(target_raster):
                logger.info("Rasterizing river network vector to: %s", target_raster)
                try:
                    convert_vector_to_global_raster(
                        sFilename_flowline_out,
                        target_raster,
                        dResolution_x_in,
                        dResolution_y_in,
                    )
                    sFilename_river_network_raster = target_raster
                except Exception as e:
                    logger.warning("Could not auto-rasterize river network: %s", e)

        # -------------------------------------------------------------------
        # 2. Coastline & Land-Ocean Mask Preprocessing
        # -------------------------------------------------------------------
        sFilename_vector_coastline_merged = os.path.join(
            sWorkspace_coastline_output, "land_ocean_mask_wo_island_merged.geojson"
        )
        sFilename_tif_wo_island = os.path.join(
            sWorkspace_coastline_output, "land_ocean_mask_wo_island.tif"
        )

        bProcess_coastline = bool(preprocessing.get("process_coastline", False))
        if bProcess_coastline:
            mask_source = preprocessing.get("source", "naturalearth")
            logger.info(
                "Processing coastline and land-ocean mask (source=%s)...", mask_source
            )
            if mask_source == "naturalearth":
                sFilename_tif_wo_island, sFilename_vector_coastline = (
                    create_land_ocean_mask_from_naturalearth(
                        sWorkspace_coastline_output=sWorkspace_coastline_output,
                        dResolution_x_in=dResolution_x_in,
                        dResolution_y_in=dResolution_y_in,
                        dThreshold_area_island=dThreshold_area_island,
                        dResolution_coastline_buffer=dResolution_coastline_buffer,
                        iRaster_buffer_pixel=iRaster_buffer_pixel,
                    )
                )
            elif mask_source in ("basins", "hydrobasins"):
                sFilename_tif_wo_island, sFilename_vector_coastline = (
                    create_land_ocean_mask_from_hydrobasins(
                        sWorkspace_output=sWorkspace_coastline_output,
                        sFolder_boundary_in=paths.get("basins_boundary"),
                        dResolution_x=dResolution_x_in,
                        dResolution_y=dResolution_y_in,
                        dThreshold_area_island=dThreshold_area_island,
                        dResolution_coastline_buffer=dResolution_coastline_buffer,
                    )
                )
            else:  # 'shoreline'
                sFilename_tif_wo_island, sFilename_vector_coastline = (
                    create_land_ocean_mask_from_shoreline(
                        sWorkspace_output=sWorkspace_coastline_output,
                        sFilename_land_vector_in=paths.get("naturalearth_land"),
                        dResolution_x=dResolution_x_in,
                        dResolution_y=dResolution_y_in,
                        dThreshold_area_island=dThreshold_area_island,
                        iRaster_buffer_pixel=iRaster_buffer_pixel,
                        iFlag_remove_internal_rings=True,
                    )
                )

            # Reconcile coastline with river flowline outlets
            aFilename_flowline = []
            for i in range(1, nOutlet_largest + 1):
                sBasin_id = f"{i:04d}"
                sFilename_flowline_basin = sFilename_flowline_out.replace(
                    ".geojson", f"_{sBasin_id}.geojson"
                )
                if not os.path.exists(sFilename_flowline_basin):
                    legacy_basin = os.path.join(
                        sWorkspace_river_network_output,
                        f"HydroRIVERS_v10_simplified_{sDistance_tolerance}_{sDrainage_area_threshold}_{sBasin_id}.geojson",
                    )
                    if os.path.exists(legacy_basin):
                        sFilename_flowline_basin = legacy_basin
                    else:
                        candidates = glob.glob(
                            os.path.join(
                                sWorkspace_river_network_output,
                                f"*_{sBasin_id}.geojson",
                            )
                        )
                        if candidates:
                            sFilename_flowline_basin = candidates[0]
                if os.path.exists(sFilename_flowline_basin):
                    aFilename_flowline.append(sFilename_flowline_basin)

            if not aFilename_flowline and os.path.exists(sFilename_flowline_out):
                aFilename_flowline.append(sFilename_flowline_out)

            if aFilename_flowline and os.path.exists(sFilename_vector_coastline):
                logger.info(
                    "Reconciling coastline with %d flowline dataset(s)...",
                    len(aFilename_flowline),
                )
                sFilename_vector_coastline_updated = os.path.join(
                    sWorkspace_coastline_output,
                    "land_ocean_mask_wo_island_fixed.geojson",
                )
                extend_land_mask_to_rivers(
                    aFilename_flowlines=aFilename_flowline,
                    sFilename_vector_land=sFilename_vector_coastline,
                    sFilename_vector_land_updated=sFilename_vector_coastline_updated,
                )
                merge_features(
                    sFilename_vector_coastline_updated,
                    sFilename_vector_coastline_merged,
                    iFlag_force=True,
                )
            else:
                sFilename_vector_coastline_merged = sFilename_vector_coastline

            # Tag attributes from region geometry if available
            sFilename_region = paths.get("region_geometry")
            if (
                sFilename_region
                and os.path.exists(sFilename_region)
                and os.path.exists(sFilename_vector_coastline_merged)
            ):
                aField, aValue = get_field_and_value(sFilename_region)
                if aField and aValue:
                    add_field_to_vector_file(
                        sFilename_vector_coastline_merged, aField, aValue
                    )

        elif os.path.exists(sFilename_vector_coastline_merged) and os.path.exists(
            sFilename_tif_wo_island
        ):
            logger.info(
                "Reusing existing coastline outputs in: %s",
                sWorkspace_coastline_output,
            )
        elif paths.get("coastline_raster"):
            sFilename_tif_wo_island = paths["coastline_raster"]
            sFilename_vector_coastline_merged = paths.get("coastline_vector", "")

        self.preprocessed_paths = {
            "river_network_raw": sFilename_flowline_raw_in,
            "river_network_vector": sFilename_flowline_out,
            "river_network_raster": sFilename_river_network_raster,
            "coastline_raster": sFilename_tif_wo_island,
            "coastline_vector": sFilename_vector_coastline_merged,
            "dam_vector": paths.get("dam_vector", ""),
            "region_geometry": paths.get("region_geometry", ""),
            "output_workspace": sWorkspace_output,
            "input_workspace": paths.get("input_workspace", ""),
            "nrow_space": int(180 / dResolution_y_in),
            "ncolumn_space": int(360 / dResolution_x_in),
        }
        return self.preprocessed_paths

    def compile_mesh_ready_config(self, save_yaml: bool = True) -> dict[str, Any]:
        """Compile a mesh-ready configuration dictionary for pyearthmesh.

        Evaluates jigsaw flags (such as iFlag_geom_river_network, iFlag_spac_river_network,
        iFlag_spac_watershed_boundary, iFlag_geom_lake_boundary, etc.) to determine
        which preprocessed or direct user-provided vectors and rasters are included
        in the mesh-ready specification.

        Parameters
        ----------
        save_yaml : bool, default=True
            If True, writes 'mesh_ready_workflow.yaml' to output_workspace.

        Returns
        -------
        dict[str, Any]
            Mesh-ready configuration dictionary compatible with pyearthmesh.
        """
        if not self.preprocessed_paths:
            self.preprocess()

        wf = self.config["workflow"]
        jigsaw_settings = dict(self.config.get("jigsaw", {}))
        jigsaw_settings["nrow_space"] = self.preprocessed_paths.get(
            "nrow_space", jigsaw_settings.get("nrow_space", 180)
        )
        jigsaw_settings["ncolumn_space"] = self.preprocessed_paths.get(
            "ncolumn_space", jigsaw_settings.get("ncolumn_space", 360)
        )

        raw_paths = self.config.get("paths", {})

        def _get_path(key: str) -> str:
            val = self.preprocessed_paths.get(key)
            if not val:
                val = raw_paths.get(key, "")
            return str(val) if val else ""

        # Master JIGSAW flags
        bGeom = bool(jigsaw_settings.get("iFlag_geom", True))
        bSpac = bool(jigsaw_settings.get("iFlag_spac", True))

        compiled_paths: dict[str, str] = {
            "input_workspace": _get_path("input_workspace"),
            "output_workspace": _get_path("output_workspace") or self.output_workspace,
        }

        # 1. River Network (preprocessed or direct)
        flag_geom_river = bool(jigsaw_settings.get("iFlag_geom_river_network", True))
        if bGeom and flag_geom_river:
            p = _get_path("river_network_vector")
            if p:
                compiled_paths["river_network_vector"] = p
            else:
                logger.warning(
                    "iFlag_geom_river_network is True but river_network_vector is not provided."
                )

        flag_spac_river = bool(jigsaw_settings.get("iFlag_spac_river_network", True))
        if bSpac and flag_spac_river:
            p = _get_path("river_network_raster")
            if p:
                compiled_paths["river_network_raster"] = p
            else:
                logger.warning(
                    "iFlag_spac_river_network is True but river_network_raster is not provided."
                )

        # 2. Coastline (preprocessed or direct)
        flag_geom_coast = bool(jigsaw_settings.get("iFlag_geom_coastline", False))
        if bGeom and flag_geom_coast:
            p = _get_path("coastline_vector")
            if p:
                compiled_paths["coastline_vector"] = p
            else:
                logger.warning(
                    "iFlag_geom_coastline is True but coastline_vector is not provided."
                )

        flag_spac_coast = bool(
            jigsaw_settings.get("iFlag_spac_coastline", False)
            or jigsaw_settings.get("iFlag_spac_ocean", False)
        )
        if bSpac and flag_spac_coast:
            p = _get_path("coastline_raster")
            if p:
                compiled_paths["coastline_raster"] = p
            else:
                logger.warning(
                    "iFlag_spac_coastline/ocean is True but coastline_raster is not provided."
                )

        # 3. Watershed boundary (direct from original config paths)
        flag_geom_wshed = bool(
            jigsaw_settings.get("iFlag_geom_watershed_boundary", False)
        )
        if bGeom and flag_geom_wshed:
            p = _get_path("watershed_boundary_vector")
            if p:
                compiled_paths["watershed_boundary_vector"] = p
            else:
                logger.warning(
                    "iFlag_geom_watershed_boundary is True but watershed_boundary_vector is not provided."
                )

        flag_spac_wshed = bool(
            jigsaw_settings.get("iFlag_spac_watershed_boundary", False)
        )
        if bSpac and flag_spac_wshed:
            p = _get_path("watershed_boundary_raster")
            if p:
                compiled_paths["watershed_boundary_raster"] = p
            else:
                logger.warning(
                    "iFlag_spac_watershed_boundary is True but watershed_boundary_raster is not provided."
                )

        # 4. Lake boundary (direct from original config paths)
        flag_geom_lake = bool(jigsaw_settings.get("iFlag_geom_lake_boundary", False))
        if bGeom and flag_geom_lake:
            p = _get_path("lake_boundary_vector")
            if p:
                compiled_paths["lake_boundary_vector"] = p
            else:
                logger.warning(
                    "iFlag_geom_lake_boundary is True but lake_boundary_vector is not provided."
                )

        flag_spac_lake = bool(jigsaw_settings.get("iFlag_spac_lake_boundary", False))
        if bSpac and flag_spac_lake:
            p = _get_path("lake_boundary_raster")
            if p:
                compiled_paths["lake_boundary_raster"] = p
            else:
                logger.warning(
                    "iFlag_spac_lake_boundary is True but lake_boundary_raster is not provided."
                )

        # 5. Dam (direct from original config paths)
        flag_geom_dam = bool(jigsaw_settings.get("iFlag_geom_dam", False))
        if (bGeom and flag_geom_dam) or (
            raw_paths.get("dam_vector")
            and jigsaw_settings.get("iFlag_geom_dam") is not False
        ):
            p = _get_path("dam_vector")
            if p:
                compiled_paths["dam_vector"] = p

        flag_spac_dam = bool(jigsaw_settings.get("iFlag_spac_dam", False))
        if bSpac and flag_spac_dam:
            p = _get_path("dam_raster")
            if p:
                compiled_paths["dam_raster"] = p

        # 6. City (direct from original config paths)
        flag_geom_city = bool(jigsaw_settings.get("iFlag_geom_city", False))
        if bGeom and flag_geom_city:
            p = _get_path("city_vector")
            if p:
                compiled_paths["city_vector"] = p

        flag_spac_city = bool(jigsaw_settings.get("iFlag_spac_city", False))
        if bSpac and flag_spac_city:
            p = _get_path("city_raster")
            if p:
                compiled_paths["city_raster"] = p

        # 7. Region geometry and optional output mesh
        if raw_paths.get("region_geometry"):
            compiled_paths["region_geometry"] = _get_path("region_geometry")
        if raw_paths.get("output_mesh"):
            compiled_paths["output_mesh"] = _get_path("output_mesh")

        mesh_ready = {
            "workflow": {
                "case_index": wf.get("case_index", 1),
                "date": wf.get("date"),
                "mesh_type": wf.get("mesh_type", "mpas"),
                "model": wf.get("model", "jigsaw"),
                "debug": wf.get("debug", True),
                "run_on_hpc": wf.get("run_on_hpc", False),
            },
            "resolutions": dict(self.config["resolutions"]),
            "paths": compiled_paths,
            "jigsaw": jigsaw_settings,
        }

        if save_yaml:
            out_file = os.path.join(
                self.preprocessed_paths["output_workspace"],
                "mesh_ready_workflow.yaml",
            )
            try:
                dump_ulrm_config(mesh_ready, out_file)
                logger.info("Saved mesh-ready configuration to: %s", out_file)
            except Exception as e:
                logger.warning("Could not save mesh_ready_workflow.yaml: %s", e)

        self._mesh_ready_config = mesh_ready
        return mesh_ready

    def create_mesh_workflow(self) -> Any:
        """Create a pyearthmesh workflow instance from the compiled mesh-ready config."""
        mesh_ready = self.compile_mesh_ready_config()
        return create_workflow(mesh_ready)

    def run(self) -> Any:
        """Execute the full end-to-end pipeline (preprocess -> mesh generation)."""
        logger.info("Starting ULRM Workflow...")
        self.preprocess()
        mesh_workflow = self.create_mesh_workflow()

        if self.config["workflow"].get("run_on_hpc", False):
            hours = self.config.get("jigsaw", {}).get("hours", 5)
            logger.info("Generating HPC SLURM job scripts (walltime=%dh)...", hours)
            mesh_workflow.create_hpc_job(hours=hours)
            logger.info(
                "HPC job scripts created in: %s", mesh_workflow.output_workspace
            )
            return mesh_workflow
        else:
            logger.info("Executing mesh generation via pyearthmesh...")
            mesh_workflow.setup()
            mesh_workflow.run()
            logger.info("Mesh generation completed successfully!")
            return mesh_workflow


# Alias for backward compatibility
MPASLandMeshWorkflow = ULRMWorkflow


def run_workflow(config_or_path: str | Path | dict[str, Any]) -> Any:
    """Convenience function to instantiate and run a ULRMWorkflow."""
    workflow = ULRMWorkflow(config_or_path)
    return workflow.run()
