#refer to the e3sm confluence page for more details of the workflow

#each step may require different input datasets
#however, all the output will be saved in the same output directory for easy access

#you can change this to your preferred output directory
import os
import argparse
import json
import logging
#using standalone dependency for this workflow, which will reduce the dependency to pyearth and pyflowline.
from mpas_land_mesh.utilities.vector import get_field_and_value, merge_features, add_field_to_vector_file
from mpas_land_mesh.utilities.raster import convert_vector_to_global_raster
from mpas_land_mesh.utilities.constants import KM2_TO_M2, ISLAND_AREA_MULTIPLIER, DRAINAGE_AREA_MULTIPLIER

from mpas_land_mesh.preprocessing.river_network import simplify_hydrorivers_network
from mpas_land_mesh.preprocessing.coastline import create_land_ocean_mask_from_naturalearth, fix_naturalearth_hydrosheds_incompatibility

from mpas_land_mesh.utilities.config import load_workflow_config, create_jigsaw_case

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('config', help='Path to the workflow JSON configuration')
args = parser.parse_args()
config = load_workflow_config(args.config)
workflow = config['workflow']
resolutions = config['resolutions']
paths = config['paths']
jigsaw = config['jigsaw']
sDate_today = workflow['date']
sMesh_type = workflow['mesh_type']
iCase_index = workflow['case_index']
sModel = workflow['model']

#resolution settings
#for rivers and watershed
dResolution_ocean = resolutions['ocean_km']
dResolution_land = resolutions['land_km']
dResolution_river_network = resolutions['river_network_km']
dResolution_coastline = resolutions['coastline_km']

#for coastline
dThreshold_area_island = dResolution_ocean * dResolution_ocean * ISLAND_AREA_MULTIPLIER * KM2_TO_M2  #unit m2, this one may need to be adjusted based on the resolution
dResolution_coastline_buffer = dResolution_coastline * 1.0E3  #buffer zone for coastline line
#small island removal threshold
dDrainage_area_threshold = dResolution_land * dResolution_land * DRAINAGE_AREA_MULTIPLIER * KM2_TO_M2  #at least 100 grid cells of drainage area, this may be adjusted as well

#setup flags for debugging
iFlag_simplify_hydrosheds_river_network = int(workflow['simplify_hydrosheds_river_network'])
iFlag_process_coastline = int(workflow['process_coastline'])

#number of largest outlet to be processed
nOutlet_largest = workflow['largest_outlets']

#thing may not need to be changed
sWorkspace_input = paths['input_workspace']
sWorkspace_output = paths['output_workspace']


#define global output directory
sWorkspace_river_network_output = os.path.join(sWorkspace_output, 'river_network')
if os.path.exists(sWorkspace_river_network_output) is False:
    os.makedirs(sWorkspace_river_network_output)


sWorkspace_coastline_output = os.path.join(sWorkspace_output, 'coastline')
if os.path.exists(sWorkspace_coastline_output) is False:
    os.makedirs(sWorkspace_coastline_output)

#for jigsaw resolution control
dResolution_x_in = 30.0/3600 * dResolution_coastline
dResolution_y_in = dResolution_x_in
nrow = int(180 / dResolution_y_in)
ncolumn = int(360 / dResolution_x_in)

#string format for file names
#river first
dDistance_tolerance = dResolution_river_network * 1.0E3 #how far away two river need to be for mesh generation

sDistance_tolerance = "{:.2E}".format(dDistance_tolerance)
sDrainage_area_threshold = "{:.2E}".format(dDrainage_area_threshold) # m2

#coastline second
sCoastline_buffer = "{:.1E}".format(dResolution_coastline_buffer  ) # to m
sThreshold_area_island = "{:.1E}".format(dThreshold_area_island ) # to m2


#add the threshol into the output folder
sWorkspace_river_network_output = os.path.join(sWorkspace_river_network_output,  sDistance_tolerance + '_' + sDrainage_area_threshold)
if os.path.exists(sWorkspace_river_network_output) is False:
    os.makedirs(sWorkspace_river_network_output)

sWorkspace_coastline_output = os.path.join(sWorkspace_coastline_output,  sCoastline_buffer + '_' + sThreshold_area_island )
if os.path.exists(sWorkspace_coastline_output) is False:
    os.makedirs(sWorkspace_coastline_output)


#Step 1
#prepare the river network and coastline line dataset

sFilename_flowline_hydrosheds_in = paths['hydrosheds_rivers']
sFilename_flowline_hydroshed_tmp = 'HydroRIVERS_v10_simplified_' + sDistance_tolerance + '_' + sDrainage_area_threshold + '.geojson'
sFilename_flowline_hydrosheds_out = os.path.join(sWorkspace_river_network_output, sFilename_flowline_hydroshed_tmp)
sFilename_geojson_geometery_feature = paths['region_geometry']

#step 1: record attribute from the MPAS tools
aField, aValue = get_field_and_value(sFilename_geojson_geometery_feature)

sFilename_dam = paths['dam_vector']
sFilename_river_network_raster = os.path.join(sWorkspace_river_network_output, 'river_network_raster.tif')
if iFlag_simplify_hydrosheds_river_network == 1:
    simplify_hydrorivers_network(sFilename_flowline_hydrosheds_in,
                       sFilename_flowline_hydrosheds_out,
                       dDistance_tolerance,
                        dDrainage_area_threshold,
                        nOutlet_largest=nOutlet_largest)
    convert_vector_to_global_raster(sFilename_flowline_hydrosheds_out, sFilename_river_network_raster,
                                         dResolution_x_in, dResolution_y_in )
else:
    #reuse the existing simplified river network for debug purpose, you can also change this to the original hydrosheds river network for testing
    pass

sFilename_vector_coastline_merged = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island_merged.geojson')
sFilename_tif_wo_island = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island.tif')
if iFlag_process_coastline == 1:
    sFilename_tif_wo_island, sFilename_vector_coastline = create_land_ocean_mask_from_naturalearth(sWorkspace_coastline_output,
                                                                             dResolution_x_in, dResolution_y_in,
                                                                             dThreshold_area_island,
                                                                             dResolution_coastline_buffer,
                                                                             iRaster_buffer_pixel = 2)

    ##we need to fix the incompatibilty between hydrosheds and naturalearth
    aFilename_flowline = list()
    for i in range(1, nOutlet_largest+1):
        sBasin_id = '{:04d}'.format(i)
        sFilename_flowline_simplified_basin = os.path.join(sWorkspace_river_network_output, 'HydroRIVERS_v10_simplified_' + sDistance_tolerance + '_' + sDrainage_area_threshold +'_'+ sBasin_id + '.geojson')
        aFilename_flowline.append(sFilename_flowline_simplified_basin)

    sFilename_vector_coastline_updated = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island_fixed.geojson')
    fix_naturalearth_hydrosheds_incompatibility(aFilename_flowline, sFilename_vector_coastline, sFilename_vector_coastline_updated )
    #should be merged into one single function
    merge_features(sFilename_vector_coastline_updated, sFilename_vector_coastline_merged, iFlag_force= True)
    add_field_to_vector_file(sFilename_vector_coastline_merged, aField, aValue)
else:
    #reuse
    pass

#Step 2 - 4
#run the hexwatershed model, this step include three steps merged together.
#for debug purpose, you can also run then one by one, using the iFlag_debug flag to control
#sFilename_mpas_mesh_netcdf = '/compyfs/liao313/04model/pyhexwatershed/global/pyflowline20251122001/jigsaw/out/invert_mesh.nc'

iFlag_debug = int(workflow['debug'])
if iFlag_debug == 1:
    generated_files = {
        'river_network_vector': sFilename_flowline_hydrosheds_out,
        'river_network_raster': sFilename_river_network_raster,
        'coastline_raster': sFilename_tif_wo_island,
        'dam_vector': sFilename_dam,
    }

    # Build JIGSAW deck in memory and create case instance directly
    oJigsaw = create_jigsaw_case(
        config,
        generated_files=generated_files,
        output_workspace=sWorkspace_output,
        iFlag_create_directory_in=1,
    )

    oJigsaw._jigsaw_create_hpc_job(sSlurm_in=jigsaw['slurm'], hours_in=jigsaw['hours'])
    # now you should manually submit the job

else:
    pass

logger.info('='*80)
logger.info('Workflow completed successfully!')
logger.info('='*80)
