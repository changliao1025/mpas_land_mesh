#refer to the e3sm confluence page for more details of the workflow

#each step may require different input datasets
#however, all the output will be saved in the same output directory for easy access

#you can change this to your preferred output directory
import os
import glob
import logging
from datetime import datetime
from shutil import copy2

#using standalone dependency for this workflow, which will reduce the dependency to pyearth and pyflowline.
from mpas_land_mesh.utilities.change_json_key_value import change_json_key_value
from mpas_land_mesh.utilities.vector import get_field_and_value, merge_features, add_field_to_vector_file
from mpas_land_mesh.utilities.raster import convert_vector_to_global_raster
from mpas_land_mesh.utilities.constants import KM2_TO_M2, ISLAND_AREA_MULTIPLIER, DRAINAGE_AREA_MULTIPLIER

from mpas_land_mesh.preprocessing.river_networks import simplify_hydrorivers_networks
from mpas_land_mesh.preprocessing.coastlines import create_land_ocean_mask_from_naturalearth, fix_naturalearth_hydrosheds_incompatibility

from mpas_land_mesh.utilities.config_manager import create_jigsaw_template_configuration_file, read_jigsaw_configuration_file

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

##======================================================================
# The only thing you need to change for different runs is the settings below,
# which is used to control the resolution and other settings for the mesh generation. You can also set up different flags to turn on/off certain process for debugging purpose. The output will be saved in the same output directory, which is defined below as well.
# start of common user settings
##======================================================================
sDate_today = datetime.now().strftime('%Y%m%d')

#things that may need to be changed for different runs
sDate_today = '20260402'  #use a fixed date for easy repeatability
sMesh_type = 'mpas'  #
#index for different runs
iCase_index = 3
sModel = 'jigsaw'
#flag for component

iFlag_simplify_hydrosheds_river_network = 0
iFlag_process_coastline = 0

#resolution settings

dResolution_ocean = 30
dResolution_land = 10  #unit in km
dResolution_river_network = 10
dResolution_coastline = 10  #unit in km

sWorkspace_output = '/compyfs/liao313/04model/pyhexwatershed/global/'

sFilename_flowline_hydrosheds_in = '/compyfs/liao313/00raw/hydrology/hydrosheds/hydroriver/HydroRIVERS_v10_shp/HydroRIVERS_v10_shp/HydroRIVERS_v10.shp'
sFilename_geojson_geometery_feature = '/qfs/people/liao313/data/hexwatershed/global/vector/region.geojson'

##======================================================================
# end of common user settings
##======================================================================
##======================================================================
# start uncommon user settings
##======================================================================
iFlag_dam = 0
sFilename_dam = '/compyfs/liao313/00raw/dam/GRanD_Version_1_3/GRanD_dams_v1_3_merged.geojson' #should consider both on and snapped dams in this dataset

##======================================================================
# end of uncommon user settings
##======================================================================
#for coastline
dThreshold_area_island = dResolution_ocean * dResolution_ocean * ISLAND_AREA_MULTIPLIER * KM2_TO_M2  #unit m2, this one may need to be adjusted based on the resolution
dResolution_coastline_buffer = dResolution_coastline * 1.0E3  #buffer zone for coastline line
#small island removal threshold
dDrainage_area_threshold = dResolution_land * dResolution_land * DRAINAGE_AREA_MULTIPLIER * KM2_TO_M2  #at least 100 grid cells of drainage area, this may be adjusted as well


#number of largest outlet to be processed
nOutlet_largest = 100




#define global output directory
sWorkspace_river_network_output = '/compyfs/liao313/04model/pyhexwatershed/global/river_network'
if os.path.exists(sWorkspace_river_network_output) is False:
    os.makedirs(sWorkspace_river_network_output)


sWorkspace_coastline_output = '/compyfs/liao313/04model/pyhexwatershed/global/coastline'
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


sFilename_flowline_hydroshed_tmp = 'HydroRIVERS_v10_simplified_' + sDistance_tolerance + '_' + sDrainage_area_threshold + '.geojson'
sFilename_flowline_hydrosheds_out = os.path.join(sWorkspace_river_network_output, sFilename_flowline_hydroshed_tmp)
sWorkspace_watershed_boundary_in = '/compyfs/liao313/00raw/hydrology/hydrosheds/hydrobasin'
#step 1: record attribute from the MPAS tools
aField, aValue = get_field_and_value(sFilename_geojson_geometery_feature)

sFilename_river_network_raster = os.path.join(sWorkspace_river_network_output, 'river_network_raster.tif')
if iFlag_simplify_hydrosheds_river_network == 1:
    nOutlet_actual = simplify_hydrorivers_networks(sFilename_flowline_hydrosheds_in,
                       sFilename_flowline_hydrosheds_out,
                       dDistance_tolerance,
                        dDrainage_area_threshold,
                        nOutlet_largest=nOutlet_largest)

    logger.info(f'Processed {nOutlet_actual} largest outlet basins')
    logger.info(f'Output filename: {sFilename_flowline_hydrosheds_out}')

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

sFilename_jigsaw_configuration_json = os.path.join(sWorkspace_river_network_output, 'jigsaw_configuration.json')

try:
    create_jigsaw_template_configuration_file(sFilename_jigsaw_configuration_json)
    change_json_key_value(sFilename_jigsaw_configuration_json, "sWorkspace_output", sWorkspace_output)
    oJigsaw = read_jigsaw_configuration_file(sFilename_jigsaw_configuration_json, \
    iCase_index_in=iCase_index, sDate_in=sDate_today, iFlag_create_directory_in=1)
    sWorkspace_output_case = oJigsaw.sWorkspace_output
    sFilename_jigsaw_configuration_copy = os.path.join( sWorkspace_output_case, 'jigsaw_configuration_copy.json' )
    copy2(sFilename_jigsaw_configuration_json, sFilename_jigsaw_configuration_copy)

    #update the jigsaw configuration file below
    change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_geom", "true") # enable geometry control
    change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_geom_river_network", "true") #set the resolution
    change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_spac", "true") #enable resolution control
    change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_spac_ocean", "true")
    #change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_spac_river_network", "true")
    #change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_RRS18to6_ocean", "true")
    change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_spac_land", "true")
    change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_spac_coastline", "true") #set the resolution for coastline line
    change_json_key_value(sFilename_jigsaw_configuration_copy, "dResolution_ocean", dResolution_ocean) #set the resolution for ocean
    change_json_key_value(sFilename_jigsaw_configuration_copy, "dResolution_coastline", dResolution_coastline) #set the resolution for coastline line
    change_json_key_value(sFilename_jigsaw_configuration_copy, "dResolution_land", dResolution_land) #set the resolution for land
    change_json_key_value(sFilename_jigsaw_configuration_copy, "dResolution_river_network", dResolution_river_network) #set the resolution for river network
    change_json_key_value(sFilename_jigsaw_configuration_copy, "ncolumn_space", ncolumn) #set the resolution for x direction
    change_json_key_value(sFilename_jigsaw_configuration_copy, "nrow_space", nrow) #set the resolution for y direction
    change_json_key_value(sFilename_jigsaw_configuration_copy, "sFilename_river_network_vector", sFilename_flowline_hydrosheds_out) #set the resolution for x direction
    change_json_key_value(sFilename_jigsaw_configuration_copy, "sFilename_river_network_raster", sFilename_river_network_raster) #set the small island removal threshold
    change_json_key_value(sFilename_jigsaw_configuration_copy, "sFilename_land_ocean_mask", sFilename_vector_coastline_merged) #set the land ocean mask
    change_json_key_value(sFilename_jigsaw_configuration_copy, "sFilename_coastline_raster", sFilename_tif_wo_island) #set the resolution for y direction
    if iFlag_dam == 1:
        change_json_key_value(sFilename_jigsaw_configuration_copy, "iFlag_geom_dam", "true")
        change_json_key_value(sFilename_jigsaw_configuration_copy, "sFilename_dam_vector", sFilename_dam) #set the dam file

    #now we can set up the actual pyflowline to create the mesh
    oJigsaw = read_jigsaw_configuration_file(sFilename_jigsaw_configuration_copy,
                    iCase_index_in=iCase_index,
                    sDate_in= sDate_today)

    oJigsaw._jigsaw_create_hpc_job(sSlurm_in = 'slurm', hours_in = 5 )
    #now you should manually submit the job

except Exception as e:
    logger.error(f"An error occurred: {e}")

logger.info('='*80)
logger.info('Workflow completed successfully!')
logger.info('='*80)
