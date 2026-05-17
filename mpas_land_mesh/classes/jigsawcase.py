import os
import stat
import json
import datetime
from pathlib import Path
from shutil import copy2
from mpas_land_mesh.utilities.system import get_python_environment, get_extension_from_path
#create a simple jigsaw run case class
pDate = datetime.datetime.today()
sDate_default = (
    "{:04d}".format(pDate.year)
    + "{:02d}".format(pDate.month)
    + "{:02d}".format(pDate.day)
)

class jigsawcase:
    """Class to represent a JIGSAW run case with configuration and file paths"""
    iFlag_save_mesh = 1
    sFilename_mpas_mesh_netcdf = None
    sFilename_jigsaw_mesh_netcdf = None
    sFilename_land_ocean_mask = None
    pBoundary_wkt = None
    aConfig_jigsaw = None
    def __init__(self, aConfig_in,
        iFlag_standalone_in=None,
        iFlag_create_directory_in=1,
        sModel_in=None,
        sDate_in=None,
        sWorkspace_output_in=None):
        # flags
        if iFlag_standalone_in is not None:
            self.iFlag_standalone = iFlag_standalone_in
        else:
            if "iFlag_standalone" in aConfig_in:
                self.iFlag_standalone = int(aConfig_in["iFlag_standalone"])
            else:
                self.iFlag_standalone = 1

        if "iFlag_global" in aConfig_in:
            self.iFlag_global = int(aConfig_in["iFlag_global"])
        else:
            self.iFlag_global = 1

        if "iCase_index" in aConfig_in:
            iCase_index = int(aConfig_in["iCase_index"])
        else:
            iCase_index = 1



        sCase_index = "{:03d}".format(iCase_index)
        self.iCase_index = iCase_index

        if sWorkspace_output_in is not None:
            self.sWorkspace_output = sWorkspace_output_in
        else:
            if "sWorkspace_output" in aConfig_in:
                self.sWorkspace_output = aConfig_in["sWorkspace_output"]

        if "sFilename_mpas_mesh_netcdf" in aConfig_in:
            self.sFilename_mpas_mesh_netcdf = aConfig_in["sFilename_mpas_mesh_netcdf"]

        if "sFilename_jigsaw_mesh_netcdf" in aConfig_in:
            self.sFilename_jigsaw_mesh_netcdf = aConfig_in["sFilename_jigsaw_mesh_netcdf"]

        if "sFilename_land_ocean_mask" in aConfig_in:
            self.sFilename_land_ocean_mask = aConfig_in["sFilename_land_ocean_mask"]

        if "pBoundary_wkt" in aConfig_in:
            self.pBoundary_wkt = aConfig_in["pBoundary_wkt"]

        sDate = aConfig_in["sDate"]
        if sDate is not None:
            self.sDate = sDate
        else:
            self.sDate = sDate_default

        self.sModel = 'jigsaw'

        sCase = self.sModel + self.sDate + sCase_index
        self.sCase = sCase

        # the model can be run as part of hexwatershed or standalone
        if self.iFlag_standalone == 1:
            # in standalone case, will add case information and update output path
            sPath = str(Path(self.sWorkspace_output) / sCase)
            self.sWorkspace_output = sPath
        else:
            # use specified output path, also do not add output or input tag
            sPath = self.sWorkspace_output

        if iFlag_create_directory_in == 1:
            Path(sPath).mkdir(parents=True, exist_ok=True)

        if "sFilename_model_configuration" in aConfig_in:
            self.sFilename_model_configuration = aConfig_in["sFilename_model_configuration"]
        else:
            self.sFilename_model_configuration = None

        self.sFilename_mesh = os.path.join(
            str(Path(self.sWorkspace_output)),  "mpas.geojson"
        )

        if "iFlag_run_jigsaw" in aConfig_in:
            self.iFlag_run_jigsaw = int(aConfig_in["iFlag_run_jigsaw"])
        else:
            self.iFlag_run_jigsaw = 1

        if self.iFlag_run_jigsaw == 1:
            self.aConfig_jigsaw = aConfig_in

        return

    def jigsaw_setup(self):
        """Prepare the JIGSAW workspace directories.

        Creates the required subdirectories (tmp, out) inside
        self.sWorkspace_output so that run_jigsaw can write its
        intermediate and output files.

        Also copies the land ocean mask to a GeoJSON file for record
        if provided.
        """
        import jigsawpy

        # Ensure the main output directory exists
        Path(self.sWorkspace_output).mkdir(parents=True, exist_ok=True)

        # Create standard JIGSAW subdirectories
        for subdir in ("tmp", "out"):
            Path(os.path.join(self.sWorkspace_output, subdir)).mkdir(
                parents=True, exist_ok=True
            )

        # Copy land ocean mask to GeoJSON for record if provided
        if self.sFilename_land_ocean_mask is not None:
            self._convert_land_ocean_mask_to_geojson()

        # Verify jigsawpy is available and its binary is accessible
        try:
            _opts = jigsawpy.jigsaw_jig_t()
            print(f"JIGSAW setup complete in: {self.sWorkspace_output}")
        except Exception as e:
            print(f"Warning: jigsawpy setup check failed: {e}")

        return

    def jigsaw_run(self):
        """Run the JIGSAW mesh generation.

        Calls run_jigsaw using the configuration stored in
        self.sFilename_model_configuration. The configuration JSON
        must contain all required run_jigsaw parameters.
        """
        from mpas_land_mesh.mesh.jigsaw.run_jigsaw import run_jigsaw
        import json

        if self.sFilename_model_configuration is None:
            raise ValueError("sFilename_model_configuration is not set")

        with open(self.sFilename_model_configuration, "r") as f:
            aConfig = json.load(f)

        projector = [0.0, 0.0]
        geom, gprj, mesh, mprj = run_jigsaw(
            self.sWorkspace_output,
            projector,
            aConfig_in=aConfig,
        )
        return geom, gprj, mesh, mprj

    def _convert_land_ocean_mask_to_geojson(self):
        """Convert the land ocean mask to GeoJSON format for record.

        This function follows the pattern from pyflowline's setup function,
        copying the land ocean mask file to a GeoJSON file in the output
        directory for record keeping.
        """
        if self.sFilename_land_ocean_mask is None:
            return

        sFilename_raw = self.sFilename_land_ocean_mask
        sFilename_out = os.path.join(
            str(Path(self.sWorkspace_output)), "land_ocean_mask.geojson"
        )

        # Check whether the file exists
        if not os.path.isfile(sFilename_raw):
            print(f"The land ocean mask file does not exist: {sFilename_raw}")
            return

        # Check the file type of the input file
        sExtension = get_extension_from_path(sFilename_raw)

        if sExtension == ".geojson" or sExtension == ".json":
            # If already GeoJSON, just copy it
            copy2(sFilename_raw, sFilename_out)
            print(f"Copied land ocean mask to: {sFilename_out}")
        else:
            # For other formats (shapefile, etc.), use GDAL to convert
            try:
                from osgeo import ogr
                # Open the input file
                pDataset_in = ogr.Open(sFilename_raw, 0)

                if pDataset_in is None:
                    print(f"Could not open land ocean mask file: {sFilename_raw}")
                    return

                # Create GeoJSON output
                pDriver_json = ogr.GetDriverByName("GeoJSON")
                if os.path.exists(sFilename_out):
                    pDriver_json.DeleteDataSource(sFilename_out)

                pDataset_out = pDriver_json.CreateDataSource(sFilename_out)
                pLayer_in = pDataset_in.GetLayer(0)

                # Copy layer to GeoJSON
                pDataset_out.CopyLayer(pLayer_in, "land_ocean_mask", ["OVERWRITE=YES"])

                # Clean up
                pDataset_in = None
                pDataset_out = None

                print(f"Converted land ocean mask to GeoJSON: {sFilename_out}")
            except Exception as e:
                print(f"Error converting land ocean mask to GeoJSON: {e}")

        return

    def _jigsaw_create_hpc_job(self, sSlurm_in=None, hours_in=10):
        """Create a HPC job for this JIGSAW simulation.

        Generates two files in self.sWorkspace_output:
          - run_jigsaw.py  : Python driver script that reads the configuration
                             and runs the JIGSAW mesh-generation workflow.
          - submit.job     : SLURM batch script that activates the conda
                             environment and executes run_jigsaw.py.

        Args:
            sSlurm_in (str, optional): SLURM partition name. Defaults to 'slurm'.
            hours_in  (int, optional): Wall-clock time limit in hours. Defaults to 10.
        """
        os.chdir(self.sWorkspace_output)
        sConda_env_path, sConda_env_name, env_type = get_python_environment()

        # ------------------------------------------------------------------
        # Part 1 – Python driver script (run_jigsaw.py)
        # ------------------------------------------------------------------
        sFilename_jigsaw = os.path.join(
            str(Path(self.sWorkspace_output)), "run_jigsaw.py"
        )
        ofs_jigsaw = open(sFilename_jigsaw, "w")

        sLine = (
            "#!/qfs/people/liao313/.conda/envs/"
            + sConda_env_name
            + "/bin/"
            + "python3"
            + "\n"
        )
        ofs_jigsaw.write(sLine)
        sLine = "import os" + "\n"
        ofs_jigsaw.write(sLine)
        sLine = (
            "os.environ['PROJ_LIB']=" + '"' + sConda_env_path + "/share/proj" + '"' + "\n"
        )
        ofs_jigsaw.write(sLine)
        sLine = (
            "os.environ['LD_LIBRARY_PATH']="
            + '"'
            + sConda_env_path
            + '/lib:${LD_LIBRARY_PATH}"'
            + "\n"
        )
        ofs_jigsaw.write(sLine)
        sLine = (
            "from mpas_land_mesh.utilities.config_manager import read_jigsaw_configuration_file"
            + "\n"
        )
        ofs_jigsaw.write(sLine)
        sLine = (
            "from mpas_land_mesh.mesh.mpas.create_mpas_mesh import create_mpas_mesh" + "\n"
        )
        ofs_jigsaw.write(sLine)
        sLine = (
            "sFilename_configuration_in = "
            + '"'
            + self.sFilename_model_configuration
            + '"\n'
        )
        ofs_jigsaw.write(sLine)
        sLine = (
            "oJigsaw = read_jigsaw_configuration_file(sFilename_configuration_in,"
            + "iCase_index_in="
            + str(self.iCase_index)
            + ","
            + 'sDate_in="'
            + str(self.sDate)
            + '"'
            + ")"
            + "\n"
        )
        ofs_jigsaw.write(sLine)
        sLine = "oJigsaw.jigsaw_setup()" + "\n"
        ofs_jigsaw.write(sLine)

        sLine = (
            "aMpas = create_mpas_mesh("
            + "oJigsaw.sFilename_mesh,"
            + " iFlag_global_in=oJigsaw.iFlag_global,"
            + " iFlag_save_mesh_in=oJigsaw.iFlag_save_mesh,"
            + " iFlag_run_jigsaw_in=1,"
            + " sWorkspace_jigsaw_in=oJigsaw.sWorkspace_output,"
            + " sFilename_mpas_mesh_netcdf_in=oJigsaw.sFilename_mpas_mesh_netcdf,"
            + " sFilename_jigsaw_mesh_netcdf_in=oJigsaw.sFilename_jigsaw_mesh_netcdf,"
            + " sFilename_land_ocean_mask_in=oJigsaw.sFilename_land_ocean_mask,"
            + " aConfig_jigsaw_in=oJigsaw.aConfig_jigsaw,"
            + ")"
            + "\n"
        )
        ofs_jigsaw.write(sLine)
        sLine = "print('Finished')" + "\n"
        ofs_jigsaw.write(sLine)
        ofs_jigsaw.close()
        os.chmod(sFilename_jigsaw, stat.S_IREAD | stat.S_IWRITE | stat.S_IXUSR)

        # ------------------------------------------------------------------
        # Part 2 – SLURM batch script (submit.job)
        # ------------------------------------------------------------------
        sFilename_job = os.path.join(str(Path(self.sWorkspace_output)), "submit.job")
        ofs = open(sFilename_job, "w")
        sLine = "#!/bin/bash\n"
        ofs.write(sLine)
        sLine = "#SBATCH -A E3SM\n"
        ofs.write(sLine)
        sLine = "#SBATCH --job-name=" + self.sCase + "\n"
        ofs.write(sLine)
        sHour = "{:02d}".format(hours_in)
        sLine = "#SBATCH -t " + sHour + ":00:00" + "\n"
        ofs.write(sLine)
        sLine = "#SBATCH --nodes=1" + "\n"
        ofs.write(sLine)
        sLine = "#SBATCH --ntasks-per-node=1" + "\n"
        ofs.write(sLine)
        if sSlurm_in is not None:
            sSlurm = sSlurm_in
        else:
            sSlurm = "slurm"
        sLine = "#SBATCH --partition=" + sSlurm + "\n"
        ofs.write(sLine)
        sLine = "#SBATCH -o stdout.out\n"
        ofs.write(sLine)
        sLine = "#SBATCH -e stderr.err\n"
        ofs.write(sLine)
        sLine = "module purge\n"
        ofs.write(sLine)
        sLine = "module load gcc/8.1.0" + "\n"
        ofs.write(sLine)
        sLine = "module load python/miniconda2024May29 " + "\n"
        ofs.write(sLine)
        sLine = "source /share/apps/python/miniconda2024May29/etc/profile.d/conda.sh" + "\n"
        ofs.write(sLine)
        sLine = "conda activate " + sConda_env_name + "\n"
        ofs.write(sLine)
        sLine = "cd $SLURM_SUBMIT_DIR\n"
        ofs.write(sLine)
        sLine = "JOB_DIRECTORY=" + self.sWorkspace_output + "\n"
        ofs.write(sLine)
        sLine = "cd $JOB_DIRECTORY" + "\n"
        ofs.write(sLine)
        sLine = "python3 run_jigsaw.py" + "\n"
        ofs.write(sLine)
        sLine = "conda deactivate" + "\n"
        ofs.write(sLine)
        ofs.close()
        return

