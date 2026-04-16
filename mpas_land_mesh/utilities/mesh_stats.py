import os
import numpy as np
import netCDF4 as nc
def get_mesh_stats(sFilename_mpas_mesh):
    """
    Get statistics about the mesh.

    Parameters
    ----------
    sFilename_mpas_mesh : str
        The path to the mesh file.

    Returns
    -------
    dict
        A dictionary containing the mesh statistics.
    """

    #check if the file exists
    if not os.path.isfile(sFilename_mpas_mesh):
        raise FileNotFoundError(f'Mesh file {sFilename_mpas_mesh} not found.')

    #open the mesh file
    pDataset = nc.Dataset(sFilename_mpas_mesh)
    #read the mesh varaibles and compute the statistics
    aVertex_longititude = pDataset.variables['lonVertex'][:] * (180.0 / np.pi)
    #fix range of longitudes if necessary from [0, 360] to [-180, 180]
    aVertex_longititude = np.where(aVertex_longititude > 180, aVertex_longititude - 360, aVertex_longititude)
    aVertex_latitude = pDataset.variables['latVertex'][:] * (180.0 / np.pi)
    aConnectivity = pDataset.variables['verticesOnCell'][:]
    aCellID = pDataset.variables['indexToCellID'][:]

    aCellArea = pDataset.variables['areaCell'][:]
    dArea_min  = aCellArea.min().item()
    dArea_max  = aCellArea.max().item()
    dArea_mean = aCellArea.mean().item()
    nCell = aCellArea.size

    dLength_min  = np.sqrt(dArea_min)
    dLength_max  = np.sqrt(dArea_max)
    dLength_mean = np.sqrt(dArea_mean)

    stats = {}
    stats['nCells'] = nCell

    stats['minCellArea']   = dArea_min
    stats['maxCellArea']   = dArea_max
    stats['avgCellArea']   = dArea_mean

    stats['minEdgeLength'] = dLength_min
    stats['maxEdgeLength'] = dLength_max
    stats['avgEdgeLength'] = dLength_mean

    return stats