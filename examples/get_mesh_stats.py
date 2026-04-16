from mpas_land_mesh.utilities.mesh_stats import get_mesh_stats

sFilename_mpas_mesh_in = '/compyfs/liao313/04model/pyhexwatershed/global/pyflowline20260302004/jigsaw/out/invert_mesh.nc'

stats = get_mesh_stats(sFilename_mpas_mesh_in)

print('Number of cells:      ', stats['nCells'])
print('Min cell area:        ', stats['minCellArea'])
print('Max cell area:        ', stats['maxCellArea'])
print('Average cell area:    ', stats['avgCellArea'])
print('Min edge length:      ', stats['minEdgeLength'])
print('Max edge length:      ', stats['maxEdgeLength'])
print('Average edge length:  ', stats['avgEdgeLength'])