import mbuild as mb
from mbuild.lib.bulk_materials import AmorphousSilica

tile_x = 1
tile_y = 1
surface_id = 2
silica_interface = mb.SilicaInterface(bulk_silica=AmorphousSilica(), tile_x=tile_x, tile_y=tile_y, thickness=1.2)
traj = silica_interface.to_trajectory()
traj.save('silica_surface-{}x{}_{}.pdb'.format(tile_x,tile_y,surface_id))
