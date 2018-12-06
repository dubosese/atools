import inspect
import os

import mbuild as mb

from atools.lib import surface_cache
from atools.recipes.silica_interface_carve import SilicaInterfaceCarve


class SilicaInterface(SilicaInterfaceCarve):
    """Recipe for generating an amorphous silica interface.

    This is essentially a wrapper around the `SilicaInterfaceCarve` recipe, which is an extension of that included with
    the mBuild package (https://github.com/mosdef-hub/mbuild). Please refer to the
    mBuild documentation for further details into how the interface is carved.

    A surface cache located in `atools/lib` is first checked to see if a silica
    interface has already been generated with the specified `tile_x`, `tile_y`,
    `thickness`, and `seed` parameters. If so, the surface is read from a PDB file,
    otherwise, the `SilicaInterfaceCarve` recipe is called and the newly generated
    surface is added to the surface cache.

    Parameters
    ----------
    tile_x : int, optional, default=1
        Number of times to replicate the surface in the x dimension. The default
        length in the x dimension is 5nm, so increasing tile_x will lead to higher
        multiples of this number.
    tile_y : int, optional, default=1
        Number of times to replicate the surface in the y dimension. The default
        length in the y dimension is 5nm, so increasing tile_y will lead to higher
        multiples of this number.
    thickness : float, optional, default=1.0
        Desired thickness of the surface (in nm; not including oxygen layers on the
        top and bottom of the surface)
    seed : int, optional, default=12345
        Seed for the random number generator used in bridging surface oxygens
    add_to_cache : bool, optional, True
        Whether or not the surface should be added to the surface cache in the
        `atools` library
    read_from_cache : bool, optional, True
        Whether or not the surface should be loaded from the `atools` library if
        a surface matching the specified parameters is found.
    """
    def __init__(self, tile_x=1, tile_y=1, thickness=1.0, seed=12345,
                 add_to_cache=True, read_from_cache=True):
        thickness = float(thickness)

        cache_dir = os.path.dirname(inspect.getfile(surface_cache))
        filename ='silica-interface-tx{}-ty{}-thickness{}-seed{}.mol2'.format(tile_x,
            tile_y, thickness, seed)
        if (read_from_cache and 
                any(file == filename for file in os.listdir(cache_dir))):
            mb.Compound.__init__(self)

            mb.load(os.path.join(cache_dir, filename), compound=self)
            self._add_ports()
            self.periodicity = [5. * tile_x, 5. * tile_y, 0]

        else:
            from mbuild.lib.bulk_materials import AmorphousSilica

            super(SilicaInterface, self).__init__(bulk_silica=AmorphousSilica(),
                tile_x=tile_x, tile_y=tile_y, thickness=thickness, seed=seed)
            if add_to_cache:
                self.save(os.path.join(cache_dir, filename))

    def _add_ports(self):
        """Add ports above surface sites """
        for atom in self.particles_by_name('OS'):
            port = mb.Port(anchor=atom, orientation=[0, 0, 1],
                           separation=0.1)
            self.add(port, "port_{}".format(len(self.referenced_ports())))

if __name__ == "__main__":
    silica_interface = SilicaInterface(thickness=1.2, seed=10)
