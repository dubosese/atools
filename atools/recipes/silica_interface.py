import inspect
import os

import mbuild as mb

from atools.lib import surface_cache

class SilicaInterface(mb.SilicaInterface):

    def __init__(self, tile_x=1, tile_y=1, thickness=1.0, seed=12345):
        thickness = float(thickness)

        cache_dir = os.path.dirname(inspect.getfile(surface_cache))
        filename ='silica-interface-tx{}-ty{}-thickness{}-seed{}.mol2'.format(tile_x, tile_y, thickness, seed)
        if any(file == filename for file in os.listdir(cache_dir)):
            mb.Compound.__init__(self)

            mb.load(os.path.join(cache_dir, filename), compound=self)
            self._add_ports()
            self.periodicity = [5. * tile_x, 5. * tile_y, 0]

        else:
            from mbuild.lib.bulk_materials import AmorphousSilica

            super(SilicaInterface, self).__init__(bulk_silica=AmorphousSilica(),
                tile_x=tile_x, tile_y=tile_y, thickness=thickness, seed=seed)
            self.save(os.path.join(cache_dir, filename))

    def _add_ports(self):
        """ Add ports above surface sites """
        for atom in self.particles_by_name('OS'):
            port = mb.Port(anchor=atom, orientation=[0, 0, 1],
                           separation=0.1)
            self.add(port, "port_{}".format(len(self.referenced_ports())))

if __name__ == "__main__":
    silica_interface = SilicaInterface()
