from __future__ import division

import os
import mbuild as mb
import numpy as np


class SilicaSurface(mb.Compound):
    """ """
    def __init__(self, surface_id=1):
        super(SilicaSurface, self).__init__()

        mb.load('silica_surface-1x1_{}.pdb'.format(surface_id),
                compound=self,
                relative_to_module=self.__module__)

        self.periodicity = np.array([5.0, 5.0, 0.0])

        for atom in self.particles():
            if atom.name == 'OS':
                port = mb.Port(anchor=atom)
                mb.rotate_around_x(port, np.pi/2)
                mb.translate(port, atom.pos + np.array([0.0, 0.0, 0.1]))
                self.add(port, "port_{}".format(len(self.referenced_ports())))

if __name__ == "__main__":
    surface = SilicaSurface()
