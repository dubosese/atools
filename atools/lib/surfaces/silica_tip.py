from __future__ import division

import os
import mbuild as mb
import numpy as np


class SilicaTip(mb.Compound):
    """ """
    def __init__(self, tip_radius=2):
        super(SilicaTip, self).__init__()

        mb.load('silica_tip-{}nm.pdb'.format(int(tip_radius)),
                compound=self,
                relative_to_module=self.__module__)
        self.periodicity = np.array([0.0, 0.0, 0.0])
        mb.spin_x(self,np.pi)

        for atom in self.particles():
            if atom.name == 'OS':
                port = mb.Port(anchor=atom)
                mb.rotate_around_x(port, np.pi/2)
                mb.translate(port, atom.pos + np.array([0.0, 0.0, 0.1]))
                self.add(port, "port_{}".format(len(self.referenced_ports())))

if __name__ == "__main__":
    tip = SilicaTip(2) 
