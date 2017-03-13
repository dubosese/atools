import mbuild as mb
import numpy as np


class Aldehyde(mb.Compound):
    """ A aldehyde group. """
    def __init__(self):
        super(Aldehyde, self).__init__()

        mb.load('aldehyde.pdb', compound=self, relative_to_module=self.__module__)
        #mb.spin_x(self, 90.0 * (np.pi / 180.0))
        mb.translate(self, -self[0].pos)

        self.add(mb.Port(anchor=self[0]), 'down')
        mb.translate(self['down'], [0, -0.07, 0])

if __name__ == '__main__':
    aldehyde = Aldehyde()
    aldehyde.save('aldehyde-test.mol2', overwrite=True)
