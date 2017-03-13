import mbuild as mb
import numpy as np


class Methyl(mb.Compound):
    """ A methyl group. """
    def __init__(self):
        super(Methyl, self).__init__()

        mb.load('methyl.pdb', compound=self, relative_to_module=self.__module__)
        #mb.spin_x(self, 90.0 * (np.pi / 180.0))
        mb.translate(self, -self[0].pos)

        self.add(mb.Port(anchor=self[0]), 'down')
        mb.translate(self['down'], [0, -0.07, 0])

if __name__ == '__main__':
    methyl = Methyl()
    methyl.save('methyl-test.mol2', overwrite=True)
