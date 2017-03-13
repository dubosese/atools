import mbuild as mb
import numpy as np


class Ethylene(mb.Compound):
    """ A ethylene group. """
    def __init__(self):
        super(Ethylene, self).__init__()

        mb.load('ethylene.pdb', compound=self, relative_to_module=self.__module__)
        #mb.spin_x(self, 90.0 * (np.pi / 180.0))
        mb.translate(self, -self[0].pos)

        self.add(mb.Port(anchor=self[0]), 'down')
        mb.translate(self['down'], [0, -0.07, 0])

if __name__ == '__main__':
    ethylene = Ethylene()
    ethylene.save('ethylene-test.mol2', overwrite=True)
