import mbuild as mb
import numpy as np


class Nitro(mb.Compound):
    """ A nitro group. """
    def __init__(self):
        super(Nitro, self).__init__()

        mb.load('nitro.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)

        self.add(mb.Port(anchor=self[0], orientation=[0, -1, 0], separation=0.07),
            'down')

if __name__ == '__main__':
    nitro = Nitro()
    nitro.save('nitro-test.mol2', overwrite=True)
