import numpy as np

import mbuild as mb

from atools.lib.moieties.multiple_ports.carbonyl import Carbonyl
from atools.lib.moieties.multiple_ports.nh import NH


class Amide(mb.Compound):
    """A useful description """
    def __init__(self):
        super(Amide, self).__init__()

        carbonyl = Carbonyl()
        nh = NH()
        self.add(carbonyl, 'carbonyl')
        self.add(nh, 'nh')
        mb.force_overlap(carbonyl, self['carbonyl']['right'], self['nh']['up'])
        self['nh'][1].translate(2 * (self['nh'][0].pos - self['nh'][1].pos))

        self.add(carbonyl['left'], 'up', containment=False)
        self.add(nh['down'], 'down', containment=False)

if __name__ == '__main__':
    amide = Amide()
    amide.save('test.mol2', overwrite=True)
