from __future__ import division

import numpy as np

import mbuild as mb
from mbuild.examples.alkane.alkane import Alkane

from atools.lib.moieties.multiple_ports.carbonyl import Carbonyl
from atools.lib.moieties.multiple_ports.ch import CH
from atools.lib.moieties.multiple_ports.nh import NH
from atools.lib.moieties.one_port.dihydroxyphenyl import DihydroxyPhenyl

class APD(mb.Compound):
    """ Aminopropyl Dopamine terminated with an amide """
    def __init__(self):
        super(APD, self).__init__()

        carbonyl = Carbonyl()
        ch = CH()
        ch['left'].spin(np.pi/2, [1, 0, 0])
        dihydroxyphenyl = DihydroxyPhenyl()
        nh = NH()
        propyl = Alkane(3, cap_end=False, cap_front=False)

        amide_nh = mb.clone(nh)
        amide = mb.Compound()
        amide.add(amide_nh, 'nh')
        amide.add(carbonyl, 'carbonyl')
        mb.force_overlap(amide['nh'], amide['nh']['up'], amide['carbonyl']['left'])

        amide.add(amide_nh['down'], 'down', containment=False)
        amide.add(carbonyl['right'], 'up', containment=False)

        self.add(propyl, 'propyl')
        self.add(nh, 'nh')
        mb.force_overlap(self['nh'], self['nh']['down'], self['propyl']['up'])

        self.add(ch, 'ch')
        mb.force_overlap(self['ch'], self['ch']['down'], self['nh']['up'])

        self.add(dihydroxyphenyl, 'dhp')
        mb.force_overlap(self['dhp'], self['dhp']['down'], self['ch']['left'])
    
        self.add(amide, 'amide')
        mb.force_overlap(self['amide'], self['amide']['down'], self['ch']['up'])

        self.add(amide['up'], 'up', containment=False)
        self.add(propyl['down'], 'down', containment=False)
        self['up'].spin(-60.0 * (np.pi/180), [0, 0, 1])

if __name__ == '__main__':
    molecule = APD()
