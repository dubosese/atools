import numpy as np

import mbuild as mb
from mbuild.examples.alkane.alkane import Alkane
from mbuild.lib.moieties import Silane

from atools.lib.moieties.multiple_ports.aminopropyldopamine import APD

class nC_APD(mb.Compound):
    """Aminopropyl dopamine functionalized with an alkane chain. """
    def __init__(self, chain_length):
        super(nC_APD, self).__init__()

        alkane = Alkane(chain_length, cap_end=False)
        self.add(alkane, 'alkane')
        apd = APD()
        self.add(apd, 'apd')
        mb.force_overlap(self['alkane'], self['alkane']['down'], 
                         self['apd']['up'])

        silane = Silane()
        self.add(silane, 'silane')
        mb.force_overlap(self['silane'], self['silane']['up'], 
                         self['apd']['down'])

        # Hoist silane port to AlkylSilane level.
        self.add(silane['down'], 'down', containment=False)
        self.spin(np.pi/2, [1, 0, 0])

if __name__ == "__main__":
    octadecyl_apd = nC_APD(18)
    octadecyl_apd.save('c18-APD.mol2', overwrite=True)
