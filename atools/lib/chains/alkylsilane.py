import mbuild as mb
from mbuild.examples.alkane.alkane import Alkane
from mbuild.lib.moieties import Silane

import atools


class Alkylsilane(mb.Compound):
    """A terminal-functionalized alkylsilane chain.

    An alkylsilane chain featuring a user-specified functional group at one
    terminus and a silane group (featuring an open port to attach to a surface)
    at the other terminus.

    Parameters
    ----------
    chain_length : int
        Length of the chain (number of carbons)
    terminal_group : str
        Functional group to attach to the chain terminus. Valid options are:
        'acetyl', 'amino', 'carboxyl', 'cyano', 'cyclopropyl', 'ethylene',
        'fluorophenyl', 'hydroxyl', 'isopropyl', 'methoxy', 'methyl', 'nitro',
        'nitrophenyl', 'perfluoromethyl', 'phenyl', or 'pyrrole'
    """
    def __init__(self, chain_length, terminal_group):
        super(Alkylsilane, self).__init__()

        module = __import__('atools.lib.moieties.one_port.'+terminal_group)
        class_ = getattr(module.lib.moieties.one_port, terminal_group.title())
        tgroup = class_()

        alkane = Alkane(chain_length, cap_front=False, cap_end=False)
        self.add(alkane, 'alkane')
        self.add(tgroup, 'terminal_group')
        mb.force_overlap(self['alkane'], self['alkane']['up'], 
                         self['terminal_group']['down'])
        silane = Silane()
        self.add(silane, 'silane')
        mb.force_overlap(self['silane'], self['silane']['up'], self['alkane']['down'])

        self.add(silane['down'], 'down', containment=False)

if __name__ == "__main__":
    chain = Alkylsilane(chain_length=8, terminal_group='cyclopropyl')
    chain.save('chain.mol2', overwrite=True)
