import mbuild as mb

from atools.lib.moieties.multiple_ports.united_atom.ch2 import CH2
from atools.lib.moieties.one_port.united_atom.ch3 import CH3


class Alkane(mb.Compound):
    """ A united-atom alkane chain, optionally terminated by CH3 particles"""
    def __init__(self, n, cap_front=True, cap_end=True):
        """Initialize an united atom Alkane Compound.

        Parameters
        ----------
        n: int
            Number of backbone particles.
        cap_front: bool, optional, default=True
            Add methyl group to beginning of chain ('down' port).
        cap_end: bool, optional, default=True
            Add methyl group to end of chain ('up' port).

        """

        if n < 2:
            raise Exception('n must be 1 or more')
        super(Alkane, self).__init__()

        # Adjust length of Polmyer for absence of methyl terminations.
        if not cap_front:
            n += 1
        if not cap_end:
            n += 1
        chain = mb.Polymer(CH2(), n=n-2, port_labels=('up', 'down'))
        self.add(chain, 'chain')

        if cap_front:
            self.add(CH3(), "methyl_front")
            mb.force_overlap(self['chain'], self['chain']['up'], 
                             self['methyl_front']['up'])
        else:
            # Hoist port label to UA_alkane level.
            self.add(chain['up'], 'up', containment=False)

        if cap_end:
            self.add(CH3(), 'methyl_end')
            mb.force_overlap(self['methyl_end'], self['methyl_end']['up'], 
                             self['chain']['down'])
        else:
            # Hoist port label to UA_alkane level.
            self.add(chain['down'], 'down', containment=False)

if __name__ == "__main__":
    n = 8
    chain = Alkane(n=n, cap_front=False, cap_end=True)
    chain.save('octane.mol2', overwrite=True)
