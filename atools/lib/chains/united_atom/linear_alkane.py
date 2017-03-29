import mbuild as mb

from atools.lib.moieties.multiple_ports.united_atom.ch2 import CH2
from atools.lib.moieties.one_port.united_atom.ch3 import CH3
from atools.lib.moieties.no_ports.united_atom.ch4 import CH4


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
        super(Alkane, self).__init__()

        if n == 1:
            if cap_front and cap_end:
                ch4 = CH4()
                self.add(ch4)
            else:
                ch3 = CH3()
                self.add(ch3, 'ch3')
                self.add(self['ch3']['up'], 'up', containment=False)
        elif n == 2:
            if cap_front:
                ua1 = CH3()
            else:
                ua1 = CH2()
            if cap_end:
                ua2 = CH3()
            else:
                ua2 = CH2()
            self.add(ua1, 'ua1')
            self.add(ua2, 'ua2')
            mb.force_overlap(self['ua2'], self['ua2']['up'], self['ua1']['up'])
        else:
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
                self.add(chain['up'], 'up', containment=False)

            if cap_end:
                self.add(CH3(), 'methyl_end')
                mb.force_overlap(self['methyl_end'], self['methyl_end']['up'], 
                                 self['chain']['down'])
            else:
                self.add(chain['down'], 'down', containment=False)

if __name__ == "__main__":
    n = 3
    chain = Alkane(n=n, cap_front=True, cap_end=True)
    chain.save('chain.mol2', overwrite=True)
