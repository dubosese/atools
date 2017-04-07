from warnings import warn

import mbuild as mb
from mbuild.lib.moieties import CH2
from atools.lib.moieties.one_port import CH3


class Alkane(mb.Compound):
    """An alkane which may optionally end with a hydrogen or a Port."""
    def __init__(self, n=3, cap_front=True, cap_end=True):
        """Initialize an Alkane Compound.

        Args:
            n: Number of carbon atoms.
            cap_front: Add methyl group to beginning of chain ('down' port).
            cap_end: Add methyl group to end of chain ('up' port).
        """
        if n < 1:
            raise ValueError('n must be 1 or more')
        elif n == 1:
            if cap_front or cap_end:
                warn('Returning CH3 group')
            else:
                warn('Returning CH2 group')
        super(Alkane, self).__init__()

        if n > 1:
            if cap_front:
                n -= 1
                ch3_front = CH3()
                self.add(ch3_front, 'methyl_front')
            if cap_end:
                n -= 1
                ch3_end = clone(ch3_front)
                self.add(ch3_end, 'methyl_end')
            try:
                chain = mb.Polymer(CH2(), n=n, port_labels=('up', 'down'))
                self.add(chain, 'chain')
                if cap_end:
                    mb.force_overlap(self['methyl_end'], self['methyl_end']['down'],
                                     self['chain']['down'])
                else:
                    self.add(chain['down'], 'down', containment=False)
                if cap_front:
                    mb.force_overlap(self['chain'], self['chain']['up'],
                                     self['methyl_front']['down'])
                else:
                    self.add(chain['up'], 'up', containment=False)
            except ValueError:
                mb.force_overlap(self['methyl_end'], self['methyl_end']['down'],
                                 self['methyl_front']['down'])
        else:
            if cap_end or cap_front:
                ch3 = CH3()
                self.add(ch3, 'methyl')
                self.add(ch3['down'], 'down', containment=False)
            else:
                ch2 = CH2()
                self.add(ch2, 'methylene')
                self.add(ch2['up'], 'up', containment=False)
                self.add(ch2['down'], 'down', containment=False)
