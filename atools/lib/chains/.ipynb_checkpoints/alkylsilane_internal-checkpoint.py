import mbuild as mb

from mbuild.lib.atoms import H
from mbuild.lib.moieties import Silane
import atools
from atools.lib.chains.alkane import Alkane
import math


class Alkylsilane(mb.Compound):
    """A functionalized alkylsilane chain. """
    def __init__(self, chain_length, internal_group, locations, terminal_group):
        super(Alkylsilane, self).__init__()

        module = __import__('atools.lib.moieties.one_port.'+terminal_group)
        class_ = getattr(module.lib.moieties.one_port, terminal_group.title())
        tgroup = class_()
 
        hmodule = __import__('atools.lib.moieties.multiple_ports.'+internal_group)
        hclass_ = getattr(hmodule.lib.moieties.multiple_ports, internal_group.title())
        hgroup = hclass_()

        # Determine alkane segments
        if isinstance(locations, int):
            locations = [locations]
        locations.sort()

        silane = Silane()
        self.add(silane, 'silane')
        self.add(silane['down'], 'down', containment=False)

        if 0 in locations:
            '''
            self.add(hgroup, 'hgroup')
            mb.force_overlap(self['silane'], self['silane']['up'], 
                             self['hgroup']['down'])
            '''
            current_segment = silane
        else:
            first_length = locations[0]
            first_segment = Alkane(first_length, cap_front=False, cap_end=False)
            self.add(first_segment, 'bottom_chain')
            mb.force_overlap(self['silane'], self['silane']['up'], 
                             self['bottom_chain']['down'])
            current_segment = first_segment

        c_remove = 0
        if internal_group in ['amide', 'hemiacetal']:
            c_remove += 1

        for i, loc in enumerate(locations[1:]):
            hgroup_clone = mb.clone(hgroup)
            self.add(hgroup_clone, 'hgroup{}'.format(i+1))
            mb.force_overlap(self['hgroup{}'.format(i+1)], 
                             self['hgroup{}'.format(i+1)]['down'],
                             current_segment['up'])
            current_segment = hgroup_clone
            length = loc - locations[i] - 1 - c_remove
            if length > 0:
                segment = Alkane(length, cap_front=False, cap_end=False)
                self.add(segment, 'internal_chain{}'.format(i+1))
                current_segment = segment
                mb.force_overlap(self['internal_chain{}'.format(i+1)],
                                 self['internal_chain{}'.format(i+1)]['down'],
                                 self['hgroup{}'.format(i+1)]['up'])

        self.add(hgroup, 'hgroup')
        mb.force_overlap(self['hgroup'], 
                         self['hgroup']['down'],
                         current_segment['up'])

        last_length = chain_length - locations[-1] - 1 - c_remove
        if last_length:
            last_segment = Alkane(last_length, cap_front=False, cap_end=False)
            self.add(last_segment, 'top_chain')
            mb.force_overlap(self['top_chain'], self['top_chain']['down'],
                             self['hgroup']['up'])
        else:
            hydrogen = H()
            self.add(hydrogen, 'H-cap')
            mb.force_overlap(self['H-cap'], self['H-cap']['up'],
                             self['hgroup']['up'])


        self.add(tgroup, 'terminal_group')
        mb.force_overlap(self['terminal_group'], 
                         self['terminal_group']['down'],
                         self['top_chain']['up'])
                    
if __name__ == "__main__":
    chain = Alkylsilane(chain_length=18, internal_group='phenyl',
                        locations=17, terminal_group='phenyl')
    chain.save('chain-hbond.mol2', overwrite=True)

