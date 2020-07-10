import mbuild as mb
import atools
import math
from mbuild.lib.moieties import Silane
        
name = 'o'
module = __import__('compound_with_ports.'+name)
O = getattr(module, name.title())
        
name = 'ch2'
module = __import__('compound_with_ports.'+name)
CH2 = getattr(module, name.title())

name = 'h'
module = __import__('compound_with_ports.'+name)
H = getattr(module, name.title())
        
    
class Pegsilane(mb.Compound):
    def __init__(self, length=1, terminal_group=None):
        super(Pegsilane, self).__init__()
        
        module = __import__('atools.lib.moieties.one_port.'+terminal_group)
        class_ = getattr(module.lib.moieties.one_port, terminal_group.title())
        t_group = class_()
        
        self.add(t_group, 'terminal_group')
        
        chain = mb.recipes.Polymer(monomers=[CH2(), O()], sequence='AAB',
                                n=math.floor(length/3), port_labels=('up', 'down'))
        
        mb.force_overlap(self['terminal_group'],
                        self['terminal_group']['down'],
                        to_positions=chain['up'])
        
        self.add(chain)
        
        silane = Silane()
        self.add(silane)
        mb.force_overlap(move_this=silane,
                        from_positions=silane.all_ports()[1],
                        to_positions=chain['down'])
        
        self.add(silane['down'], 'down', containment=False)
        
if __name__ == "__main__":
    pegsilane = Pegsilane(length=12, terminal_group='methyl')
