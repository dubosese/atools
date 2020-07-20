import mbuild as mb

from atools.lib.moieties.one_port.nitro import Nitro
from atools.lib.moieties.one_port.phenyl import Phenyl


class Nitrophenyl(mb.Compound):
    """ A nitrophenyl group. """
    def __init__(self):
        super(Nitrophenyl, self).__init__()

        phenyl = Phenyl()
        nitro = Nitro()

        phenyl.remove(phenyl[6])
        
        self.add(phenyl, 'phenyl')
        self.add(nitro, 'nitro')
        mb.force_overlap(self['nitro'], self['nitro']['down'], self['phenyl'].all_ports()[0])
        
        self.add(phenyl['down'], 'down', containment=False)

if __name__ == '__main__':
    nitrophenyl = Nitrophenyl()
    nitrophenyl.save('nitrophenyl-test.mol2', overwrite=True)
