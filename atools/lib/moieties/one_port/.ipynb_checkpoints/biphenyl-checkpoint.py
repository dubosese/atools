import mbuild as mb
        
class Phenyl(mb.Compound):
    """ A phenyl group. """
    def __init__(self):
        super(Phenyl, self).__init__()

        mb.load('phenyl.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)

        self.add(mb.Port(anchor=self[0], orientation=[0, -1, 0], separation=0.07),
            'down')    

class Biphenyl(mb.Compound):
    """ A biphenyl group. """
    def __init__(self):
        super(Biphenyl, self).__init__()

        phenyl1 = Phenyl()
        phenyl2 = Phenyl()

        phenyl1.remove(phenyl1[6])

        self.add(phenyl1, 'phenyl1')
        self.add(phenyl2, 'phenyl2')
        mb.force_overlap(self['phenyl2'], self['phenyl2']['down'], 
                         self['phenyl1'].all_ports()[0])

        self.add(phenyl1['down'], 'down', containment=False)

if __name__ == '__main__':
    biphenyl = Biphenyl()
    biphenyl.save('biphenyl-test.mol2', overwrite=True)
