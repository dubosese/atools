import mbuild as mb

class Phenyl(mb.Compound):
    """ A phenyl group. """
    def __init__(self):
        super(Phenyl, self).__init__()

        mb.load('phenyl.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)

        self.add(mb.Port(anchor=self[0], orientation=[0, -1, 0], separation=0.07),
            'down')

if __name__ == '__main__':
    phenyl = Phenyl()
    phenyl.save('phenyl-test.mol2', overwrite=True)
