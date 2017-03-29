import mbuild as mb


class Amine(mb.Compound):
    """A secondary amine group. """
    def __init__(self):
        super(Amine, self).__init__()

        mb.load('nh.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)

        self.add(mb.Port(anchor=self[0],
                         orientation=[0, 1, 0],
                         separation=0.075), 'up')

        self.add(mb.Port(anchor=self[0],
                         orientation=[0, -1, 0],
                         separation=0.075), 'down')

if __name__ == '__main__':
    amine = Amine()
