import mbuild as mb


class Hemiacetal(mb.Compound):
    """A nice description. """
    def __init__(self):
        super(Hemiacetal, self).__init__()

        mb.load('hemiacetal.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)

        self.add(mb.Port(anchor=self[0],
                         orientation=[0, 1, 0],
                         separation=0.075), 'up')

        self.add(mb.Port(anchor=self[1],
                         orientation=[0, -1, 0],
                         separation=0.075), 'down')

if __name__ == '__main__':
    hemiacetal = Hemiacetal()
