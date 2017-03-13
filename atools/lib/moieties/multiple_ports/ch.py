import mbuild as mb


class CH(mb.Compound):
    """A carbon with a hydrogen and three open ports. """
    def __init__(self):
        super(CH, self).__init__()

        mb.load('ch.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)  # Move carbon to origin.

        self.add(mb.Port(anchor=self[0],
                         orientation=[0, 1, 0],
                         separation=0.075), 'up')

        self.add(mb.Port(anchor=self[0],
                         orientation=[0, -1, 0],
                         separation=0.075), 'down')

        self.add(mb.Port(anchor=self[0],
                         orientation=[-1, 0, 0],
                         separation=0.075), 'left')

if __name__ == '__main__':
    ch = CH()
