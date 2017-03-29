import mbuild as mb


class Carbonyl(mb.Compound):
    """A carbonyl group with two ports 120 degrees apart. """
    def __init__(self):
        super(Carbonyl, self).__init__()

        mb.load('carbonyl.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)  # Move carbon to origin.

        self.add(mb.Port(anchor=self[0],
                         orientation=[1, 0, 0],
                         separation=0.075), 'right')

        self.add(mb.Port(anchor=self[0],
                         orientation=[-1, 0, 0],
                         separation=0.075), 'left')

if __name__ == '__main__':
    carbonyl = Carbonyl()
