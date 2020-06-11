import mbuild as mb


class Ch3(mb.Compound):
    """A methyl group. """
    def __init__(self):
        super(Ch3, self).__init__()

        mb.load('ch3.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)  # Move carbon to origin.

        self.add(mb.Port(anchor=self[0]), 'down')
        self['down'].translate([0, -0.07, 0])

if __name__ == '__main__':
    m = Ch3()
    m.visualize(show_ports=True)

