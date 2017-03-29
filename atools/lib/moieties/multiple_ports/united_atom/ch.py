import mbuild as mb


class CH(mb.Compound):
    def __init__(self):
        super(CH, self).__init__()

        self.add(mb.Particle(name='_CH'))

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
