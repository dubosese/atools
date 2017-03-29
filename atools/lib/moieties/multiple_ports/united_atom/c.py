import mbuild as mb


class C(mb.Compound):
    def __init__(self):
        super(C, self).__init__()

        self.add(mb.Particle(name='C'))

        self.add(mb.Port(anchor=self[0],
                         orientation=[0, 1, 0],
                         separation=0.075), 'up')

        self.add(mb.Port(anchor=self[0],
                         orientation=[0, -1, 0],
                         separation=0.075), 'down')

        self.add(mb.Port(anchor=self[0],
                         orientation=[-1, 0, 0],
                         separation=0.075), 'left')

        self.add(mb.Port(anchor=self[0],
                         orientation=[1, 0, 0],
                         separation=0.075), 'right')

if __name__ == '__main__':
    carbon = C()
