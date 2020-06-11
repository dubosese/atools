import mbuild as mb


class Oxygen(mb.Compound):
    def __init__(self):
        super(Oxygen, self).__init__()
                
        self.add(mb.Particle(name='O'))
        up_port = mb.Port(anchor=self[0], orientation=[0, 0, 1], separation=.075)
        down_port = mb.Port(anchor=self[0], orientation=[0, 0, -1], separation=.075)
        self.add(up_port, label='up')
        self.add(down_port, label='down')


if __name__ == '__main__':
    oxygen = Oxygen()
