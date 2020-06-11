import mbuild as mb

class Cf2(mb.Compound):
    def __init__(self):
        super(Cf2, self).__init__()
        
        mb.load('cf2.pdb', compound=self, relative_to_module=self.__module__)
        
        #cf2 = mb.load('./cf2.pdb', compound=self)
        carbon = list(self.particles_by_name('C'))[0]
        up_port = mb.Port(anchor=carbon, orientation=[0, 0, 1],
                         separation=.075)
        down_port = mb.Port(anchor=carbon, orientation=[0, 0, -1],
                           separation=.075)
        self.add(up_port, label='up')
        self.add(down_port, label='down')


if __name__ == '__main__':
    cf2 = Cf2()
