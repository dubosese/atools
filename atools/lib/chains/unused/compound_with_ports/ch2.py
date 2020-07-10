import mbuild as mb

class Ch2(mb.Compound):
    def __init__(self):
        super(Ch2, self).__init__()
        
        mb.load('ch2.pdb', compound=self, relative_to_module=self.__module__)
        carbon = list(self.particles_by_name('C'))[0]
        self.add(mb.Port(anchor=carbon, orientation=[0, 0, 1], separation=.075), label='up')
        self.add(mb.Port(anchor=carbon, orientation=[0, 0, -1], separation=0.075), label='down')

if __name__ == '__main__':
    ch2 = Ch2()
    