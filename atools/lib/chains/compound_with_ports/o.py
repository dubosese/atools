import mbuild as mb

class O(mb.Compound):
    def __init__(self):
        super(O, self).__init__()
                
        mb.load('o.pdb', compound=self, relative_to_module=self.__module__)
        self.add(mb.Port(anchor=self[0], orientation=[0, 0, 1], 
                         separation=.075), label='up')
        self.add(mb.Port(anchor=self[0], orientation=[0, 0, -1], 
                         separation=.075), label='down')
        
        
if __name__ == '__main__':
    o = O()
    