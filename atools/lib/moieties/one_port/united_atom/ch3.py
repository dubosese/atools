import numpy as np
import mbuild as mb

class CH3(mb.Compound):
    """ """
    def __init__(self):
        super(CH3, self).__init__()
        self.add(mb.Particle(name='_CH3'))

        self.add(mb.Port(anchor=self[0], separation=0.075), 'up')

if __name__ == '__main__':
    ua = CH3()
