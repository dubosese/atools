import numpy as np
import mbuild as mb

class CH4(mb.Compound):
    """ """
    def __init__(self):
        super(CH4, self).__init__()
        self.add(mb.Particle(name='_CH4'))

if __name__ == '__main__':
    ua = CH4()
