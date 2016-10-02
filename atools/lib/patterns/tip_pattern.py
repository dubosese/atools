import sys
import os
import mbuild as mb
import numpy as np


class TipPattern(mb.Pattern):
    def __init__(self, tip_radius=2, pattern_num=1, orientations=None):

        script_path = os.path.realpath(sys.modules[self.__module__].__file__)
        file_dir = os.path.dirname(script_path)
        filename = os.path.join(file_dir,'tip_{}nm_pattern{}.xyz'.format(int(tip_radius),int(pattern_num)))
        
        points = np.loadtxt(filename,skiprows=2)[:,1:]

        super(TipPattern, self).__init__(points=points, orientations=None)
