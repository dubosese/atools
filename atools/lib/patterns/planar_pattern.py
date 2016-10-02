import sys
import os
import mbuild as mb
import numpy as np


class PlanarPattern(mb.Pattern):
    def __init__(self, tile_x, tile_y, surface_id=1, pattern_id=1, orientations=None):

        script_path = os.path.realpath(sys.modules[self.__module__].__file__)
        file_dir = os.path.dirname(script_path)
        filename = os.path.join(file_dir,'planar_pattern-{}x{}_{}_{}.xyz'.format(int(tile_x),int(tile_y),int(surface_id),int(pattern_id)))
        
        points = np.loadtxt(filename,skiprows=2)[:,1:]

        super(PlanarPattern, self).__init__(points=points, orientations=None)
