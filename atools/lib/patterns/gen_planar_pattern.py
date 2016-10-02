from __future__ import division

import random
import mbuild as mb
import numpy as np


class PlanarPattern(mb.Pattern):
    def __init__(self, surface, chain_density):
        center = [surface.periodicity[0]/2, surface.periodicity[1]/2]

        points = []
        for atom in surface.particles():
            if atom.name == 'OS':
                points.append(atom.pos)

        # Areas in square nanometers
        length = surface.periodicity[0]
        width = surface.periodicity[1]
        area = length * width

        chains = int(round(area * chain_density))

        points = random.sample(points,chains)
        points = np.asarray(points) * 10.0
    
        super(PlanarPattern, self).__init__(points=points, orientations=None)

if __name__ == "__main__":
    from atools.lib.surfaces.silica_surface import SilicaSurface

    tile_x = 1
    tile_y = 1
    surface_id = 5
    pattern_num = 1
    chain_density = 4.0 # chains per square nm

    surface = SilicaSurface(surface_id=surface_id)
    pattern = PlanarPattern(surface,chain_density)

    from groupy.mdio import write_xyz
    write_xyz(np.array(pattern.points), np.zeros(len(pattern.points)), 'planar_pattern-{}x{}_{}_{}.xyz'.format(tile_x,tile_y,surface_id,pattern_num))
