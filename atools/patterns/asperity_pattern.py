from __future__ import division

from math import cos, sin
import numpy as np

import mbuild as mb

from atools.patterns.tip_pattern import TipPattern


def _check_circle(xy, center, radius):
    check = sum((xy[i]-coord)**2 for i,coord in enumerate(center))
    if check <= radius ** 2:
        return True
    else:
        return False

class AsperityPattern(mb.Pattern):
    def __init__(self, density_asperity, asperity_radius, density_surface,
                 xlo, xhi, ylo, yhi, seed=12345):

        center = ((xlo + xhi) / 2, (ylo + yhi) / 2)

        tip_pattern = TipPattern(density_asperity, asperity_radius, seed=seed)

        points = np.stack((x, y, z), axis=1)

        super(RandomHemispherePattern, self).__init__(points=points, **kwargs)

if __name__ == "__main__":
    pattern = RandomHemispherePattern(n=1000, radius=1.0)

    lj_proto = mb.Compound(name='LJ')
    lj_box = mb.Compound()
    for pos in pattern:
        lj_particle = mb.clone(lj_proto)
        lj_particle.translate(pos)
        lj_box.add(lj_particle)
    lj_box.save('hemisphere.xyz', overwrite=True)
