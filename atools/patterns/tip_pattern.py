from __future__ import division

import numpy as np

import mbuild as mb

from atools.patterns.random_hemisphere_pattern import RandomHemispherePattern


class TipPattern(RandomHemispherePattern):
    def __init__(self, chain_density, tip_radius, seed=12345):

        area = 2 * np.pi * tip_radius**2
        n = int(round(area * chain_density))

        super(TipPattern, self).__init__(n=n, seed=seed, scale=tip_radius)

if __name__ == "__main__":
    pattern = TipPattern(chain_density=2.0, tip_radius=2.0)

    lj_proto = mb.Compound(name='LJ')
    lj_box = mb.Compound()
    for pos in pattern:
        lj_particle = mb.clone(lj_proto)
        lj_particle.translate(pos)
        lj_box.add(lj_particle)
    lj_box.save('tip-pattern.xyz', overwrite=True)
