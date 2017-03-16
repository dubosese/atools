import os
from pkg_resources import resource_filename

import numpy as np

import mbuild as mb

class DualSurface(mb.Compound):
    """ A recipe for creating a system with two opposing surfaces.
    """

    def __init__(self, bottom, top=None, separation=0.8):
        super(DualSurface, self).__init__()

        if top is None:
            top = mb.clone(bottom)
        
        top.spin(np.pi, [0, 1, 0])
        top_box = top.boundingbox
        bot_box = bottom.boundingbox

        top_of_bot = bot_box.maxs[2]
        bot_of_top = top_box.mins[2]

        top.translate([0, 0, top_of_bot - bot_of_top + separation])

        self.add(bottom)
        self.add(top)

if __name__ == "__main__":
    from atools.lib.chains import Alkylsilane
    from atools.recipes.silica_interface import SilicaInterface
    from atools.recipes.surface_monolayer import SurfaceMonolayer

    from mbuild.lib.atoms import H

    fractions = [0.75, 0.25]

    chain_a = Alkylsilane(chain_length=6, terminal_group='amino')
    chain_b = Alkylsilane(chain_length=18, terminal_group='carboxyl')
    hydrogen = H()

    seed = 12345

    planar_surface = SilicaInterface(seed=seed)
    planar_monolayer = SurfaceMonolayer(surface=planar_surface,
        chains=[chain_a, chain_b], n_chains=100, seed=seed, fractions=fractions,
        backfill=hydrogen)
    dual_monolayer = DualSurface(planar_monolayer)
    forcefield_dir = resource_filename('atools', 'forcefields')
    box = dual_monolayer.boundingbox
    dual_monolayer.periodicity += np.array([0, 0, 5. * box.lengths[2]])
    dual_monolayer.save('dual_mono.gro',
        forcefield_files=os.path.join(forcefield_dir, 'oplsaa-silica.xml'),
        overwrite=True)
    dual_monolayer.save('dual_mono.top',
        forcefield_files=os.path.join(forcefield_dir, 'oplsaa-silica.xml'),
        overwrite=True)
