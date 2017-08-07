import os
from pkg_resources import resource_filename

import mbuild as mb
import numpy as np


class DualSurface(mb.Compound):
    """ A recipe for creating a system with two opposing surfaces.

    Parameters
    ----------
    bottom : mb.Compound
        The bottom surface (typically a monolayer-functionalized silica surface).
    top : mb.Compound, optional, default=None
        The top surface. If not specified, `bottom` is cloned, rotated, and shifted
        to yield a system of two identical, opposing surfaces.
    separation : float, optional, default=0.8
        The separation (in nm) between the top and bottom surfaces.
    shift : bool, optional, default=False
        Option to shift the top surface so that it has an equivalent center of mass
        in the XY plane to the bottom surface.
    """

    def __init__(self, bottom, top=None, separation=0.8, shift=False):
        super(DualSurface, self).__init__()

        if top is None:
            top = mb.clone(bottom)
        
        top.spin(np.pi, [0, 1, 0])
        top_box = top.boundingbox
        bot_box = bottom.boundingbox

        top_of_bot = bot_box.maxs[2]
        bot_of_top = top_box.mins[2]

        top.translate([0, 0, top_of_bot - bot_of_top + separation])
        if shift:
            top.translate([bottom.pos[0] - top.pos[0],
                           bottom.pos[1] - top.pos[1],
                           0])

        self.add(bottom, 'bottom')
        self.add(top, 'top')

if __name__ == "__main__":
    from atools.lib.chains import Alkylsilane
    from atools.recipes.silica_interface import SilicaInterface
    from atools.recipes.silica_tip import SilicaTip
    from atools.recipes.surface_monolayer import SurfaceMonolayer

    from mbuild.lib.atoms import H

    chain = Alkylsilane(chain_length=12, terminal_group='methyl')
    hydrogen = H()

    seed = 12345

    planar_surface = SilicaInterface(seed=seed, tile_x=2, tile_y=3)
    planar_monolayer = SurfaceMonolayer(surface=planar_surface,
        chains=chain, n_chains=100, seed=seed, backfill=hydrogen)
    top_surface = SilicaTip(tip_radius=2.0, seed=seed)
    tip = SurfaceMonolayer(surface=top_surface, chains=None, n_chains=0,
        seed=seed, backfill=hydrogen)
    dual_monolayer = DualSurface(planar_monolayer, top=tip, shift=True)
    forcefield_dir = resource_filename('atools', 'forcefields')
    box = dual_monolayer.boundingbox
    dual_monolayer.periodicity += np.array([0, 0, 5. * box.lengths[2]])
    dual_monolayer.save('dual_mono.gro',
        forcefield_files=os.path.join(forcefield_dir, 'oplsaa-silica.xml'),
        overwrite=True)
    dual_monolayer.save('dual_mono.top',
        forcefield_files=os.path.join(forcefield_dir, 'oplsaa-silica.xml'),
        overwrite=True)
