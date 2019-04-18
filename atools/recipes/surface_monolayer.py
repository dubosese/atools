import os
from pkg_resources import resource_filename

import mbuild as mb

from atools.patterns import RandomHemispherePattern
from atools.recipes import Monolayer

class SurfaceMonolayer(mb.Compound):

    def __init__(self, surface, chains, n_chains, seed, fractions=None,
                 backfill=None, **kwargs):
        super(SurfaceMonolayer, self).__init__()

        if surface.name == 'SilicaInterface' and n_chains > 0:
            pattern = mb.Random2DPattern(n_chains, seed=seed)
        elif surface.name == 'SilicaTip':
            shift = 0.25
            surface.translate_to([0, 0, 0])
            surface.translate([0, 0, -1 * min(surface.xyz[:,2]) - shift])
            radius = max(surface.xyz[:,2]) - min(surface.xyz[:,2])
            pattern = RandomHemispherePattern(n=n_chains, seed=seed)
            pattern.scale(radius)
        elif surface.name == 'SilicaAsperity':
            pass
        else:
            pattern = mb.Random2DPattern(n_chains, seed=seed)

        if chains and n_chains > 0:
            monolayer = Monolayer(surface=surface, chains=chains, pattern=pattern,
                                  fractions=fractions, backfill=backfill, **kwargs)
        else:
            monolayer = Monolayer(surface=surface, chains=backfill, guest_port_name='up', **kwargs)

        self.add(monolayer)


if __name__ == "__main__":
    from atools.lib.chains import Alkylsilane
    from atools.recipes.silica_interface import SilicaInterface
    from atools.recipes.silica_asperity import SilicaAsperity
    from atools.recipes.silica_tip import SilicaTip

    from mbuild.lib.atoms import H

    hydrogen = H()
    seed = 12345

    terminal_groups = ['nitrophenyl', 'perfluoromethyl', 'phenyl', 'pyrrole']

    for terminal_group in terminal_groups:
        chain = Alkylsilane(chain_length=6, terminal_group=terminal_group)

        planar_surface = SilicaInterface(seed=seed)
        planar_monolayer = SurfaceMonolayer(surface=planar_surface,
            chains=chain, n_chains=100, seed=seed, backfill=hydrogen)
        forcefield_dir = resource_filename('atools', 'forcefields')
        planar_monolayer.save('planar-{}.gro'.format(terminal_group), overwrite=True)
        planar_monolayer.save('planar-{}.top'.format(terminal_group),
            forcefield_files=os.path.join(forcefield_dir, 'oplsaa-silica.xml'),
            overwrite=True)

    '''
    tip = SilicaTip(tip_radius=2.0, seed=seed)
    tip_monolayer = SurfaceMonolayer(surface=tip, chains=[chain_a, chain_b],
        n_chains=50, seed=2, fractions=fractions, backfill=hydrogen)
    tip_monolayer.save('tip.mol2', overwrite=True)
    '''

    '''
    asperity = SilicaAsperity(tile_x=3, tile_y=3, asperity_radius=2.0, seed=seed)
    asperity_monolayer = SurfaceMonolayer(surface=asperity,
        chains=[chain_a, chain_b], n_chains=100, seed=seed, fractions=fractions,
        backfill=hydrogen)
    asperity_monolayer.save('asperity.mol2', overwrite=True)
    '''
