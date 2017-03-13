import mbuild as mb


class SurfaceMonolayer(mb.Compound):

    def __init__(self, surface, chains, n_chains, seed, fractions=None,
                 backfill=None):
        super(SurfaceMonolayer, self).__init__()

        if surface.name == 'SilicaInterface':
            pattern = mb.Random2DPattern(n_chains, seed=seed)
        elif surface.name == 'SilicaTip':
            pass
        elif surface.name == 'SilicaAsperity'
            pass
        
        monolayer = mb.Monolayer(surface=surface, chains=chains, pattern=pattern,
                                 fractions=fractions, backfill=backfill)

        self.add(monolayer)

if __name__ == "__main__":
    from atools.recipes.silica_interface import SilicaInterface
    from atools.recipes.silica_asperity import SilicaAsperity
    from atools.recipes.silica_tip import SilicaTip

    from mbuild.examples.alkane_monolayer.alkylsilane import AlkylSilane
    from mbuild.lib.atoms import H

    fractions = [0.75, 0.25]

    chain_a = AlkylSilane(6)
    chain_b = AlkylSilane(18)
    hydrogen = H()

    seed = 12345

    '''
    planar_surface = SilicaInterface(seed=seed)
    planar_monolayer = SurfaceMonolayer(surface=planar_surface, 
        chains=[chain_a, chain_b], n_chains=100, seed=seed, fractions=fractions,
        backfill=hydrogen)
    planar_monolayer.save('planar.mol2', overwrite=True)
    '''

    tip = SilicaTip(tip_radius=2.0, seed=seed)
    tip_monolayer = SurfaceMonolayer(surface=tip, chains=[chain_a, chain_b], 
        n_chains=50, seed=seed, fractions=fractions, backfill=hydrogen)
    tip_monolayer.save('tip.mol2', overwrite=True)

    '''
    asperity = SilicaAsperity(tile_x=3, tile_y=3, asperity_radius=2.0, seed=seed)
    asperity_monolayer = SurfaceMonolayer(surface=asperity, 
        chains=[chain_a, chain_b], n_chains=100, seed=seed, fractions=fractions,
        backfill=hydrogen)
    asperity_monolayer.save('asperity.mol2', overwrite=True)
    '''
