import mbuild as mb
from mbuild.lib.bulk_materials import AmorphousSilica

from atools.lib.surfaces.silica_surface import SilicaSurface

class Planar(mb.Compound):
    """ An optionally functionalized, planar amorphous
        silica surface.
    """

    def __init__(self, chains, n_chains, seed=12345, tile_x=1, tile_y=1, **kwargs):
        super(Planar, self).__init__()

        surface = mb.SilicaInterface(bulk_silica=AmorphousSilica(), tile_x=tile_x, 
                                     tile_y=tile_y, thickness=1.2, seed=seed)
        pattern = mb.Random2DPattern(n_chains, seed=seed)

        monolayer = mb.Monolayer(surface=surface, chains=chains, pattern=pattern,
                                 **kwargs)
        self.add(monolayer)

if __name__ == "__main__":
    from mbuild.examples.alkane_monolayer.alkylsilane import AlkylSilane
    from mbuild.lib.atoms import H

    fractions = [0.75,0.25]

    chain_a = AlkylSilane(6)
    chain_b = AlkylSilane(18)
    hydrogen = H()
    monolayer = Planar(chains=[chain_a,chain_b], n_chains=100, fractions=fractions,
                       backfill=hydrogen, tile_x=1, tile_y=1, seed=12345)
    monolayer.save('test.mol2')
