import mbuild as mb
from atools.lib.surfaces.silica_surface import SilicaSurface

class Planar(mb.Compound):
    """ An optionally functionalized, planar amorphous
        silica surface.
    """

    def __init__(self, chains, fractions=None, backfill=None, pattern=None,
                 tile_x=1, tile_y=1, surface_id=1, **kwargs):
        super(Planar, self).__init__()

        surface = SilicaSurface(surface_id)

        if pattern:
            pattern.points /= 10.0
            pattern.points -= surface.boundingbox.mins
            pattern.points /= surface.boundingbox.lengths

        monolayer = mb.Monolayer(surface,
                                 chains,
                                 fractions=fractions,
                                 backfill=backfill,
                                 pattern=pattern,
                                 tile_x=tile_x,
                                 tile_y=tile_y,
                                 **kwargs)
        self.add(monolayer)

if __name__ == "__main__":
    from mbuild.examples.alkane_monolayer.alkylsilane import AlkylSilane
    from mbuild.lib.atoms import H

    pattern = mb.Grid2DPattern(30,30)
    fractions = [0.75,0.25]

    chain_a = AlkylSilane(6)
    chain_b = AlkylSilane(18)
    hydrogen = H()
    monolayer = Planar([chain_a,chain_b],fractions=fractions,backfill=hydrogen,
                       pattern=pattern,tile_x=3,tile_y=3)
    monolayer.save('test.mol2')
