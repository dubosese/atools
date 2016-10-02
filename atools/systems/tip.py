import mbuild as mb
from atools.lib.surfaces.silica_tip import SilicaTip

class Tip(mb.Compound):

    def __init__(self, tip_radius, chains, fractions=None, backfill=None, pattern=None, **kwargs):
        super(Tip, self).__init__()

        surface = SilicaTip(tip_radius)

        pattern.points /= 10.0
        pattern.points -= surface.boundingbox.mins
        pattern.points /= surface.boundingbox.lengths

        monolayer = mb.Monolayer(surface,
                                 chains,
                                 fractions=fractions,
                                 backfill=backfill,
                                 pattern=pattern,
                                 tile_x=1,
                                 tile_y=1,
                                 **kwargs)
        self.add(monolayer)

if __name__ == "__main__":
    from mbuild.examples.alkane_monolayer.alkylsilane import AlkylSilane
    from mbuild.lib.atoms import H
    from atools.lib.patterns.tip_pattern import TipPattern

    pattern = TipPattern(2)
    fractions = [0.75,0.25]

    chain_a = AlkylSilane(6)
    chain_b = AlkylSilane(18)
    hydrogen = H()
    monolayer = Tip(2,[chain_a,chain_b],fractions=fractions,backfill=hydrogen,pattern=pattern)
    monolayer.save('test.mol2')
