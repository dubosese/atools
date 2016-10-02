import numpy as np
import mbuild as mb

class DualSurface(mb.Compound):
    """ A recipe for creating a system with two opposing surfaces.
    """

    def __init__(self, top, bottom, top_shift=False):
        super(DualSurface, self).__init__()
        
        mb.spin_y(top,np.pi)
        top_box = top.boundingbox
        bot_box = bottom.boundingbox

        top_of_bot = bot_box.maxs[2]
        bot_of_top = top_box.mins[2]
        if top_shift:
            mb.translate(top, [top_box.lengths[0]/2., 0, top_of_bot-bot_of_top+0.5])
        else:
            mb.translate(top, [0, 0, top_of_bot-bot_of_top+0.5])

        self.add(bottom)
        self.add(top)
