import numpy as np
from math import atan2, cos, sin
from atools.math_azs import anint

class NeighborList(object):

    def __init__(self, coords, box, cutoff, build="full"):
        self.coords = coords
        self.box = box
        self.cutoff = cutoff
        if build == "full":
            self.neighbors, self.neighbor_coords, self.bod = self.build_full()
        '''
        else:
            assert build == "half",
                   "Build must be either 'full' or 'half'"
            self.neighbors, self.neighbor_coords, self.bod = build.half()
        '''


    def build_full(self):
        neighbors = [[] for atom in self.coords]
        neighbor_coords = [[] for atom in self.coords]
        bod = []
        for i,atom1 in enumerate(self.coords):
            for j,atom2 in enumerate(self.coords):
                if i != j:
                    dr = np.array([(atom1[k]-atom2[k]) for k in range(3)])
                    dr -= [self.box.lengths[k]*anint(dr[k]/self.box.lengths[k])]
                    r = np.linalg.norm(np.zeros(3) - dr)

                    if r < self.cutoff:
                        neighbors[i].append(j)
                        neighbors[j].append(i)
                        neighbor_coords[i].append(dr)
                        neighbor_coords[j].append(dr*-1.)

                        bod.append(dr/(r**2))
                        bod.append((dr*-1.)/(r**2))

        return neighbors, neighbor_coords, bod

    def calc_hex_order(self):
        c = complex(0.0,0.0)
        for coord in self.bod:
            theta = atan2(coord[1],coord[0]) + np.pi
            c += complex(cos(6.*theta), -sin(6.*theta))
        c /= len(self.bod)
        c = np.abs(c)
        return np.real(c)
