from __future__ import division

import mbuild as mb
import networkx as nx
import numpy as np

from random import randint, choice
from math import ceil

def _check_sphere(xyz, center, radius):
    check = sum((xyz[i]-coord)**2 for i,coord in enumerate(center))
    if check <= radius ** 2:
        return True
    else:
        return False

class SilicaTip(mb.Compound):
    """ A recipe for creating a hemispherical tip from bulk silica,

    Carves silica tip from bulk, adjusts to desired surface
    hydroxyl density by creating Si-O-Si bridges, and yields a 2:1
    Si:O ratio (excluding surface binding sites)

    Parameters
    ----------
    bulk_silica : mb.Compound
        Bulk silica from which to cleave an interface
    tip_radius :  float, optional, default=1.0
        Radius of asperity (in nm)

    """

    def __init__(self, bulk_silica, tip_radius=1.0):
        super(SilicaTip, self).__init__()

        _oh_density = 5.0
        _O_buffer = 0.275

        self._cleave_tip(bulk_silica, _O_buffer, tip_radius)
        self.generate_bonds(name_a='Si', name_b='O', dmin=0.0, dmax=0.20419)
        self._strip_stray_atoms()
        self._bridge_dangling_Os(_oh_density, _O_buffer, tip_radius)
        self._identify_surface_sites(_O_buffer, tip_radius)
        import pdb
        pdb.set_trace()
        self._adjust_stoichiometry(_O_buffer)

    def _cleave_tip(self, bulk_silica, O_buffer, tip_radius):
        """ Carve tip from bulk silica, include a buffer of O's around the
            tip to ensure the tip is coated.
        """

        tile = [int(ceil((2*O_buffer + tip_radius) / dim)) for dim in bulk_silica.periodicity]
        bulk = mb.TiledCompound(bulk_silica, n_tiles=tile)

        tip = mb.Compound(periodicity=np.zeros(3))
        center = np.array([(np.max(bulk.xyz[:,dim]) + np.min(bulk.xyz[:,dim]))/2 for dim in range(3)])

        for i, particle in enumerate(bulk.particles()):
            if ((particle.name == 'Si' and _check_sphere(particle.pos, center, tip_radius) and particle.pos[2] < center[2]) or (particle.name == 'O' and _check_sphere(particle.pos, center, tip_radius+O_buffer) and particle.pos[2] < center[2]+O_buffer)):
                tip_particle = mb.Compound(name=particle.name, pos=particle.pos)
                tip.add(tip_particle, particle.name + "_{}".format(i))
        self.add(tip)

    def _strip_stray_atoms(self):
        """ Remove stray atoms and surface pieces """

        components = self.bond_graph.connected_components()
        major_component = max(components, key=len)
        for atom in list(self.particles()):
            if atom not in major_component:
                self.remove(atom)

    def _bridge_dangling_Os(self, oh_density, O_buffer, tip_radius):
        """ Create Si-O-Si bridges on the surface to yield the desired
            density of reactive surface sites
        """

        area_curved = 2 * np.pi * tip_radius**2
        area = area_curved
        target = int(oh_density * area)

        dangling_Os = []
        for atom in list(self.particles()):
            if atom.name == 'O' and atom.pos[2] < np.max(self.xyz[:,2])-O_buffer and len(self.bond_graph.neighbors(atom)) == 1:
                dangling_Os.append(atom)

        n_bridges = int((len(dangling_Os) - target) / 2)

        for _ in range(n_bridges):
            bridged = False
            while not bridged:
                O1 = choice(dangling_Os)
                Si1 = self.bond_graph.neighbors(O1)[0]
                for O2 in dangling_Os:
                    if O2 == O1:
                        continue
                    Si2 = self.bond_graph.neighbors(O2)[0]
                    if Si1 == Si2:
                        continue
                    if any(neigh in self.bond_graph.neighbors(Si2) for neigh in self.bond_graph.neighbors(Si1)):
                        continue
                    r = self.min_periodic_distance(Si1.pos, Si2.pos)
                    if r < 0.45:
                        bridged = True
                        self.add_bond((O1, Si2))
                        dangling_Os.remove(O1)
                        dangling_Os.remove(O2)
                        self.remove(O2)
                        break

    def _identify_surface_sites(self, O_buffer, tip_radius):
        """ Label surface sites and add ports below them """

        center = np.array([(np.max(self.xyz[:,0])+np.min(self.xyz[:,0]))/2,
                           (np.max(self.xyz[:,1])+np.min(self.xyz[:,1]))/2,
                           np.max(self.xyz[:,2])-O_buffer])

        for atom in list(self.particles()):
            if len(self.bond_graph.neighbors(atom)) == 1:
                if atom.name == 'O' and atom.pos[2] < np.max(self.xyz[:,2])-O_buffer and not _check_sphere(atom.pos,center,tip_radius-O_buffer):
                    atom.name = 'OS'
                    port = mb.Port(anchor=atom)
                    #mb.rotate_around_x(port, np.pi/2)
                    mb.translate(port, atom.pos + np.array([0.0, 0.0, -0.1]))
                    self.add(port, "port_{}".format(len(self.referenced_ports())))

    def _adjust_stoichiometry(self, O_buffer):
        """ Remove O's from the top of the tip to yield a 2:1 Si:O ratio """

        num_O = len(list(self.particles_by_name('O')))
        num_Si = len(list(self.particles_by_name('Si')))
        n_deletions = num_O - 2*num_Si

        top_Os = []
        for atom in list(self.particles()):
            if atom.name == 'O' and atom.pos[2] > np.max(self.xyz[:,2])-O_buffer and len(self.bond_graph.neighbors(atom)) == 1:
                top_Os.append(atom)

        for _ in range(n_deletions):
            O1 = choice(top_Os)
            top_Os.remove(O1)
            self.remove(O1)

if __name__ == "__main__":
    from mbuild.lib.bulk_materials import AmorphousSilica
    tip = SilicaTip(bulk_silica=AmorphousSilica(), tip_radius = 2.0)
    traj = tip.to_trajectory()
    traj.save('silica_tip-2nm.pdb')
    #tip.save('silica_tip-2nm.mol2')
    #tip.save('silica_tip-2nm.pdb')
