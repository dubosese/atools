from __future__ import division

import inspect
import math
import numpy as np
import os
import random

import mbuild as mb
from mbuild.lib.bulk_materials import AmorphousSilica

from atools.lib import surface_cache

def _check_sphere(xyz, center, radius):
    check = sum((xyz[i]-coord)**2 for i,coord in enumerate(center))
    if check <= radius ** 2:
        return True
    else:
        return False

class SilicaAsperity(mb.Compound):
    """ A recipe for creating an interface from bulk silica,
        containing a hemispherical asperity

    Carves silica interface from bulk, adjusts to desired surface
    hydoxyl density by creating Si-O-Si bridges, and yields a 2:1
    Si:O ratio (excluding surface binding sites)

    Parameters
    ----------
    tile_x : int, optional, default=1
        Number of times to replicate bulk silica in x-direction
    tile_y : int, optional, default=1
        Number of times to replicate bulk silica in y-direction
    thickness : float, optional, default=1.0
        Thickness of the interface (in nm)
    asperity_radius :  float, optional, default=1.0
        Radius of asperity (in nm)

    """

    def __init__(self, tile_x=1, tile_y=1, thickness=1.0, asperity_radius=1.0, 
                 seed=12345):
        super(SilicaAsperity, self).__init__()

        thickness = float(thickness)
        asperity_radius = float(asperity_radius)

        self._oh_density = 5.0
        self._O_buffer = 0.275
        self._asperity_radius = asperity_radius
        self._thickness = thickness

        cache_dir = os.path.dirname(inspect.getfile(surface_cache))
        for file in os.listdir(cache_dir):
            if file == 'silica-asperity-tx{}-ty{}-thickness{}-r{}-seed{}.mol2'.format(tile_x, tile_y, thickness, asperity_radius, seed):
                mb.load(os.path.join(cache_dir, file), compound=self)
                self._add_ports()

        if not self.children:
            random.seed(seed)

            self._cleave_interface(AmorphousSilica(), tile_x, tile_y)
            self.generate_bonds(name_a='Si', name_b='O', dmin=0.0, dmax=0.20419)
            self._strip_stray_atoms()
            self._bridge_dangling_Os()
            self._identify_surface_sites()
            self._add_ports()
            self._adjust_stoichiometry()

            filename = 'silica-asperity-tx{}-ty{}-thickness{}-r{}-seed{}.mol2'.format(tile_x, tile_y, thickness, asperity_radius, seed)
            self.save(os.path.join(cache_dir, filename))

    def _cleave_interface(self, bulk_silica, tile_x, tile_y):
        """ Carve interface from bulk silica, include a buffer of O's above and
            below the surface to ensure the interface is coated.
        """

        asperity_radius = self._asperity_radius
        O_buffer = self._O_buffer
        thickness = self._thickness        

        block_width = thickness + 2*O_buffer + asperity_radius
        tile_z = int(math.ceil(block_width / bulk_silica.periodicity[2]))
        bulk = mb.TiledCompound(bulk_silica, n_tiles=(tile_x, tile_y, tile_z))

        interface = mb.Compound(periodicity=(bulk.periodicity[0],
            bulk.periodicity[1], 0.0))
        center = np.array([(np.max(bulk.xyz[:,0]) + np.min(bulk.xyz[:,0])) / 2,
                           (np.max(bulk.xyz[:,1]) + np.min(bulk.xyz[:,1])) / 2,
                           thickness + O_buffer])

        for i, particle in enumerate(bulk.particles()):
            if ((particle.name == 'Si' and O_buffer < particle.pos[2] and 
                    (particle.pos[2] < (thickness + O_buffer) or 
                    _check_sphere(particle.pos, center, asperity_radius))) or 
                    (particle.name == 'O' and 
                    (particle.pos[2] < (thickness + 2*O_buffer) or 
                    _check_sphere(particle.pos, center, asperity_radius+O_buffer)))):
                interface_particle = mb.Compound(name=particle.name, pos=particle.pos)
                interface.add(interface_particle, particle.name + "_{}".format(i))
        self.add(interface)

    def _strip_stray_atoms(self):
        """ Remove stray atoms and surface pieces """
        components = self.bond_graph.connected_components()
        major_component = max(components, key=len)
        for atom in list(self.particles()):
            if atom not in major_component:
                self.remove(atom)

    def _bridge_dangling_Os(self):
        """ Create Si-O-Si bridges on the surface to yield the desired
            density of reactive surface sites
        """

        asperity_radius = self._asperity_radius
        oh_density = self._oh_density
        O_buffer = self._O_buffer
        thickness = self._thickness

        area_flat = self.periodicity[0] * self.periodicity[1]
        area_removed = np.pi * asperity_radius**2
        area_asperity = 2 * np.pi * asperity_radius**2
        area = area_flat - area_removed + area_asperity
        target = int(oh_density * area)

        center = np.array([(np.max(self.xyz[:,0]) + np.min(self.xyz[:,0])) / 2,
                           (np.max(self.xyz[:,1]) + np.min(self.xyz[:,1])) / 2,
                           thickness + O_buffer])

        dangling_Os = [atom for atom in self.particles()
                       if atom.name == 'O' and
                       thickness < atom.pos[2] and
                       not _check_sphere(atom.pos, center, asperity_radius-0.15) and
                       len(self.bond_graph.neighbors(atom)) == 1]

        n_bridges = int((len(dangling_Os) - target) / 2)

        for _ in range(n_bridges):
            bridged = False
            while not bridged:
                O1 = random.choice(dangling_Os)
                Si1 = self.bond_graph.neighbors(O1)[0]
                for O2 in dangling_Os:
                    if O2 == O1:
                        continue
                    Si2 = self.bond_graph.neighbors(O2)[0]
                    if Si1 == Si2:
                        continue
                    if any(neigh in self.bond_graph.neighbors(Si2) 
                            for neigh in self.bond_graph.neighbors(Si1)):
                        continue
                    r = self.min_periodic_distance(Si1.pos, Si2.pos)
                    if r < 0.45:
                        bridged = True
                        self.add_bond((O1, Si2))
                        dangling_Os.remove(O1)
                        dangling_Os.remove(O2)
                        self.remove(O2)
                        break

    def _identify_surface_sites(self):
        """ Label surface sites """

        asperity_radius = self._asperity_radius
        O_buffer = self._O_buffer
        thickness = self._thickness

        center = np.array([(np.max(self.xyz[:,0]) + np.min(self.xyz[:,0])) / 2,
                           (np.max(self.xyz[:,1]) + np.min(self.xyz[:,1])) / 2,
                           thickness + O_buffer])

        for atom in list(self.particles()):
            if len(self.bond_graph.neighbors(atom)) == 1:
                if (atom.name == 'O' and atom.pos[2] > thickness and 
                        not _check_sphere(atom.pos, center, asperity_radius-0.15)):
                    atom.name = 'OS'

    def _add_ports(self):
        """ Add ports above surface sites. """
        for atom in self.particles_by_name('OS'):
            port = mb.Port(anchor=atom, orientation=[0, 0, 1],
                           separation=0.1)
            self.add(port, "port_{}".format(len(self.referenced_ports())))

    def _adjust_stoichiometry(self):
        """ Remove O's from the underside of the surface to yield a 2:1 Si:O ratio """

        O_buffer = self._O_buffer

        num_O = len(list(self.particles_by_name('O')))
        num_Si = len(list(self.particles_by_name('Si')))
        n_deletions = num_O - 2*num_Si

        bottom_Os = [atom for atom in self.particles()
                     if atom.name == 'O' and atom.pos[2] < O_buffer and
                        len(self.bond_graph.neighbors(atom)) == 1]

        for _ in range(n_deletions):
            O1 = random.choice(bottom_Os)
            bottom_Os.remove(O1)
            self.remove(O1)

if __name__ == "__main__":
    from mbuild.lib.bulk_materials import AmorphousSilica
    silica_interface = SilicaAsperity(bulk_silica=AmorphousSilica(), tile_x=3, tile_y=3, thickness=1.2, asperity_radius = 2.0)
    silica_interface.save('test.pdb', overwrite=True)
