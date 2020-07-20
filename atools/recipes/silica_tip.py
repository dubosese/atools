from __future__ import division

import inspect
import math
import numpy as np
import os
import random

import mbuild as mb
from mbuild.lib.bulk_materials import AmorphousSilicaBulk

from atools.lib import surface_cache

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
    tip_radius :  float, optional, default=1.0
        Radius of asperity (in nm)

    """

    def __init__(self, tip_radius=1.0, seed=12345):
        super(SilicaTip, self).__init__()

        tip_radius = float(tip_radius)

        self._oh_density = 5.0
        self._O_buffer = 0.275
        self._radius = tip_radius

        cache_dir = os.path.dirname(inspect.getfile(surface_cache))
        for file in os.listdir(cache_dir):
            if file == 'silica-tip-r{}-seed{}.mol2'.format(tip_radius, seed):
                mb.load(os.path.join(cache_dir, file), compound=self)
                self._add_ports()

        if not self.children:
            random.seed(seed)

            self._cleave_tip(AmorphousSilica())
            self.generate_bonds(name_a='Si', name_b='O', dmin=0.0, dmax=0.20419)
            self._strip_stray_atoms()
            self._bridge_dangling_Os()
            self._identify_surface_sites()
            self._add_ports()
            self._adjust_stoichiometry()

            filename = 'silica-tip-r{}-seed{}.mol2'.format(tip_radius, seed)
            self.save(os.path.join(cache_dir, filename))

    def _cleave_tip(self, bulk_silica):
        """ Carve tip from bulk silica, include a buffer of O's around the
            tip to ensure the tip is coated.
        """

        O_buffer = self._O_buffer
        tip_radius = self._radius
        tile = [int(math.ceil((2*O_buffer + 2*tip_radius) / dim)) 
                for dim in bulk_silica.periodicity]
        bulk = mb.TiledCompound(bulk_silica, n_tiles=tile)

        tip = mb.Compound(periodicity=np.zeros(3))
        center = bulk.periodicity / 2

        check_radius = {'O': tip_radius + O_buffer, 'Si': tip_radius}
        check_z = {'O': center[2] - O_buffer, 'Si': center[2]}
        for i, particle in enumerate(bulk.particles()):
            if (_check_sphere(particle.pos, center, check_radius[particle.name]) 
                    and particle.pos[2] > check_z[particle.name]):
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

    def _bridge_dangling_Os(self):
        """ Create Si-O-Si bridges on the surface to yield the desired
            density of reactive surface sites
        """

        O_buffer = self._O_buffer
        oh_density = self._oh_density
        tip_radius = self._radius

        area = 2 * np.pi * tip_radius**2
        target = int(oh_density * area)

        dangling_Os = [atom for atom in self.particles()
                       if atom.name == 'O' and
                       atom.pos[2] > np.min(self.xyz[:,2]) + O_buffer and
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

        O_buffer = self._O_buffer
        tip_radius = self._radius

        center = np.array([(np.max(self.xyz[:,0]) + np.min(self.xyz[:,0])) / 2,
                           (np.max(self.xyz[:,1]) + np.min(self.xyz[:,1])) / 2,
                           (np.min(self.xyz[:,2]) + O_buffer)])

        for atom in list(self.particles()):
            if len(self.bond_graph.neighbors(atom)) == 1:
                if (atom.name == 'O' and 
                        atom.pos[2] > np.min(self.xyz[:,2]) + O_buffer and
                        not _check_sphere(atom.pos, center, tip_radius - O_buffer)):
                    atom.name = 'OS'

    def _add_ports(self):
        """ Add ports above surface sites. """

        O_buffer = self._O_buffer

        center = np.array([(np.max(self.xyz[:,0]) + np.min(self.xyz[:,0])) / 2,
                           (np.max(self.xyz[:,1]) + np.min(self.xyz[:,1])) / 2,
                           (np.min(self.xyz[:,2]) + O_buffer)])

        for atom in self.particles_by_name('OS'):
            port = mb.Port(anchor=atom, orientation=atom.pos - center,
                           separation=0.1)
            self.add(port, "port_{}".format(len(self.referenced_ports())))

    def _adjust_stoichiometry(self):
        """ Remove O's from the bottom of the tip to yield a 2:1 Si:O ratio """

        O_buffer = self._O_buffer

        num_O = len(list(self.particles_by_name('O')))
        num_Si = len(list(self.particles_by_name('Si')))
        n_deletions = num_O - 2*num_Si

        top_Os = [atom for atom in self.particles()
                  if atom.name == 'O' and
                     atom.pos[2] < np.min(self.xyz[:,2]) + O_buffer and
                     len(self.bond_graph.neighbors(atom)) == 1]

        for _ in range(n_deletions):
            O1 = random.choice(top_Os)
            top_Os.remove(O1)
            self.remove(O1)

if __name__ == "__main__":
    tip = SilicaTip(tip_radius=2.0)
    tip.save('tip.mol2')
