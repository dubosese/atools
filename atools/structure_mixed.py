from __future__ import division

import mdtraj as md
from mdtraj.geometry.order import _compute_director
import numpy as np
from scipy.integrate import simps
from scipy.stats import binned_statistic_2d

from atools.fileio import read_ndx


def calc_nematic_order(traj_filename, top_filename, output_filename,
                       ndx_filename, n_chains):
    """Calculate the nematic order of both monolayers in a two monolayer system

    Returns the nematic order of each monolayer at each frame of a trajectory
    in a dual monolayer system.

    Parameters
    ----------
    traj_filename : str
        Name of trajectory file
    top_filename : str
        Name of topology file
    output_filename : str
        Name of output file
    ndx_filename : str
        Name of Gromacs .ndx file which specifies atom groups
    n_chains : int
        Number of chains per monolayer

    """
    groups = read_ndx(ndx_filename)
    bottom_chains = [id-1 for id in groups['bottom_chains']]
    bottom_chains = [chain.tolist() for chain in
        np.array_split(bottom_chains, n_chains)]
    top_chains = [id-1 for id in groups['top_chains']]
    top_chains = [chain.tolist() for chain in np.array_split(top_chains, n_chains)]

    traj = md.load(traj_filename, top=top_filename)
    S2_bottom = md.compute_nematic_order(traj, indices=bottom_chains)
    S2_top = md.compute_nematic_order(traj, indices=top_chains)
    np.savetxt(output_filename, np.column_stack((traj.time, S2_bottom, S2_top)),
        header='Time\tBottom\tTop')

def _tilt_angle(v1, v2):
    cosang = np.dot(v1, v2)
    sinang = np.linalg.norm(np.cross(v1, v2))
    angle = np.arctan2(sinang, cosang) * (180.0 / np.pi)
    if angle > 90 and angle < 180:
        angle -= 180
    return abs(angle)

def calc_avg_tilt_angle(traj_filename, top_filename, output_filename,
                        ndx_filename, n_chains):
    """Calculate the avg. tilt angle of both monolayers in a two monolayer system

    Returns the average tilt angle of each monolayer at each frame of a trajectory
    in a dual monolayer system.

    Parameters
    ----------
    traj_filename : str
        Name of trajectory file
    top_filename : str
        Name of topology file
    output_filename : str
        Name of output file
    ndx_filename : str
        Name of Gromacs .ndx file which specifies atom groups
    n_chains : int
        Number of chains per monolayer

    """
    groups = read_ndx(ndx_filename)
    bottom_chains = [id-1 for id in groups['bottom_chains']]
    bottom_chains = [chain.tolist() for chain in
        np.array_split(bottom_chains, n_chains)]
    top_chains = [id-1 for id in groups['top_chains']]
    top_chains = [chain.tolist() for chain in np.array_split(top_chains, n_chains)]

    traj = md.load(traj_filename, top=top_filename)
    
    tilts = []
    tilts_err = []
    for monolayer in [bottom_chains, top_chains]:
        monolayer_tilt = []
        for i, chain in enumerate(monolayer):
            chain_slice = traj.atom_slice(chain)
            directors = _compute_director(chain_slice)
            monolayer_tilt.append([_tilt_angle(director, [0.0, 0.0, 1.0])
                for director in directors])
        monolayer_tilt = np.asarray(monolayer_tilt)
        monolayer_tilt = np.transpose(monolayer_tilt)
        tilts_by_frame = np.mean(monolayer_tilt, axis=1)
        tilts_err_by_frame = np.std(monolayer_tilt, axis=1)
        tilts.append(tilts_by_frame)
        tilts_err.append(tilts_err_by_frame)

    tilts = np.asarray(tilts)
    tilts = np.transpose(tilts)

    np.savetxt(output_filename, np.column_stack((traj.time, tilts)),
        header='Time\tBottom\tTop')

def calc_interdigitation(traj_filename, top_filename, output_filename,
                         ndx_filename, bin_size=0.025):
    groups = read_ndx(ndx_filename)
    bottom_monolayer = [id-1 for id in groups['bottom_chains']]
    top_monolayer = [id-1 for id in groups['top_chains']]

    traj = md.load(traj_filename, top=top_filename)

    top_z = traj.xyz[:, top_monolayer, 2]
    bottom_z = traj.xyz[:, bottom_monolayer, 2]
    masses = np.array([atom.element.mass for atom in traj.top.atoms])
    top_masses = masses[top_monolayer]
    bottom_masses = masses[bottom_monolayer]
    top_max = np.max(top_z, axis=1)
    bottom_min = np.min(bottom_z, axis=1)
    bins = [int((zmax-zmin)/bin_size) for zmax, zmin in zip(top_max, bottom_min)]
    interdigitation = []
    for i, (top_frame_z, bottom_frame_z) in enumerate(zip(top_z, bottom_z)):
        hist_top, bin_edges = np.histogram(top_frame_z, bins=bins[i],
            range=(bottom_min[i], top_max[i]), normed=False, weights=top_masses)
        hist_bottom, bin_edges = np.histogram(bottom_frame_z, bins=bins[i],
            range=(bottom_min[i], top_max[i]), normed=False, weights=bottom_masses)
        hist_overlap = []
        bin_middles = [(bin_edges[j+1] + edge)/2 for j, edge in enumerate(bin_edges[:-1])]
        for count_top, count_bottom in zip(hist_top, hist_bottom):
            count_prod = count_top * count_bottom
            count_sum = count_top + count_bottom
            overlap = 0
            if count_sum != 0:
                overlap = (4 * count_prod) / (count_sum ** 2)
            hist_overlap.append(overlap)
        hist_overlap = np.asarray(hist_overlap)
        lam = simps(hist_overlap, bin_middles)
        interdigitation.append(lam)
    np.savetxt(output_filename, np.column_stack((traj.time, interdigitation)))

def calc_friction(trr_filename, output_filename, ndx_filename):
    import MDAnalysis as mda

    groups = read_ndx(ndx_filename)
    bottom = groups['bottom']
    bottom = np.array(bottom)
    bottom -= 1

    fric = []
    trr = mda.coordinates.TRR.TRRReader(trr_filename)
    for frame in trr:
        forces_on_bottom = frame.forces[bottom]
        fric.append([frame.time, np.sum(forces_on_bottom[:,0]) * 0.0166])
    np.savetxt(output_filename, fric)

def count_hydrogen_bonds(traj_filename, top_filename, output_filename,
                         mol2_filename, ndx_filename, top_group,
                         bottom_group):
    """Count the number of inter- or intra-monolayer hydrogen bonds

    The number of inter- and intra-monolayer hydrogen bonds is determined for each
    frame in a trajectory using the Wernet-Nilsson method implemented in MDTraj.

    Parameters
    ----------
    traj_filename : str
        Name of trajectory file
    top_filename : str
        Name of topology file
    output_filename : str
        Name of output file
    mol2_filename : str
        Name of mol2 file used to read bond information
    ndx_filename : str
        Name of Gromacs .ndx file which specifies atom groups
    top_group : str
        Only atom indices from this group (read from a .ndx file) will be considered
        as part of the top monolayer.
    bottom_group : str
        Only atom indices from this group (read from a .ndx file) will be considered
        as part of the bottom monolayer.

    """
    top = md.load(top_filename).topology
    top_atoms = [atom for atom in top.atoms]

    groups = read_ndx(ndx_filename)

    top_monolayer = np.array(groups[top_group]) - 1
    bottom_monolayer = np.array(groups[bottom_group]) - 1
    monolayers = np.hstack((bottom_monolayer, top_monolayer))

    h_bonds_top = []
    h_bonds_bottom = []
    h_bonds_interface = []
    time_traj = []
    for traj_chunk in md.iterload(traj_filename, top=mol2_filename, chunk=10):
        for i, atom in enumerate(traj_chunk.top.atoms):
            atom.element = top_atoms[i].element

        h_bonds_total = md.wernet_nilsson(traj_chunk)
        for frame in h_bonds_total:
            h_bonds_top_frame = [tuple(bond) for bond in frame
                                 if all(atom_id in top_monolayer for atom_id in bond)]
            h_bonds_bottom_frame = [tuple(bond) for bond in frame
                                    if all(atom_id in bottom_monolayer
                                    for atom_id in bond)]
            h_bonds_interface_frame = [tuple(bond) for bond in frame
                                       if(all(atom_id in monolayers for atom_id in bond)
                                       and tuple(bond) not in h_bonds_top_frame
                                       and tuple(bond) not in h_bonds_bottom_frame)]
            h_bonds_top.append(len(h_bonds_top_frame))
            h_bonds_bottom.append(len(h_bonds_bottom_frame))
            h_bonds_interface.append(len(h_bonds_interface_frame))
        time_traj.append(traj_chunk.time)

    time_traj = np.concatenate(time_traj).ravel()
    h_bonds_top = np.asarray(h_bonds_top)
    h_bonds_bottom = np.asarray(h_bonds_bottom)
    h_bonds_interface = np.asarray(h_bonds_interface)
    np.savetxt(output_filename,
        np.column_stack((time_traj, h_bonds_interface, h_bonds_top, h_bonds_bottom)),
        header='Time\tInter-\tIntra-top\tIntra-bottom')


def interfacial_fraction(traj_filename, top_filename, output_filename,
                         ndx_filename, n_chains):
    """Estimate the fraction of chains where termini are at the interface

    Some description of what is returned

    Parameters
    ----------
    traj_filename : str
        Name of trajectory file
    top_filename : str
        Name of topology file
    output_filename : str
        Name of output file
    ndx_filename : str
        Name of Gromacs .ndx file which specifies atom groups
    n_chains : int
        Number of chains per monolayer

    """
    traj = md.load(traj_filename, top=top_filename)

    interfacial_frac_bottom = []
    interfacial_frac_top = []

    groups = read_ndx(ndx_filename)
    for group in ['bottom_termini', 'top_termini']:
        termini = np.array(ndx[group]) - 1
        termini = np.array_split(termini, n_chains)
        xyz = []
        for chain in termini:
            chain_slice = traj.atom_slice(chain)
            chain_com = md.compute_center_of_mass(chain_slice)
            xyz.append(chain_com[0])
        termini_z = np.array([pos[2] for pos in bottom_xyz])
        if group == 'bottom_termini':
            cut = np.median(termini_z) - 0.5
            termini_xy = np.array([coords[:2] for i, coords in enumerate(xyz)
                                   if termini_z[i] > cut])
            interfacial_frac_bottom.append(len(termini_xy) / 100)
        else:
            cut = np.median(termini_z) + 0.5
            termini_xy = np.array([coords[:2] for i, coords in enumerate(xyz)
                                   if termini_z[i] < cut])
            interfacial_frac_top.append(len(termini_xy) / 100)

    print('Under construction!')
    '''
    np.savetxt(output_filename,
        np.column_stack((traj.time, interfacial_frac_top, interfacial_frac_bottom)),
        header='Time\tTop\tBottom')
    '''

def voronoi_termini(traj_filename, top_filename, output_filename,
                    ndx_filename, n_chains):
    """Perform a Voronoi tessellation to estimate coordination number and area

    Some description of what is returned

    Parameters
    ----------
    traj_filename : str
        Name of trajectory file
    top_filename : str
        Name of topology file
    output_filename : str
        Name of output file
    ndx_filename : str
        Name of Gromacs .ndx file which specifies atom groups
    n_chains : int
        Number of chains per monolayer

    """
    traj = md.load(traj_filename, top=top_filename)
    x_length.unitcell_lengths[0,0]
    y_length.unitcell_lengths[0,1]

    coordination_number = []
    coordination_err = []
    interfacial_frac = []

    groups = read_ndx(ndx_filename)
    for group in ['bottom_termini', 'top_termini']:
        termini = np.array(ndx[group]) - 1
        termini = np.array_split(termini, n_chains)
        xyz = []
        for chain in termini:
            chain_slice = traj.atom_slice(chain)
            chain_com = md.compute_center_of_mass(chain_slice)
            xyz.append(chain_com[0])
        termini_z = np.array([pos[2] for pos in bottom_xyz])
        if group == 'bottom_termini':
            cut = np.median(termini_z) - 0.5
            termini_xy = np.array([coords[:2] for i, coords in enumerate(xyz)
                                   if termini_z[i] > cut])
        else:
            cut = np.median(termini_z) + 0.5
            termini_xy = np.array([coords[:2] for i, coords in enumerate(xyz)
                                   if termini_z[i] < cut])
        termini_xy = [arr.tolist() for arr in termini_xy]
        cells = pyvoro.compute_2d_voronoi(points=termini_xy,
            limits=[[0., x_length], [0., y_length]], dispersion=2,
            periodic=[False, False])
        coordination_number.append(np.mean([float(len(cell['adjacency']))
                                            for cell in cells]))
        coordination_err.append(np.std([float(len(cell['adjacency']))
                                        for cell in cells]))

    print('Under construction!')
    '''
    np.savetxt(output_filename,
        np.column_stack((time_traj, h_bonds_interface, h_bonds_top, h_bonds_bottom)),
        header='Time\tInter-\tIntra-top\tIntra-bottom')
    '''

def identify_rigid_groups(monolayer, terminal_group=None, freeze_thickness=5):
    bounding_box = monolayer.boundingbox
    bot_of_box = bounding_box.mins[2]
    top_of_box = bounding_box.maxs[2]

    full_system = []
    bot_rigid = []
    top_rigid = []
    if terminal_group:
        nonterminal = []
        terminal = []
    for i, particle in enumerate(monolayer.particles()):
        full_system.append(i + 1)
        if particle.pos[2] < bot_of_box + freeze_thickness:
            bot_rigid.append(i + 1)
        elif particle.pos[2] > top_of_box - freeze_thickness:
            top_rigid.append(i + 1)
        if terminal_group:
            if terminal_group.title() in [parent.name for parent in particle.ancestors()]:
                terminal.append(i + 1)
            else:
                nonterminal.append(i + 1)

    rigid_groups = {'System': full_system,
                    'bot': bot_rigid,
                    'top': top_rigid}
    if terminal_group:
        rigid_groups['nonterminal'] = nonterminal
        rigid_groups['terminal'] = terminal
    return rigid_groups

def identify_rigid_groups_single(monolayer, terminal_group=None, freeze_thickness=5):
    bounding_box = monolayer.boundingbox
    bot_of_box = bounding_box.mins[2]

    full_system = []
    bot_rigid = []
    chains = []
    if terminal_group:
        nonterminal = []
        terminal = []
    for i, particle in enumerate(monolayer.particles()):
        full_system.append(i + 1)
        if particle.pos[2] < bot_of_box + freeze_thickness:
            bot_rigid.append(i + 1)
        if 'SilicaInterface' in [parent.name for parent in particle.ancestors()]:
            chains.append(i + 1)
        if terminal_group:
            if terminal_group.title() in [parent.name for parent in particle.ancestors()]:
                terminal.append(i + 1)
            else:
                nonterminal.append(i + 1)

    rigid_groups = {'System': full_system,
                    'bot': bot_rigid,
                    'chains': chains}
    if terminal_group:
        rigid_groups['nonterminal'] = nonterminal
        rigid_groups['terminal'] = terminal
    return rigid_groups
