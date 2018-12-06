from __future__ import division

import mdtraj as md
from mdtraj.core.element import get_by_symbol
from mdtraj.geometry.order import _compute_director
import numpy as np
from scipy.integrate import simps
from scipy.stats import binned_statistic_2d

from atools.fileio import read_ndx

def calc_nematic_order(traj_filename, top_filename, output_filename,
                       ndx_filename, chainlength):
    """Calculate the nematic order of a monolayer

    Returns the average nematic order at each frame of a trajectory
    between the top and bottom monolayers in a dual monolayer system.

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
    chainlength : int
        Number of carbons per chain

    """
    topology = md.load(top_filename).topology
    atoms = np.array(list(topology.atoms))
    atom_names = [atom.name for atom in atoms]

    groups = read_ndx(ndx_filename)
    bottom_chains = [id-1 for id in groups['bottom_chains']
                     if atom_names[id-1] == 'C']
    n_chains = int(len(bottom_chains) / chainlength)
    bottom_chains = [chain.tolist() for chain in
                     np.array_split(bottom_chains, n_chains)]
    top_chains = [id-1 for id in groups['top_chains']
                  if atom_names[id-1] == 'C']
    top_chains = [chain.tolist() for chain in
                  np.array_split(top_chains, n_chains)]
    bottom_chemisorbed = [id-1 for id in groups['bottom_chemisorbed']
                          if atom_names[id-1] == 'C']
    n_chemisorbed = int(len(bottom_chemisorbed) / chainlength)
    bottom_chemisorbed = [chain.tolist() for chain in
                          np.array_split(bottom_chemisorbed, n_chemisorbed)]
    top_chemisorbed = [id-1 for id in groups['top_chemisorbed']
                       if atom_names[id-1] == 'C']
    top_chemisorbed = [chain.tolist() for chain in
                       np.array_split(top_chemisorbed, n_chemisorbed)]
    n_crosslinked = n_chains - n_chemisorbed
    if n_crosslinked > 0:
        bottom_crosslinked = [id-1 for id in groups['bottom_crosslinked']
                              if atom_names[id-1] == 'C']
        bottom_crosslinked = [chain.tolist() for chain in
                              np.array_split(bottom_crosslinked, n_crosslinked)]
        top_crosslinked = [id-1 for id in groups['top_crosslinked']
                           if atom_names[id-1] == 'C']
        top_crosslinked = [chain.tolist() for chain in
                           np.array_split(top_crosslinked, n_crosslinked)]

    traj = md.load(traj_filename, top=top_filename)
    S2_bottom = md.compute_nematic_order(traj, indices=bottom_chains)
    S2_top = md.compute_nematic_order(traj, indices=top_chains)
    S2_bottom_chemisorbed = md.compute_nematic_order(traj, indices=bottom_chemisorbed)
    S2_top_chemisorbed = md.compute_nematic_order(traj, indices=top_chemisorbed)
    if n_crosslinked > 0:
        S2_bottom_crosslinked = md.compute_nematic_order(traj, indices=bottom_crosslinked)
        S2_top_crosslinked = md.compute_nematic_order(traj, indices=top_crosslinked)
        S2_mean_crosslinked = np.mean([S2_bottom_crosslinked, S2_top_crosslinked], axis=0)

    S2_mean_chains = np.mean([S2_bottom, S2_top], axis=0)
    S2_mean_chemisorbed = np.mean([S2_bottom_chemisorbed, S2_top_chemisorbed], axis=0)

    if n_crosslinked > 0:
        np.savetxt(output_filename, np.column_stack((traj.time, S2_bottom, S2_top,
            S2_mean_chains, S2_bottom_chemisorbed, S2_top_chemisorbed,
            S2_mean_chemisorbed, S2_bottom_crosslinked, S2_top_crosslinked,
            S2_mean_crosslinked)), header='Time\tBottom\tTop\tAll-mean\t'+
            'Bottom-chemisorbed\tTop-chemisorbed\tChemisorbed-mean\t'+
            'Bottom-crosslinked\tTop-crosslinked\tCrosslinked-mean')
    else:
        np.savetxt(output_filename, np.column_stack((traj.time, S2_bottom, S2_top,
            S2_mean_chains, S2_bottom_chemisorbed, S2_top_chemisorbed,
            S2_mean_chemisorbed)), header='Time\tBottom\tTop\tAll-mean\t'+
            'Bottom-chemisorbed\tTop-chemisorbed\tChemisorbed-mean')

def _tilt_angle(v1, v2):
    cosang = np.dot(v1, v2) 
    sinang = np.linalg.norm(np.cross(v1, v2))
    angle = np.arctan2(sinang, cosang) * (180.0 / np.pi)
    if angle > 90 and angle < 180:
        angle -= 180 
    return abs(angle)

def calc_avg_tilt_angle(traj_filename, top_filename, output_filename,
                        ndx_filename, chainlength):
    """Calculate the average tilt angle of a monolayer

    Returns the average tilt angle at each frame of a trajectory,
    averaged between the top and bottom monolayers in a dual monolayer 
    system.

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
    chainlength : int
        Number of carbons per chain

    """
    topology = md.load(top_filename).topology
    atoms = np.array(list(topology.atoms))
    atom_names = [atom.name for atom in atoms]

    groups = read_ndx(ndx_filename)
    bottom_chains = [id-1 for id in groups['bottom_chains']
                     if atom_names[id-1] == 'C']
    n_chains = int(len(bottom_chains) / chainlength)
    bottom_chains = [chain.tolist() for chain in
                     np.array_split(bottom_chains, n_chains)]
    top_chains = [id-1 for id in groups['top_chains']
                  if atom_names[id-1] == 'C']
    top_chains = [chain.tolist() for chain in
                  np.array_split(top_chains, n_chains)]
    bottom_chemisorbed = [id-1 for id in groups['bottom_chemisorbed']
                          if atom_names[id-1] == 'C']
    n_chemisorbed = int(len(bottom_chemisorbed) / chainlength)
    bottom_chemisorbed = [chain.tolist() for chain in
                          np.array_split(bottom_chemisorbed, n_chemisorbed)]
    top_chemisorbed = [id-1 for id in groups['top_chemisorbed']
                       if atom_names[id-1] == 'C']
    top_chemisorbed = [chain.tolist() for chain in
                       np.array_split(top_chemisorbed, n_chemisorbed)]
    n_crosslinked = n_chains - n_chemisorbed
    if n_crosslinked > 0:
        bottom_crosslinked = [id-1 for id in groups['bottom_crosslinked']
                              if atom_names[id-1] == 'C']
        bottom_crosslinked = [chain.tolist() for chain in
                              np.array_split(bottom_crosslinked, n_crosslinked)]
        top_crosslinked = [id-1 for id in groups['top_crosslinked']
                           if atom_names[id-1] == 'C']
        top_crosslinked = [chain.tolist() for chain in
                           np.array_split(top_crosslinked, n_crosslinked)]

    traj = md.load(traj_filename, top=top_filename)
    if n_crosslinked > 0:
        tilt_groups = [bottom_chains, top_chains, bottom_chemisorbed,
                       top_chemisorbed, bottom_crosslinked, top_crosslinked]
    else:
        tilt_groups = [bottom_chains, top_chains, bottom_chemisorbed,
                       top_chemisorbed]

    tilts = []
    tilts_err = []
    for group in tilt_groups:
        group_tilt_angles = []
        for i, chain in enumerate(group):
            chain_slice = traj.atom_slice(chain)
            directors = _compute_director(chain_slice)
            group_tilt_angles.append([_tilt_angle(director, [0.0, 0.0, 1.0])
                for director in directors])
        group_tilt_angles = np.asarray(group_tilt_angles)
        group_tilt_angles = np.transpose(group_tilt_angles)
        tilts_by_frame = np.mean(group_tilt_angles, axis=1)
        tilts_err_by_frame = np.std(group_tilt_angles, axis=1)
        tilts.append(tilts_by_frame)
        tilts_err.append(tilts_err_by_frame)

    tilts_mean_chains = np.mean([tilts[0], tilts[1]], axis=0)
    tilts_mean_chemisorbed = np.mean([tilts[2], tilts[3]], axis=0)
    tilts.insert(2, tilts_mean_chains)
    tilts.insert(5, tilts_mean_chemisorbed)
    if n_crosslinked > 0:
        tilts_mean_crosslinked = np.mean([tilts[6], tilts[7]], axis=0)
        tilts.append(tilts_mean_crosslinked)

    tilts = np.asarray(tilts)
    tilts = np.transpose(tilts)

    if n_crosslinked > 0:
        np.savetxt(output_filename, np.column_stack((traj.time, tilts)),
            header='Time\tBottom\tTop\tAll-mean\tBottom-chemisorbed\t'+
            'Top-chemisorbed\tChemisorbed-mean\tBottom-crosslinked\t'+
            'Top-crosslinked\tCrosslinked-mean')
    else:
        np.savetxt(output_filename, np.column_stack((traj.time, tilts)),
            header='Time\tBottom\tTop\tAll-mean\tBottom-chemisorbed\t'+
            'Top-chemisorbed\tChemisorbed-mean')

def _film_level(carbon_z, monolayer):
    n_bins = int((np.max(carbon_z)-np.min(carbon_z))/0.01)
    assert n_bins > 0 
    if n_bins == 0:
        return 0.0 
    else:
        hist, bin_edges = np.histogram(carbon_z,
            bins=n_bins,
            normed=False)
        cdist = np.cumsum(hist)/np.sum(hist)
        if monolayer == 'bottom':
            return bin_edges[np.argmax(cdist > 0.95)]
        elif monolayer == 'top':
            return bin_edges[np.argmin(cdist < 0.05)]

def _film_level_bottom(carbon_z):
    return _film_level(carbon_z, monolayer='bottom')

def _film_level_top(carbon_z):
    return _film_level(carbon_z, monolayer='top')

def calc_roughness(traj_filename, top_filename, output_filename):
    """Calculate the roughness of a monolayer

    Parameters
    ----------
    traj_filename : str
        Name of trajectory file
    top_filename : str
        Name of topology file
    output_filename : str
        Name of output file

    """
    topology = md.load(top_filename).topology
    atoms = np.array(list(topology.atoms))
    half_atoms = int(len(atoms)/2)
    bottom_monolayer = [atom.index for atom in atoms if (atom.name == 'C'
        and atom.index < half_atoms)]
    top_monolayer = [atom.index for atom in atoms if (atom.name == 'C'
        and atom.index >= half_atoms)]

    traj = md.load(traj_filename, top=top_filename)
    xdim = traj.unitcell_lengths[0,0]
    ydim = traj.unitcell_lengths[0,1]
    x_carbons_bottom = traj.xyz[:,bottom_monolayer,0]
    y_carbons_bottom = traj.xyz[:,bottom_monolayer,1]
    z_carbons_bottom = traj.xyz[:,bottom_monolayer,2]
    x_carbons_top = traj.xyz[:,top_monolayer,0]
    y_carbons_top = traj.xyz[:,top_monolayer,1]
    z_carbons_top = traj.xyz[:,top_monolayer,2]
    roughness_bottom = []
    roughness_top = []
    for xbot, ybot, zbot, xtop, ytop, ztop in zip(x_carbons_bottom,
            y_carbons_bottom, z_carbons_bottom, x_carbons_top, y_carbons_top,
            z_carbons_top):
        bottom_levels = binned_statistic_2d(x=xbot, y=ybot, values=zbot,
            statistic=_film_level_bottom, bins=5, range=[[0.0, xdim],[0.0, ydim]])
        top_levels = binned_statistic_2d(x=xtop, y=ytop, values=ztop,
            statistic=_film_level_top, bins=5, range=[[0.0, xdim],[0.0, ydim]])
        roughness_bottom.append(np.std(bottom_levels[0]))
        roughness_top.append(np.std(top_levels[0]))
    np.savetxt(output_filename, np.column_stack((traj.time, roughness_bottom,
        roughness_top)), header='Time\tBottom\tTop')

def _count_gauche(dihedral_angles):
    if len(dihedral_angles) == 0:
        return 0
    else:
        positive_angles = np.array([angle if angle > 0 else np.pi*2 + angle
            for angle in dihedral_angles])
        small_dihedrals = np.where((positive_angles < np.pi/2))[0]
        large_dihedrals = np.where((positive_angles > np.pi*(3/2)))[0]
        return len(small_dihedrals) + len(large_dihedrals)

def calc_gauche_defects(traj_filename, top_filename, output_filename,
                        ndx_filename, chainlength):
    """Calculate the roughness of a monolayer

    Parameters
    ----------
    traj_filename : str
        Name of (unwrapped) trajectory file
    top_filename : str
        Name of topology file
    output_filename : str
        Name of output file

    """
    topology = md.load(top_filename).topology
    atoms = np.array(list(topology.atoms))
    atom_names = [atom.name for atom in atoms]

    groups = read_ndx(ndx_filename)
    bottom_chains = [id for id in groups['bottom_chains'] if atom_names[id-1] == 'C']
    n_chains = int(len(bottom_chains) / chainlength)
    bottom_chains = [chain.tolist() for chain in
        np.array_split(bottom_chains, n_chains)]
    top_chains = [id for id in groups['top_chains'] if atom_names[id-1] == 'C']
    top_chains = [chain.tolist() for chain in
        np.array_split(top_chains, n_chains)]
    bottom_chemisorbed = [id for id in groups['bottom_chemisorbed']
        if atom_names[id-1] == 'C']
    n_chemisorbed = int(len(bottom_chemisorbed) / chainlength)
    bottom_chemisorbed = [chain.tolist() for chain in
        np.array_split(bottom_chemisorbed, n_chemisorbed)]
    top_chemisorbed = [id for id in groups['top_chemisorbed']
        if atom_names[id-1] == 'C']
    top_chemisorbed = [chain.tolist() for chain in
        np.array_split(top_chemisorbed, n_chemisorbed)]
    n_crosslinked = n_chains - n_chemisorbed
    if n_crosslinked > 0:
        bottom_crosslinked = [id for id in groups['bottom_crosslinked']
            if atom_names[id-1] == 'C']
        bottom_crosslinked = [chain.tolist() for chain in
            np.array_split(bottom_crosslinked, n_crosslinked)]
        top_crosslinked = [id for id in groups['top_crosslinked']
            if atom_names[id-1] == 'C']
        top_crosslinked = [chain.tolist() for chain in
            np.array_split(top_crosslinked, n_crosslinked)]

    if n_crosslinked > 0:
        groups = [bottom_chains, top_chains, bottom_chemisorbed,
                  top_chemisorbed, bottom_crosslinked, top_crosslinked]
        dihedral_groups = [[],[],[],[],[],[]]
    else:
        groups = [bottom_chains, top_chains, bottom_chemisorbed,
                  top_chemisorbed]
        dihedral_groups = [[],[],[],[]]

    for i, group in enumerate(groups):
        for chain in group:
            chain_dihedrals = [chain[j:j+4] for j in range(len(chain)-3)]
            for dihedral in chain_dihedrals:
                dihedral_groups[i].append(dihedral)
    
    traj = md.load(traj_filename, top=top_filename)

    gauche_defects = []
    for dihedrals in dihedral_groups:
        dihedral_angles = md.compute_dihedrals(traj, dihedrals, periodic=False)
        gauche_defects_group = []
        for angles in dihedral_angles:
            gauche_defects_group.append(_count_gauche(angles))
        gauche_defects.append(gauche_defects_group)

    gauche_mean_chains = np.mean([gauche_defects[0], gauche_defects[1]], axis=0)
    gauche_mean_chemisorbed = np.mean([gauche_defects[2], gauche_defects[3]], axis=0)
    gauche_defects.insert(2, gauche_mean_chains)
    gauche_defects.insert(5, gauche_mean_chemisorbed)
    if n_crosslinked > 0:
        gauche_mean_crosslinked = np.mean([gauche_defects[6], gauche_defects[7]],
            axis=0)
        gauche_defects.append(gauche_mean_crosslinked)

    gauche_defects = np.asarray(gauche_defects)
    gauche_defects = np.transpose(gauche_defects)

    if n_crosslinked > 0:
        np.savetxt(output_filename, np.column_stack((traj.time, gauche_defects)),
            header='Time\tBottom\tTop\tAll-mean\tBottom-chemisorbed\t'+
            'Top-chemisorbed\tChemisorbed-mean\tBottom-crosslinked\t'+
            'Top-crosslinked\tCrosslinked-mean')
    else:
        np.savetxt(output_filename, np.column_stack((traj.time, gauche_defects)),
            header='Time\tBottom\tTop\tAll-mean\tBottom-chemisorbed\t'+
            'Top-chemisorbed\tChemisorbed-mean')

def calc_gauche_position(traj_filename, top_filename, output_filename,
                         ndx_filename, chainlength):
    """Calculate the roughness of a monolayer

    Parameters
    ----------
    traj_filename : str
        Name of (unwrapped) trajectory file
    top_filename : str
        Name of topology file
    output_filename : str
        Name of output file

    """
    topology = md.load(top_filename).topology
    atoms = np.array(list(topology.atoms))
    atom_names = [atom.name for atom in atoms]

    groups = read_ndx(ndx_filename)
    bottom_chains = [id for id in groups['bottom_chains'] if atom_names[id-1] == 'C']
    n_chains = int(len(bottom_chains) / chainlength)
    bottom_chains = [chain.tolist() for chain in
        np.array_split(bottom_chains, n_chains)]
    top_chains = [id for id in groups['top_chains'] if atom_names[id-1] == 'C']
    top_chains = [chain.tolist() for chain in
        np.array_split(top_chains, n_chains)]
    bottom_chemisorbed = [id for id in groups['bottom_chemisorbed']
        if atom_names[id-1] == 'C']
    n_chemisorbed = int(len(bottom_chemisorbed) / chainlength)
    bottom_chemisorbed = [chain.tolist() for chain in
        np.array_split(bottom_chemisorbed, n_chemisorbed)]
    top_chemisorbed = [id for id in groups['top_chemisorbed']
        if atom_names[id-1] == 'C']
    top_chemisorbed = [chain.tolist() for chain in
        np.array_split(top_chemisorbed, n_chemisorbed)]
    n_crosslinked = n_chains - n_chemisorbed
    if n_crosslinked > 0:
        bottom_crosslinked = [id for id in groups['bottom_crosslinked']
            if atom_names[id-1] == 'C']
        bottom_crosslinked = [chain.tolist() for chain in
            np.array_split(bottom_crosslinked, n_crosslinked)]
        top_crosslinked = [id for id in groups['top_crosslinked']
            if atom_names[id-1] == 'C']
        top_crosslinked = [chain.tolist() for chain in
            np.array_split(top_crosslinked, n_crosslinked)]

    if n_crosslinked > 0:
        groups = [bottom_chains, top_chains, bottom_chemisorbed,
                  top_chemisorbed, bottom_crosslinked, top_crosslinked]
    else:
        groups = [bottom_chains, top_chains, bottom_chemisorbed,
                  top_chemisorbed]

    dihedral_groups = []
    for group in groups:
        group_dihedrals = [[] for i in range(chainlength-3)]
        for chain in group:
            chain_dihedrals = [chain[i:i+4] for i in range(len(chain)-3)]
            for i, dihedral in enumerate(chain_dihedrals):
                group_dihedrals[i].append(dihedral)
        dihedral_groups.append(group_dihedrals)

    traj = md.load(traj_filename, top=top_filename)

    gauche_defects = []
    for group in dihedral_groups:
        gauche_defects_group = []
        for position in group:
            dihedral_angles = md.compute_dihedrals(traj, position, periodic=False)
            gauche_defects_position = []
            for angles in dihedral_angles:
                gauche_defects_position.append(_count_gauche(angles))
            gauche_defects_group.append(gauche_defects_position)
        gauche_defects.append(gauche_defects_group)

    gauche_sum_chains = np.sum([gauche_defects[0], gauche_defects[1]], axis=0)
    gauche_sum_chemisorbed = np.sum([gauche_defects[2], gauche_defects[3]], axis=0)
    gauche_defects.insert(2, gauche_sum_chains)
    gauche_defects.insert(5, gauche_sum_chemisorbed)
    if n_crosslinked > 0:
        gauche_sum_crosslinked = np.sum([gauche_defects[6], gauche_defects[7]],
            axis=0)
        gauche_defects.append(gauche_sum_crosslinked)
        filename_modifiers = ['bottom', 'top', 'all', 'bottom-chemisorbed',
                              'top-chemisorbed', 'chemisorbed', 'bottom-crosslinked',
                              'top-crosslinked', 'crosslinked']
    else:
        filename_modifiers = ['bottom', 'top', 'all', 'bottom-chemisorbed',
                              'top-chemisorbed', 'chemisorbed']

    gauche_defects = np.asarray(gauche_defects)
    for i, group in enumerate(gauche_defects):
        header = 'Time\t' + '\t'.join(str(num) for num in range(len(group)))
        new_filename = output_filename.split('.')[0] + \
            '-{}.'.format(filename_modifiers[i]) + output_filename.split('.')[1]
        np.savetxt(new_filename,
                   np.column_stack((traj.time, np.transpose(group))),
                   header=header)

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

def calc_hexagonal_order_freud(traj_filename, top_filename, output_filename,
                               ndx_filename, chainlength):
    from freud.box import Box 
    from freud.order import HexOrderParameter

    topology = md.load(top_filename).topology
    atoms = np.array(list(topology.atoms))
    atom_names = [atom.name for atom in atoms]

    groups = read_ndx(ndx_filename)
    bottom_chains = [id-1 for id in groups['bottom_chains']
                     if atom_names[id-1] == 'C']
    n_chains = int(len(bottom_chains) / chainlength)
    bottom_chains = np.array_split(bottom_chains, n_chains)
    top_chains = [id-1 for id in groups['top_chains'] if atom_names[id-1] == 'C']
    top_chains = np.array_split(top_chains, n_chains)

    traj = md.load(traj_filename, top=top_filename)

    for i, group in enumerate([bottom_chains, top_chains]):
        group_COM = []
        for chain in group:
            chain_slice = traj.atom_slice(chain)
            chain_com = md.compute_center_of_mass(chain_slice)
            for com in chain_com:
                for j in range(2):
                    if com[j] < 0:
                        com[j] += traj.unitcell_lengths[0,j]
                    elif com[j] > traj.unitcell_lengths[0,j]:
                        com[j] -= traj.unitcell_lengths[0,j]
            group_COM.append(chain_com)
        group_COM = np.array(group_COM).transpose(1,0,2)

        rmax = 0.75
        box = Box(Lx=traj.unitcell_lengths[0,0], Ly=traj.unitcell_lengths[0,1],
                  Lz=traj.unitcell_lengths[0,2])

        order_parameter = []
        for xyz in group_COM:
            hex_order = HexOrderParameter(rmax, 6)
            hex_order.compute(box, xyz.astype(np.float32))
            order_parameter.append(abs(hex_order.getPsi()).mean())
        print(np.mean(order_parameter[int(len(order_parameter)/2):]))


def _angle(v1, v2):
    cosang = np.dot(v1, v2) 
    sinang = np.linalg.norm(np.cross(v1, v2))
    angle = np.arctan2(sinang, cosang)
    '''
    if angle > 90 and angle < 180:
        angle -= 180 
    return abs(angle)
    '''
    return angle


def _create_com_traj(xyz, time, unitcell_lengths, unitcell_angles, masses=None):
    top = md.Topology()
    chain = top.add_chain()
    for i, com in enumerate(xyz[0]):
        com_res = top.add_residue('COM{}'.format(i), chain)
        top.add_atom('COM', get_by_symbol('C'), residue=com_res)
    return md.Trajectory(xyz, topology=top, time=time,
                         unitcell_lengths=unitcell_lengths,
                         unitcell_angles=unitcell_angles)


def calc_OCF(traj_filename, top_filename, output_filetag, ndx_filename,
             chainlength):
    """Calculate the orientational correlation function of monolayer chains

    Parameters
    ----------
    traj_filename : str
        Name of (unwrapped) trajectory file
    top_filename : str
        Name of topology file
    output_filename : str
        Name of output file

    """
    topology = md.load(top_filename).topology
    atoms = np.array(list(topology.atoms))
    atom_names = [atom.name for atom in atoms]

    groups = read_ndx(ndx_filename)
    bottom_chains = [id-1 for id in groups['bottom_chains']
                     if atom_names[id-1] == 'C']
    n_chains = int(len(bottom_chains) / chainlength)
    bottom_chains = [chain.tolist() for chain in
        np.array_split(bottom_chains, n_chains)]
    top_chains = [id-1 for id in groups['top_chains'] if atom_names[id-1] == 'C']
    top_chains = [chain.tolist() for chain in
        np.array_split(top_chains, n_chains)]

    traj = md.load(traj_filename, top=top_filename)
    directors = []
    locations = []
    for i, chain in enumerate(bottom_chains):
        chain_slice = traj.atom_slice(chain)
        chain_director = _compute_director(chain_slice)
        for director in chain_director:
            if director[-1] < 0:
                director *= -1.
        directors.append(chain_director)
        com = md.compute_center_of_mass(chain_slice)
        for com_frame in com:
            for j in range(2):
                if com_frame[j] < 0:
                    com_frame[j] += traj.unitcell_lengths[0,j]
                elif com_frame[j] > traj.unitcell_lengths[0,j]:
                    com_frame[j] -= traj.unitcell_lengths[0,j]
        locations.append(com)
    directors = np.array(directors).transpose(1,0,2)
    film_directors = np.mean(directors, axis=1)
    film_directors /= np.linalg.norm(film_directors, axis=1).reshape(len(film_directors),1)
    locations = np.array(locations).transpose(1,0,2)
    com_traj = _create_com_traj(locations, traj.time, traj.unitcell_lengths,
                                traj.unitcell_angles)
    pairs = com_traj.top.select_pairs('all', 'all')
    distances = md.compute_distances(com_traj, pairs, periodic=True)
    cof = []
    cof_6 = []
    for i, frame in enumerate(directors):
        frame_angles = np.array([_angle(film_directors[i], director) for director in frame])
        dist_hist = np.histogram(distances[i], bins=np.arange(0.25, 2.5, 0.1))
        correlations = np.cos(np.array([2.*(frame_angles[pairs[j,0]] - frame_angles[pairs[j,1]]) for j, dist in enumerate(distances[i])]))
        correlation_hist = np.histogram(distances[i], bins=np.arange(0.25, 2.5, 0.1),
                                        weights=correlations)
        correlations_6 = np.cos(np.array([6.*(frame_angles[pairs[j,0]] - frame_angles[pairs[j,1]]) for j, dist in enumerate(distances[i])]))
        correlation_hist_6 = np.histogram(distances[i], bins=np.arange(0.25, 2.5, 0.1), weights=correlations_6)
        cof.append(correlation_hist[0] / dist_hist[0])
        cof_6.append(correlation_hist_6[0] / dist_hist[0])
    cof = np.array(cof)
    cof = np.mean(cof, axis=0)
    cof_6 = np.array(cof_6)
    cof_6 = np.mean(cof_6, axis=0)
    r_vals = np.arange(0.25, 2.5, 0.1)
    r_vals = np.array([(val + r_vals[i+1])/2 for i, val in enumerate(r_vals[:-1])])
    np.savetxt(output_filetag + '2.txt', np.vstack((r_vals, cof)))
    np.savetxt(output_filetag + '6.txt', np.vstack((r_vals, cof_6)))


def calc_com_2Dmsd(traj_filename, top_filename, output_filetag, ndx_filename,
                 chainlength):
    topology = md.load(top_filename).topology
    atoms = np.array(list(topology.atoms))
    atom_names = [atom.name for atom in atoms]

    groups = read_ndx(ndx_filename)
    bottom_chains = [id-1 for id in groups['bottom_chains']
                     if atom_names[id-1] == 'C']
    n_chains = int(len(bottom_chains) / chainlength)
    bottom_chains = [chain.tolist() for chain in
        np.array_split(bottom_chains, n_chains)]
    top_chains = [id-1 for id in groups['top_chains'] if atom_names[id-1] == 'C']
    top_chains = [chain.tolist() for chain in
        np.array_split(top_chains, n_chains)]

    traj = md.load(traj_filename, top=top_filename)
    locations = []
    for i, chain in enumerate(bottom_chains):
        chain_slice = traj.atom_slice(chain)
        com = md.compute_center_of_mass(chain_slice)
        for com_frame in com:
            for j in range(2):
                if com_frame[j] < 0:
                    com_frame[j] += traj.unitcell_lengths[0,j]
                elif com_frame[j] > traj.unitcell_lengths[0,j]:
                    com_frame[j] -= traj.unitcell_lengths[0,j]
        locations.append(com)
    locations = np.array(locations).transpose(1,0,2)
    com_traj = _create_com_traj(locations, traj.time, traj.unitcell_lengths,
                                traj.unitcell_angles)

    msd = []
    for frame in com_traj:
        msd_frame = 0
        for i, atom in enumerate(com_traj.top.atoms):
            disp = frame.xyz[0,i,:2] - com_traj.xyz[0,i,:2]
            for j, dim in enumerate(disp):
                if dim > com_traj.unitcell_lengths[0,j] / 2:
                    dim -= traj.unitcell_lengths[0,j]
                elif dim < -com_traj.unitcell_lengths[0,j] / 2:
                    dim += traj.unitcell_lengths[0,j]
            msd_frame += np.linalg.norm(dim)**2
        msd_frame /= com_traj.top.n_atoms
        msd.append(msd_frame)

    np.savetxt(output_filetag + '.txt', np.column_stack((traj.time, msd)))
