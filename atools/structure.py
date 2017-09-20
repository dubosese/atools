from __future__ import division

import mdtraj as md
from mdtraj.geometry.order import _compute_director
import numpy as np

from atools.fileio import read_ndx


def calc_nematic_order(traj_filename, top_filename, output_filename,
                       n_chains, average=True):
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
    n_chains : int
        Number of monolayer chains per surface
    average : bool, optional, default=True
        Average nematic order values for the top and bottom monolayers

    Notes
    -----
    Assumes a unique chain prototype.
    Assumes identical top and bottom monolayers.

    """
    topology = md.load(top_filename).topology
    atoms = np.array(list(topology.atoms))
    atom_names = [atom.name for atom in atoms]
    monolayer_begins = atom_names.index('C')
    half_atoms = int(len(atoms)/2)
    backfill_sites = int(atom_names.count('OS')/2 - n_chains)
    bottom_monolayer = [atom.index for atom in atoms[monolayer_begins:half_atoms - backfill_sites]]
    top_monolayer = [atom.index for atom in atoms[monolayer_begins + half_atoms:-backfill_sites]]
    bottom_chains = np.array_split(bottom_monolayer, n_chains)
    top_chains = np.array_split(top_monolayer, n_chains)
    bottom_chains_list = [chain.tolist() for chain in bottom_chains]
    top_chains_list = [chain.tolist() for chain in top_chains]

    traj = md.load(traj_filename, top=top_filename)
    S2_bottom = md.compute_nematic_order(traj, indices=bottom_chains_list)
    S2_top = md.compute_nematic_order(traj, indices=top_chains_list)
    if average:
        S2_mean = np.mean([S2_bottom, S2_top], axis=0)
        np.savetxt(output_filename, np.column_stack((traj.time, S2_mean)))
    else:
        bottom_filename = output_filename.split('.')[0] + '-bottom.' + \
            output_filename.split('.')[1]
        np.savetxt(bottom_filename, np.column_stack((traj.time, S2_bottom)))
        top_filename = output_filename.split('.')[0] + '-top.' + \
            output_filename.split('.')[1]
        np.savetxt(top_filename, np.column_stack((traj.time, S2_top)))

def _tilt_angle(v1, v2):
    cosang = np.dot(v1, v2)
    sinang = np.linalg.norm(np.cross(v1, v2))
    angle = np.arctan2(sinang, cosang) * (180.0 / np.pi)
    if angle > 90 and angle < 180:
        angle -= 180
    return abs(angle)

def calc_avg_tilt_angle(traj_filename, top_filename, output_filename,
                        n_chains, average=True):
    """Calculate the average tilt angle of a monolayer

    Returns the average tilt angle at each frame of a trajectory,
    averaged between the top and bottom monolayers in a dual monolayer 
    system.

    NOTE: There is a lot of duplicate code between this function and the
          nematic order calculation. Should probably use a helper function.

    Notes
    -----
    Assumes a unique chain prototype.
    Assumes identical top and bottom monolayers.

    """
    topology = md.load(top_filename).topology
    atoms = np.array(list(topology.atoms))
    atom_names = [atom.name for atom in atoms]
    monolayer_begins = atom_names.index('C')
    half_atoms = int(len(atoms)/2)
    backfill_sites = int(atom_names.count('OS')/2 - n_chains)
    bottom_monolayer = [atom.index for atom in atoms[monolayer_begins:half_atoms - backfill_sites]]
    top_monolayer = [atom.index for atom in atoms[monolayer_begins + half_atoms:-backfill_sites]]
    bottom_chains = np.array_split(bottom_monolayer, n_chains)
    top_chains = np.array_split(top_monolayer, n_chains)

    traj = md.load(traj_filename, top=top_filename)
    chain_tilts_bottom = []
    for i, chain in enumerate(bottom_chains):
        chain_slice = traj.atom_slice(chain)
        directors = _compute_director(chain_slice)
        chain_tilts_bottom.append([_tilt_angle(director, [0.0, 0.0, 1.0])
            for director in directors])
    chain_tilts_bottom = np.asarray(chain_tilts_bottom)
    chain_tilts_bottom = np.transpose(chain_tilts_bottom)
    bottom_tilts_by_frame = np.mean(chain_tilts_bottom, axis=1)
    bottom_err_by_frame = np.std(chain_tilts_bottom, axis=1)

    chain_tilts_top = []
    for i, chain in enumerate(top_chains):
        chain_slice = traj.atom_slice(chain)
        directors = _compute_director(chain_slice)
        chain_tilts_top.append([_tilt_angle(director, [0.0, 0.0, -1.0])
            for director in directors])
    chain_tilts_top = np.asarray(chain_tilts_top)
    chain_tilts_top = np.transpose(chain_tilts_top)
    top_tilts_by_frame = np.mean(chain_tilts_top, axis=1)
    top_err_by_frame = np.std(chain_tilts_top, axis=1)

    if average:
        tilt_mean = np.mean([bottom_tilts_by_frame, top_tilts_by_frame], axis=0)
        np.savetxt(output_filename, np.column_stack((traj.time, tilt_mean)))
    else:
        bottom_filename = output_filename.split('.')[0] + '-bottom.' + \
            output_filename.split('.')[1]
        np.savetxt(bottom_filename, np.column_stack((traj.time, bottom_tilts_by_frame,
            bottom_err_by_frame)))
        top_filename = output_filename.split('.')[0] + '-top.' + \
            output_filename.split('.')[1]
        np.savetxt(top_filename, np.column_stack((traj.time, top_tilts_by_frame,
            top_err_by_frame)))

def count_hydrogen_bonds(traj_filename, top_filename, output_filename,
                         mol2_filename, ndx_filename=None, top_group=None,
                         bottom_group=None):
    """Count the number of inter- or intra-monolayer hydrogen bonds

    The number of inter- and intra-monolayer hydrogen bonds is determined for each
    frame in a trajectory using the Wernet-Nilsson method implemented in MDTraj. By
    default, atom indices in the first half of the system are considered to be part
    of the bottom monolayer and those in the second half are part of the top monolayer.
    This behavior can be overridden by providing names for `top_group` and
    `bottom_group`, where these groups will be read from an .ndx file and only indices
    in these groups are considered.

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
    ndx_filename : str, optional, default=None
        Name of Gromacs .ndx file which specifies atom groups. Required if `top_group`
        or `bottom_group` is not None.
    top_group : str, optional, default=None
        Only atom indices from this group (read from a .ndx file) will be considered
        as part of the top monolayer.
    bottom_group : str, optional, default=None
        Only atom indices from this group (read from a .ndx file) will be considered
        as part of the bottom monolayer.

    """
    top = md.load(top_filename).topology
    top_atoms = [atom for atom in top.atoms]
    half_atoms = int(len(top_atoms)/2)

    if ndx_filename:
        groups = read_ndx(ndx_filename)

    if top_group:
        top_monolayer = np.array(groups[top_group]) - 1
    else:
        top_monolayer = np.arange(half_atoms, 2*half_atoms)

    if bottom_group:
        bottom_monolayer = np.array(groups[bottom_group]) - 1
    else:
        bottom_monolayer = np.arange(half_atoms)
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
