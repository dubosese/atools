from __future__ import division

import mdtraj as md
from mdtraj.geometry.order import _compute_director
import numpy as np


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
                        n_chains):
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

    time = np.array([traj[frame].time[0] for frame in range(traj.n_frames)])
    bottom_filename = output_filename.split('.')[0] + '-bottom.' + \
        output_filename.split('.')[1]
    np.savetxt(bottom_filename, np.column_stack((time, bottom_tilts_by_frame,
        bottom_err_by_frame)))
    top_filename = output_filename.split('.')[0] + '-top.' + \
        output_filename.split('.')[1]
    np.savetxt(top_filename, np.column_stack((time, top_tilts_by_frame,
        top_err_by_frame)))

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
