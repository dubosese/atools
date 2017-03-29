from __future__ import division

import numpy as np

import mdtraj as md

def calc_nematic_order(traj_filename, top_filename, output_filename,
                       n_chains):
    """Calculate the nematic order of a monolayer

    Returns the average nematic order at each frame of a trajectory
    between the top and bottom monolayers in a dual monolayer system.

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
    S2_mean = np.mean([S2_bottom, S2_top], axis=0)
    np.savetxt(output_filename, np.column_stack((traj.time, S2_mean)))

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
