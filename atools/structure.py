from __future__ import division

import numpy as np

import mdtraj as md

def calc_nematic_order(traj_filename, top_filename, output_filename):
    first = md.load_frame(traj_filename, 0, top=top_filename)
    atoms = list(first.top.atoms)

    count = 0
    chain_indices_bot = []
    chain_indices_top = []
    for i, atom in enumerate(atoms):
        if atom.name == ""


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
