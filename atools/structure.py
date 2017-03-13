import networkx as nx

def find_angles(bond_graph):
    angles = []
    dihedrals = []
    for i,atom in enumerate(bond_graph.nodes()):
        for j,atom2 in enumerate(bond_graph.nodes()):
            if i < j:
                paths = nx.all_simple_paths(bond_graph,source=atom,target=atom2,cutoff=3)
                for path in paths:
                    if len(path) == 3:
                        angles.append(path)
                    elif len(path) == 4:
                        dihedrals.append(path)

    return angles,dihedrals

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
    for i, particle in monolayer.particles():
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

if __name__ == '__main__':
    from mbuild.examples.alkane.alkane import Alkane

    alkane = Alkane(n=3)
    G = alkane.bond_graph
    angles,dihedrals = find_angles(G)
    print len(angles),len(dihedrals)
