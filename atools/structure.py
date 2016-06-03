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

if __name__ == '__main__':
    from mbuild.examples.alkane.alkane import Alkane

    alkane = Alkane(n=3)
    G = alkane.bond_graph
    angles,dihedrals = find_angles(G)
    print len(angles),len(dihedrals)
