import textwrap

def write_monolayer_ndx(rigid_groups, filename):
    with open(filename, 'w') as f:
        for name, indices in rigid_groups.items():
            f.write('[ {name} ]\n'.format(name=name))
            atoms = '{}\n'.format(' '.join(str(x) for x in indices))
            f.write(textwrap.fill(atoms, 80))
            f.write('\n')
