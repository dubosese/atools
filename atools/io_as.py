import os
import numpy as np
import mdtraj as md

def xtc_unwrap_chains(xtc_file,top_file,chain_indices):
    traj = md.load(xtc_file,top=top_file)
    for frame in range(traj.n_frames):
        for chain in chain_indices:
            for atom in chain:
                print atom
                for dim in range(3):
                    box_length = traj[frame].unitcell_lengths[0,dim]*10.
                    traj[frame].xyz[0,atom,dim] *= (10.*box_length)
                    r = traj[frame].xyz[0,atom,dim] - traj[frame].xyz[0,chain[0],dim]
                    if r > box_length/2.:
                        traj[frame].xyz[0,atom,dim] -= box_length
                    elif r < -box_length/2.:
                        traj[frame].xyz[0,atom,dim] += box_length
    traj.save_xtc("test.xtc")

def write_lammps_data(filename,xyz,box,types=None,charges=None,mol_id=None,mass_dict=None,velocities=None,image_flags=None,bonds=None,bond_types=None,angles=None,angle_types=None,dihedrals=None,dihedral_types=None):
    if types is not None:
        types = np.asarray(types)
    else:
        types = np.full(len(xyz),1)
    if charges is not None:
        charges = np.asarray(charges)
    else:
        charges = np.zeros(len(xyz))
    if mol_id is not None:
        mol_id = np.asarray(mol_id)
    else:
        mol_id = np.zeros(len(xyz))
    if mass_dict is None:
        mass_dict = {1:1}
    if image_flags is not None:
        image_flags = np.asarray(image_flags)
    else:
        image_flags = np.zeros((len(xyz),3))
    out = open(filename,'w')
    out.write(filename+'\n\n')
    out.write(str(len(xyz))+' atoms\n')
    if bonds is not None:
        out.write(str(len(bonds))+' bonds\n')
    if angles is not None:
        out.write(str(len(angles))+' angles\n')
    if dihedrals is not None:
        out.write(str(len(dihedrals))+' dihedrals\n\n')
    out.write('{:d} atom types\n'.format(len(np.unique(types))))
    if bonds is not None:
        out.write('{:d} bond types\n'.format(len(np.unique(bond_types))))
    if angles is not None:
        out.write('{:d} angle types\n'.format(len(np.unique(angle_types))))
    if dihedrals is not None:
        out.write('{:d} dihedral types\n\n'.format(len(np.unique(dihedral_types))))
    for i,dim in enumerate(['x','y','z']):
        out.write('{0:.4f} {1:.4f} {2}lo {2}hi\n'.format(box.mins[i],box.maxs[i],dim))
    out.write('\nMasses\n\n')
    for i in range(np.max(types)):
        out.write('{:d}\t{:.4f}\n'.format(i+1,mass_dict[i+1]))
    out.write('\nAtoms\n\n')
    for i,coords in enumerate(xyz):
        out.write('{:d}\t{:.0f}\t{:d}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.0f}\t{:.0f}\t{:.0f}\n'.format(i+1,mol_id[i],types[i],charges[i],coords[0],coords[1],coords[2],image_flags[i,0],image_flags[i,1],image_flags[i,2]))
    if velocities is not None:
        out.write('\nVelocities\n\n')
        for i,vel in enumerate(velocities):
            out.write('{:d}\t{:.6f}\t{:.6f}\t{:.6f}\n'.format(i+1,vel[0],vel[1],vel[2]))
    if bonds is not None:
        out.write('\nBonds\n\n')
        for i,bond in enumerate(bonds):
            out.write('{:d}\t{:.0f}\t{:d}\t{:d}\n'.format(i+1,bond_types[i],bond[0],bond[1]))
    if angles is not None:
        out.write('\nAngles\n\n')
        for i,angle in enumerate(angles):
            out.write('{:d}\t{:.0f}\t{:d}\t{:d}\t{:d}\n'.format(i+1,angle_types[i],angle[0],angle[1],angle[2]))
    if dihedrals is not None:
        out.write('\nDihedrals\n\n')
        for i,dihedral in enumerate(dihedrals):
            out.write('{:d}\t{:.0f}\t{:d}\t{:d}\t{:d}\t{:d}\n'.format(i+1,dihedral_types[i],dihedral[0],dihedral[1],dihedral[2],dihedral[3]))

def read_topology(filename):
    """
    Currently supports the following directives:
        bonds

    Args:
        filename (str): name of TOP file to read in
    Returns:
        top (dict):
            'bonds': bonds (numpy.ndarray)
    """
    bonds = []
    with open(filename, 'r') as f:
        data_lines = f.readlines()

    directives = re.compile(r"""
        ((?P<bonds>\s*bonds))
        """, re.VERBOSE)

    i = 0
    while i < len(data_lines):
        match = directives.match(data_lines[i])
        if match:
            if match.group('bonds'):
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = map(int, data_lines.pop(i).split())
                    bond = [fields[0], fields[1]]
                    bonds.append(bond)
            else:
                i += 1
        else:
            i += 1

    top = {'bonds': np.asarray(bonds)}

    return top

def save_hoomdxml(compound, box, filename, traj, forcefield, force_overwrite=False):
    """ """
    from foyer.forcefield import apply_forcefield

    # Create separate file paths for .gro and .top
    filepath, filename = os.path.split(filename)
    basename = os.path.splitext(filename)[0]
    top_filename = os.path.join(filepath, basename + '.top')
    gro_filename = os.path.join(filepath, basename + '.gro')

    structure = compound.to_parmed()
    if forcefield:
        structure = apply_forcefield(structure, forcefield=forcefield)
    xml_writer(structure, box, filename)

def xml_writer(structure, box, filename):
    '''
    atoms: structure.atoms
    bonds: structure.bonds
    bond_types: structure.bond_types
    angles: structure.angles
    angle_types: structure.angle_types
    dihedrals: structure.dihedrals
    dihedral_types: structure.dihedral_types
    '''

    f = open(filename, 'w')
    f.write('<?xml version="1.3" encoding="UTF-8"?>\n')
    f.write('<hoomd xml>\n')
    f.write('<configuration time_step="0">\n')
    f.write('<box units="sigma"  Lx="{}" Ly="{}" Lz="{}"/>\n'.format(*box.lengths*10.))
    f.write('<position units="sigma" num="{}">\n'.format(len(structure.atoms)))
    for coord in structure.coordinates:
        f.write('{}\t{}\t{}\n'.format(*coord))
    f.write('</position>\n')

    f.write('<type>\n')
    for atom in structure.atoms:
        f.write('{}\n'.format(atom.type))
    f.write('</type>\n')

    f.write('<mass>\n')
    for atom in structure.atoms:
        f.write('{}\n'.format(atom.mass))
    f.write('</mass>\n')

    f.write('<charge>\n')
    for atom in structure.atoms:
        f.write('{}\n'.format(atom.charge))
    f.write('</charge>\n')

    f.write('<bond>\n')
    for bond in structure.bonds:
        f.write('bond{}\t{}\t{}\n'.format(bond.type.idx,bond.atom1.idx,bond.atom2.idx))
    f.write('</bond>\n')

    f.write('<angle>\n')
    for angle in structure.angles:
        f.write('angle{}\t{}\t{}\t{}\n'.format(angle.type.idx,angle.atom1.idx,angle.atom2.idx,angle.atom3.idx))
    f.write('</angle>\n')

    f.write('<dihedral>\n')
    for dihedral in structure.dihedrals:
        f.write('dihedral{}\t{}\t{}\t{}\t{}\n'.format(dihedral.type.idx,dihedral.atom1.idx,dihedral.atom2.idx,dihedral.atom3.idx,dihedral.atom4.idx))
    f.write('</dihedral>\n')

    f.write('</configuration>\n')
    f.write('</hoomd_xml>\n')

def read_lammps_data(data_file, verbose=False):
    """Reads a LAMMPS data file

    *** Only works for directives delimited by blank lines ***

    Currently supports the following directives:
        Masses
        Pair Coeffs (must be mix geometry)
        Bond Coeffs (must be harmonic)
        Angle Coeffs (must be harmonic)
        Dihedral Coeffs (must be OPLS)
        Atoms
        Bonds
        Angles
        Dihedrals

    Args:
        data_file (str): name of LAMMPS data file to read in
    Returns:
        lmp_data (dict):
            'xyz': xyz (numpy.ndarray)
            'types': types (numpy.ndarray)
            'masses': masses (numpy.ndarray)
            'bonds': bonds (numpy.ndarray)
            'angles': angles (numpy.ndarray)
            'dihedrals': dihedrals (numpy.ndarray)
            'pair_types': pair_types (dict)
            'bond_types': bond_types (dict)
            'angle_types': angle_types (dict)
            'dihedral_types': dihedral_type (dict)

        box (numpy.ndarray): box dimensions
    """
    bonds = np.empty(shape=(0, 3), dtype='int')
    angles = np.empty(shape=(0, 4), dtype='int')
    dihedrals = np.empty(shape=(0, 5), dtype='int')

    pair_types = dict()
    bond_types = dict()
    angle_types = dict()
    dihedral_types = dict()

    print "Reading '" + data_file + "'"
    with open(data_file, 'r') as f:
        data_lines = f.readlines()

    # TODO: improve robustness of xlo regex
    directives = re.compile(r"""
        ((?P<n_atoms>\s*\d+\s+atoms)
        |
        (?P<n_bonds>\s*\d+\s+bonds)
        |
        (?P<n_angles>\s*\d+\s+angles)
        |
        (?P<n_dihedrals>\s*\d+\s+dihedrals)
        |
        (?P<box>.+xlo)
        |
        (?P<Masses>\s*Masses)
        |
        (?P<PairCoeffs>\s*Pair\sCoeffs)
        |
        (?P<BondCoeffs>\s*Bond\sCoeffs)
        |
        (?P<AngleCoeffs>\s*Angle\sCoeffs)
        |
        (?P<DihedralCoeffs>\s*Dihedral\sCoeffs)
        |
        (?P<Atoms>\s*Atoms)
        |
        (?P<Velocities>\s*Velocities)
        |
        (?P<Bonds>\s*Bonds)
        |
        (?P<Angles>\s*Angles)
        |
        (?P<Dihedrals>\s*Dihedrals))
        """, re.VERBOSE)

    i = 0
    while i < len(data_lines):
        match = directives.match(data_lines[i])
        if match:
            if verbose:
                print(match.groups())

            elif match.group('n_atoms'):
                fields = data_lines.pop(i).split()
                n_atoms = int(fields[0])
                xyz = np.empty(shape=(n_atoms, 3))
                types = np.empty(shape=(n_atoms), dtype='int')
                masses = np.empty(shape=(n_atoms))

            elif match.group('n_bonds'):
                fields = data_lines.pop(i).split()
                bonds = np.empty(shape=(float(fields[0]), 3), dtype='int')

            elif match.group('n_angles'):
                fields = data_lines.pop(i).split()
                angles = np.empty(shape=(float(fields[0]), 4), dtype='int')

            elif match.group('n_dihedrals'):
                fields = data_lines.pop(i).split()
                dihedrals = np.empty(shape=(float(fields[0]), 5), dtype='int')

            elif match.group('box'):
                dims = np.zeros(shape=(3, 2))
                for j in range(3):
                    fields = map(float, data_lines.pop(i).split()[:2])
                    dims[j, 0] = fields[0]
                    dims[j, 1] = fields[1]
                box = Box(dims[:, 1], dims[:, 0])

            elif match.group('Masses'):
                if verbose:
                    print 'Parsing Masses...'
                data_lines.pop(i)
                data_lines.pop(i)

                mass_dict = dict()  # type:mass
                #     not end of file         not blank line
                while i < len(data_lines) and data_lines[i].strip():
                    fields = data_lines.pop(i).split()
                    mass_dict[int(fields[0])] = float(fields[1])

            elif match.group('Atoms'):
                if verbose:
                    print 'Parsing Atoms...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = data_lines.pop(i).split()
                    if len(fields) == 7:
                        a_id = int(fields[0])
                        types[a_id - 1] = int(fields[2])
                        masses[a_id - 1] = mass_dict[int(fields[2])]
                        charges[a_id - 1] = float(fields[3])
                        xyz[a_id - 1] = np.array([float(fields[4]),
                                             float(fields[5]),
                                             float(fields[6])])

                    if len(fields) == 10:
                        a_id = int(fields[0])
                        types[a_id - 1] = int(fields[2])
                        masses[a_id - 1] = mass_dict[int(fields[2])]
                        charges[a_id - 1] = float(fields[3])
                        xyz[a_id - 1] = np.array([float(fields[4]),
                                             float(fields[5]),
                                             float(fields[6])])
                        image_flags[a_id - 1] = np.array([float(fields[7]),
                                                    float(fields[8]),
                                                    float(fields[9])])

                    # non-official file format
                    if len(fields) == 8:
                        a_id = int(fields[0])
                        types[a_id - 1] = int(fields[1])
                        masses[a_id - 1] = mass_dict[int(fields[1])]
                        charges[a_id - 1] = 0.0
                        xyz[a_id - 1] = np.array([float(fields[2]),
                                             float(fields[3]),
                                             float(fields[4])])

            elif match.group('Velocities'):
                if verbose:
                    print 'Parsing Velocities...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = data_lines.pop(i).split()
                    a_id = int(fields[0])
                    velocities[a_id - 1] = np.array([float(fields[1]),
                                            float(fields[2]),
                                            float(fields[3])])

            elif match.group('PairCoeffs'):
                if verbose:
                    print 'Parsing Pair Coeffs...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = data_lines.pop(i).split()
                    pair_types[int(fields[0])] = (float(fields[1]),
                                                  float(fields[2]))
            elif match.group('BondCoeffs'):
                if verbose:
                    print 'Parsing Bond Coeffs...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = map(float, data_lines.pop(i).split())
                    bond_types[int(fields[0])] = fields[1:]

            elif match.group('AngleCoeffs'):
                if verbose:
                    print 'Parsing Angle Coeffs...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = map(float, data_lines.pop(i).split())
                    angle_types[int(fields[0])] = fields[1:]

            elif match.group('DihedralCoeffs'):
                if verbose:
                    print 'Parsing Dihedral Coeffs...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = map(float, data_lines.pop(i).split())
                    dihedral_types[int(fields[0])] = fields[1:]

            elif match.group('ImproperCoeffs'):
                if verbose:
                    print 'Parsing Improper Coeffs...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = map(float, data_lines.pop(i).split())
                    improper_types[int(fields[0])] = fields[1:]

            elif match.group('Bonds'):
                if verbose:
                    print 'Parsing Bonds...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = map(int, data_lines.pop(i).split())
                    bonds[fields[0] - 1] = fields[1:]

            elif match.group('Angles'):
                if verbose:
                    print 'Parsing Angles...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = map(int, data_lines.pop(i).split())
                    angles[fields[0] - 1] = fields[1:]

            elif match.group('Dihedrals'):
                if verbose:
                    print 'Parsing Dihedrals...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = map(int, data_lines.pop(i).split())
                    dihedrals[fields[0] - 1] = fields[1:]

            elif match.group('Impropers'):
                if verbose:
                    print 'Parsing Impropers...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = map(int, data_lines.pop(i).split())
                    impropers[fields[0] - 1] = fields[1:]

            else:
                i += 1
        else:
            i += 1

    lmp_data = {'xyz': xyz,
                'types': types,
                'masses': masses,
                'charges': charges,
                'velocities': velocities,
                'image_flags': image_flags,
                'bonds': bonds,
                'angles': angles,
                'dihedrals': dihedrals,
                'impropers': impropers,
                'pair_types': pair_types,
                'bond_types': bond_types,
                'angle_types': angle_types,
                'dihedral_types': dihedral_types,
                'improper_types': improper_types
                }
    return lmp_data, box
