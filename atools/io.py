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
